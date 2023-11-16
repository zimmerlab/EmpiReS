package empires.rnaseq.reads;

import empires.rnaseq.GenomeSequenceExtractor;
import lmu.utils.*;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;

import java.io.File;
import java.io.PrintWriter;
import java.util.*;
import java.util.function.Consumer;

import static lmu.utils.ObjectGetter.*;

public class BAMMappingCheck {

    static HashSet<String> getChromosomes(File bamfile) {
        SAMFileReader sam_reader  = new SAMFileReader(bamfile, false);
        return mapToSet(sam_reader.getFileHeader().getSequenceDictionary().getSequences(),
                (_s) -> _s.getSequenceName());
    }

    public static void readBAM(File bamfile, Consumer<UPair<SAMRecord>> handler, String selectedChr) {
        long t1 = System.currentTimeMillis();
        Logger log = LogConfig.getLogger();
        SAMFileReader sam_reader  = new SAMFileReader(bamfile, new File(bamfile + ".bai"), false);
        sam_reader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
        Iterator<SAMRecord> bamIt = null;
        if(selectedChr != null) {
            int length = sam_reader.getFileHeader().getSequence(selectedChr).getSequenceLength();
            log.info("%s got length: %d", selectedChr, length);
            bamIt = sam_reader.query(selectedChr, 0, length, true);
        } else {
            bamIt = sam_reader.iterator();
        }



        HashMap<Integer, SAMRecord> open = new HashMap<>();
        String lastchr = null;

        while(bamIt.hasNext()) {
            SAMRecord sr = bamIt.next();

            if(sr.getReadUnmappedFlag() || sr.getMateUnmappedFlag() || sr.getNotPrimaryAlignmentFlag())
                continue;

            if(lastchr == null ||  !sr.getReferenceName().equals(lastchr)) {

                if(lastchr != null) {
                    log.info("chr %s dropped %d partly mapped reads", lastchr, open.size());
                }
                open.clear();
                lastchr = sr.getReferenceName();
            }

            int readId = Integer.parseInt(sr.getReadName());
            if(!open.containsKey(readId)) {
                open.put(readId, sr);
                continue;
            }

            SAMRecord other = open.remove(readId);
            handler.accept(UPair.createU(other.getFirstOfPairFlag() ? other : sr, other.getFirstOfPairFlag() ? sr : other));
        }

        sam_reader.close();
        long t2 = System.currentTimeMillis();
        log.info("parsing through %s took %.2f sec", bamfile.getName(), (t2 - t1) / 1000.0);
    }

    static String getGenomeSeq(empires.rnaseq.GenomeSequenceExtractor gex, String chr, RegionVector rv) {
        StringBuilder sb = new StringBuilder();
        for(Region1D r : rv.getRegions()) {
            sb.append(gex.getGenomeSeq(chr, r.getX1(), r.getX2()));
        }
        return sb.toString();
    }

    static void writeBitsStringSorted(File of, BitSet bs) {
        Vector<String> strings = new Vector<>();
        for(int i=0; i<bs.size(); i++) {
            if(!bs.get(i))
                continue;

            strings.add(""+i);
        }
        Collections.sort(strings);
        PrintWriter pw = FileUtils.getWriter(of);
        apply(strings, (_s) -> pw.println(_s));
        pw.close();
    }
    public static void main(String[] args) {
        SimpleOptionParser cmd = new SimpleOptionParser("info", "bam", "od", "firstN", "total", "genome", "genomeidx", "prefix", "detailerrorperchr");
        cmd.setFile("info", "bam", "genome", "genomeidx");
        cmd.setDir("od");
        cmd.setInt("firstN");
        cmd.setInt("total", "detailerrorperchr");
        cmd.setDefault("firstN", "-1");
        cmd.setDefault("detailerrorperchr", "1000");

        cmd.setDefault("total", "-1");
        cmd.setOptional("genome", "genomeidx", "prefix");
        if(!OptionParser.parseParams(args, true, false, true, cmd))
            return;

        HashMap<Integer, Pair<String, UPair<RegionVector>>> reference = new HashMap<>();

        empires.rnaseq.GenomeSequenceExtractor gex = (cmd.isOptionSet("genome") && cmd.isOptionSet("genomeidx")) ? new GenomeSequenceExtractor(cmd.getFile("genome"), cmd.getFile("genomeidx")) : null;
        Logger log = LogConfig.getLogger();
        File info = cmd.getFile("info");
        int numSimulated = cmd.getInt("total");
        if(numSimulated < 0) {
            log.info("num simulated not provided, read info file");
            HashSet<Integer> simulated = new HashSet<>();
            apply(ReadReferenceInfo.getIterator(info), (_rri) -> simulated.add(_rri.readid));
            log.info("got %d simulated reads", simulated.size());
            numSimulated = simulated.size();
        }

        int NUM_DETAIL_ERRORS_PER_CHR = cmd.getInt("detailerrorperchr");
        File bamfile = cmd.getFile("bam");

        HashSet<String> chromosomes = getChromosomes(bamfile);

        log.info("got chromosomes: %s\n", chromosomes);
        HashSet<Integer> mapped = new HashSet<>();

        String prefix = cmd.getOptionalValue("prefix", null);
        if(prefix == null) {
            prefix = bamfile.getName().split("\\.bam")[0];
        }
        File od = cmd.getFile("od");

        log.info("will write error/missing/partly mappings with prefix: %s", prefix);


        class Statistic {
            BitSet simulated;
            BitSet fullOK;
            BitSet mapped;
            BitSet wrongChrMapped;
            BitSet errorMapping;
            BitSet partlyMapping;

            Statistic(int numSimulated) {
                simulated = new BitSet(numSimulated);
                fullOK = new BitSet(numSimulated);
                mapped = new BitSet(numSimulated);
                wrongChrMapped = new BitSet(numSimulated);
                errorMapping = new BitSet(numSimulated);
                partlyMapping = new BitSet(numSimulated);
            }
            public String getInfo(String prefix) {
                double simul_norm = 100.0 / simulated.cardinality();

                BitSet missing =  SetInfo.minus(simulated, mapped);
                BitSet correctchr =  SetInfo.intersect(simulated, mapped);
                double mapped_norm = 100.0 / correctchr.cardinality();
                return String.format("%s got %d/%d mapped (%.2f%%) missing: %d (%.2f%%) mapped full ok: %d (%.2f%%) partly mapped %d (%.2f%%) wrong chr: %d (%.2f%%) err: %d: (%.2f%%)",
                        prefix,
                        correctchr.cardinality(),  simulated.cardinality(),
                        correctchr.cardinality() * simul_norm,
                        missing.cardinality(), missing.cardinality() * simul_norm,
                        fullOK.cardinality(), fullOK.cardinality() * mapped_norm,
                        partlyMapping.cardinality(), partlyMapping.cardinality() * mapped_norm,
                        wrongChrMapped.cardinality(),  wrongChrMapped.cardinality() * (100.0 / mapped.size()),
                        errorMapping.cardinality(),  errorMapping.cardinality() * mapped_norm);

            }

            public void update(Statistic statistic) {
                simulated = SetInfo.union(simulated, statistic.simulated);
                fullOK = SetInfo.union(fullOK, statistic.fullOK);
                mapped = SetInfo.union(mapped, statistic.mapped);
                wrongChrMapped = SetInfo.union(wrongChrMapped, statistic.wrongChrMapped);
                errorMapping = SetInfo.union(errorMapping, statistic.errorMapping);
                partlyMapping = SetInfo.union(partlyMapping, statistic.partlyMapping);
            }
        }

        Statistic globalStatistic = new Statistic(numSimulated);

        for(String chr : chromosomes) {

            log.info("check chromosome: %s", chr);
            File chrref = new File(info+"."+chr);
            if(!chrref.exists()) {
                log.info("chr-splitted ref missing: %s", chrref.getAbsolutePath());
                continue;
            }
            HashMap<Integer, ReadReferenceInfo> chrSimulated = new HashMap<>();
            filter_and_apply(ReadReferenceInfo.getIterator(chrref), (_rri) -> _rri.chr.equals(chr), (_rri) -> chrSimulated.put(_rri.readid, _rri));
            log.info("chr: %s got %d simulated reads", chr, chrSimulated.size());


            Statistic chrStatistic = new Statistic(numSimulated);
            apply(chrSimulated.keySet(), (_s) -> chrStatistic.simulated.set(_s));

            HashMap<UPair<RegionVector>, Integer> wrongExamples = new HashMap<>();
            HashMap<UPair<RegionVector>, UPair<RegionVector>> wrongExampleDetails = new HashMap<>();


            readBAM(bamfile, (_p) -> {
                        int readId = Integer.parseInt(_p.getFirst().getReadName());
                        chrStatistic.mapped.set(readId);
                        ReadReferenceInfo ref = chrSimulated.get(readId);
                        if(ref == null) {
                            chrStatistic.wrongChrMapped.set(readId);
                            return;
                        }

                        RegionVector fw = BAMIterator.toRV(_p.getFirst());
                        boolean ok = true;
                        boolean sub = ok;
                        if(!ref.fw.equals(fw)) {
                            if(!ref.fw.isSubVector(fw)) {
                                sub = false;
                            }


                            ok = false;
                        }
                        RegionVector rw = BAMIterator.toRV(_p.getSecond());
                        if(!ref.rw.equals(rw)) {
                            if(!ref.rw.isSubVector(rw)) {
                                sub = false;
                            }
                            ok = false;
                        }
                        if(ok) {
                            chrStatistic.fullOK.set(readId);
                            return;
                        }


                        if(sub) {
                            chrStatistic.partlyMapping.set(readId);
                            return;
                        }
                        chrStatistic.errorMapping.set(readId);

                        UPair<RegionVector> key = UPair.createU(ref.fw, ref.rw);
                        if(wrongExamples.containsKey(key)) {
                            MapBuilder.update(wrongExamples, key, 1);
                            return;
                        }
                        if(wrongExamples.size() >= NUM_DETAIL_ERRORS_PER_CHR)
                            return;

                        wrongExamples.put(key, 1);
                        wrongExampleDetails.put(key, UPair.createU(fw, rw));
                        System.out.printf("chr:%s err: %d %d FW ref: %s check: %s sub: %s\n", chr, wrongExamples.size(), readId, ref.fw, fw, ref.fw.isSubVector(fw));
                        System.out.printf("chr:%s err: %d %d RW ref: %s check: %s sub: %s\n", chr, wrongExamples.size(), readId, ref.rw, fw, ref.rw.isSubVector(rw));
                        if(gex  != null) {
                            System.out.printf("simulated fw: %s\nmapped fw: %s\n", getGenomeSeq(gex, chr, ref.fw), getGenomeSeq(gex, chr, fw));
                            System.out.printf("simulated rw: %s\nmapped rw: %s\n", getGenomeSeq(gex, chr, ref.rw), getGenomeSeq(gex, chr, rw));
                        }
                    },
                    chr);

            log.info(chrStatistic.getInfo("chr: " + chr));
            globalStatistic.update(chrStatistic);
        }

        log.info(globalStatistic.getInfo("GLOBAL:" ));

        BitSet allsimul = new BitSet(numSimulated);
        for(int i=0; i<numSimulated; i++) {
            allsimul.set(i);
        }

        BitSet unmapped = SetInfo.minus(allsimul, globalStatistic.mapped);

        writeBitsStringSorted(new File(od, prefix+".unmapped.readids"), unmapped);
        writeBitsStringSorted(new File(od, prefix+".wrong.readids"), globalStatistic.errorMapping);
        writeBitsStringSorted(new File(od, prefix+".partially_mapped.readids"), globalStatistic.partlyMapping);
        writeBitsStringSorted(new File(od, prefix+".wrongchr.readids"), globalStatistic.wrongChrMapped);

    }

}
