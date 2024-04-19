package nlEmpiRe.test.rnaseq;

import lmu.utils.*;
import nlEmpiRe.rnaseq.*;
import nlEmpiRe.rnaseq.simulation.SplicingSimulation;


import java.io.File;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Vector;

import static lmu.utils.IteratorUtils.rangev;
import static lmu.utils.ObjectGetter.*;

public class FastQGenerator {
    GenomeSequenceExtractor gex;
    IsoformRegionGetter isoformRegionGetter;
    File fastqOutDir = null;
    SplicingSimulation.ReadInReplicateConsumer readInReplicateConsumer;
    HashMap<Pair<String, Integer>, PrintWriter[]> fastqAndInfoMap = new HashMap<>();
    HashMap<Pair<String, Integer>, BAMGenerator> bamInfoMap = new HashMap<>();


    HashMap<Pair<String, Integer>, Integer> fastqCount = new HashMap<>();

    boolean writeInfoFile = true;

    double mutation_rate = 0.01;

    final static String QUALSTRING = StringUtils.multiply("IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
                50);


    static Vector<Character> DEFAULT_MISSING = toVector('N');
    static HashMap<Character, Vector<Character>> mutmap = new HashMap<>();

    static {
        Vector<Character> tcga = ObjectGetter.toVector('A', 'T', 'C', 'G');
        HashSet<Character> h = toSet(tcga);
        for (Character c : tcga) {
            mutmap.put(c, ObjectGetter.filter(tcga, (_x) -> _x != c));
        }

    }

    public void setWriteInfoFile(boolean b) {
        writeInfoFile = b;
    }

    public FastQGenerator(IsoformRegionGetter isoformRegionGetter, File genome_fasta, File genome_idx, double mutation_rate) {
        this.isoformRegionGetter = isoformRegionGetter;
        this.mutation_rate = mutation_rate;
        setGenomeSequenceInfo(genome_fasta, genome_idx);
    }

    public void setGenomeSequenceInfo(File genome_fasta, File genome_idx) {
        gex = null;
        if (genome_fasta == null || genome_idx == null)
            return;

        gex = new GenomeSequenceExtractor(genome_fasta, genome_idx);
    }

    BAMGenerator getBamWriter(String cond, int replicateId) {
        if(!generate_bams)
            return null;

        Pair<String, Integer> key = Pair.create(cond, replicateId);
        BAMGenerator generator = bamInfoMap.get(key);
        if(generator != null)
            return generator;

        File bamOutFile = new File(fastqOutDir, String.format("%s_%02d.bam", cond, replicateId));
        generator = new BAMGenerator(bamfileReference, gex, bamOutFile);
        bamInfoMap.put(key, generator);
        return generator;
    }

    public static String getReplicatePrefix(String condition, int replicateId) {
        return String.format("%s_%02d", condition, replicateId);
    }

    public static String getFastqFileName(String condition, int replicateId, boolean fw) {
        return  String.format("%s_%d.fastq.gz", getReplicatePrefix(condition, replicateId), (fw) ? 1 : 2);
    }

    synchronized PrintWriter[] getFastQAndInfoWriter(String cond, int replicateId) {
        Pair<String, Integer> key = Pair.create(cond, replicateId);
        PrintWriter[] pw = fastqAndInfoMap.get(key);
        if (pw != null)
            return pw;

        File fw_outfile = new File(fastqOutDir, getFastqFileName(cond, replicateId, true));
        File rw_outfile = new File(fastqOutDir, getFastqFileName(cond, replicateId, false));
        File infoFile = new File(fastqOutDir, getReplicatePrefix(cond, replicateId)+".info");

        pw = new PrintWriter[]{
                FileUtils.getWriter(fw_outfile, false, true),
                FileUtils.getWriter(rw_outfile, false, true),
                (!writeInfoFile) ? null :initMappingInfo(infoFile)
        };

        fastqAndInfoMap.put(key, pw);
        return pw;
    }


    public void close() {
        for (PrintWriter[] pws : fastqAndInfoMap.values()) {
            for (int i = 0; i < pws.length; i++) {
                PrintWriter pw = pws[i];
                if (pw == null)
                    continue;
                pw.close();
            }
        }
        apply(bamInfoMap.values(), (_b) -> _b.close());
    }

    boolean checkSeq(String gene, String tr, String readseq) {
        if(0 == filteredSize(rangev(readseq.length()), (_i) -> !mutmap.containsKey(readseq.charAt(_i))))
            return true;

        throw new FRuntimeException("invalid character for chr: %s %s:%s readseq: %s\n", isoformRegionGetter.getRegionById(gene).chr, gene, tr, readseq);

    }

    boolean generate_bams = false;
    File bamfileReference;

    public void setGenerateBams(boolean b) {
        generate_bams = b;
    }

    public void setBamfileReference(File bamfileReference) {
        this.bamfileReference = bamfileReference;
        if(bamfileReference != null) {
            setGenerateBams(true);
        }
    }

    public void setFastQOutDir(File fastQOutDir) {
        readInReplicateConsumer = null;

        if (fastQOutDir == null)
            return;

        this.fastqOutDir = fastQOutDir;
        readInReplicateConsumer = (_cond, _repid, _read) ->
        {
            PrintWriter[] pws = getFastQAndInfoWriter(_cond, _repid);
            String trseq = getTrSeq(_read.gene, _read.transcriptId);



            Region1D fw_reg = _read.tr_fw.getRegion(0);
            Region1D rw_reg = _read.tr_rw.getRegion(0);

            String fw_read = trseq.substring(fw_reg.getX1(), fw_reg.getX2());
            String rw_read = GenomicUtils.reverse_complement(trseq.substring(rw_reg.getX1(), rw_reg.getX2()));


            BAMGenerator bamGenerator  = getBamWriter(_cond, _repid);

            Pair<String, Integer> key = Pair.create(_cond, _repid);
            int readId = fastqCount.getOrDefault(key, 0) + 1;
            if(bamGenerator != null) {
                MultiIsoformRegion gene = isoformRegionGetter.getRegionById(_read.gene);
                bamGenerator.writeRecords(""+readId, fw_read, rw_read, QUALSTRING.substring(0 , fw_read.length()), gene.chr, gene.strand, _read.fw, _read.rw);
            }
            fastqCount.put(key, readId);


            writeRecord(pws[0], readId, fw_read);
            writeRecord(pws[1], readId, rw_read);
            writeMappingInfo(pws[2], readId, _read);

        };

    }

    HashMap<String, String> transcriptSeqCache = new HashMap<>();

    String getTrSeq(String gene, String transcript) {
        String cached = transcriptSeqCache.get(transcript);
        if (cached != null)
            return cached;

        MultiIsoformRegion mir = isoformRegionGetter.getRegionById(gene);
        RegionVector trv = mir.isoforms.get(transcript);
        try {
            String tseq = GenomicUtils.getSplicedSeq(mir.chr, mir.strand, trv, gex, true);

            //checkSeq(mir.id, transcript, tseq);

            transcriptSeqCache.put(transcript, tseq);
            return tseq;
        } catch (Exception e) {
            throw new FRuntimeException("error:%s at getting transcript seq for %s:%s chr: %s strand: %s trv: %s gex: %s", gene, e, e.getMessage(), transcript, mir.chr, mir.strand, gex);
        }




    }

    static String mutate(String s, double rate) {
        Vector<Double> props = ObjectGetter.mapIndex(s.length(), (_i) -> Math.random());
        Vector<Integer> mutations = ObjectGetter.filter(rangev(props), (_i) -> props.get(_i) < rate);
        if (mutations.size() == 0)
            return s;

        StringBuffer sb = new StringBuffer();
        sb.append(s);

        for (int i : mutations) {
            Vector<Character> vc = mutmap.getOrDefault(s.charAt(i), DEFAULT_MISSING);
            if(vc == null) {
                return null;
            }
            int ridx = (int) (Math.random() * vc.size());
            sb.setCharAt(i, vc.get(ridx));
        }
        return sb.toString();
    }


    void writeRecord(PrintWriter pw, int rid, String seq) {
        String mutated = mutate(seq, mutation_rate);
        if(mutated == null)
            throw new FRuntimeException("mutate null for %s", seq);
        pw.printf("@%d\n%s\n+%d\n%s\n", rid, mutated, rid, QUALSTRING.substring(0, mutated.length()));
    }


    public PrintWriter initMappingInfo(File infoFile) {
        PrintWriter mapping = FileUtils.getWriter(infoFile);
        mapping.printf("readid\tchr\tgene\ttranscript\tt_fw_regvec\tt_rw_regvec\tfw_regvec\trw_regvec\n");
        return mapping;
    }

    public void writeMappingInfo(PrintWriter mapping, int readId, nlEmpiRe.rnaseq.simulation.SimulatedRead read) {
        if(mapping == null)
            return;

        MultiIsoformRegion mir = isoformRegionGetter.getRegionById(read.gene);

        mapping.printf("%d\t%s\t%s\t%s", readId, mir.chr, mir.id, read.transcriptId);
        mapping.printf("\t%s\t%s", read.tr_fw.getSimpleRepresentation(), read.tr_rw.getSimpleRepresentation());
        mapping.printf("\t%s\t%s\n", read.fw.getSimpleRepresentation(), read.rw.getSimpleRepresentation());
    }

}