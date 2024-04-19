package nlEmpiRe.rnaseq.mapping;

import lmu.utils.*;
import nlEmpiRe.rnaseq.*;
import nlEmpiRe.rnaseq.reads.*;

import java.io.*;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;

import static lmu.utils.ObjectGetter.*;
import org.apache.logging.log4j.Logger;

public class TranscriptInfoBasedGenomicMapper {

    Logger log = LogConfig.getLogger();
    List<TranscriptInfo> transcriptInfoList;
    SuffixArray transcriptomeSA;
    int[] startposlist;
    String saInput;

    HashMap<String, TranscriptInfo> trid2info;
    boolean verbose = false;

    HashMap<String, SuffixArray> tr2suffixArray = new HashMap<>();
    int maxMismatches = 10;
    static int MIN_MMP_LENGTH = 10;

    public TranscriptInfoBasedGenomicMapper(List<TranscriptInfo> transcriptInfoList, StringBuilder sb) {
        log.info("start building SA from %.1f Mio positions", sb.length() / 1_000_000.0);
        long t1 = System.currentTimeMillis();
        transcriptomeSA = new SuffixArray(saInput = sb.toString());
        transcriptomeSA.calcMidPoints();

        long t2 = System.currentTimeMillis();
        log.info("SA build and LCP took %.2f sec", (t2 - t1) / 1000.0);

        this.transcriptInfoList = transcriptInfoList;
        startposlist = new int[transcriptInfoList.size()];
        for(int i=0; i<transcriptInfoList.size(); i++) {
            startposlist[i] = transcriptInfoList.get(i).strStart;

        }

        trid2info = buildReverseMap(transcriptInfoList, (_t) -> _t.transcriptId);
    }

    Integer num_hashes_calced = 0;

    public TranscriptInfoBasedGenomicMapper(File f) {
        try
        {
            long t1 = System.currentTimeMillis();
            DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(f)));
            this.transcriptomeSA = SuffixArray.readSuffixArray(dis);
            saInput = FileUtils.readString(dis);
            int numTrInfos = dis.readInt();
            transcriptInfoList = new ArrayList<>();
            startposlist = new int[numTrInfos];
            for(int i=0; i<numTrInfos; i++) {
                TranscriptInfo tinfo = TranscriptInfo.read(dis);
                transcriptInfoList.add(tinfo);
                startposlist[i] = tinfo.strStart;
            }
            dis.close();
            trid2info = buildReverseMap(transcriptInfoList, (_t) -> _t.transcriptId);
            long t2 = System.currentTimeMillis();
            log.info("read SA and %d infos from %s took: %.2f", transcriptInfoList.size(), f.getName(), (t2 - t1) / 1_000.0);

            //apply(transcriptInfoList, (_t) -> _t.calcHashes(saInput.substring(_t.strStart, _t.strStart + _t.length)));
            //long t3 = System.currentTimeMillis();
            //log.info("calc hashes per transcript took: %.2f", (t3 - t2) / 1_000.0);
        }catch (IOException ie) {
            throw new FRuntimeException("i/o error: %s while reading transcript based genomic info from: %s", ie, ie.getMessage(), f.getAbsolutePath());
        }


    }

    public String getTranscriptSequence(String trid) {
        TranscriptInfo tinfo = trid2info.get(trid);
        return saInput.substring(tinfo.strStart, tinfo.strStart + tinfo.length);
    }

    public SuffixArray getTranscriptSuffixArray(String trid) {
        synchronized (tr2suffixArray) {
            SuffixArray sa = tr2suffixArray.get(trid);
            if(sa != null)
                return sa;

            tr2suffixArray.put(trid, sa = new SuffixArray(getTranscriptSequence(trid)));
            sa.calcMidPoints();
            return sa;
        }
    }
    public void serialize(File f) {


        try {

            BufferedOutputStream bos = new BufferedOutputStream(new FileOutputStream(f.getAbsolutePath()));
            DataOutputStream dos = new DataOutputStream(bos);

            transcriptomeSA.serialize(dos);
            FileUtils.writeString(saInput, dos);
            dos.writeInt(transcriptInfoList.size());
            for(TranscriptInfo ti : transcriptInfoList) {
                ti.write(dos);
            }
            dos.close();
            bos.close();

        } catch(IOException ie) {
            throw new FRuntimeException("i/o error: %s while writing suffix array file: %s", ie, ie.getMessage(), f.getAbsolutePath());
        }

    }



    HashMap<String, Vector<Integer>> getPartialHits(String seq, boolean fw_read, boolean reversecomplement) {

        List<int[]> longest_hits = new ArrayList<>();

        int hit_start = 0;
        int maxlength = 0;
        while(hit_start + MIN_MMP_LENGTH < seq.length()) {
            int[] hit = transcriptomeSA.search_longest_prefix(seq, hit_start, false);
            maxlength = Math.max(maxlength, hit[2]);
            longest_hits.add(new int[]{hit_start, hit[2]});
            hit_start += hit[2] + 1;
        }

        int MINLENGTH = Math.max(MIN_MMP_LENGTH, (int)(0.8 * maxlength));

        HashMap<Pair<String, Integer>, Integer> pos2matches = new HashMap<>();

        int maxmatches = 0;

        for(int[] hit : longest_hits) {
            if(hit[1] < MINLENGTH)
                continue;

            int start = hit[0];
            Vector<Integer> hitposlist = transcriptomeSA.getHitPoslist(seq.substring(start, start + hit[1]));

            for(int hp : hitposlist) {
                int hitIdx = Arrays.binarySearch(startposlist, hp);
                if(hitIdx < 0) {
                    hitIdx = - hitIdx - 1;
                    hitIdx -= 1;
                }

                TranscriptInfo thit = transcriptInfoList.get(hitIdx);

                int hitstart = hp - thit.strStart;

                if(start > hitstart)
                    continue; // can not be part

                int trstart = hitstart - start;
                if(trstart + seq.length() > thit.length)
                    continue;

                Pair<String, Integer> k = Pair.create(thit.transcriptId, trstart);

                MapBuilder.update(pos2matches, k, hit[1]);

                maxmatches = Math.max(maxmatches, pos2matches.get(k));
            }
        }

        HashMap<String, Vector<Integer>> tr2poslist = new HashMap<>();

        int accepted_maxmatch = (int)(0.8 * maxmatches);
        for(Map.Entry<Pair<String, Integer>, Integer> e : pos2matches.entrySet()) {
            if(e.getValue() < accepted_maxmatch)
                continue;

            MapBuilder.updateV(tr2poslist, e.getKey().getFirst(), e.getKey().getSecond());
            //hits.add(new PartialHit(trid2info.get(e.getKey().getFirst()),  0, e.getKey().getSecond(), fw_read, reversecomplement));
        }

        return tr2poslist;
    }

    String toRevComp(String s) {
        return (s == null) ? s : GenomicUtils.reverse_complement(s);
    }


    static class ReadMatchInfo {
        String readseq;
        int startPos;
        int length;
        List<Integer> mismatches;

        ReadMatchInfo(int startPos, int length) {
            this.startPos = startPos; this.length = length;
        }

        public int getNumMismatches(){
            return (mismatches == null) ? 0 : mismatches.size();
        }

        void setMismatch(int readPos) {
            if(mismatches == null) {
                mismatches = new ArrayList<>();
            }
            mismatches.add(readPos);
        }

        public String toString() {
            return String.format("%d-%d mm:%d%s",
                    startPos, startPos + length, getNumMismatches(), (mismatches == null) ? "" : " " + mismatches);
        }
    }

    ReadMatchInfo getReadMatchInfo(String read, TranscriptInfo tinfo, int hitstart) {
        ReadMatchInfo rmi = new ReadMatchInfo(tinfo.strStart + hitstart, read.length());
        rmi.readseq = read;
        String tr_fragment = saInput.substring(rmi.startPos, rmi.startPos + read.length());

        for(int i=0; i<read.length(); i++) {
            if(read.charAt(i) == tr_fragment.charAt(i))
                continue;
            rmi.setMismatch(i);
        }
        return rmi;
    }


    public Collection<TranscriptHit> getTranscriptHits(String readid, HashMap<UPair<Boolean>, String> readseqs , Boolean strandness) {

        List<TranscriptHit> transcriptHits = new ArrayList<>();

        for(int strandness_fw = 0; strandness_fw < 2; strandness_fw++) {

            if(strandness != null) {
                if(strandness != (strandness_fw == 0))
                    continue;
            }

            String fw_read = (strandness_fw == 0) ? readseqs.get(UPair.createU(true, false)) : readseqs.get(UPair.createU(false, false));
            String rw_read = (strandness_fw == 0) ? readseqs.get(UPair.createU(false, true)) : readseqs.get(UPair.createU(true, true));

            //Evi: inserted check if read is null for fw_hits and rw_hits
            HashMap<String, Vector<Integer>> fw_hits = fw_read==null ? new HashMap<>() : getPartialHits(fw_read, true, false);
            HashMap<String, Vector<Integer>> rw_hits = rw_read==null ? new HashMap<>() : getPartialHits(rw_read, false, false);

            Set<String> smaller = (fw_hits.size() <= rw_hits.size()) ? fw_hits.keySet() : rw_hits.keySet();
            Set<String> bigger = (smaller == fw_hits.keySet()) ? rw_hits.keySet() : fw_hits.keySet();

            if (fw_read == null ||  rw_read == null) {
                //single read
                String s_read = (fw_read == null) ? rw_read : fw_read;
                HashMap<String, Vector<Integer>> single_hits = (fw_hits.size() > 0) ? fw_hits : rw_hits;
                for (Map.Entry<String, Vector<Integer> > e : single_hits.entrySet()) {
                    TranscriptInfo tinfo = trid2info.get(e.getKey());
                    TranscriptHit bestHit = new TranscriptHit(null, 0, 1, maxMismatches + 1);

                    for(int pos: e.getValue()) {
                        ReadMatchInfo rmi = getReadMatchInfo(s_read, tinfo, pos);
                        if( rmi.getNumMismatches() >= bestHit.numMismatches)
                            continue;

                        bestHit = new TranscriptHit(tinfo, rmi.startPos, rmi.startPos + rmi.length, rmi.getNumMismatches());
                    }
                    if(bestHit.transcriptInfo == null)
                        continue;

                    transcriptHits.add(bestHit);
                }

                return transcriptHits;
            }

            for(String trid : smaller) {
                if(!bigger.contains(trid))
                    continue;

                TranscriptInfo tinfo = trid2info.get(trid);
                TranscriptHit bestHit = new TranscriptHit(null, 0, 1, maxMismatches + 1);
                TreeMap<Integer, ReadMatchInfo> fw2hit = new TreeMap<>();
                for(int pos : fw_hits.get(trid)) {
                    fw2hit.put(pos, getReadMatchInfo(fw_read, tinfo, pos));
                }

                for(int pos : rw_hits.get(trid)) {
                    ReadMatchInfo rmi = getReadMatchInfo(rw_read, tinfo, pos);
                    for(int fw_pos : fw2hit.keySet()) {
                        if(fw_pos > pos)
                            break;

                        ReadMatchInfo fw_rmi = fw2hit.get(fw_pos);
                        int sumMM = rmi.getNumMismatches() + fw_rmi.getNumMismatches();
                        if(sumMM >= bestHit.numMismatches)
                            continue;

                        bestHit = new TranscriptHit(tinfo, fw_pos, pos + rw_read.length(), sumMM);
                    }
                }

                if(bestHit.transcriptInfo == null)
                    continue;

                transcriptHits.add(bestHit);
            }

        }

        return transcriptHits;
    }

    static class MapperJob extends Thread{
        String label;
        HashMap<String, HashMap<Tuple, Integer>> gene2eqclasscounts = new HashMap<>();
        HashMap<String, Integer> genecounts = new HashMap<>();
        TranscriptInfoBasedGenomicMapper mapper;
        File tmpfile;
        PrintWriter ambigPW;
        int nprocessed = 0;
        CountUpdater hitInfos = new CountUpdater("hitinfos");
        CountUpdater validationInfos = new CountUpdater("validation");
        CountUpdater mismatchInfos = new CountUpdater("mismatches");

        Boolean strandness;
        long t1 = System.currentTimeMillis();
        Logger log = LogConfig.getLogger();
        FastqPairReader src;
        int PACKSIZE = 1000;
        boolean ready = false;
        int numGeneAmbigReads = 0;

        MapperJob(TranscriptInfoBasedGenomicMapper mapper, FastqPairReader src) {
            this.label = src.label;
            this.mapper = mapper;
            this.ambigPW = FileUtils.getWriter(tmpfile = new File(label + StringUtils.generateRandomString(20)));
            this.strandness = src.strandness;
            this.src = src;
            tmpfile.deleteOnExit();
        }

        public void run() {
            while(src.hasNext()) {
                Vector<ReadPairInfo> pack = src.readNextPack(PACKSIZE + (int)(Math.random() * PACKSIZE));
                process(pack);
            }
            ready = true;
            ambigPW.close();
        }

        public void resolveAmbig(HashMap<String, Integer> total_gene2count) {
            CountUpdater resolveInfo = new CountUpdater();

            for(String[] sp : FileUtils.getFieldSetsIterable(tmpfile, "\t")) {
                Vector<TranscriptHit.TrHitShort> hits = map(sp, (_s) -> new TranscriptHit.TrHitShort(_s));
                HashSet<String> genes = mapToSet(hits, (_h) -> _h.gene);
                Vector<Pair<String, Integer>> v = NumUtils.sort(map(genes, (_s) -> Pair.create(_s, total_gene2count.getOrDefault(_s, 0))), (_p) -> _p.getSecond(), true);

                boolean firstOk = genes.size() == 1 || v.get(0).getSecond() > 4 * v.get(1).getSecond();
                String selectedGene = (firstOk) ? v.get(0).getFirst() : null;

                if(selectedGene == null) {
                    resolveInfo.update("unresolved", 1);
                    continue;
                }

                int logratio = (v.get(1).getSecond() == 0) ? 100 : (int)NumUtils.logN(v.get(0).getSecond() / (0.0 + v.get(1).getSecond()), 2.0);
                resolveInfo.update("resolved logratio="+logratio, 1);

                Tuple trEQClass = Tuple.tupleFromCollection(toSortedVector(mapToSet(filter(hits, (_h) -> _h.gene.equals(selectedGene)), (_h) -> _h.trid), true));
                MapBuilder.updateMCount(gene2eqclasscounts, selectedGene, trEQClass, 1);

            }
            resolveInfo.printBestN(-1);
        }



        Region1D toGlobal(TranscriptHit th) {
            int x1 = th.transcriptInfo.rv.toGlobalPosition(th.fragment.getX1(), th.transcriptInfo.strand);
            int x2 = th.transcriptInfo.rv.toGlobalPosition(th.fragment.getX2(), th.transcriptInfo.strand);
            return new Region1D(Math.min(x1, x2), Math.max(x1, x2));
        }

        void process(Vector<ReadPairInfo> jobs) {

            for(ReadPairInfo rpi: jobs) {

                nprocessed++;
                if(mapper.verbose) {
                    System.out.printf("next: %s %s %s\n", (rpi.ref == null) ? "unknown" : rpi.ref.transcript, rpi.fw, rpi.rw);

                }
                HashMap<UPair<Boolean>, String> readseqs = new HashMap<>();
                for(int i=0; i<2; i++) {
                    String read = (i == 0) ? rpi.fw : rpi.rw;
                    if(read == null)
                        continue;

                    for(int j=0; j<2; j++) {
                        String seq = (j == 0) ? read : GenomicUtils.reverse_complement(read);
                        readseqs.put(UPair.createU(i == 0, j == 1), seq);

                    }
                }

                Collection<TranscriptHit> trhits = mapper.getTranscriptHits(rpi.readId, readseqs, strandness);

                if(mapper.verbose) {
                    log.info("read mapped to %d transcripts strandness: %s\n", trhits.size(), strandness);
                }

                HashSet<String> trcandidates = new HashSet<>();
                if(trhits.size() == 0) {
                    //check what kind of mapping would be needed....
                    hitInfos.update("unmapped", 1);

                    if(hitInfos.getCount("unmapped") % 10_000 == 0) {
                       // hitInfos.printBestN(-1);

                        long t2 = System.currentTimeMillis();
                        double secs = (t2 - t1) / 1_000.0;
                        double mioreads = nprocessed / 1_000_000.0;
                        log.info("%s: %s ready with %.2f Mio reads took %.2f per Mio reads", src.label, getName(), mioreads, secs / mioreads);

                    }
                    if(rpi.ref != null) {
                        if(trcandidates.contains(rpi.ref.transcript)) {

                            validationInfos.update("unmapped-gotcandidate", 1);
                        } else {
                            validationInfos.update("unmapped-NOCANDIDATE", 1);
                        }
                    } else {
                        validationInfos.update("unmapped", 1);
                    }

                    mismatchInfos.update("unmapped", 1);

                } else {

                    TranscriptHit bestHit = NumUtils.minObj(trhits, (_t) -> _t.numMismatches).getSecond();

                    mismatchInfos.update("mm="+bestHit.numMismatches, 1);
                    Set<String> besthitGenes = mapToSet(filter(trhits, (_t) -> _t.numMismatches == bestHit.numMismatches), (_t) -> _t.transcriptInfo.gene);

                    int MAXMISMATCHES = bestHit.numMismatches + 2;

                    boolean gotcorrecttr = false;
                    boolean gotcorrecttr_poscorrect = false;
                    boolean ambig = false;
                    if(besthitGenes.size() == 1) {
                        String gene = first(besthitGenes);
                        MapBuilder.update(genecounts, gene);

                        Vector<TranscriptHit> geneTrHits = filter(trhits, (_t) -> gene.equals(_t.transcriptInfo.gene) && _t.numMismatches <= MAXMISMATCHES);
                        Tuple trEQClass = Tuple.tupleFromCollection(toSortedVector(mapToSet(geneTrHits, (_t) -> _t.transcriptInfo.transcriptId), true));
                        if(rpi.ref != null) {
                            TranscriptHit sensRealHit = filterOne(trhits, (_t) -> _t.transcriptInfo.transcriptId.equals(rpi.ref.transcript));
                            TranscriptHit realHit = filterOne(geneTrHits, (_t) -> _t.transcriptInfo.transcriptId.equals(rpi.ref.transcript));

                            if(gene.equals(rpi.ref.gene) && sensRealHit != null && realHit == null) {

                                System.out.printf("ref: %s sens real hit: %s best hit: %s\nL1: %d L2: %d\n%s vs %s\n", rpi.ref, sensRealHit, bestHit,
                                        sensRealHit.fragment.getLength(), bestHit.fragment.getLength(),
                                        toGlobal(sensRealHit), toGlobal(bestHit));
                            }
                            gotcorrecttr = (realHit != null);

                            gotcorrecttr_poscorrect = (gotcorrecttr &&
                                            realHit.fragment.getX1() == rpi.ref.t_fw.getX1() &&
                                            realHit.fragment.getX2() == rpi.ref.t_rw.getX2());

                            if(realHit != null && !gotcorrecttr_poscorrect) {
                                System.out.printf("got realhit but wrong choord: ref: %s real hit: %s\n", rpi.ref, realHit);
                            }
                            if(false)
                                System.out.printf("mm-gene-uniq geneok: %s gotcorrect-tr: %s trs: %d ref: %s realhit: %s\n", gene.equals(rpi.ref.gene), gotcorrecttr,
                                        geneTrHits.size(),
                                        rpi.ref, realHit);
                        }
                        MapBuilder.updateMCount(gene2eqclasscounts, gene, trEQClass, 1);
                    } else {

                        ambig = true;
                        TranscriptHit refHit = (rpi.ref == null) ? null : filterOne(trhits, (_t) -> _t.transcriptInfo.transcriptId.equals(rpi.ref.transcript));

                        Vector<TranscriptHit> ambig_hits = filter(trhits, (_t) -> besthitGenes.contains(_t.transcriptInfo.gene) && _t.numMismatches <= MAXMISMATCHES);

                        if(rpi.ref != null) {
                            TranscriptHit realHit = filterOne(ambig_hits, (_h) -> _h.transcriptInfo.transcriptId.equals(rpi.ref.transcript));
                            gotcorrecttr = realHit != null;
                            gotcorrecttr_poscorrect = gotcorrecttr && realHit.fragment.getX1() == rpi.ref.t_fw.getX1() &&
                                    realHit.fragment.getX2() == rpi.ref.t_rw.getX2();
                        }
                        Vector<String> infos = filter_and_map(trhits, (_t) -> besthitGenes.contains(_t.transcriptInfo.gene) && _t.numMismatches <= MAXMISMATCHES, (_t) -> _t.getShortInfo());
                        if(rpi.ref != null) {

                            if(refHit != null) {
                                infos.insertElementAt(refHit.getRefShortInfo(), 0);
                            }
                        }


                        this.ambigPW.printf("%s\n", StringUtils.joinObjects("\t", infos))   ;
                        numGeneAmbigReads++;

                    }
                    boolean bestMismatchAmbigGene = besthitGenes.size() > 1;
                    boolean gene_correct = rpi.ref != null && besthitGenes.contains(rpi.ref.gene);

                    if(gotcorrecttr) {
                        String label =  String.format("%s-%s-correct-%sTR", "regular",
                                (ambig) ? "ambig": "regular", (gotcorrecttr_poscorrect) ? "Pos" : "");

                        validationInfos.update(label, 1);

                    } else {
                        validationInfos.update( "regular" + "-gene_correct", 1);
                    }


                    if(rpi.ref != null) {
                        if(!gene_correct) {
                            HashMap<String, Integer> gene2min = new HashMap<>();
                            for(TranscriptHit th : trhits){
                                gene2min.put(th.transcriptInfo.gene, Math.min(th.numMismatches, gene2min.getOrDefault(th.transcriptInfo.gene, 100)));
                            }
                            //System.out.printf("gene incorrect valid would be: %s minmm: %d hits: %s\n", info.gene, gene2min.get(info.gene), Pair.convert_reverse_sorted(gene2min, true));
                        }

                    }
                    HashSet<String> genehits = mapToSet(trhits, (_t) -> _t.transcriptInfo.gene);


                    hitInfos.updateMulti(1,
                            toVector("tr-" + ((trhits.size() == 1) ? "uniq" : "ambig"),
                                    "gene-" + ((genehits.size() == 1) ? "uniq" : "ambig"),
                                    "bestMM-" + ((bestMismatchAmbigGene) ? "ambig" : "uniq")
                            )
                    );

                    if(nprocessed % 1_000_000 != 0)
                        continue;

                    long t2 = System.currentTimeMillis();
                    double secs = (t2 - t1) / 1_000.0;
                    double mioreads = nprocessed / 1_000_000.0;
                    log.info("%s: %s ready with %.2f Mio reads took %.2f per Mio reads", src.label, getName(), mioreads, secs / mioreads);
                    //hitInfos.printBestN(-1);
                    //mismatchInfos.printBestN(15);
                    if(rpi.ref != null) {
                        validationInfos.printBestN(-1);
                    }

                }
            }
        }
    }




    static class ReadPairInfo {
        int readIdx;
        String readId;
        String fw;
        String rw;
        ReadReferenceInfo ref;
    }

    public static class FastqPairReader {

        Logger log = LogConfig.getLogger();
        String label;
        FastQReader fw;
        FastQReader rw;
        Iterator<ReadReferenceInfo> refInfo;
        boolean closed = false;
        int nreads;
        FastQRecord fwRecord;
        FastQRecord rwRecord;
        long t1 = System.currentTimeMillis();

        Boolean strandness;

        FastqPairReader(String label, File fwq, File rwq, File info, Boolean strandness) {
            this.label = label;
            fw = new FastQReader(fwq);
            rw = (rwq == null) ? null : new FastQReader(rwq);
            this.strandness = strandness;
            refInfo = (info == null) ? null : ReadReferenceInfo.getIterator(info);
        }


        boolean hasNext() {
            return !closed;
        }

        public synchronized Vector<ReadPairInfo> readNextPack(int numRequested) {

            Vector<ReadPairInfo> rv = new Vector<>();

            if(closed)
                return rv;

            while (null != (fwRecord = fw.next())) {
                rwRecord = (rw == null) ? null : rw.next();
                ReadPairInfo rpi = new ReadPairInfo();
                rpi.readId = fwRecord.header.toString();
                rpi.readIdx = nreads;
                rpi.fw = fwRecord.readseq.toString();
                rpi.rw = (rwRecord == null) ? null : rwRecord.readseq.toString();
                rpi.ref = (refInfo == null) ? null : refInfo.next();
                nreads++;

                rv.add(rpi);

                if(nreads % 1_000_000 == 0) {
                    long t2 = System.currentTimeMillis();
                    double mio = nreads / 1_000_000.;
                    double sec_per_mio = (t2 - t1) / (mio * 1000.0);
                    log.info("%s read : %.2f Mio reads %.2f sec per mio", label, mio, sec_per_mio);
                }
                if (rv.size() == numRequested)
                    break;

            }
            if (fwRecord == null) {
                closed = true;
            }
            return rv;
        }
    }

    public static void main(String[] args) {
        //gtf, genome, extract all (or only certain type of ) transcripts (from an optionally limited set of chromosomes */

        SimpleOptionParser cmd = new SimpleOptionParser("index", "table", "v", "o", "nthreads", "basedir");
        cmd.setFile("table", "index");
        cmd.setDir("basedir");
        cmd.setInt("nthreads");
        cmd.setOptional("basedir");
        cmd.setDefault("nthreads", "10");
        cmd.setSwitches("v");

        if (!OptionParser.parseParams(args, true, true, true, cmd))
            return;

        long t1 = System.currentTimeMillis();
        File basedir = cmd.getOptionalFile("basedir");


        Logger log = LogConfig.getLogger();
        //put it into  asuffix array
        TranscriptInfoBasedGenomicMapper transcriptInfoBasedGenomicMapper = new TranscriptInfoBasedGenomicMapper(cmd.getFile("index"));

        transcriptInfoBasedGenomicMapper.verbose = cmd.isSet("v");

        HashMap<String, Vector<String>> cond2reps = new HashMap<>();
        DataTable dt = DataTable.readFile(cmd.getFile("table"), "\t");

        HashMap<String, String> rep2cond = new HashMap<>();
        Vector<FastqPairReader> inputs = new Vector<>();

        boolean gotinfo = null != filterOne(dt.header.getHeaderNames(), (_n) -> "info".equals(_n));
        Vector<String> labelvec = new Vector<>();
        for (int i = 0; i < dt.getSize(); i++) {

            ObjectGetter og = dt.getDataAt(i);
            String cond = og.getString("condition");
            String label = og.getString("label");
            labelvec.add(label);
            inputs.add(new FastqPairReader(label, og.getPath("fw", basedir), !og.gotHeader("rw") ? null : og.getPath("rw", basedir), (!gotinfo) ? null : og.getPath("info"), og.getBoolean("strandness", null)));
            MapBuilder.updateV(cond2reps, cond, label);
            if (rep2cond.get(label) != null) {
                System.err.printf("label not unique: %s!\n", label);
                return;
            }
            rep2cond.put(label, cond);
        }


        LinkedList<MapperJob> mapperjobs = new LinkedList<>();

        int nthreads = cmd.getInt("nthreads");


        while (mapperjobs.size() < nthreads) {
            for (FastqPairReader src : inputs) {
                mapperjobs.add(new MapperJob(transcriptInfoBasedGenomicMapper, src));
            }
        }

        List<MapperJob> readyjobs = new LinkedList<>();
        LinkedList<MapperJob> startedJobs = new LinkedList<>();

        while (true) {
            UPair<Integer> preStatus = UPair.createU(startedJobs.size(), readyjobs.size());

            while (startedJobs.size() < nthreads && mapperjobs.size() > 0) {
                MapperJob mj = mapperjobs.poll();
                mj.start();
                startedJobs.add(mj);
            }

            List<MapperJob> stillRunning = new ArrayList<>();
            while (startedJobs.size() > 0) {
                MapperJob mj = startedJobs.poll();
                if (!mj.ready) {
                    stillRunning.add(mj);
                } else {
                    readyjobs.add(mj);
                }
            }
            startedJobs.addAll(stillRunning);

            UPair<Integer> postStatus = UPair.createU(startedJobs.size(), readyjobs.size());

            if(!preStatus.equals(postStatus)) {
                log.info("pre status (started/ready jobs) %s post status: %s", preStatus, postStatus);
            }
            if (startedJobs.size() == 0)
                break;

            try {
                Thread.sleep(500);
            } catch (InterruptedException iex) {
                break;
            }

        }

        //create a condition based gene2count
        HashMap<String, HashMap<String, Integer>> cond2gene2count = new HashMap<>();

        for (MapperJob mj : readyjobs) {
            String cond = rep2cond.get(mj.label);

            for (Map.Entry<String, Integer> e : mj.genecounts.entrySet()) {
                MapBuilder.updateMCount(cond2gene2count, cond, e.getKey(), e.getValue());
            }
        }
        HashSet<String> allgenes = new HashSet<>();

        for (MapperJob mj : readyjobs) {
            String cond = rep2cond.get(mj.label);
            mj.resolveAmbig(cond2gene2count.get(cond));
            log.info("mapper job: %s cond: %s got %d genes with counts\n",
                    mj, cond, (cond2gene2count.get(cond) == null) ? -1 : cond2gene2count.get(cond).size());
            allgenes.addAll(cond2gene2count.get(cond).keySet());
        }

        log.info("write eq-class counts on %d genes", allgenes.size());


        HashMap<String, HashMap<String, HashMap<Tuple, Integer>>> label2gene2eqclasscounts = new HashMap<>();
        for (MapperJob mj : readyjobs) {
            HashMap<String, HashMap<Tuple, Integer>> m = label2gene2eqclasscounts.get(mj.label);
            if(m == null) {
                label2gene2eqclasscounts.put(mj.label, m = new HashMap<>());
            }

            for(String gene : mj.gene2eqclasscounts.keySet()) {
                for(Map.Entry<Tuple, Integer>  e : mj.gene2eqclasscounts.get(gene).entrySet()) {
                    MapBuilder.updateMCount(m, gene, e.getKey(), e.getValue());
                }
            }

        }

        log.info("labels keys: %s labelvec: %s", label2gene2eqclasscounts.keySet(), labelvec);


        //and write out eqclass output
        PrintWriter pw = cmd.getWriter("o");
        for(String gene : allgenes) {


            Vector<HashMap<Tuple, Integer>> data = map(labelvec, (_l) -> label2gene2eqclasscounts.get(_l).getOrDefault(gene, new HashMap<>()));

            Vector<Integer> total = map(data, (_m) -> NumUtils.sum(_m.values(), 0));
            if(NumUtils.sum(total) == 0)
                continue;

            pw.printf(">%s\tTOTAL\n", gene);

            HashSet<Tuple> merged = new HashSet<>();
            apply(data, (_d) -> merged.addAll(_d.keySet()));

            for(String type : toVector("reads")) {

                pw.printf("%s\t%s\n", type, StringUtils.joinObjects("\t", total, (_t) -> String.format("%.3f", _t + 0.0)));
            }

            for(Tuple t : merged) {
                pw.printf(">%s\t%s\n", gene, StringUtils.joinObjects(",", mapIndex(t.cardinality(), (_i) -> t.getAsString(_i))));
                Vector<Integer> counts = map(data, (_m) -> _m.getOrDefault(t, 0));
                for(String type : toVector("reads")) {
                    pw.printf("%s\t%s\n", type, StringUtils.joinObjects("\t", counts, (_t) -> String.format("%.3f", _t + 0.0)));
                }
            }
        }
        pw.close();
        long t2 = System.currentTimeMillis();
        log.info("total processing time (wallclock time): %.2f sec", (t2 - t1) / 1_000.0);
    }
}
