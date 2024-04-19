package nlEmpiRe.rnaseq.reads;

import nlEmpiRe.rnaseq.*;
import lmu.utils.*;
import lmu.utils.tuple.Tuple3;
import net.sf.samtools.AlignmentBlock;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

import java.io.File;
import java.io.PrintWriter;
import java.util.*;
import java.util.function.BiConsumer;

import static lmu.utils.ObjectGetter.*;
import static lmu.utils.ObjectGetter.apply;
import org.apache.logging.log4j.Logger;

public class BAMIterator<A>
{
    static boolean NO_PCR_DOWNSAMPLING = false;

    Logger log = LogConfig.getLogger();
    IsoformRegionGetter annot;

    public BAMIterator(IsoformRegionGetter annot)
    {
        this.annot = annot;
    }


    public static class CSAMRecord implements Comparable<CSAMRecord>
    {
        SAMRecord sr;

        CSAMRecord(SAMRecord sr)
        {
            this.sr = sr;
        }

        public int compareTo(CSAMRecord csr)
        {int d = csr.sr.getReferenceIndex().compareTo(sr.getReferenceIndex());
            return (d != 0) ? d : ( 0 != (d = csr.sr.getAlignmentStart() - sr.getAlignmentStart())) ? d : csr.sr.getAlignmentEnd() - sr.getAlignmentEnd();
        }

        public String toString()
        {
            String mateinfo = (sr == null || !sr.getReadPairedFlag()) ? "" : (sr.getMateUnmappedFlag()) ? " mate unmapped" : String.format(" mate: %s:%d-", sr.getMateReferenceName(), sr.getMateAlignmentStart());
            return (sr == null) ? "<null>" : String.format("rid: %s %s:%d-%d%s", correctReadId(sr), sr.getReferenceName(), sr.getAlignmentStart(), sr.getAlignmentEnd(), mateinfo);
        }
    }
    static class ClosableSamIterator<A>  implements  Comparable<ClosableSamIterator>
    {
        Logger log = LogConfig.getLogger();
        Boolean strandness;
        File bam;
        File bami;

        A obj;
        SAMRecordIterator it;
        SAMFileReader sam_reader;
        PeakableIterator<CSAMRecord> pit;

        public ClosableSamIterator(A obj, File bam, File bami, Boolean strandness)
        {
            this.obj = obj;
            this.bam = bam; this.bami = bami; this.strandness = strandness;

            sam_reader  = new SAMFileReader(bam, false);
            sam_reader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
            it = sam_reader.iterator();
            pit = new PeakableIterator<>(getTransformedIterator(it, (_sr) -> new CSAMRecord(_sr)));
        }

        public void query(String chr, int start, int end)
        {
            it.close();
            it = sam_reader.query(chr, start, end, false);
            log.debug("query %s for %s %d-%d", bam.getAbsolutePath(), chr, start, end);
            pit = new PeakableIterator<>(getTransformedIterator(it, (_sr) -> new CSAMRecord(_sr)));
        }



        public int compareTo(ClosableSamIterator csi)
        {
            CSAMRecord cs1 = (CSAMRecord)csi.pit.peak();
            CSAMRecord cs2 = pit.peak();
            if(cs1 == null)
                return 1;
            if(cs2 == null)
                return -1;

            return cs1.compareTo(cs2);
        }

        public void close()
        {
            it.close();
            sam_reader.close();
        }

        public String toString()
        {
            return String.format("(%s: %s)", obj, pit.peak());
        }
    }


    PriorityQueue<ClosableSamIterator<A>> bams = new PriorityQueue<ClosableSamIterator<A>>();

    public void addBAM(A obj,File bam, File bami, Boolean strandness)
    {
        bams.add(new ClosableSamIterator<>(obj, bam, bami, strandness));
    }



    public static String getAliBlockInfo(SAMRecord sr)
    {
        Vector<String> blocks = new Vector<>();
        for (AlignmentBlock ab : sr.getAlignmentBlocks())
        {
            int s = ab.getReferenceStart();
            int rs = ab.getReadStart();
            blocks.add(String.format("ref: (%d-%d) read: (%d-%d)", s, s + ab.getLength(), rs, rs + ab.getLength()));
        }
        return "" + blocks;

    }

    public static class SRInfo
    {
        TreeSet<Region1D> splits = new TreeSet<>();
        TreeSet<Region1D> deletions = new TreeSet<>();
        SAMRecord sr;
        RegionVector rv;

        SRInfo(SAMRecord sr)
        {
            this.sr = sr;
            rv = (sr == null) ? null : toRV(sr, splits, deletions, 0);
        }

        public String toString()
        {
            return (sr == null) ? "sr: null": String.format("%s%c %s splits: %s dels: %s", sr.getReferenceName(), (sr.getReadNegativeStrandFlag()) ? '-' : '+',  rv, splits, deletions);
        }
    }

    public static RegionVector toRV(SAMRecord sr, SAMRecord mate, TreeSet<Region1D> splits, TreeSet<Region1D> deletions, int minsplitsize)
    {
        RegionVector mv = toRV(sr, splits, deletions, minsplitsize);
        return (mate == null) ? mv : RegionVector.merge(mv, toRV(mate, splits, deletions, minsplitsize));

    }

    public static RegionVector toRV(SAMRecord sr) {
        return toRV(sr, new TreeSet<>(), new TreeSet<>(), 0);
    }

    public static RegionVector toRV(SAMRecord sr, TreeSet<Region1D> splits, TreeSet<Region1D> deletions, int minsplitsize)
    {

        int last = -1;
        Vector<Region1D> regs = new Vector<>();
        for (AlignmentBlock ab : sr.getAlignmentBlocks())
        {
            final int L = regs.size() -1;
            int s = ab.getReferenceStart();

            if (last >= 0)
            {
                Region1D skipped = new Region1D(last, s);

                if (s - last < minsplitsize)
                {
                    deletions.add(skipped);
                    //extend last region
                    last = s + ab.getLength();
                    regs.get(L).setX2(last);
                    continue;
                }
                if (last < s)
                {
                    splits.add(new Region1D(last, s));
                }


            }

            if(L >= 0 && s == regs.get(L).getX2())
            {
                regs.get(L).setX2(last = s + ab.getLength());
                continue;

            }


            regs.add(new Region1D(s, last = s + ab.getLength()));
        }
        return new RegionVector(regs, true);
    }

    public static UPair<Region1D> isInconsistent3(RegionVector merged, TreeSet<Region1D> splits)
    {
        if(splits.size() == 0)
            return null;

        TreeSet<Region1D> implied = new TreeSet<>();
        implied.addAll(merged.getInverseRegionVector().getRegions());
        for(Region1D r  : splits)
        {
            if(!implied.contains(r))
                return UPair.createU(r, null);
        }

        return null;
    }

    public static UPair<Region1D> isInconsistent2(RegionVector merged, TreeSet<Region1D> splits)
    {
        if(splits.size() == 0)
            return null;

        for(Region1D r : splits) {
            Pair<Boolean, Integer> p1 = VectorUtils.binarySearch(merged.getRegions(), r.getX1(), (_r1, _p) -> Integer.compare(_r1.getX2(), _p));
            if (!p1.getFirst())
                return UPair.createU(r, null);

            Pair<Boolean, Integer> p2 = VectorUtils.binarySearch(merged.getRegions(), r.getX2(), (_r1, _p) -> Integer.compare(_r1.getX1(), _p));

            if (!p2.getFirst())
                return UPair.createU(r, null);

            if (p2.getSecond() - p1.getSecond() != 1)
                return UPair.createU(r, null);
        }

        return null;

    }
    public static UPair<Region1D> isInconsistent(RegionVector merged, TreeSet<Region1D> splits)
    {
        if (splits.size() ==0)
            return null;

        for (Region1D sp: splits)
        {
            for (Region1D r : merged.getRegions())
            {
                if (r.getOverlap(sp) > 0)
                    return UPair.createU(r, sp);


                if (sp.getX1() < r.getX1() && sp.getX2() == r.getX1())
                    continue; //consinstent

                if (sp.getX2() > r.getX2() && sp.getX1() == r.getX2())
                    continue; //consinstent

                if (!r.contains(sp.getX1()) && !r.contains(sp.getX2()))
                    continue;

                return UPair.createU(r, sp);
            }
        }
        return null;
    }


    public enum READTYPE
    {
        INCONSISTENT,
        GENE_OVERLAPPING,
        INTRONIC,
        MERGED_TR,
        TRANSCRIPTOMIC
        ;
    }

    /** if set to true the internal regions, splits etc do not matter
     *  -> relevant for DNA-sequencing e.g. ATAC-seq
     */
    boolean ignore_fragment_internal_structure = false;

    public void setIgnoreFragmentInternalStructure(boolean b)
    {
        ignore_fragment_internal_structure = b;
    }

    public static class ReadAnnotInfo<A>
    {
        public String chr;
        AnnotKey key;
        public RegionVector merged;
        RegionVector splits;
        public Boolean strand;
        public Vector<MultiIsoformRegion> contained_genes;
        public Vector<MultiIsoformRegion> in_genes;
        public HashSet<MultiIsoformRegion> MERGED_GENES = new HashSet<>();
        public Vector<MultiIsoformRegion> INTRONIC_GENES = new Vector<>();
        public Vector<MultiIsoformRegion> ANTISENSE_GENES = new Vector<>();
        public HashMap<MultiIsoformRegion, HashSet<String>> TR_GENES = new HashMap<>();

        Pair<Integer, Region1D<String>> closest_strandless = null;
        Pair<Integer, Region1D<String>> closest_strandspec = null;

        int gdist = -1;

        public HashMap<A, Integer> pcr_index = new HashMap<>();


        READTYPE rtype = READTYPE.GENE_OVERLAPPING;

        ReadAnnotInfo(String chr, RegionVector fw_rv, RegionVector rw_rv, RegionVector merged, RegionVector splits, Boolean strand, IsoformRegionGetter annot)
        {
            this.chr = chr;
            this.merged = merged; this.strand = strand; this.splits = splits;
            key = new AnnotKey(merged, strand, splits);
            Vector<MultiIsoformRegion> genes = toVector(annot.getRegions(chr, merged.getX1(), merged.getX2()));

            //ANTISENSE_GENES = filter(genes, (_g) -> strand != null && _g.strand != strand);

            //System.out.printf("%s %d-%d get genes: %d\n", chr, merged.getX1(), merged.getX2(), genes.size());
            contained_genes = filter(genes, (_g) -> (strand == null || strand == _g.strand) && _g.start >= merged.getX1() && _g.end  <= merged.getX2());
            in_genes= filter(genes, (_g) -> (strand == null || strand == _g.strand) && _g.start <= merged.getX1() && _g.end  >= merged.getX2());
            ANTISENSE_GENES = filter(genes, (_g) -> strand != null && _g.strand != strand && _g.start <= merged.getX1() && _g.end  >= merged.getX2());


            if(in_genes.size() == 0 && contained_genes.size() > 0)
                return;




            if(in_genes.size() == 0)
            {
                rtype = READTYPE.INTRONIC;
                //calc distance
                closest_strandless = annot.getClosestRegion(chr, null, merged.getX1(), merged.getX2());
                closest_strandspec =  (strand == null) ? closest_strandless : annot.getClosestRegion(chr, strand, merged.getX1(), merged.getX2());
                //System.out.printf("%s dist: %s\n", merged, closest);

                gdist = (closest_strandspec == null) ? Integer.MIN_VALUE : closest_strandspec.getFirst();
                Vector<MultiIsoformRegion> dgenes = toVector(annot.getRegionsAround(chr, merged.getX1(), merged.getX2()));

                //System.out.printf("genes around: %s\n", map(dgenes, (_g) -> _g.getLocationInfo()));
                return;
            }


            for(MultiIsoformRegion mir : in_genes)
            {
                if(!mir.getMergedTranscript().containsRV2(merged)) {
                    INTRONIC_GENES.add(mir);
                    continue; //not in merged
                }

                HashSet<String> trset = new HashSet<>();
                for(Map.Entry<String, RegionVector> e : mir.isoforms.entrySet())
                {
                    RegionVector tr = e.getValue();
                    //boolean old = tr.isSubVector(merged);
                    boolean n = (merged.getNumRegions() == 1) ? tr.isSubVector3(merged, false) : tr.isSubVector3(fw_rv, false)
                            && (rw_rv == null || tr.isSubVector3(rw_rv, false));
                    if(!n)
                        continue;

                    trset.add(e.getKey());

                }
                if(trset.size() == 0)
                {
                    MERGED_GENES.add(mir);
                    continue;
                }
                TR_GENES.put(mir, trset);

            }

            if(TR_GENES.size()  > 0)
            {
                rtype = READTYPE.TRANSCRIPTOMIC;
                return;
            }
            if(MERGED_GENES.size() > 0)
            {
                rtype = READTYPE.MERGED_TR;
                return;
            }

            rtype = READTYPE.INTRONIC;

        }

        public HashMap<A, Integer> getPCR_index()
        {
            return this.pcr_index;
        }

        public RegionVector getMerged()
        {
            return this.merged;
        }
    }

    static boolean gobi_mode = true;
    static boolean gobi_verbose = false;

    static HashMap<String, String> chr2chr = new HashMap<>();
    static Logger chrModLog = LogConfig.getLogger();

    static String correctChr(String chr, IsoformRegionGetter annot) {
        String mapped = chr2chr.get(chr);
        if(mapped != null)
            return mapped;

        if(annot.getChrLength(chr) != null) {
            chr2chr.put(chr, chr);
            return chr;
        }

        //check if it is due to prefix
        String modchr = (chr.startsWith("chr")) ? chr.substring(3) : "chr"+chr;
        if(annot.getChrLength(modchr) == null) {
            //maybe MT -> M or the vice versa`
            if(chr.equals("MT") || chr.equals("M")) {
                modchr = (chr.equals("MT")) ? "M" : "MT";
                if(annot.getChrLength(modchr) != null) {
                    chrModLog.info("input chr: %s change it to :%s", chr, modchr);
                    chr2chr.put(chr, modchr);
                    return modchr;
                }
            }
        } else {

            chrModLog.info("input chr: %s change it to :%s", chr, modchr);
            chr2chr.put(chr, modchr);
            return modchr;
        }

        chr2chr.put(chr, chr);
        return chr;
    }

    static String correctReadId(SAMRecord sr) {
        String readId = sr.getReadName();
        if(readId.endsWith("/1") || readId.endsWith("/2")) {
            return readId.substring(readId.length() - 2);
        }
        return readId;

    }
    public static class AnnotatedRead<A>
    {
        public A obj;
        public String rid;
        public String chr;
        public SAMRecord fw;
        public SAMRecord rw;

        int second_start = -1;
        Tuple3<Integer, Integer, Integer> fw_infos;
        Tuple3<Integer, Integer, Integer> rw_infos;

        public RegionVector fw_rv;
        public RegionVector rw_rv;
        public RegionVector merged;
        Boolean strand;

        UPair<Region1D> incons = null;
        TreeSet<Region1D> splits = new TreeSet<>();
        TreeSet<Region1D> deletions = new TreeSet<>();



        int gdist = -1;

        Region1D partly_defined_region = null;
        boolean gene_overlapping_read = false;

        boolean fr_strand;
        public int pcr_index = 0;
        public ReadAnnotInfo<A> info = null;
        IsoformRegionGetter annot;

        AnnotatedRead(A obj, SAMRecord sr1, int mate_start,  Boolean strand, IsoformRegionGetter annot)
        {
            this.obj = obj; this.annot = annot;

            this.strand = strand;

            rid = correctReadId(sr1);
            chr = correctChr(sr1.getReferenceName(), annot);

            int minpos = Math.min(sr1.getAlignmentStart(), mate_start);
            int maxpos = Math.max(sr1.getAlignmentStart(), mate_start);
            partly_defined_region = (minpos >= maxpos) ? null : new Region1D(minpos, maxpos); //Math.max(sr1.getAlignmentEnd(), mate_start));
            gene_overlapping_read = true;
        }

        public boolean getReadStrand(boolean first)
        {
            return !((first) ? fw : rw).getReadNegativeStrandFlag();
        }

        public READTYPE getType()
        {
            return (gene_overlapping_read) ? READTYPE.GENE_OVERLAPPING : (incons != null) ? READTYPE.INCONSISTENT : (info == null) ? READTYPE.GENE_OVERLAPPING : info.rtype;
        }

        public String toString()
        {
            return String.format("anread<%s>: %s gov: %s", obj, rid, gene_overlapping_read);
        }

        AnnotatedRead(A obj, IsoformRegionGetter annot, Cache<A> cache, Boolean strandness, SAMRecord sr1, SAMRecord sr2)
        {
            this.obj = obj;
            this.annot = annot;
            second_start = (sr2 == null) ? sr1.getAlignmentStart() : Math.max(sr1.getAlignmentStart(), sr2.getAlignmentStart());



            rid = correctReadId(sr1);


            fw = (sr2 == null || sr1.getFirstOfPairFlag()) ? sr1 : sr2;
            rw =  (sr2 == null) ? null :  (sr1.getFirstOfPairFlag()) ? sr2 : sr1;

            fr_strand = !fw.getReadNegativeStrandFlag();
            strand = (strandness == null) ? null : (strandness) ? !fw.getReadNegativeStrandFlag() : fw.getReadNegativeStrandFlag();
            fw_infos = getClipAndMM(fw);
            rw_infos = (sr2 == null) ? null : getClipAndMM(rw);

            chr = correctChr(fw.getReferenceName(), annot);


            fw_rv = toRV(fw, splits, deletions, 0);
            rw_rv = (sr2 == null) ? null : toRV(rw, splits, deletions, 0);
            merged = (sr2 == null) ? fw_rv : RegionVector.merge(fw_rv, rw_rv);

            if(cache != null && cache.ignore_fragment_internal_structure) {
                merged = new RegionVector(merged.getX1(), merged.getX2());
                splits.clear();
            }


            incons = (sr2 == null) ? null : isInconsistent2(merged, splits);

            if(incons != null)
            {
                if(gobi_mode) {
                    if(gobi_verbose) {
                        info = new ReadAnnotInfo<A>(chr, fw_rv, rw_rv, merged, null, strand, annot);
                    }
                    gene_overlapping_read = checkSkipDueContainment(annot, chr, merged.getX1(), merged.getX2(), strand);
                }

                return;
            }




            RegionVector splits_rv = (splits == null || splits.size() == 0) ? null : new RegionVector(toSortedVector(splits, true), true);

            info = (cache == null) ? null : cache.rv2ar.get(new AnnotKey(merged, strand, splits_rv));
            if(incons == null && info != null)
            {

                pcr_index = 1 + info.pcr_index.getOrDefault(obj, -1);
                info.pcr_index.put(obj, pcr_index);
                return;

            }


            info = new ReadAnnotInfo<A>(chr, fw_rv, rw_rv, merged, splits_rv, strand, annot);
            gene_overlapping_read = (info.in_genes.size() == 0 && info.contained_genes.size() > 0);

            if(incons != null)
                return;

            info.pcr_index.put(obj, 0);

            if(cache != null) {
                cache.add(chr, info);
            } else {

            }

        }

        public void check(String chinfo)
        {
            PrintWriter pw = new PrintWriter(System.out);
            check(pw, chinfo);
            pw.flush();
        }
        public void check(PrintWriter pw, String chinfo)
        {
            pw.printf("\n---------------------------------------------\n");

            pw.printf("read: %s", rid);
            if(chinfo != null)
            {
                pw.printf(" check: %s", chinfo);
            }
            pw.printf("\nr1: %s r2: %s merged: %s splits: %s\nREADPAIR strand: %c\n", new SRInfo(fw), new SRInfo(rw), merged, splits, GenomicUtils.getStrand(strand));


            if(info != null)
            {
                pw.printf("gene_overlapping_read: %s\n\tin genes: %d  (%s)\n\toverlapping genes: %d (%s)\n",gene_overlapping_read,
                        (info.in_genes == null) ? -1 : info.in_genes.size(), (info.in_genes == null) ? "()" : info.in_genes,
                        (info.contained_genes == null) ? -1 : info.contained_genes.size(), (info.contained_genes == null) ? "()" : info.contained_genes);
            }
            else
            {
                pw.printf("info null gene_overlapping_read: %s!\n", gene_overlapping_read);
            }

            if(partly_defined_region != null)
            {
                Vector<Tuple3<String, Boolean, Region1D<String>>> reginfos = annot.getRegionInfos(chr, partly_defined_region.getX1(), partly_defined_region.getX2());

                Vector<MultiIsoformRegion> genes = toVector(annot.getRegions(chr, partly_defined_region.getX1(), partly_defined_region.getX2()));

                Vector<MultiIsoformRegion> contained_genes = filter(genes, (_g) -> (strand == null || strand == _g.strand) && _g.start >= partly_defined_region.getX1() && _g.end  <= partly_defined_region.getX2());
                Vector<MultiIsoformRegion> in_genes= filter(genes, (_g) -> (strand == null || strand == _g.strand) && _g.start <= partly_defined_region.getX1() && _g.end  >= partly_defined_region.getX2());

                pw.printf("partly defined region: %s\n", partly_defined_region);
                pw.printf("contained genes (%d): \n%s\nin genes(%d): %s\n",
                        contained_genes.size(), StringUtils.joinObjects("\n\t", contained_genes, (g) -> "" +g),
                        in_genes.size(), StringUtils.joinObjects("\n\t", in_genes, (g) -> "" +g)
                );


            }
            pw.printf("pcr: %d gdist: %d strand: %s incons: %s trgenes: %s merged: %s intronic: %s\n", pcr_index, gdist, strand, incons,
                    (info == null) ? "null" : info.TR_GENES,
                    (info == null) ? "null" : info.MERGED_GENES,
                    (info == null) ? "null" : info.INTRONIC_GENES);

            pw.printf("antisense: %s genes: %s\n", (info != null && info.ANTISENSE_GENES != null && info.ANTISENSE_GENES.size() > 0), (info == null) ? "info empty" : info.ANTISENSE_GENES);
            pw.printf("fw clip/ mm info: %s\n", (fw == null) ? null : getClipAndMM(fw));
            pw.printf("rw clip/ mm info: %s\n", (rw == null) ? null : getClipAndMM(rw));
            if(merged != null) {
                LIntervalTree<Tuple3<String, Boolean, Region1D<String>>> tree = annot.region_tree.getTree(chr);
                if(tree == null)
                {
                    pw.printf("TREE IS NULL FOR >%s<\n", chr);
                }
                else {
                    pw.printf("left neightbours: %s\n", tree.getNeighBours(merged.getX1(), merged.getX2(), true));
                    pw.printf("right neightbours: %s\n", tree.getNeighBours(merged.getX1(), merged.getX2(), false));
                }

            }

            pw.printf("closest strandspec: %s\n", (info == null) ? "null" : info.closest_strandspec);
            pw.printf("closest strandless: %s\n", (info == null) ? "null" : info.closest_strandless);
            pw.printf("\n\n");
        }

        Tuple3<Integer, Integer, Integer> getClipAndMM(SAMRecord sr)
        {
            boolean strand = !sr.getReadNegativeStrandFlag();
            int _fw_clip = sr.getAlignmentStart() - sr.getUnclippedStart();
            int _rw_clip = sr.getUnclippedEnd() - sr .getAlignmentEnd();
            int fw_clip = (strand) ? _fw_clip : _rw_clip;
            int rw_clip = (strand) ? _rw_clip : _fw_clip;

            Integer nm = (Integer)sr.getAttribute("NM");
            nm = (nm != null) ? nm : (Integer)sr.getAttribute("nM");
            nm = (nm != null) ? nm : (Integer)sr.getAttribute("XM");


            return Tuple3.create(fw_clip, rw_clip, nm);



        }


        void writePCR(PrintWriter pw)
        {
            pw.printf("\tpcrindex: %d\n", pcr_index);
        }
        void write(PrintWriter pw)
        {

            if(incons != null)
            {
                pw.printf("%s\tsplit-inconsistent:true\n", rid);
                return;
            }
            Tuple3<Integer, Integer, Integer> t_fw = getClipAndMM(fw);
            Tuple3<Integer, Integer, Integer> t_rw = (rw == null) ? null : getClipAndMM(rw);
            int mm = t_fw.get2() + ((t_rw == null) ? 0 : t_rw.get2());
            int clipping = t_fw.get0() + t_fw.get1() + ((t_rw == null) ? 0 : (t_rw.get0() +  t_rw.get1()));

            pw.printf("%s\tmm:%d\tclipping:%d\tnsplit:%d", rid, mm, clipping, splits.size());



            if(info.in_genes == null || info.in_genes.size() == 0)
            {
                pw.printf("\tgcount:0\tgdist:%d\tantisense:%s", (info.closest_strandspec == null) ? Integer.MIN_VALUE : info.closest_strandspec.getFirst(),
                        (info.ANTISENSE_GENES != null && info.ANTISENSE_GENES.size() > 0));
                writePCR(pw);
                return;

            }
            if(info.TR_GENES.size() > 0)
            {
                pw.printf("\tgcount:%d\t%s", info.TR_GENES.size(),
                        StringUtils.joinObjects("|",

                                info.TR_GENES.entrySet(), (e) -> String.format("%s,%s:%s",
                                        e.getKey().id, e.getKey().biotype, StringUtils.join(",", toSortedVector(e.getValue(), true)))
                        ));

                writePCR(pw);
                return;
            }
            if(info.MERGED_GENES.size() > 0)
            {
                pw.printf("\tgcount:%d\t%s", info.MERGED_GENES.size(), StringUtils.joinObjects("|", info.MERGED_GENES, (g) -> String.format("%s,%s:MERGED", g.id, g.biotype)));
                writePCR(pw);
                return;
            }

            pw.printf("\tgcount:%d\t%s",
                    info.INTRONIC_GENES.size(),
                    StringUtils.joinObjects("|", info.INTRONIC_GENES, (g) -> String.format("%s,%s:INTRON", g.id, g.biotype)));
            writePCR(pw);




        }
    }


    static class Statistics
    {
        CountUpdater skips = new CountUpdater("skipped reads");
        int nprimary;
        int non_primary;
        int ntranscriptomic = 0;

        public void write(PrintWriter pw, String prefix)
        {
            pw.printf("%sprimary: %d\n%snon-primary: %d\n%stranscriptomic: %d\n", prefix, nprimary, prefix, non_primary, prefix, ntranscriptomic);
            skips.printBestN(pw, -1);
        }
    }

    static class CStat
    {
        Tuple3<Statistics, Statistics, Statistics> stat = Tuple3.create(new Statistics(), new Statistics(), new Statistics());
        Statistics pairstat = new Statistics();

        Statistics getRStat(SAMRecord sr)
        {
            int sidx = (!sr.getReadPairedFlag()) ? 0 : (sr.getFirstOfPairFlag()) ? 1 : 2;

            return (Statistics) stat.get(sidx);
        }

        public String toString()
        {
            BufferPrintWriter bpw = new BufferPrintWriter();
            bpw.printf("unpaired stat\n");
            stat.get0().write(bpw.getWriter(), "\t");
            bpw.printf("fw stat\n");
            stat.get1().write(bpw.getWriter(), "\t");
            bpw.printf("rw stat\n");
            stat.get2().write(bpw.getWriter(), "\t");

            bpw.printf("pair stat\n");
            pairstat.write(bpw.getWriter(), "\t");

            return bpw.getBuffer();
        }
    }

    public static class AnnotKey extends Tuple3<RegionVector, Boolean, RegionVector>
    {
        AnnotKey(RegionVector read_rv, Boolean strandness, RegionVector split_rv)
        {
            super(read_rv, strandness, split_rv);
        }
    }

    BiConsumer<AnnotKey, ReadAnnotInfo<A>> pcr_handler = (_a, _b) -> {};



    static class Cache<A>
    {
        boolean ignore_fragment_internal_structure = false;
        BiConsumer<AnnotKey, ReadAnnotInfo<A>> pcr_handler = null;
        HashSet<Integer> known = new HashSet<>();
        PriorityQueue<Integer> count =new PriorityQueue<>();
        //HashMap<Integer, Integer> count = new HashMap<>();
        HashMap<Integer, HashSet<AnnotKey>> st2rv = new HashMap<>();
        TreeMap<AnnotKey, ReadAnnotInfo<A>> rv2ar = new TreeMap<>();


        int nremoved = 0;

        long ntimes_removed = 0;
        Logger log = LogConfig.getLogger();

        long in_remove = 0;
        void addPos(String chr, int pos)
        {
            if(count.size() > 0 && count.peek() == pos)
                return;

            if(known.contains(pos))
                return;

            //System.out.printf("will add %s:%d\n", chr, pos);
            count.add(pos);
            known.add(pos);
        }
        void clear()
        {
            nremoved = 0;
            known.clear();
            remove(null, -1);
            count.clear();
            st2rv.clear();
            rv2ar.clear();
        }

        void add(String chr, ReadAnnotInfo ar)
        {
            addPos(chr, ar.merged.getX2());
            MapBuilder.update(st2rv, ar.merged.getX2(), ar.key);
            rv2ar.put(ar.key, ar);
        }


        void remove(String chr, int startpos)
        {
            if(count.size() >0 && (startpos >= 0 && count.peek() >= startpos))
                return;

            long t1 = System.currentTimeMillis();

            while(count.size() > 0 && (startpos < 0 || count.peek() < startpos))
            {
                int pos = count.poll();

                known.remove(pos);

                HashSet<AnnotKey> h = st2rv.remove(pos);
                if(h == null)
                    return;

                apply(h, (_rv) -> pcr_handler.accept(_rv, rv2ar.remove(_rv)));
            }
            long t2 = System.currentTimeMillis();
            in_remove += (t2 - t1);
            nremoved = 0;
            if(++ntimes_removed % 1_000_000 == 0)
            {
                log.info("cache clearing invoked %d times took: %2.2f sec\n", ntimes_removed, (in_remove / 1_000.0));
            }
        }



    }


    class AnnotIterator implements Iterator<AnnotatedRead>
    {

        HashMap<A, CStat> stats = new HashMap<>();

        Cache cache = (NO_PCR_DOWNSAMPLING) ? null : new Cache();


        String lastchr = null;

        HashSet<Pair<String, A>> ignored = new HashSet<>();

        HashMap<Pair<String, A>, SAMRecord> id2read = new HashMap<>();

        int peak_open = 0;
        int peak_cache = 0;

        long t1 = System.currentTimeMillis();
        long nproc = 0;

        int nar = 0;
        AnnotatedRead<A> nextRead = null;

        AnnotIterator()
        {
            this(null, null, null);
        }

        AnnotIterator(String chr, Integer start, Integer end)
        {
            if(!NO_PCR_DOWNSAMPLING) {
                cache.pcr_handler = pcr_handler;
                cache.clear();
                cache.ignore_fragment_internal_structure = ignore_fragment_internal_structure;
            }

            log.info("start reading bam iterator on %d bams", bams.size());
            apply(bams, (_b) -> stats.put(_b.obj, new CStat()));


            if(chr != null && start != null && end != null)
            {
                //retrive reads from a region
                apply(bams, (_csi) -> _csi.query(chr, start, end));
            }
            nextRead = readNext();

            invokePCRHandler();

        }

        void invokePCRHandler() {
            if(nextRead != null && nextRead.info != null && cache == null && pcr_handler != null) {
                pcr_handler.accept(nextRead.info.key, nextRead.info);
            }

        }

        public boolean hasNext()
        {
            return nextRead != null;
        }

        public AnnotatedRead next()
        {
            if(!hasNext())
                return null;

            AnnotatedRead rv = nextRead;

            nextRead = readNext();
            invokePCRHandler();
            return rv;
        }
        public AnnotatedRead readNext() {

            nar++;

            while(bams.size() > 0)
            {
                //log.info("peaks : %s", bams);
                ClosableSamIterator<A> csi = bams.poll();
                if(!csi.pit.hasNext())
                    continue;

                SAMRecord sr = csi.pit.next().sr;
                String readId = correctReadId(sr);

                CStat stat = stats.get(csi.obj);
                decide(csi.it.hasNext(), () -> bams.add(csi), () -> csi.close());

                /*if(sr.getReadName().equals("225303"))
                {
                    System.out.printf("CHECK sr: %s first: %s primary: %s strand: %s\n", sr.getReadName(), sr.getFirstOfPairFlag(), !sr.getNotPrimaryAlignmentFlag(), !sr.getReadNegativeStrandFlag());

                    TreeSet<Region1D> ts = new TreeSet<>();
                    TreeSet<Region1D> td = new TreeSet<>();

                    RegionVector chrv = toRV(sr, ts, td, 0);

                    System.out.printf("read info: %s splits; %s\n", chrv, ts);
                }
                */
                nproc++;
                //log.info("next read: obj: %s: %s: fw: %s %s %d-%d", csi.obj, sr.getReadName(), (!sr.getReadPairedFlag() || sr.getFirstOfPairFlag()), sr.getReferenceName(),  sr.getAlignmentStart(), sr.getAlignmentEnd());
                Statistics st = stat.getRStat(sr);

                if (sr.getNotPrimaryAlignmentFlag()) {

                    st.non_primary++;
                    continue;
                }
                st.nprimary++;
                final String chr = sr.getReferenceName();

                if(lastchr == null || !chr.equals(lastchr))
                {
                    log.info("will switch to next chr: %s peak: %s", chr, bams);
                    lastchr = chr;
                    if(cache != null) {
                        cache.clear();
                    }

                    id2read.clear();

                }
                if(false) {
                    if(nproc % 10_000_000 == 0) {
                        long t2 = System.currentTimeMillis();

                        log.info("after: %2.2f Mio records took: %2.2f sec ignored: %d peak open: %d current: %d peak cache: %d current: %d chr: %s lastchr: %s\n",
                                nproc / 1_000_000.0, (t2 - t1) / 1_000.0, ignored.size(),
                                peak_open, id2read.size(),
                                peak_cache, (cache == null) ? 0 : cache.rv2ar.size(),
                                sr.getReferenceName(), lastchr);
                    }

                    continue;
                }


                if(!sr.getReadPairedFlag())
                {
                    AnnotatedRead ar = new AnnotatedRead(csi.obj, annot, cache, csi.strandness, sr, null);


                    if(ar != null && ar.info != null && ar.info.TR_GENES.size() > 0) {
                        st.ntranscriptomic++;
                    }
                    //log.info("next read: %s %s key: %d", sr.getReadName(), ar.merged, ar.merged.getX2());
                    if(cache != null) {
                        cache.remove(sr.getReferenceName(), ar.second_start);
                        peak_cache = Math.max(peak_cache, cache.rv2ar.size());
                    }

                    return ar;
                }

                if (sr.getMateUnmappedFlag()) {

                    st.skips.update("mate unmapped", 1);
                    continue;
                }

                st = stat.pairstat;

                if (sr.getMateNegativeStrandFlag() == sr.getReadNegativeStrandFlag()) {
                    st.skips.update("mate on same strand", 1);
                    continue;
                }

                if (sr.getReferenceIndex() - sr.getMateReferenceIndex() != 0) {

                    st.skips.update("mate on diff chr", 1);
                    continue;
                }

                final Pair<String, A> rid = Pair.create(readId, csi.obj);


                if (ignored.contains(rid))
                    continue;

                SAMRecord mate = id2read.remove(rid);



                if (mate != null)
                {
                    //ReadRecord rr = new ReadRecord(mate);
                    st.nprimary++;

                    AnnotatedRead ar = new AnnotatedRead(csi.obj, annot, cache, csi.strandness, mate, sr);

                    if(ar != null && ar.info != null && ar.info.TR_GENES.size() > 0) {
                        st.ntranscriptomic++;
                    }

                    //log.info("next read: %s %s key: %d", sr.getReadName(), ar.merged, ar.merged.getX2());
                    if(cache != null) {
                        cache.remove(sr.getReferenceName(), ar.second_start);
                        peak_cache = Math.max(peak_cache, cache.rv2ar.size());

                    }
                    return ar;

                }


                int rstart = Math.min(sr.getAlignmentStart(), sr.getMateAlignmentStart());
                int rend = Math.max(sr.getAlignmentStart(), sr.getMateAlignmentStart());

                //MapBuilder.update(cache.count, rstart);



                //cache.addPos(chr, rend);


                boolean fw_strand = (sr.getFirstOfPairFlag()) ? !sr.getReadNegativeStrandFlag() : !sr.getMateNegativeStrandFlag();


                Boolean strand = (csi.strandness == null) ? null : (csi.strandness) ? fw_strand : !fw_strand;


                if (checkSkipDueContainment(annot, chr, rstart, rend, strand)) {

                    st.skips.update("contains a whole gene", 1);
                    ignored.add(rid);
                    //log.info("last anread: %s new is gene overlapping %s: %s:%d-%s:%d", nextRead, rid, sr.getReferenceName(), sr.getAlignmentStart(), sr.getMateReferenceName(), sr.getMateAlignmentStart());
                    //SAMRecord f = (sr.getAlignmentStart() < sr.getMateAlignmentStart())
                    return new AnnotatedRead(csi.obj, sr, sr.getMateAlignmentStart(), strand, annot);
                }
                id2read.put(rid, sr);

                peak_open = Math.max(peak_open, id2read.size());


                if (nproc % 1_000_000 != 0)
                    continue;

                long t2 = System.currentTimeMillis();

                log.info("after: %2.2f Mio records took: %2.2f sec ignored: %d peak open: %d current: %d peak cache: %d current: %d current chr: %s -> %s\n",
                        nproc / 1_000_000.0, (t2 - t1) / 1_000.0, ignored.size(),
                        peak_open, id2read.size(),
                        peak_cache, (cache == null) ? 0 : cache.rv2ar.size(),
                        chr, correctChr(chr, annot)
                );

                System.out.printf("%s\n", stat.toString());


            }


            if(cache != null)
                cache.clear();

            long t2 = System.currentTimeMillis();
            System.out.printf("finished: total of %2.2f Mio records took: %2.2f sec ignored: %d peak open: %d current: %d peak cache: %d current: %d\n",
                    nproc / 1_000_000.0, (t2 - t1) / 1_000.0, ignored.size(),
                    peak_open, id2read.size(),
                    peak_cache, (cache == null) ? 0 : cache.rv2ar.size());
            return null;
        }

    }

    static boolean checkSkipDueContainment(IsoformRegionGetter annot, String chr, int rstart, int rend, Boolean strand)
    {
        UPair<Integer> up = checkGeneContainment(annot, chr, rstart, rend, strand);

        return (up.getSecond() == 0 && up.getFirst() > 0);
    }

    static UPair<Integer> checkGeneContainment(IsoformRegionGetter annot, String chr, int rstart, int rend, Boolean strand)
    {
        Vector<Tuple3<String, Boolean, Region1D<String>>> reginfos = annot.getRegionInfos(chr, rstart, rend);

        int n_containsgenes = filteredSize(reginfos, (_t) -> (strand == null || strand == _t.get1()) && _t.get2().getX1() >= rstart && _t.get2().getX2() <= rend);
        int n_genes = filteredSize(reginfos, (_t) -> (strand == null || strand == _t.get1()) && _t.get2().getX1() <= rstart && _t.get2().getX2() >= rend);

        return UPair.createU(n_containsgenes, n_genes);

    }
    /*
    public Iterator<AnnotatedRead> getAnnotIterator(Iterator<SAMRecord> it)
    {
        return new AnnotIterator(it);
    }
    */

    public Iterator<AnnotatedRead> getMultiBamIterator()
    {
        return new AnnotIterator();
    }

    public Iterator<AnnotatedRead> getMultiBamIterator(String chr, int start, int end)
    {
        return new AnnotIterator(chr, start, end);
    }
    static class GeneCountInfo<A>
    {
        MultiIsoformRegion gene;
        Vector<ReadAnnotInfo<A>> counts = new Vector<>();
        Vector<ReadAnnotInfo<A>> ov_counts = new Vector<>();

        public GeneCountInfo(MultiIsoformRegion g)
        {
            this.gene = gene;
        }

        public void add(ReadAnnotInfo<A> rai, boolean multi)
        {
            ((multi) ? ov_counts : counts).add(rai);
        }
    }

    static class GeneCountInfos<A>
    {
        PriorityQueue<Integer> gene_pos_lookup =new PriorityQueue<>();
        HashMap<Integer, HashSet<GeneCountInfo<A>>> gene_counts = new HashMap<>();
        HashMap<MultiIsoformRegion, GeneCountInfo<A>> g2count = new HashMap<>();

        String lastchr = null;

        void clear(int pos)
        {
            while(pos >= 0 || gene_pos_lookup.size() > 0 && gene_pos_lookup.peek() < pos)
            {
                HashSet<GeneCountInfo<A>> h = gene_counts.remove(gene_pos_lookup.poll());
                //handle the genes
            }

        }
        void add(AnnotatedRead<A> read)
        {
            final int start = read.merged.getX1();

            clear(read.merged.getX1());

            Collection<MultiIsoformRegion> genes = null;
            switch(read.getType())
            {
                case TRANSCRIPTOMIC:
                    genes = read.info.TR_GENES.keySet();
                    break;
                case MERGED_TR:
                    genes = read.info.MERGED_GENES;
                    break;
                case INTRONIC:
                    genes = read.info.INTRONIC_GENES;
                    break;

            }
            if(genes == null)
                return;

            final boolean multi = genes.size() > 1;
            apply(genes, (g) ->
            {
                GeneCountInfo<A> gci = g2count.get(g);
                if(gci == null)
                {
                    g2count.put(g, gci = new GeneCountInfo<A>(g));
                    HashSet<GeneCountInfo<A>> h = gene_counts.get(g.end);
                    if( h == null)
                    {
                        gene_counts.put(g.end, h = new HashSet<>());
                        gene_pos_lookup.add(g.end);
                    }
                    h.add(gci);

                }
                gci.add(read.info, multi);
            });


        }

    }

    static class CountInfo<A>
    {
        int nprocessed = 0;
        A idx;
        HashMap<READTYPE, Integer> counts = new HashMap<>();


        HashMap<Integer, Integer> pcr2freq = new HashMap<>();

        HashMap<Pair<String, READTYPE>, Integer> g2count = new HashMap<>();
        HashMap<Pair<String, READTYPE>, Vector<Integer>> g2readid = new HashMap<>();

        CountInfo(A idx)
        {
            this.idx = idx;
        }

        public String toString()
        {
            return String.format("%s %s PCR: %s", idx, counts, pcr2freq);
        }
    }


    public void setPcrHandler(BiConsumer<AnnotKey, ReadAnnotInfo<A>> handler)
    {
        pcr_handler = handler;
    }

    static<A> HashMap<A, CountInfo> process(BAMIterator<A> bam, HashMap<A, PrintWriter> pwmap, boolean verbose)
    {
        return process(bam, pwmap, false, verbose);
    }

    static<A> HashMap<A, CountInfo> process(BAMIterator<A> bam, HashMap<A, PrintWriter> pwmap, boolean write_readids, boolean verbose)
    {

        HashMap<A, CountInfo> counts = buildMap(bam.bams, (_c) -> _c.obj, (_c) -> new CountInfo(_c.obj));

        bam.pcr_handler = (_rv, _ri) -> apply(_ri.pcr_index.entrySet(), (Map.Entry<A, Integer> _e) -> MapBuilder.update(counts.get(_e.getKey()).pcr2freq, _e.getValue()));
        int[] nproc = new int[4];

        Iterator<AnnotatedRead> it = bam.getMultiBamIterator();

        Logger log = LogConfig.getLogger();

        PrintWriter vspw = (!verbose) ? null : new PrintWriter(System.out);
        while (it.hasNext()) {

            AnnotatedRead ar = it.next();

            if(verbose)
            {
                ar.check(vspw, null);
            }
            CountInfo ci = counts.get(ar.obj);

            PrintWriter pw = (pwmap == null) ? null : pwmap.get(ar.obj);

            //log.info("next ar %s %s%c %s type: %s", ar, ar.chr, GenomicUtils.getStrand(ar.strand), ar.merged, ar.getType());

            //System.err.printf("errtest: %s\n", ar);

            MapBuilder.update(ci.counts, ar.getType());

            if(ar.info != null) {
                Collection<MultiIsoformRegion> genes = ar.info.TR_GENES.keySet();

                if (genes.size() == 0) {
                    genes = ar.info.MERGED_GENES;
                }

                if (genes.size() == 0) {
                    genes = ar.info.INTRONIC_GENES;
                }

                for (MultiIsoformRegion mir : genes) {
                    MapBuilder.update(ci.g2count, Pair.create(mir.id, ar.getType()));
                    if(write_readids)
                    {
                        int rid = Integer.parseInt(ar.fw.getReadName());
                        MapBuilder.updateV(ci.g2readid, Pair.create(mir.id, ar.getType()), rid);
                    }
                }
            }


            if (++ci.nprocessed % 1_000_000 == 0) {
                log.info("%s processed %2.2f Mio : %s", ci.idx, ci.nprocessed / 1_000_000.0, ci.counts);
            }
            if (pw == null)
                continue;

            switch (ar.getType()) {
                case GENE_OVERLAPPING:
                    break;
                default:
                    ar.write(pw);
                    break;
            }
        }

        if(vspw != null)
        {
            vspw.flush();
        }
        return counts;
    }

    public enum STRANDNESS
    {
        POS("true"),
        NEG("false"),
        UNSPEC("")
        ;

        String name;
        STRANDNESS(String n)
        {
            this.name = n;
            EnumGetter.add(n, this);
        }

        public Boolean toBool()
        {
            switch (this)
            {
                case POS:
                    return true;
                case NEG:
                    return false;
                case UNSPEC:
                    return null;
            }
            return null;
        }

        public static STRANDNESS get(String s)
        {
            return EnumGetter.get(s, false);
        }
    }
    public static void main(String[] args)
    {
        SimpleOptionParser cmd = new SimpleOptionParser("gtf", "bam", "o", "frstrand", "check", "onlycount", "table", "getstrandness", "nreadtocheck", "verbose",
                "writereadids", "getrrnareads");
        cmd.setFile("gtf", "bam", "check");
        cmd.setOptional("frstrand");
        cmd.setInt("nreadtocheck");
        cmd.getOption("frstrand").setHelp(" true or false");
        cmd.setOptional("check", "bam", "frstrand", "table");
        cmd.setSwitches("onlycount", "getstrandness", "verbose", "writereadids", "getrrnareads");
        cmd.setDefault("nreadtocheck", "100000");

        if(!OptionParser.parseParams(args, true, false, true, true, cmd))
            return;

        boolean writereadids = cmd.isSet("writereadids");
        boolean verbose = cmd.isSet("verbose");
        File cf = cmd.getOptionalFile("check");
        HashMap<String, String> check = (cf == null) ? null : buildMap(FileUtils.getLineIterator(cf), (l) -> l.substring(0, l.indexOf('\t')), (l) -> l);

        String strand = cmd.getOptionalValue("frstrand", null);

        if(cmd.isSet("getrrnareads")) {
            BAMIterator<Integer> bam = new BAMIterator<>(new GFFBasedIsoformRegionGetter(cmd.getFile("gtf"), null, null));
            bam.addBAM(0, cmd.getFile("bam"), null, null);

            int [] counts = new int[3];
            Iterator<AnnotatedRead> it = bam.getMultiBamIterator();
            CountUpdater cu = new CountUpdater();
            int REPORT = 10_000;

            while(it.hasNext())
            {
                AnnotatedRead<Integer> ar = it.next();

                counts[0]++;
                switch(ar.getType())
                {
                    case TRANSCRIPTOMIC:
                        break;
                    default:
                        continue;
                }
                if(ar.info.TR_GENES.size() > 1)
                    continue;


                MultiIsoformRegion mir = first(ar.info.TR_GENES.keySet());

                cu.update(mir.biotype, 1);


                if(++counts[1] % REPORT != 0)
                    continue;

                System.out.printf("after %.2f/%.2f (%2.2f%%) Mio tr reads got: \n", counts[1] / 1_000_000.0, counts[0] / 1_000_000.0,
                        (100.0 * counts[1]) / counts[0]);

                cu.printBestN(10);

            }

            System.out.printf("READY after %.2f/%.2f (%2.2f%%) Mio tr reads got: \n", counts[1] / 1_000_000.0, counts[0] / 1_000_000.0,
                    (100.0 * counts[1]) / counts[0]);

            cu.printBestN(10);

            return;
        }
        if(cmd.isSet("getstrandness"))
        {
            int tocheck = cmd.getInt("nreadtocheck");
            BAMIterator<Integer> bam = new BAMIterator<>(new GFFBasedIsoformRegionGetter(cmd.getFile("gtf"), null, null));
            bam.addBAM(0, cmd.getFile("bam"), null, null);



            int [] counts = new int[3];
            Iterator<AnnotatedRead> it = bam.getMultiBamIterator();
            while(it.hasNext())
            {
                AnnotatedRead<Integer> ar = it.next();

                switch(ar.getType())
                {
                    case TRANSCRIPTOMIC:
                        break;
                    default:
                        continue;
                }
                if(ar.info.TR_GENES.size() > 1)
                    continue;


                MultiIsoformRegion mir = first(ar.info.TR_GENES.keySet());
                counts[(mir.strand == ar.fr_strand) ? 1 : 2]++;

                if(++counts[0] >= tocheck)
                    break;


            }
            double diff = NumUtils.logN((counts[1] + 1.0) / (counts[2] + 1.0), 2.0);

            String label = (diff > 1.5) ? STRANDNESS.POS.name : (diff < -1.5) ? STRANDNESS.NEG.name : STRANDNESS.UNSPEC.name;

            System.out.printf("after %d reads got +: %d -: %d : %2.3f %s\n", counts[0], counts[1], counts[2], diff, label);
            PrintWriter pw = cmd.getWriter("o");
            pw.println(label);
            pw.close();

            return;
        }


        boolean onlycount = cmd.isSet("onlycount");


        Vector<CountInfo> counts = new Vector<>();
        if(cmd.isOptionSet("bam")) {
            BAMIterator<Integer> bam = new BAMIterator<>(new GFFBasedIsoformRegionGetter(cmd.getFile("gtf"), null, null));
            bam.gobi_verbose = verbose;
            HashMap<Integer, PrintWriter> outmap = (onlycount) ? null : new HashMap<>();
            bam.addBAM(0, cmd.getFile("bam"), null, (strand == null) ? null : Boolean.parseBoolean(strand));
            if(!onlycount)
            {
                outmap.put(0, cmd.getWriter("o"));
            }
            counts.addAll(process(bam, outmap, writereadids, verbose).values());
            if(outmap != null)
            {
                apply(outmap.values(), (_pw) -> _pw.close());
            }
            if(onlycount)
            {

                CountInfo ci = first(counts);
                PrintWriter pw = cmd.getWriter("o");
                if(writereadids)
                {

                    apply(toSortedVector(ci.g2readid.keySet(), true), (Pair<String, READTYPE> _k) -> pw.printf("%s\t%s\t%s\n",
                            _k.getFirst(), _k.getSecond(), StringUtils.joinObjects(",", (Vector<Integer>)ci.g2readid.get(_k))));
                }
                else {
                    apply(ci.g2count.entrySet(), (Map.Entry<Pair<String, READTYPE>, Integer> _e) -> pw.printf("%s\t%s\t%d\n", _e.getKey().getFirst(), _e.getKey().getSecond(), _e.getValue()));
                    apply(ci.counts.entrySet(), (Map.Entry<READTYPE, Integer> _e) -> pw.printf("_\t%s\t%d\n", _e.getKey(), _e.getValue()));
                }
                pw.close();

            }
            return;
        }

        File table = cmd.getOptionalFile("table");
        if(table != null)
        {
            BAMIterator<String> bam = new BAMIterator<>(new GFFBasedIsoformRegionGetter(cmd.getFile("gtf"), null, null));
            DataTable dt = DataTable.readFile(table, "\t");

            HashMap<String, PrintWriter> outmap = (onlycount) ? null : new HashMap<>();
            File od = cmd.getFile("o");

            if(!od.exists() || ! od.isDirectory())
            {
                System.err.printf("-o must be a directory if you set -table!\n");
                return;
            }

            for(int i=0; i<dt.getSize(); i++)
            {
                ObjectGetter og = dt.getDataAt(i);
                String id = og.getString("id");
                bam.addBAM(id, og.getPath("bam"), null, (og.getString("strandness").trim().length() == 0) ? null : og.getBoolean("strandness"));
                if(!onlycount)
                {
                    outmap.put(id, FileUtils.getWriter(od, "%s.annot",  id));
                }

            }
            counts.addAll(process(bam, outmap, false).values());
        }


        apply(counts, (_c) -> System.out.println(_c));


    }
}

