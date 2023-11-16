package empires.rnaseq.reads;


import empires.rnaseq.MultiIsoformRegion;
import lmu.utils.*;
import lmu.utils.tuple.Tuple3;

import java.io.File;
import java.io.PrintWriter;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;

import static lmu.utils.ObjectGetter.*;

public class TranscriptEQClassWriter {

    Logger log = LogConfig.getLogger();
    HashMap<Integer, Integer> ambig2count = new HashMap<>();

    Vector<String> ids;
    HashMap<String, Integer> id2idx;

    HashMap<Tuple, String> eqClass2Gene = new HashMap<>();
    HashMap<Tuple, Vector<Tuple>> multiGene2Parts = new HashMap<>();

    HashMap<Tuple, DownSampleVector> counts = new HashMap<>();

    empires.rnaseq.IsoformRegionGetter irg;
    int nonTranscript = 0;
    int nTranscript = 0;
    ObjectGetter.MapGetter<Integer, Double> downsampler = (_d) -> (_d == null) ? 0 : NumUtils.logN(_d + 1, 2.0);
    PrintWriter pw;

    DownSampleVector EMPTY;

    String lastchr = null;
    TranscriptEQClassWriter(empires.rnaseq.IsoformRegionGetter irg, Vector<String> ids, PrintWriter pw) {
        this.pw = pw;
        this.ids = ids;
        id2idx = buildIndexMap(this.ids);
        EMPTY = new DownSampleVector(ids.size());
        this.irg = irg;
    }

    String downsample2str(DownSampleVector dsv) {
        Vector<Integer> reads = mapIndex(ids.size(), (_i) -> (int)dsv.reads[_i]);
        return String.format("total: %d  :: %s", NumUtils.sum(reads), reads);
    }

    double corr(DownSampleVector dsv1, DownSampleVector dsv2) {
        Vector<Double> reads1 = mapIndex(ids.size(), (_i) -> 1.0 + dsv1.reads[_i]);
        Vector<Double> reads2 = mapIndex(ids.size(), (_i) -> 1.0 + dsv2.reads[_i]);
        return NumUtils.corr(reads1, reads2);
    }

    static Vector<Pair<String, Function<DownSampleVector, double[]>>> DSvariants =
            toVector(
                    Pair.create("reads", (_dsv) -> _dsv.reads)
                    //,Pair.create("frags", (_dsv) -> _dsv.frags)
                    //,Pair.create("ds", (_dsv) -> _dsv.downsampled)
            );

    static Function<double[], String> fmt = (_v) -> StringUtils.joinObjects("\t", mapIndex(_v.length, (_i) -> String.format("%.2f",_v[_i]) ));


    static final int MIN_READS_FOR_AMBIG = 10;
    int total_ambig_reads = 0;
    int[] resolved = new int[2];

    public void writeChromosome() {
        if(lastchr == null)
            return;

        log.info("chr: %s got %d eq-classes, %d ambigs so far %d/%d non-transcript uniq reads", lastchr, counts.size(), multiGene2Parts.size(), nonTranscript, nTranscript);
        if(multiGene2Parts.size() > 0)
            log.info("ambig size distrib: %s\n", NumUtils.getNumInfo(multiGene2Parts.values(), (_v) -> filteredSize(_v, (_t) -> counts.containsKey(_t))).getInfoWithQ());

        for(Map.Entry<Tuple, Vector<Tuple>> e : multiGene2Parts.entrySet()) {
            DownSampleVector ambig = counts.get(e.getKey());
            int ntotal = NumUtils.sum(mapIndex(ids.size(), (_i) -> (int)ambig.reads[_i]));
            if(ntotal < MIN_READS_FOR_AMBIG)
                continue;

            Vector<Tuple> valids = filter(e.getValue(), (_t) -> counts.containsKey(_t));
            if(valids.size() == 0)
                continue;

            if(valids.size() == 1) {
                DownSampleVector other = counts.get(valids.get(0));
                other.update(ambig);
                continue;
            }


            total_ambig_reads += ntotal;

            System.out.printf("ambig of %s (%s)\n\t%s\n", valids.size(), map(valids, (_t) -> eqClass2Gene.get(_t)), downsample2str(ambig));
            Vector<Tuple3<Double, Integer, Tuple>> candidates = new Vector<>();
            for(Tuple t : valids) {
                DownSampleVector other = counts.get(t);
                int cand_total = NumUtils.sum(mapIndex(ids.size(), (_i) -> (int)other.reads[_i]));
                double corr = corr(ambig, other);
                candidates.add(Tuple3.create(corr, cand_total, t));
                System.out.printf("\tcheck %s: corr: %.2f %s\n", eqClass2Gene.get(t), corr, downsample2str(other));
            }

            NumUtils.sort(candidates, (_t) -> _t.get0(), true);

            if(candidates.get(0).get0() > 0.5 && candidates.get(0).get0() - candidates.get(1).get0() > 0.2) {
                //select first candidate
                resolved[0] += ntotal;
                counts.get(candidates.get(0).get2()).update(ambig);
                continue;
            }
            NumUtils.sort(candidates, (_t) -> _t.get1(), true);
            if(candidates.get(0).get1() > 3 * candidates.get(1).get1()) {
                //select first candidate
                resolved[1] += ntotal;
                counts.get(candidates.get(0).get2()).update(ambig);
                continue;
            }

        }

        HashMap<String, Vector<Tuple>> gene2tuple = new HashMap<>();
        for(Tuple t : counts.keySet()) {
            String gene = eqClass2Gene.get(t);
            if(gene == null)
                continue;
            MapBuilder.updateV(gene2tuple, gene, t);
        }

        log.info("will write %d genes for chr %s\n", gene2tuple.size(), lastchr);
        for(String gene : gene2tuple.keySet()) {

            //Armin
            RegionVector merged = RegionVector.merge(irg.getRegionById(gene).isoforms.values());

            Vector<Tuple> eqClasses = gene2tuple.get(gene);
            DownSampleVector Tdsv = new DownSampleVector(ids.size());
            apply(eqClasses, (_t) -> Tdsv.update(counts.get(_t)));
            pw.printf(">%s\tTOTAL\n", gene);

            apply(DSvariants, (_p) ->
            {
                double[] d = _p.getSecond().apply(Tdsv);
                pw.printf("%s\t%s\n", _p.getFirst(), fmt.apply(d));

            });

            for(Tuple t : eqClasses) {
                DownSampleVector dsv = counts.get(t);
//                pw.printf(">%s\t%s\n", gene, StringUtils.joinObjects(",", mapIndex(t.cardinality(), (_i) -> t.getAsString(_i))));

                //Armin
                List<String> tr_names = Arrays.stream(t.values()).map(Object::toString).collect(Collectors.toList());
                Set<String> other_names = new HashSet<>(irg.getRegionById(gene).isoforms.keySet());
                other_names.removeAll(tr_names);
                List<RegionVector> test = tr_names.stream().map(_tr -> irg.getRegionById(gene).isoforms.get(_tr)).collect(Collectors.toList());
                RegionVector intersect = test.get(0);
                for (RegionVector rv : test) {
                    intersect = intersect.intersect(rv);
                }

                RegionVector regions = new RegionVector(intersect);

//                RegionVector regions = RegionVector.merge(test);
                RegionVector to_remove = RegionVector.merge(other_names.stream().map(_tr -> irg.getRegionById(gene).isoforms.get(_tr)).collect(Collectors.toList()));

                if (other_names.size() > 0) {
                    regions = regions.substract(to_remove);
//                    regions = RegionVector.convertToInner(merged, irg.getRegionById(gene).strand, regions.substract(to_remove));
                }
//                else {
//                    regions = RegionVector.convertToInner(merged, irg.getRegionById(gene).strand, regions);
//                }
                pw.printf(">%s\t%s;%s\n", gene, StringUtils.joinObjects(",", mapIndex(t.cardinality(), (_i) -> t.getAsString(_i))), StringUtils.joinObjects("", regions.getRegions()));

//                if (gene.equals("ENSG00000143842")) {
//                    System.out.printf(">%s\t%s;%s\n", gene, StringUtils.joinObjects(",", mapIndex(t.cardinality(), (_i) -> t.getAsString(_i))), StringUtils.joinObjects("", regions.getRegions()));;
//                }
                apply(DSvariants, (_p) ->
                {
                    double[] d = _p.getSecond().apply(dsv);
                    pw.printf("%s\t%s\n", _p.getFirst(), fmt.apply(d));

                });

            }

        }
        pw.flush();

        log.info("after %s total ambig reads: %.4f Mio resolved by corr: %.4f Mio by size: %.4f Mio ", lastchr,
                total_ambig_reads / 1_000_000.0, resolved[0] / 1_000_000.0, resolved[1] / 1_000_000.0);
        eqClass2Gene.clear();
        multiGene2Parts.clear();
        counts.clear();
    }

    public void close() {
        writeChromosome();
        pw.close();
        lastchr = null;
    }

    boolean use_gene_ambigous = false;
    public void process(BAMIterator.AnnotKey akey, BAMIterator.ReadAnnotInfo<String> annot) {
        if(annot.TR_GENES.size() == 0) {
            nonTranscript++;
            return;
        }

        nTranscript++;

        if(lastchr == null) {
            lastchr = annot.chr;
        }
        if(!lastchr.equals(annot.chr)) {
            writeChromosome();
        }

        lastchr = annot.chr;

        Tuple dataKey = null;

        if(!use_gene_ambigous && annot.TR_GENES.size() > 1)
            return;


        Vector<Tuple> keys = new Vector<>();
        for(Map.Entry<MultiIsoformRegion, HashSet<String>> e :  annot.TR_GENES.entrySet()) {
            Vector<String> trKey = toSortedVector(e.getValue(), true);
            Tuple key = Tuple.tupleFromCollection(trKey);
            eqClass2Gene.put(key, e.getKey().id);
            keys.add(key);
        }
        if(keys.size() > 1) {
            Vector<Tuple> sortedMultiKey = toSortedVector(keys, true);
            dataKey = Tuple.tupleFromCollection(sortedMultiKey);
            multiGene2Parts.put(dataKey, sortedMultiKey);

        } else {
            dataKey = keys.get(0);
        }

        DownSampleVector dsv = counts.get(dataKey);
        if(dsv == null) {
            counts.put(dataKey, dsv = new DownSampleVector(ids.size()));
        }
        DownSampleVector.update(dsv, map(ids, (_s) -> annot.pcr_index.getOrDefault(_s, -1) + 1), downsampler);
        //System.out.printf("%s pcr index: %s\n",annot.TR_GENES, annot.pcr_index);
        MapBuilder.update(ambig2count, annot.TR_GENES.size());
    }

    public static void main(String[] args) {
        SimpleOptionParser cmd = new SimpleOptionParser("gtf", "o",  "table", "basedir");
        cmd.setFile("gtf", "table");
        cmd.setDir("basedir");
        cmd.setOptional("basedir");
        cmd.setOutFile("table");

        if(!OptionParser.parseParams(args, true, false, true, true, cmd))
            return;

        Logger log = LogConfig.getLogger();
        long t1 = System.currentTimeMillis();
        empires.rnaseq.IsoformRegionGetter irg = new empires.rnaseq.GFFBasedIsoformRegionGetter(cmd.getFile("gtf"), null, null);
        BAMIterator.NO_PCR_DOWNSAMPLING = true;
        BAMIterator<String> bam = new BAMIterator<>(irg);
        DataTable dt = DataTable.readFile(cmd.getFile("table"), "\t");

        Set<String> headers = toSet(dt.header.getHeaderNames());
        Vector<String> idalternatives = toVector("id", "label", "ID", "Id", "name");
        String idheader = filterOne(idalternatives, (_s) -> headers.contains(_s));
        if(idheader == null) {
            System.err.printf("did not find an id header in: %s (tested header names: %s got headers: %s)", cmd.getFile("table").getAbsolutePath(),
                    idalternatives, headers);
            return;
        }

        File basedir = cmd.getOptionalFile("basedir");
        HashSet<String> chromosomes = mapToSet(irg.getRegions(), (_r) -> _r.chr);

        Vector<String> ids = new Vector<>();
        for(int i=0; i<dt.getSize(); i++)
        {
            ObjectGetter og = dt.getDataAt(i);
            String id = og.getString(idheader);
            ids.add(id);
            String strandinfo = og.getString("strandness");
            bam.addBAM(id, og.getPath("bam", basedir) , null, (strandinfo == null ||  strandinfo.trim().length() == 0) ? null : og.getBoolean("strandness"));

        }
        TranscriptEQClassWriter teqW = new TranscriptEQClassWriter(irg, ids, cmd.getWriter("o"));
        bam.setPcrHandler((_k, _a) -> teqW.process(_k, _a));


        Iterator<BAMIterator.AnnotatedRead> it =  bam.getMultiBamIterator();
        while(it.hasNext()) {
            BAMIterator.AnnotatedRead ar = it.next();
            //System.out.printf("%s\n", ar);
            //System.out.printf("read: %s %s%c %s: %s\n", ar.rid, ar.chr, GenomicUtils.getStrand(ar.strand), ar.merged, ar.info.TR_GENES);

        }

        teqW.close();
        long t2 = System.currentTimeMillis();
        log.info("took %.2f sek ambig: %s", (t2 - t1) / 1000.0, teqW.ambig2count);
    }
}
