package empires.test.microarray;

import empires.EmpiRe;
import empires.input.MicroArrayInputRaw;
import empires.input.ReplicateSetInfo;
import empires.plotting.DiffExpTable;
import lmu.utils.*;
import lmu.utils.fdr.BenjaminiHochberg;
import lmu.utils.fdr.PerformanceResult;
import lmu.utils.fdr.RocInfo;
import lmu.utils.plotting.CachedPlotCreator;
import lmu.utils.plotting.PlotCreator;
import empires.input.GeneralOptions;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.awt.image.BufferedImage;
import java.io.File;
import java.util.*;
import java.util.function.Function;

import static lmu.utils.IteratorUtils.rangev;
import static lmu.utils.ObjectGetter.*;

public class Simulation {

    static public int getIndexForFoldChange(int baseIndex, double foldchange, Vector<Pair<empires.ReplicatedMeasurement, String>> data, boolean[] masked) {
        double base_meanval = data.get(baseIndex).getFirst().mean;
        int step  = (foldchange > 0) ? -1 : 1;
        Pair<Boolean,Integer> hit = VectorUtils.binarySearch(data, base_meanval + foldchange,  (_p, _d) -> Double.compare(_p.getFirst().mean, _d));
        int hitIdx = hit.getSecond();
        while(hitIdx >=0 && hitIdx < masked.length && masked[hitIdx]) {
            hitIdx += step;
        }
        if(hitIdx < 0 || hitIdx >= masked.length)
            return -1;

        masked[hitIdx] = true;
        return hitIdx;
    }

    static Vector<String> getSortedSubFeatures(empires.NormalizedReplicateSet nrs, String gene) {
        Vector<String> ids = new Vector<>();
        ids.addAll(nrs.getInData().getSubFeatures(gene));
        return NumUtils.sort(ids, (_id) -> nrs.getNormed(_id).mean, false);
    }

    static final String SUMMARIZED_SUFFIX = "_summarized.mat";
    public static HashMap<String, empires.NormalizedReplicateSet> readSummarized(File sumdir, File labels, File mappings, String condition) {
        HashMap<String, empires.NormalizedReplicateSet> rv = new HashMap<>();

        for(File f : FileUtils.getFiles(sumdir, ".+" + SUMMARIZED_SUFFIX)) {
            String name = f.getName().split(SUMMARIZED_SUFFIX)[0];
            empires.input.MicroArrayInputRaw raw = new empires.input.MicroArrayInputRaw(f, labels, mappings, true, new empires.AutoBackGroundContextProvider());
            rv.put(name, raw.getNormalized(condition));
        }
        return rv;
    }

    static class SwapPairs {
        String gene1;
        String gene2;
        Vector<String> oligos1;
        Vector<String> oligos2;
        int num_oligos;
        double fc;

        public SwapPairs(String g1, String g2, int num_oligos, double fc) {
            this.gene1 = g1; this.gene2 = g2; this.num_oligos = num_oligos;
            this.fc = fc;
        }

        public String toString() {
            return String.format("%s <=> %s (%d oligos, fc: %.2f)", gene1, gene2, num_oligos, fc);
        }
    }

    static class SummarizedOligo {
        String oligo;
        Vector<Double> c1;
        Vector<Double> c2;
        double pval;
        double fdr;
        double fc;

        public SummarizedOligo(String name, Vector<Double> c1, Vector<Double> c2) {
            this.oligo = name;
            this.c1 = c1; this.c2 = c2;
            fc = NumUtils.mean(c1) - NumUtils.mean(c2);
            pval = NumUtils.getTTestPvalue(c1, c2);
        }
    }
    static class SummarizedGene {
        String gene;

        Vector<SummarizedOligo> oligo_vals;

        SummarizedOligo selected;

        public SummarizedGene(String gene, Vector<SummarizedOligo> oligos) {
            this.gene = gene;
            this.oligo_vals = oligos;

            selected = NumUtils.minObj(oligos, (_o) -> _o.pval).getSecond();
        }
    }


    static public Vector<SummarizedGene> simulateSummarized(empires.NormalizedReplicateSet rsi, Vector<String> cond1_replicates, Vector<String> cond2_replicates, Vector<SwapPairs> swappairs, Set<String> selected_trues) {
        HashMap<String, Integer> repname2index = buildIndexMap(map(rsi.getInData().getReplicateNames(), (_s) -> _s.split(".CEL")[0]));
        Vector<Integer> cond1indeces = map(cond1_replicates, (_r) -> repname2index.get(_r));
        Vector<Integer> cond2indeces = map(cond2_replicates, (_r) -> repname2index.get(_r));

        Function<String, SummarizedOligo> getReplicateData = (_n) ->
        {
            Vector<Double> ndata = rsi.getNormed(_n).replicates;
            return new SummarizedOligo(_n,
                    map_and_filter(cond1indeces, (_i) -> ndata.get(_i), (_d) -> !Double.isNaN(_d)),
                    map_and_filter(cond2indeces, (_i) -> ndata.get(_i), (_d) -> !Double.isNaN(_d))
            );
        };


        Vector<SummarizedGene> genes = new Vector<>();
        HashMap<Integer, Integer> diffs = new HashMap<>();
        for(SwapPairs sp : swappairs) {
            Vector<String> subfeatures1 = NumUtils.sort(rsi.getInData().getSubFeatures(sp.gene1), (_n) -> rsi.getNormed(_n).mean, true);
            Vector<String> subfeatures2 = NumUtils.sort(rsi.getInData().getSubFeatures(sp.gene2), (_n) -> rsi.getNormed(_n).mean, true);

            int numsub = Math.min(subfeatures1.size(), subfeatures2.size());

            MapBuilder.update(diffs, subfeatures1.size() - subfeatures2.size());
            //System.out.printf("subfeatures: %d vs %d\n", subfeatures1.size(), subfeatures2.size());

            Vector<SummarizedOligo> oligs1 = new Vector<>();
            Vector<SummarizedOligo> oligs2 = new Vector<>();

            for(int j=0; j< numsub; j++) {
                SummarizedOligo s1 = getReplicateData.apply(subfeatures1.get(j));
                SummarizedOligo s2 = getReplicateData.apply(subfeatures2.get(j));

                oligs1.add(new SummarizedOligo(subfeatures1.get(j), s1.c1, s2.c2));
                oligs2.add(new SummarizedOligo(subfeatures2.get(j), s2.c1, s1.c2));

            }

            genes.add(new SummarizedGene(sp.gene1, oligs1));
            genes.add(new SummarizedGene(sp.gene2, oligs2));
        }

        for(String feature : rsi.getInData().getFeatureNamesToTest()) {
            if(selected_trues.contains(feature))
                continue;

            genes.add(new SummarizedGene(feature, map(rsi.getInData().getSubFeatures(feature), (_n) -> getReplicateData.apply(_n))));
        }
        System.out.printf("size diffs: %s\n", diffs);

        BenjaminiHochberg.adjust_pvalues(genes, (_g) -> _g.selected.pval, (_p) -> _p.getFirst().selected.fdr = _p.getSecond());
        return genes;
    }


    public static UPair<Vector<String>> selectReplicates(empires.NormalizedReplicateSet nrs, int sample_per_replicate) {

        Vector<String> replicateNames = nrs.getInData().getReplicateNames();
        Vector<Vector<Double>> data = nrs.getInData().getLog2Data();
        Vector<Integer> selectedIdx = empires.Normalization.getBestSubSample(data, 2 * sample_per_replicate);


        Vector<Vector<Double>> seldata = map(selectedIdx, (_i) -> data.get(_i));

        empires.PairScore[][] errMatrix = empires.Normalization.getPairWiseErrors(seldata);

        HashSet<Integer> oneSide = new HashSet<>();
        empires.PairScore maxPS = null;
        for(empires.PairScore[] psa : errMatrix) {
            for(empires.PairScore ps : psa) {
                if(ps == null || ps.err == null)
                    continue;
                if(maxPS != null && maxPS.err.err >= ps.err.err)
                    continue;

                maxPS = ps;
            }
        }
        oneSide.addAll(maxPS.sg1.getMembers());
        oneSide.addAll(maxPS.sg2.getMembers());

        while(oneSide.size() < sample_per_replicate) {
            int nextPick = -1;
            double nextPickScore = -1.0;

            for(int i=0; i<seldata.size(); i++) {
                if(oneSide.contains(i))
                    continue;

                double score = 0.0;

                for(int src : oneSide) {
                    score += errMatrix[src][i].err.err;
                }
                if(nextPickScore >= score)
                    continue;

                nextPickScore = score;
                nextPick = i;
            }

            oneSide.add(nextPick);
        }
        Vector<String> replicatesOne = filter_and_map(rangev(selectedIdx.size()), (_i) -> oneSide.contains(_i), (_i) -> replicateNames.get(selectedIdx.get(_i)));
        Vector<String> replicatesTwo = filter_and_map(rangev(selectedIdx.size()), (_i) -> !oneSide.contains(_i), (_i) -> replicateNames.get(selectedIdx.get(_i)));

        return UPair.createU(replicatesOne, replicatesTwo);
    }


    static Vector<Double> applyPseudo(Vector<Double> vals, double pseudo) {
        return (pseudo == 0.0) ? vals : map(vals, (_d) -> Double.isNaN(_d) ? _d : NumUtils.logN(Math.pow(2.0, _d) + pseudo, 2.0));
    }

    static Double applyPseudo(Double val, double pseudo) {
        return (Double.isNaN(val) || pseudo == 0.0) ? val : NumUtils.logN(Math.pow(2.0, val) + pseudo, 2.0);
    }

    public static void main(String[] args) {
        GeneralOptions generalOptions = new GeneralOptions();
        SimpleOptionParser cmd = new SimpleOptionParser("labels", "rawdata", "mappings", "summarizeddir", "od",  "nreps", "percentdiff", "minfc", "check", "calcdistribs",
                "pseudo", "plot", "nofilter");
        cmd.setFile("labels", "rawdata", "mappings");
        cmd.setDir("summarizeddir", "od");
        cmd.setInt("nreps");
        cmd.setDouble("percentdiff", "minfc", "pseudo");
        cmd.setDefault("nreps", "3");
        cmd.setDefault("percentdiff", "10");
        cmd.setDefault("pseudo", "0.0");
        cmd.setDefault("minfc", "0.6");
        cmd.setSwitches("check", "calcdistribs", "plot", "nofilter");

        if(!OptionParser.parseParams(args, true, false, true, true, cmd, generalOptions))
            return;

        generalOptions.apply();
        empires.SingleFeatureDiffExp.USE_NO_COMBINED_DISTRIBS = !cmd.isSet("calcdistribs");

        empires.input.MicroArrayInputRaw raw = new MicroArrayInputRaw(cmd.getFile("rawdata"), cmd.getFile("labels"), cmd.getFile("mappings"), false, new empires.AutoBackGroundContextProvider());
        int maxNumRep = NumUtils.max(raw.getCondition2NumReplicates().values());
        int nreps = cmd.getInt("nreps");
        String selcond = filterOne(raw.getCondition2NumReplicates().entrySet(), (_e) -> _e.getValue() == maxNumRep).getKey();

        if(maxNumRep < 2 * nreps) {
            System.err.printf("the maximum number of replicates: %d (condition: %s) less then 2 x %d (= requested number of replicates), quitting.\n", maxNumRep, selcond, nreps);
            return;
        }


        empires.NormalizedReplicateSet nrs = raw.getNormalized(selcond);

        UPair<Vector<String>> selected_reps = selectReplicates(nrs, nreps);
        //Vector<String> selected_reps = shuffleN(nrs.getInData().getReplicateNames(), 2 * nreps);
        Vector<String> features = nrs.getInData().getFeatureNames();
        Vector<String> cond1_replicates = selected_reps.getFirst();
        Vector<String> cond2_replicates = selected_reps.getSecond();

        Vector<String> joined = VectorUtils.join(cond1_replicates, cond2_replicates);
        HashMap<String, Integer> imap = buildIndexMap(joined);
        Vector<Integer> joinedidx = map(joined, (_j) -> imap.get(_j));


        Vector<Pair<String, Double>> f2sd = new Vector<>();
        for(String f : features) {
            Vector<Double> reps = nrs.getNormed(f).replicates;
            Vector<Double> vals = map(joinedidx, (_i) -> reps.get(_i));
            f2sd.add(Pair.create(f, NumUtils.getNumInfo(vals).SD));
        }

        Vector<Double> quants = NumUtils.getQuantiles(map(f2sd, (_p) -> _p.getSecond()), toVector(0.75, 0.85, 0.9, 0.95, 0.99), false);

        double FILTER_SD = cmd.isSet("nofilter") ? Double.POSITIVE_INFINITY : quants.get(quants.size() - 2);
        
        HashSet<String> REMOVE_FEATURES = toSet(filter_and_map(f2sd, (_p) -> _p.getSecond() > FILTER_SD, (_p) -> _p.getFirst()));

        System.out.printf("quants: %s\n", quants);


        System.out.printf("selected %s with %d replicates: %s \n", selcond, maxNumRep, selected_reps);

        double pseudo = cmd.getDouble("pseudo");
        if(cmd.isSet("plot")) {
            HashSet<String> c1set = toSet(cond1_replicates);


            PlotCreator  pc = CachedPlotCreator.getPlotCreator();

            File od = cmd.getFile("od");
            int PIDX = 0;
            for(UPair<String> up : getPairs(joined, true)) {
                boolean c1 = c1set.contains(up.getFirst());
                boolean c2 = c1set.contains(up.getSecond());

                boolean within_same = c1 == c2;

                int IDX1 = imap.get(up.getFirst());
                int IDX2 = imap.get(up.getSecond());
                Vector<Double> v1 = map(features, (_f) -> nrs.getNormed(_f).replicates.get(IDX1));
                Vector<Double> v2 = map(features, (_f) -> nrs.getNormed(_f).replicates.get(IDX2));
                Vector<Integer> valids = filter(rangev(features.size()), (_i) -> !Double.isNaN(v1.get(_i)) && !Double.isNaN(v2.get(_i)));


                pc.setTitle("%s vs %s (%s, %s)", up.getFirst(), up.getSecond(), c1, c2);
                pc.setLabels(up.getFirst(), up.getSecond(), null);
                pc.scatter("", valids, (_i) -> v1.get(_i), (_i) -> v2.get(_i));
                pc.setLog(true, true);
                pc.abline("", null, null, 0.0, 1.0);

                BufferedImage bim1 = pc.getImage();

                pc.setTitle("%s vs %s (%s, %s)", up.getFirst(), up.getSecond(), c1, c2);
                pc.setLabels(up.getFirst(), up.getSecond(), null);
                pc.dens_scatter2("", valids, (_i) -> v1.get(_i), (_i) -> v2.get(_i));
                pc.setLog(true, true);
                pc.abline("", null, null, 0.0, 1.0);

                BufferedImage bim2 = pc.getImage();
                //pc.writeImage(new ßFile(od, String.format("%s_%d.png", (within_same) ? "same" : "diff", PIDX++)));
                ImageUtils.saveImage(ImageUtils.concat(bim1, bim2), new File(od, String.format("%s_%d.png", (within_same) ? "same" : "diff", PIDX++)));

            }
        }

        HashMap<String, Vector<String>> feature2subfeature = new HashMap<>();
        HashMap<String, Integer> gene2numoligo = new HashMap<>();
        HashMap<Integer, Vector<String>> numoligo2genes = new HashMap<>();

        int removed_genes = 0;
        for(Map.Entry<String, Vector<String>> e : raw.getGene2Oligos().entrySet()) {
            Vector<String> filtered = filter(e.getValue(), (_e) -> !REMOVE_FEATURES.contains(_e));
            if(filtered.size() == 0) {
                removed_genes++;
                continue;
            }
            gene2numoligo.put(e.getKey(), filtered.size());
            feature2subfeature.put(e.getKey(), filtered);
            MapBuilder.updateV(numoligo2genes, filtered.size(), e.getKey());
        }

        System.out.printf("got %d/%d removed genes\n", removed_genes, raw.getGene2Oligos().size());


        if(cmd.isSet("check")) {
            empires.input.ReplicateSetInfo rsi = nrs.getInData();
            HashMap<String, Integer> repname2index = buildIndexMap(rsi.getReplicateNames());
            Vector<Integer> cond1indeces = map(cond1_replicates, (_r) -> repname2index.get(_r));
            Vector<Integer> cond2indeces = map(cond2_replicates, (_r) -> repname2index.get(_r));

            Vector<String> out_features = new Vector<>();
            Vector<Vector<Double>> out_vals_cond1 = map(cond1_replicates, (_r) -> new Vector<>());
            Vector<Vector<Double>> out_vals_cond2 = map(cond2_replicates, (_r) -> new Vector<>());

            for(String feature : nrs.getInData().getFeatureNamesToTest()) {
                for(String subf : rsi.getSubFeatures(feature)) {
                    if(REMOVE_FEATURES.contains(subf))
                        continue;

                    out_features.add(subf);
                    Vector<Double> subf_data = rsi.getReplicateData(subf);
                    applyIndex(cond1indeces.size(), (_i) -> out_vals_cond1.get(_i).add(applyPseudo(subf_data.get(cond1indeces.get(_i)), pseudo)));
                    applyIndex(cond2indeces.size(), (_i) -> out_vals_cond2.get(_i).add(applyPseudo(subf_data.get(cond2indeces.get(_i)), pseudo)));
                }

            }

            PlotCreator pc = CachedPlotCreator.getPlotCreator();
            for(int i=0; i<2; i++) {
               Vector<String> names = (i == 0) ? cond1_replicates : cond2_replicates;
               Vector<Vector<Double>> vals = (i == 0) ? out_vals_cond1 : out_vals_cond2;
               pc.cumhist(names.get(i), vals.get(i), vals.get(i).size(), false, false);
            }
            ImageUtils.showImage("distribs", pc.getImage());

            empires.input.ReplicateSetInfo out_rs1 = new empires.input.ReplicateSetInfo("cond1", cond1_replicates, out_features);
            empires.input.ReplicateSetInfo out_rs2 = new empires.input.ReplicateSetInfo("cond2", cond2_replicates, out_features);

            for(int i=0; i<out_vals_cond1.size(); i++) {
                out_rs1.setLog2Data(i, out_vals_cond1.get(i));
            }
            for(int i=0; i<out_vals_cond2.size(); i++) {
                out_rs2.setLog2Data(i, out_vals_cond2.get(i));
            }

            out_rs1.setCombinedFeatures(feature2subfeature);
            out_rs2.setCombinedFeatures(feature2subfeature);

            empires.NormalizedReplicateSet out_nrs1 = new empires.NormalizedReplicateSet(out_rs1);
            empires.NormalizedReplicateSet out_nrs2 = new empires.NormalizedReplicateSet(out_rs2);
            Vector<empires.DiffExpResult> res = new empires.EmpiRe().getDifferentialResults(out_nrs1, out_nrs2);

            empires.plotting.DiffExpTable diffExpTable = new empires.plotting.DiffExpTable(out_nrs1, out_nrs2, res);

            diffExpTable.showInteractiveTable();


        }
        HashMap<String, empires.NormalizedReplicateSet> summarized = readSummarized(cmd.getFile("summarizeddir"), cmd.getFile("labels"), cmd.getFile("mappings"), selcond);

        System.out.printf("summarized: %s\n", summarized.keySet());


        Vector<Integer> usable = filter(numoligo2genes.keySet(), (_k) -> numoligo2genes.get(_k).size() > 4);

        System.out.printf("usable oligos: %s\n", usable);

        double baseFC = cmd.getDouble("minfc");

        double MIN_FC_CALL_TRESHOLD = 0.8 * baseFC;
        NormalDistribution nd = new NormalDistribution(0.8, 0.4);

        int total_usable = 0;
        Vector<SwapPairs> data = new Vector<>();
        for(int k : usable) {

            Vector<Pair<empires.ReplicatedMeasurement, String>> vals = new Vector<>();
            total_usable += numoligo2genes.get(k).size();
            for(String g : numoligo2genes.get(k)) {
                Vector<String> oligos = feature2subfeature.get(g);
                int middle = oligos.size() >> 1;
                empires.ReplicatedMeasurement strongest = NumUtils.sort(map(oligos, (_n) -> nrs.getNormed(_n)), (_n) -> _n.mean, true).get(middle);
                vals.add(Pair.create(strongest, g));
            }
            boolean[] masked = new boolean[vals.size()];

            NumUtils.sort(vals, (_p) -> _p.getFirst().mean, false);
            int toselect = (int)(0.25 * masked.length);
            for(int i=0; i<masked.length; i++) {
                if(masked[i])
                    continue;


                double targetfc = baseFC + Math.abs(nd.sample());
                int sel = getIndexForFoldChange(i, targetfc, vals, masked);
                if(sel < 0 || sel == i)
                    continue;

                double gotfc = vals.get(sel).getFirst().mean - vals.get(toselect).getFirst().mean;

                if(gotfc < baseFC)
                    continue;


                SwapPairs sp = new SwapPairs(vals.get(i).getSecond(), vals.get(sel).getSecond(), k, gotfc);
                sp.oligos1 = getSortedSubFeatures(nrs, vals.get(i).getSecond());
                sp.oligos2 = getSortedSubFeatures(nrs, vals.get(sel).getSecond());
                Vector<Double> fcs = mapIndex(sp.num_oligos, (_i) -> nrs.getNormed(sp.oligos2.get(_i)).mean - nrs.getNormed(sp.oligos1.get(_i)).mean);
                NumUtils.BasicNumberInfo<Double> bni = NumUtils.getNumInfo(fcs);
                if(Math.abs(bni.median) < baseFC)
                    continue;


                System.out.printf("%s check: %s\n", sp, bni);
                sp.fc = bni.median;
                data.add(sp);
            }

        }

        double percentdiff = cmd.getDouble("percentdiff");
        Collections.shuffle(data);
        int totake = (int)(percentdiff / 200.0 * total_usable);
        Vector<SwapPairs> selected = VectorUtils.slice(data, 0, Math.min(data.size(), totake));
        System.out.printf("got %d possible swap pairs selected: %d\n", data.size(), selected.size());

        if(selected.size() < (int)(0.8 * totake)) {
            System.out.printf("PREC: NOT ENOUGH DIFF SIMULATABLE GENE TO TAKE: %d genes got: %d\n", totake, selected.size());
            return;
        }

        HashMap<String, SwapPairs> id2swap = new HashMap<>();
        apply(selected, (_sp) -> { id2swap.put(_sp.gene1, _sp); id2swap.put(_sp.gene2, _sp);});


        for(String method : summarized.keySet()) {
            Vector<SummarizedGene>  res = simulateSummarized(summarized.get(method), cond1_replicates, cond2_replicates, selected, id2swap.keySet());
            PerformanceResult pr = new PerformanceResult(method, res, (_r) -> _r.selected.pval, (_r) -> id2swap.containsKey(_r.gene), false, (_r) -> _r.selected.fdr <= 0.05 && Math.abs(_r.selected.fc) >= MIN_FC_CALL_TRESHOLD, null);

            System.out.printf("pr: %s\n", pr);

        }

        empires.input.ReplicateSetInfo rsi = nrs.getInData();
        HashMap<String, Integer> repname2index = buildIndexMap(rsi.getReplicateNames());
        Vector<Integer> cond1indeces = map(cond1_replicates, (_r) -> repname2index.get(_r));
        Vector<Integer> cond2indeces = map(cond2_replicates, (_r) -> repname2index.get(_r));

        Vector<String> out_features = new Vector<>();
        Vector<Vector<Double>> out_vals_cond1 = map(cond1_replicates, (_r) -> new Vector<>());
        Vector<Vector<Double>> out_vals_cond2 = map(cond2_replicates, (_r) -> new Vector<>());

        for(SwapPairs sp : selected) {
            for(int i=0; i<2; i++) {
                Vector<String> f1 = (i == 0) ? sp.oligos1 : sp.oligos2;
                Vector<String> f2 = (i == 0) ? sp.oligos2 : sp.oligos1;

                for(int j=0; j< sp.num_oligos; j++) {

                    out_features.add(f1.get(j));
                    Vector<Double> src_data = rsi.getReplicateData(f1.get(j));

                    Vector<Double> target_data = rsi.getReplicateData(f2.get(j));

                    applyIndex(cond1indeces.size(), (_i) -> out_vals_cond1.get(_i).add(applyPseudo(src_data.get(cond1indeces.get(_i)), pseudo)));
                    applyIndex(cond2indeces.size(), (_i) -> out_vals_cond2.get(_i).add(applyPseudo(target_data.get(cond2indeces.get(_i)), pseudo)));
                }
            }
        }

        for(String feature : nrs.getInData().getFeatureNamesToTest()) {
            if(id2swap.containsKey(feature) || ! feature2subfeature.containsKey(feature))
                continue;

            for(String subf : feature2subfeature.get(feature)) {
                out_features.add(subf);
                Vector<Double> subf_data = rsi.getReplicateData(subf);
                applyIndex(cond1indeces.size(), (_i) -> out_vals_cond1.get(_i).add(subf_data.get(cond1indeces.get(_i))));
                applyIndex(cond2indeces.size(), (_i) -> out_vals_cond2.get(_i).add(subf_data.get(cond2indeces.get(_i))));
            }

        }

        empires.input.ReplicateSetInfo out_rs1 = new empires.input.ReplicateSetInfo("cond1", cond1_replicates, out_features);
        empires.input.ReplicateSetInfo out_rs2 = new ReplicateSetInfo("cond2", cond2_replicates, out_features);

        for(int i=0; i<out_vals_cond1.size(); i++) {
            out_rs1.setLog2Data(i, out_vals_cond1.get(i));
        }
        for(int i=0; i<out_vals_cond2.size(); i++) {
            out_rs2.setLog2Data(i, out_vals_cond2.get(i));
        }

        out_rs1.setCombinedFeatures(feature2subfeature);
        out_rs2.setCombinedFeatures(feature2subfeature);

        empires.NormalizedReplicateSet out_nrs1 = new empires.NormalizedReplicateSet(out_rs1);
        empires.NormalizedReplicateSet out_nrs2 = new empires.NormalizedReplicateSet(out_rs2);
        Vector<empires.DiffExpResult> res = new EmpiRe().getDifferentialResults(out_nrs1, out_nrs2);
        PerformanceResult pr = new PerformanceResult("nlEmpiRe", res, (_r) -> _r.pval, (_r) -> id2swap.containsKey(_r.combinedFeatureName), false, (_r) -> _r.fdr <= 0.05 && Math.abs(_r.estimatedFC) >= MIN_FC_CALL_TRESHOLD, null, false, RocInfo.PerformanceEvaluationStrategy.OPTIMISTIC);

        System.out.printf("pr: %s\n", pr);

        empires.plotting.DiffExpTable diffExpTable = new DiffExpTable(out_nrs1, out_nrs2, res);

        //diffExpTable.showInteractiveTable(Pair.create("true", (_d) -> id2swap.containsKey(_d.combinedFeatureName)));

    }
}
