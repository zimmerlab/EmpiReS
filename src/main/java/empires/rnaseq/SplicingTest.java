package empires.rnaseq;

import empires.EmpiRe;
import empires.plotting.NormalizedReplicateSetPlotting;
import lmu.utils.*;
import lmu.utils.fdr.BenjaminiHochberg;
import lmu.utils.plotting.PlotCreator;
import lmu.utils.tuple.Tuple3;
import empires.FeatureInfo;
import empires.input.RNASeqSplicingInfo;
import empires.input.ReplicateSetInfo;
import lmu.utils.plotting.CachedPlotCreator;

import java.awt.image.BufferedImage;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Vector;

import static lmu.utils.IteratorUtils.rangev;
import static lmu.utils.ObjectGetter.*;

public class SplicingTest {

    Logger log = LogConfig.getLogger();
    Vector<empires.NormalizedReplicateSet> replicateSetInfos;
    HashMap<String, empires.NormalizedReplicateSet> name2condition;

    HashMap<String, ReplicateSetInfo> name2summarizedCondition;
    HashMap<String, empires.NormalizedReplicateSet> name2summarizedNormedCondition;
    RNASeqSplicingInfo splicingInfo;
    int numThreads = 40;    //fixme

    int minReadCountTresholdInAllReplicatesInAtLeastOneCondition = 5;
    Double pseudo = 2.0;
    boolean correctDistributionsBackwards = true;
    boolean showPlots = false;
    public SplicingTest(RNASeqSplicingInfo splicingInfo, empires.BackgroundContextFuzzficationStrategyProvider backgroundContextFuzzficationStrategyProvider) {

        this.splicingInfo = splicingInfo;
        for(ReplicateSetInfo rsi : splicingInfo.getReplicateSetInfos()) {
            log.info("condition: %s num features: %s num replicates: %d", rsi.getReplicateSetName(), rsi.getNumFeatures(), rsi.getNumReplicates());
        }
        replicateSetInfos = map(splicingInfo.getReplicateSetInfos(), (_rsi) -> new empires.NormalizedReplicateSet(_rsi, backgroundContextFuzzficationStrategyProvider));
        name2condition = buildReverseMap(replicateSetInfos, (_rsi) -> _rsi.getInData().getReplicateSetName());
    }

    public void setPseudo(double pseudo) {
        this.pseudo =  (pseudo <= 0.0) ? null : pseudo;
    }

    public Double getPseudo() {
        return pseudo;
    }

    public void setShowPlots(boolean b) {
        showPlots = b;
    }
    public int getMinReadCountTresholdInAllReplicatesInAtLeastOneCondition() {
        return minReadCountTresholdInAllReplicatesInAtLeastOneCondition;
    }

    public void setMinReadCountTresholdInAllReplicatesInAtLeastOneCondition(int minReadCountTresholdInAllReplicatesInAtLeastOneCondition) {
        this.minReadCountTresholdInAllReplicatesInAtLeastOneCondition = minReadCountTresholdInAllReplicatesInAtLeastOneCondition;
    }

    public void setNumThreads(int n) {
        numThreads = n;
    }

    public Vector<String> getConditions() {
        return map(replicateSetInfos, (_r) -> _r.getInData().getReplicateSetName());
    }

    public empires.NormalizedReplicateSet getNormalized(String cond) {
        return name2condition.get(cond);
    }
/*
    public UPair<Vector<DoubleDiffResult>> getQuickTest(String cond1, String cond2, HashSet<String> trues) {
        NormalizedReplicateSet rs1 = getNormalized(cond1);
        NormalizedReplicateSet rs2 = getNormalized(cond2);
        final int nr1 = rs1.getNumReplicates();
        final int nr2 = rs1.getNumReplicates();
        final double MINLEVEL = NumUtils.logN(15, 2.0);
        DoubleDiffManager diffManager = new EmpiRe().getDoubleDiffManager(rs1, rs2);

        class TG {
            String gene;
            double pval;
            double fdr = 1.0;

            TG(String g) {
                this.gene = g;
            }
        }
        Vector<TG> tgenes = new Vector<>();
        for(String gene : splicingInfo.getGenes()) {
            Vector<String> features = splicingInfo.getFeaturesToGene(gene);
            Vector<String> filtered = filter(features, (_f) ->
            {
                return (nr1 == filteredSize(rs1.getInData().getReplicateData(_f), (_d) -> !Double.isNaN(_d) && _d >= MINLEVEL))
                    &&
                        (nr2 == filteredSize(rs2.getInData().getReplicateData(_f), (_d) -> !Double.isNaN(_d) && _d >= MINLEVEL))
                        ;
            });

            if(filtered.size() < 2)
                continue;

            Vector<SingleFeatureDiffExp> diffFeatures = NumUtils.sort(map(filtered, (_f) -> new SingleFeatureDiffExp(diffManager.getDiffExpManager(), _f)), (_sfd) -> _sfd.medianFC);

            HashMap<String, SingleFeatureDiffExp> lookup = buildReverseMap(diffFeatures, (_f) -> _f.feature);

            Double MINPVAL = null;
            double testfc = 0.0;
            for(Tuple3<String, Vector<String>, Vector<String>> t : splicingInfo.getSplicingTestFeatureCombis(gene, cond1, cond2)) {
                Vector<SingleFeatureDiffExp> fs1 = filter_and_map(t.get1(), (_f) -> lookup.containsKey(_f), (_f) -> lookup.get(_f));
                Vector<SingleFeatureDiffExp> fs2 = filter_and_map(t.get2(), (_f) -> lookup.containsKey(_f), (_f) -> lookup.get(_f));

                if(fs1.size() == 0  || fs2.size() == 0)
                    continue;

                DiffExpResult d1 = new DiffExpResult(diffManager.getDiffExpManager(), "cs1", fs1, true);
                DiffExpResult d2 = new DiffExpResult(diffManager.getDiffExpManager(), "cs2", fs2, true);

                ErrorEstimationDistribution diff = ErrorEstimationDistribution.getDifferenceJoinedErrorDistribution(d1.combinedEmpiricalFoldChangeDistrib, d2.combinedEmpiricalFoldChangeDistrib);

                double pval = diff.getCumulativeFrequencyToFoldChange(0.0);
                pval = Math.min(pval, 1.0 - pval);

                if(MINPVAL != null && MINPVAL <= pval)
                    continue;

                testfc = diff.getMostProbableFcWindowCenter();
                MINPVAL = pval;
            }

            if(MINPVAL == null)
                continue;

            //ErrorEstimationDistribution diff = ErrorEstimationDistribution.getDifferenceJoinedErrorDistribution(diffFeatures.get(0).scaledDistribution, diffFeatures.get(1).scaledDistribution);

            TG tg = new TG(gene);
            tg.pval = MINPVAL; //diff.getCumulativeFrequencyToFoldChange(0.0);
            //tg.pval = Math.min(tg.pval, 1.0 - tg.pval);
            tgenes.add(tg);

            if(MINPVAL < 0.001 && !trues.contains(gene)) {
                System.out.printf("GOES WRONG: %s effect size: %.2f\n",gene, testfc);
                for(String f : NumUtils.sort(filtered, (_f) -> lookup.get(_f).medianFC)) {
                    SingleFeatureDiffExp sdef = lookup.get(f);
                    System.out.printf("%s: %.2f %s\t%s\n", f, sdef.medianFC, rs1.getInData().getReplicateData(f), rs2.getInData().getReplicateData(f));

                }
            }
            if(tgenes.size() % 200 == 0) {
                System.out.printf("after %d tgenes\n", tgenes.size());
                BenjaminiHochberg.adjust_pvalues(tgenes, (_t) -> _t.pval, (_p) -> _p.getFirst().fdr = _p.getSecond());

                PerformanceResult pr = new PerformanceResult("quictest", tgenes, (_dr) -> _dr.pval, (_t) -> trues.contains(_t.gene),
                        false, (_dr) -> _dr.fdr <= 0.05, null, false);


                System.out.printf("%s\n", pr);

            }
            //System.out.printf("gene: %s got %d features fcs: %s\n", gene, features.size(), NumUtils.getNumInfo(diffFeatures, (_d) -> _d.medianFC).getInfoWithQ());
        }
        System.out.printf("got %d tgenes\n", tgenes.size());
        BenjaminiHochberg.adjust_pvalues(tgenes, (_t) -> _t.pval, (_p) -> _p.getFirst().fdr = _p.getSecond());

        PerformanceResult pr = new PerformanceResult("quictest", tgenes, (_dr) -> _dr.pval, (_t) -> trues.contains(_t.gene),
                false, (_dr) -> _dr.fdr <= 0.05, null, false);


        System.out.printf("%s\n", pr);

        return null;
    }

 */

    public Vector<empires.DiffExpResult> getDiffEQClasses(String cond1, String cond2) {
        Logger log = LogConfig.getLogger();

        empires.NormalizedReplicateSet rs1 = getNormalized(cond1);
        empires.NormalizedReplicateSet rs2 = getNormalized(cond2);

        if(correctDistributionsBackwards) {

            rs1.correctBackWards();
            rs2.correctBackWards();

        }

        return new empires.EmpiRe().getDifferentialResults(rs1, rs2);

    }

    public Vector<Vector<UPair<Double>>> getCountLevel2CoveragePerReplicate(String cond) {

        HashMap<Integer, HashMap<Integer, Integer>> replicate2featurecounts = new HashMap<>();

        ReplicateSetInfo conditionReplicates = filterOne(splicingInfo.getReplicateSetInfos(), (_rsi) -> _rsi.getReplicateSetName().equals(cond));
        if(conditionReplicates == null)
            throw new FRuntimeException("unknown condition: %s", cond);

        for(String g : splicingInfo.getGenes()) {

            for(String feature: splicingInfo.getFeaturesToGene(g)) {
                Vector<Integer> fdata = map(conditionReplicates.getReplicateData(feature), (_v) -> (int)((Double.isNaN(_v) ? 0.0 : Math.pow(2.0, _v))));
                int maxval = NumUtils.max(fdata);
                MapBuilder.updateMCount(replicate2featurecounts, 0, maxval);
            }
        }
        Vector<Vector<UPair<Double>>> countlevel2percentData = new Vector<>();

        int idx = 0;
        for(HashMap<Integer, Integer> repData  :replicate2featurecounts.values()) {
            idx++;
            Vector<Integer> keys = toVector(repData.keySet());
            Collections.sort(keys);
            Collections.reverse(keys);
            double sum = 0.0 + NumUtils.sum(map(keys, (_k) -> _k * repData.get(_k)));
            double covered = 0.0;
            Vector<UPair<Double>> replicateCovered = new Vector<>();
            Vector<Integer> coveragePercents = new Vector<>();
            Vector<Integer> counts = new Vector<>();

            for(int k : keys) {
                covered += k * repData.get(k);
                coveragePercents.add((int)(covered * 100 / sum));
                counts.add(k);
                //System.out.printf("%s: %d count: %d covered: %.0f sum: %.0f : %.2f%%\n", cond, idx, k, covered, sum, covered * 100 / sum);
                replicateCovered.add(UPair.createU(k + 0.0, covered * 100.0 / sum));
            }
            int hitIdx = Collections.binarySearch(coveragePercents, 99);
            if(hitIdx < 0) {
                hitIdx = - hitIdx - 1;
            }
            int minC = counts.get(hitIdx);
            double sumTest = 0.0 + NumUtils.sum(map(filter(keys, (_k) -> _k >= minC), (_k) -> _k * repData.get(_k)));
            System.out.printf("99%% coverage at count: %d check: %.2f%%\n", minC, 100.0 * sumTest / sum);
            countlevel2percentData.add(replicateCovered);
        }
        return countlevel2percentData;
    }

    public HashMap<Integer, Vector<Double>> getEqClassLevel2Count(String cond) {
        HashMap<Integer, Vector<Double>> rv = new HashMap<>();
        empires.NormalizedReplicateSet nrs = getNormalized(cond);
        for(String g : splicingInfo.getGenes()) {
            Vector<Double> vals = new Vector<>();
            for(String feature: splicingInfo.getFeaturesToGene(g)) {
                Vector<Double> fdata = map(nrs.getInData().getReplicateData(feature), (_v) -> (Double.isNaN(_v) ? 0.0 : Math.pow(2.0, _v)));
                vals.add(NumUtils.max(fdata));
            }
            Collections.sort(vals);
            Collections.reverse(vals);
            for(int i=0; i<vals.size(); i++) {
                MapBuilder.updateV(rv, i, vals.get(i));
            }
        }
        return rv;
    }

    void calcSummarizedGeneCounts() {

        if(name2summarizedCondition != null)
            return;

        name2summarizedCondition = new HashMap<>();
        name2summarizedNormedCondition = new HashMap<>();
        Vector<String> genes = toVector(splicingInfo.getGenes());
        for(ReplicateSetInfo rsi : splicingInfo.getReplicateSetInfos()) {
            ReplicateSetInfo sumRsi = new ReplicateSetInfo(rsi.getReplicateSetName(), rsi.getReplicateNames(), genes);
            Vector<Vector<Double>> sumLogData = mapIndex(rsi.getNumReplicates(), (_i) -> new Vector<>());
            for(String gene : genes) {

                Vector<String> features = splicingInfo.getFeaturesToGene(gene);
                Vector<Double> summarized = mapIndex(rsi.getNumReplicates(), (_i) -> 0.0);
                for(String feature : features) {
                    Vector<Double> fdata = map(rsi.getReplicateData(feature), (_v) -> Math.pow(2.0, _v));
                    for(int i=0; i<summarized.size(); i++) {
                        summarized.set(i, fdata.get(i) + summarized.get(i));
                    }


                }
                for(int i=0; i<rsi.getNumReplicates(); i++) {
                    double val = summarized.get(i);
                    sumLogData.get(i).add((val == 0.0) ? Double.NaN : NumUtils.logN(val, 2.0));
                }


            }
            for(int i=0; i<rsi.getNumReplicates(); i++) {
                sumRsi.setLog2Data(i, sumLogData.get(i));
            }


            name2summarizedCondition.put(rsi.getReplicateSetName(), sumRsi);
            empires.NormalizedReplicateSet normalizedReplicateSet = new empires.NormalizedReplicateSet(sumRsi);
            if(correctDistributionsBackwards) {
                normalizedReplicateSet.correctBackWards();
            }
            name2summarizedNormedCondition.put(rsi.getReplicateSetName(), normalizedReplicateSet);
        }


    }

    public HashMap<String, Double> getSummarizedMeanCounts(String cond) {
        calcSummarizedGeneCounts();
        empires.NormalizedReplicateSet nrs = name2summarizedNormedCondition.get(cond);
        HashMap<String, Double> m = new HashMap<>();
        for(String fn : nrs.getInData().getFeatureNames()) {
            m.put(fn, NumUtils.mean(map(nrs.getNormed(fn).nonNanValues, (_n) -> Math.pow(2.0, _n))));
        }
        return m;

    }
    public Vector<empires.DiffExpResult> getDiffGenes(String cond1, String cond2, boolean summarizeBefore) {
        Logger log = LogConfig.getLogger();

        if(summarizeBefore) {
            calcSummarizedGeneCounts();
            return new empires.EmpiRe().getDifferentialResults(name2summarizedNormedCondition.get(cond1), name2summarizedNormedCondition.get(cond2));
        }
        empires.NormalizedReplicateSet rs1 = getNormalized(cond1);
        empires.NormalizedReplicateSet rs2 = getNormalized(cond2);

        if(correctDistributionsBackwards) {

            rs1.correctBackWards();
            rs2.correctBackWards();
        }


        HashMap<String, Vector<String>> gene2features = new HashMap<>();
        for(String g : splicingInfo.getGenesForSplicingTest()) {
            gene2features.put(g, splicingInfo.getFeaturesToGene(g));
        }
        return new empires.EmpiRe().getDifferentialResults(rs1, rs2, gene2features);

    }


    public empires.DoubleDiffVariant doubleDiffVariant = empires.DoubleDiffVariant.QUICK_AND_DIRTY;

    public UPair<Vector<empires.DoubleDiffResult>> getDifferentialAlternativeSplicing(String cond1, String cond2) {
        Logger log = LogConfig.getLogger();

        empires.NormalizedReplicateSet rs1 = getNormalized(cond1);
        empires.NormalizedReplicateSet rs2 = getNormalized(cond2);

        PlotCreator pc = (showPlots) ? CachedPlotCreator.getPlotCreator() : null;
        BufferedImage before = null;
        if(showPlots) {
            before = ImageUtils.concat(new empires.plotting.NormalizedReplicateSetPlotting(rs1).plotBackGroundDistribs(pc), new empires.plotting.NormalizedReplicateSetPlotting(rs2).plotBackGroundDistribs(pc));
        }

        if(correctDistributionsBackwards) {

            rs1.correctBackWards();
            rs2.correctBackWards();

        }


        if(showPlots) {
            BufferedImage after = ImageUtils.concat(new empires.plotting.NormalizedReplicateSetPlotting(rs1).plotBackGroundDistribs(pc), new empires.plotting.NormalizedReplicateSetPlotting(rs2).plotBackGroundDistribs(pc));
            ImageUtils.showImage("test", ImageUtils.vconcat(before, after), false);
        }


        empires.DoubleDiffManager diffManager = new empires.EmpiRe().getDoubleDiffManager(rs1, rs2);


        //double pseudo = splicingInfo.getMinCountIntMaxCountCondition() * 0.5;
        if(pseudo != null) {
            log.info("apply pseudo : " + pseudo);
            rs1.applyPseudo(pseudo);
            rs2.applyPseudo(pseudo);

        }

        final int nr1 = rs1.getNumReplicates();
        final int nr2 = rs1.getNumReplicates();
        final double MINLEVEL = NumUtils.logN(1.0, 2.0);

        Vector<empires.DoubleDiffResult> diffResults = new Vector<>();
        Vector<empires.DoubleDiffResult> diffResults_unpaired = new Vector<>();

        log.info("got %d/%d genes with <= features", filteredSize(splicingInfo.getGenesForSplicingTest(), (_g) -> splicingInfo.getNumFeatures(_g) <= 1), splicingInfo.getGenesForSplicingTest().size());

        Vector<Vector<String>> packed_genes = VectorUtils.pack(toVector(splicingInfo.getGenesForSplicingTest()), 10);
        int[] nupdates = new int[1];
        ThreadUtils.runLambdas(packed_genes,
                (_gene_vec) ->
                {
                    Vector<empires.DoubleDiffResult> subResults = new Vector<>();
                    Vector<empires.DoubleDiffResult> subResults_up = new Vector<>();
                    for(String gene : _gene_vec){

                        Vector<String> features = splicingInfo.getFeaturesToGene(gene);

                        HashSet<String> filtered = toSet(features);
                        /*filterToSet(features, (_f) ->
                        {
                            return (nr1 == filteredSize(rs1.getInData().getReplicateData(_f), (_d) -> !Double.isNaN(_d) && _d >= MINLEVEL))
                                    &&
                                    (nr2 == filteredSize(rs2.getInData().getReplicateData(_f), (_d) -> !Double.isNaN(_d) && _d >= MINLEVEL))
                                    ;
                        });
                    */
                        if(filtered.size() < 2)
                            continue;

                        Vector<empires.DoubleDiffResult> geneDiffs = new Vector<>();

                        Vector<Tuple3<String, Vector<String>, Vector<String>>> combis = splicingInfo.getSplicingTestFeatureCombis(gene, cond1, cond2);

                        for(Tuple3<String, Vector<String>, Vector<String>> combi  : combis) {
                            Vector<String> fs1 = filter(combi.get1(), (_f) -> filtered.contains(_f));
                            Vector<String> fs2 = filter(combi.get2(), (_f) -> filtered.contains(_f));

                            Vector<FeatureInfo> merged1 = FeatureInfo.merge(fs1, rs1, rs2, minReadCountTresholdInAllReplicatesInAtLeastOneCondition, pseudo);
                            Vector<FeatureInfo> merged2 = FeatureInfo.merge(fs2, rs1, rs2, minReadCountTresholdInAllReplicatesInAtLeastOneCondition, pseudo);

                            if(fs1.size() == 0 || fs2.size() == 0)
                                continue;

                            if(merged1.size() == 0 || merged2.size() == 0)
                                continue;

                            switch (doubleDiffVariant) {
                                case QUICK_AND_DIRTY:
                                    geneDiffs.add(diffManager.getDoubleDiffResultFromFeatureInfos(gene + "." + combi.get0(), merged1, merged2));
                                    break;
                                case MISSINGVAL:
                                    geneDiffs.add(diffManager.getDoubleDiffResultFromFeatureInfosRespectMissingValues(gene + "." + combi.get0(), merged1, merged2));
                                    break;
                                case ALLPAIRS:
                                    geneDiffs.add(diffManager.getDoubleDiffResultAllPairs(gene + "." + combi.get0(), merged1, merged2));
                                    break;

                            }

                            //geneDiffs.add(diffManager.getDoubleDiffResultFromFeatureInfosByDistributionDifference(gene + "." + combi.get0(), merged1, merged2));
                        }
                        //Vector<DoubleDiffResult> geneDiffs = diffManager.getMultiSubSetResults(gene, splicingInfo.getSplicingTestFeatureCombis(gene, cond1, cond2));
                        if(geneDiffs.size() == 0)
                            continue;

                        /**
                         * here added to check for only protein coding @Constantin
                         */
                        //here to get all pairs
//                        subResults.addAll(geneDiffs);

                        empires.DoubleDiffResult best = NumUtils.minObj(geneDiffs, (_d) -> _d.pval).second;
                        if (
//                                best.testName.split("\\.")[0].equals("ENSG00000127511") || best.testName.split("\\.")[0].equals("ENSG00000119689") ||
//                            best.testName.split("\\.")[0].equals("ENSG00000129007") ||
//                            best.testName.split("\\.")[0].equals("ENSG00000123836") ||
//                            best.testName.split("\\.")[0].equals("ENSG00000148343") ||
                            best.testName.split("\\.")[0].equals("ENSG00000138028") //||
//                            best.testName.split("\\.")[0].equals("ENSG00000135447")
                        ) {
                            System.out.println("\n--------------------------\n");
                            geneDiffs.forEach(_g -> System.out.println(_g.testName + "\t" + String.format("%.2f", _g.meanFC) + "\t" + String.format("%.2f", _g.pval)));
                            System.out.println("--best--");
                            System.out.println(best + "\t" + String.format("%.2f", best.meanFC) + "\t" + String.format("%.2f", best.pval));
                            System.out.println();
                        }

                        //here to only get best pair
                        subResults.add(NumUtils.minObj(geneDiffs, (_d) -> _d.pval).second);
                        //subResults_up.add(NumUtils.minObj(geneDiffs, (_d) -> _d.unpairedPval).second);
                    }

                    synchronized (diffResults) {
                        diffResults.addAll(subResults);
                        //diffResults_unpaired.addAll(subResults_up);
                        if(++nupdates[0] % 10 == 0) {
                            log.info("processed %d/%d genes", diffResults.size(), splicingInfo.getGenesForSplicingTest().size());
                        }

                    }
                },
                numThreads, 10, true, -1, null);



        BenjaminiHochberg.adjust_pvalues(diffResults, (_d) -> _d.pval, (_p) -> _p.getFirst().fdr = _p.getSecond());
        //BenjaminiHochberg.adjust_pvalues(diffResults_unpaired, (_d) -> _d.unpairedPval, (_p) -> _p.getFirst().unpairedFDR = _p.getSecond());

        return UPair.createU(diffResults, diffResults_unpaired);
    }







    public UPair<Vector<empires.DoubleDiffResult>> getDifferentialTranscriptExpression(String cond1, String cond2) {
        Logger log = LogConfig.getLogger();

        empires.NormalizedReplicateSet rs1 = getNormalized(cond1);
        empires.NormalizedReplicateSet rs2 = getNormalized(cond2);

        PlotCreator pc = (showPlots) ? CachedPlotCreator.getPlotCreator() : null;
        BufferedImage before = null;
        if(showPlots) {
            before = ImageUtils.concat(new empires.plotting.NormalizedReplicateSetPlotting(rs1).plotBackGroundDistribs(pc), new empires.plotting.NormalizedReplicateSetPlotting(rs2).plotBackGroundDistribs(pc));
        }

        if(correctDistributionsBackwards) {
            rs1.correctBackWards();
            rs2.correctBackWards();
        }


        if(showPlots) {
            BufferedImage after = ImageUtils.concat(new empires.plotting.NormalizedReplicateSetPlotting(rs1).plotBackGroundDistribs(pc), new NormalizedReplicateSetPlotting(rs2).plotBackGroundDistribs(pc));
            ImageUtils.showImage("test", ImageUtils.vconcat(before, after), false);
        }


        empires.DoubleDiffManager diffManager = new EmpiRe().getDoubleDiffManager(rs1, rs2);


        //double pseudo = splicingInfo.getMinCountIntMaxCountCondition() * 0.5;
        if(pseudo != null) {
            log.info("apply pseudo : " + pseudo);
            rs1.applyPseudo(pseudo);
            rs2.applyPseudo(pseudo);
        }

        final int nr1 = rs1.getNumReplicates();
        final int nr2 = rs1.getNumReplicates();
        final double MINLEVEL = NumUtils.logN(1.0, 2.0);

        Vector<empires.DoubleDiffResult> diffResults = new Vector<>();
        Vector<empires.DoubleDiffResult> diffResults_unpaired = new Vector<>();

        log.info("got %d/%d genes with <= features", filteredSize(splicingInfo.getGenesForSplicingTest(), (_g) -> splicingInfo.getNumFeatures(_g) <= 1), splicingInfo.getGenesForSplicingTest().size());

        Vector<Vector<String>> packed_genes = VectorUtils.pack(toVector(splicingInfo.getGenesForSplicingTest()), 10);
        int[] nupdates = new int[1];
        ThreadUtils.runLambdas(packed_genes,
                (_gene_vec) ->
                {
                    Vector<empires.DoubleDiffResult> subResults = new Vector<>();
                    Vector<empires.DoubleDiffResult> subResults_up = new Vector<>();
                    for(String gene : _gene_vec){

                        Vector<String> features = splicingInfo.getFeaturesToGene(gene);

                        HashSet<String> filtered = toSet(features);
                        /*filterToSet(features, (_f) ->
                        {
                            return (nr1 == filteredSize(rs1.getInData().getReplicateData(_f), (_d) -> !Double.isNaN(_d) && _d >= MINLEVEL))
                                    &&
                                    (nr2 == filteredSize(rs2.getInData().getReplicateData(_f), (_d) -> !Double.isNaN(_d) && _d >= MINLEVEL))
                                    ;
                        });
                    */
                        if(filtered.size() < 2)
                            continue;

                        Vector<empires.DoubleDiffResult> geneDiffs = new Vector<>();

//                        for (Pair<Tuple, String> eq : splicingInfo.getGene2EQClassFeature().get(gene)) {
//                            if (!filtered.contains(eq.getSecond()))
//                                continue;

                            Vector<String> f1 = new Vector<>();
                            Vector<String> f2 = new Vector<>();
                            f1.addAll(filtered);
                            f2.addAll(filtered);

                            Vector<String> fs1 = filter(f1, (_f) -> filtered.contains(_f));
                            Vector<String> fs2 = filter(f2, (_f) -> filtered.contains(_f));

                            Vector<FeatureInfo> merged1 = FeatureInfo.merge(fs1, rs1, rs2, minReadCountTresholdInAllReplicatesInAtLeastOneCondition, pseudo);
                            Vector<FeatureInfo> merged2 = FeatureInfo.merge(fs2, rs1, rs2, minReadCountTresholdInAllReplicatesInAtLeastOneCondition, pseudo);

                            if(fs1.size() == 0 || fs2.size() == 0)
                                continue;

                            if(merged1.size() == 0 || merged2.size() == 0)
                                continue;

                            switch (doubleDiffVariant) {
                                case QUICK_AND_DIRTY:
                                    geneDiffs.add(diffManager.getDoubleDiffResultFromFeatureInfos(gene + "." + "all", merged1, merged2));
                                    break;
                                case MISSINGVAL:
                                    geneDiffs.add(diffManager.getDoubleDiffResultFromFeatureInfosRespectMissingValues(gene + "." + "all", merged1, merged2));
                                    break;
                                case ALLPAIRS:
                                    geneDiffs.add(diffManager.getDoubleDiffResultAllPairs(gene + "." + "all", merged1, merged2));
                                    break;
                            }
//                        }

                        if(geneDiffs.size() == 0)
                            continue;


                        //debug atm
                        empires.DoubleDiffResult best = NumUtils.minObj(geneDiffs, (_d) -> _d.pval).second;
                        if (
//                                best.testName.split("\\.")[0].equals("ENSG00000127511") || best.testName.split("\\.")[0].equals("ENSG00000119689") ||
//                            best.testName.split("\\.")[0].equals("ENSG00000129007") ||
//                            best.testName.split("\\.")[0].equals("ENSG00000123836") ||
//                            best.testName.split("\\.")[0].equals("ENSG00000148343") ||
                                best.testName.split("\\.")[0].equals("ENSG00000138028") //||
//                            best.testName.split("\\.")[0].equals("ENSG00000135447")
                        ) {
                            System.out.println("\n--------------------------\n");
                            geneDiffs.forEach(_g -> System.out.println(_g.testName + "\t" + String.format("%.2f", _g.meanFC) + "\t" + String.format("%.2f", _g.pval)));
                            System.out.println("--best--");
                            System.out.println(best.testName + "\t" + String.format("%.2f", best.meanFC) + "\t" + String.format("%.2f", best.pval));
                            System.out.println();
                        }

                        subResults.add(NumUtils.minObj(geneDiffs, (_d) -> _d.pval).second);
                    }

                    synchronized (diffResults) {
                        diffResults.addAll(subResults);
                        //diffResults_unpaired.addAll(subResults_up);
                        if(++nupdates[0] % 10 == 0) {
                            log.info("processed %d/%d genes", diffResults.size(), splicingInfo.getGenesForSplicingTest().size());
                        }

                    }
                },
                numThreads, 10, true, -1, null);

        BenjaminiHochberg.adjust_pvalues(diffResults, (_d) -> _d.pval, (_p) -> _p.getFirst().fdr = _p.getSecond());
        //BenjaminiHochberg.adjust_pvalues(diffResults_unpaired, (_d) -> _d.unpairedPval, (_p) -> _p.getFirst().unpairedFDR = _p.getSecond());

        return UPair.createU(diffResults, diffResults_unpaired);
    }

    public String getSettings() {
        return String.format("minSD to join: %g minCountInAllReps in at least one conditions: %d pseudo: %.2f, reduce to max tr: %d",
                empires.AutoBackGroundContextProvider.MERGE_BACKGROUNDS_WITH_MAX_SD_DIFF, minReadCountTresholdInAllReplicatesInAtLeastOneCondition, pseudo,
            splicingInfo.getMaximalTrNumForReduction());

    }
}
