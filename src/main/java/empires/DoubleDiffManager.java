package nlEmpiRe;


import lmu.utils.*;
import lmu.utils.tuple.Tuple3;
import lmu.utils.tuple.Tuple4;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Vector;
import java.util.function.Function;

import static lmu.utils.IteratorUtils.rangev;
import static lmu.utils.ObjectGetter.*;
import org.apache.logging.log4j.Logger;

public class DoubleDiffManager {

    Logger log = LogConfig.getLogger();

    NormalizedReplicateSet fromSet;
    NormalizedReplicateSet toSet;
    DiffExpManager diffExpManager;

    HashMap<String, Vector<Double>> fromSetWithPseudo = new HashMap<>();
    HashMap<String, Vector<Double>> toSetWithPseudo = new HashMap<>();

    protected DoubleDiffManager(NormalizedReplicateSet fromSet, NormalizedReplicateSet toSet, DiffExpManager diffExpManager) {
        this.fromSet = fromSet;
        this.toSet = toSet;
        this.diffExpManager = diffExpManager;
    }

    public DiffExpManager getDiffExpManager() {
        return diffExpManager;
    }

    public DoubleDiffResult getDoubleDiffResult(String testname, Collection<String> subMeasurementSet1, Collection<String> subMeasurementSet2) {
        return getDoubleDiffResult(testname, subMeasurementSet1, subMeasurementSet2, DoubleDiffVariant.QUICK_AND_DIRTY);
    }

    public DoubleDiffResult getDoubleDiffResult(String testname, Collection<String> subMeasurementSet1, Collection<String> subMeasurementSet2, DoubleDiffVariant doubleDiffVariant) {
        return getDoubleDiffResult(testname, subMeasurementSet1, subMeasurementSet2, false, doubleDiffVariant);
    }

    public DoubleDiffResult getDoubleDiffResult(String testname, Collection<String> subMeasurementSet1, Collection<String> subMeasurementSet2, boolean verbose, DoubleDiffVariant doubleDiffVariant) {


        return getDoubleDiffResult(testname, subMeasurementSet1, subMeasurementSet2, null, verbose, doubleDiffVariant);
    }





    public DoubleDiffResult getDoubleDiffResult(String testname, Collection<String> subMeasurementSet1, Collection<String> subMeasurementSet2,
                                                HashSet<UPair<String>> pairs2skip, boolean verbose, DoubleDiffVariant doubleDiffVariant) {

        return getDoubleDiffResult(testname, subMeasurementSet1, subMeasurementSet2, pairs2skip, (_fn) -> diffExpManager.getDiffError(_fn), verbose, doubleDiffVariant);

    }


    /**
     *
     * @param subMeasurementSet1  set of (sub) measurements for object one (e.g. isoform one)
     * @param subMeasurementSet2  set of (sub) measurements for object one (e.g. isoform two)
     * @return DoubleDiffResult containing pvalue and a rough estimate of fold changes
     */

    public Vector<DoubleDiffResult> getMultiSubSetResults(String testPrefixName, Vector<Tuple3<String, Vector<String>, Vector<String>>> combinations, DoubleDiffVariant doubleDiffVariant) {
        HashSet<String> features = new HashSet<>();
        apply(combinations, (_c) -> {
            features.addAll(_c.get1());
            features.addAll(_c.get2());
        });

        HashMap<String, SingleFeatureDiffExp> featureDiffExpLookup = null; //buildMap(features, (_f) -> new SingleFeatureDiffExp(diffExpManager, _f));

        HashMap<Tuple, DiffExpResult>  cachedDiffExps = new HashMap<>();



        Vector<DoubleDiffResult> doubleDiffResults = new Vector<>();

        for(Tuple3<String, Vector<String>, Vector<String>> combi : combinations) {
            String testname = testPrefixName + "."  + combi.get0();
            DoubleDiffResult result = getDoubleDiffResult(testname, combi.get1(), combi.get2(), null, (_fn) -> diffExpManager.getDiffError(_fn), false, doubleDiffVariant);

            doubleDiffResults.add(result);
        }
        return doubleDiffResults;
    }

    public DoubleDiffResult getDoubleDiffResult(String testname, Vector<FeatureInfo> features1, Vector<FeatureInfo> features2, DoubleDiffVariant doubleDiffVariant) {
        switch (doubleDiffVariant) {
            case QUICK_AND_DIRTY:
                return getDoubleDiffResultFromFeatureInfos(testname, features1, features2);
            case MISSINGVAL:
                return getDoubleDiffResultFromFeatureInfosRespectMissingValues(testname, features1, features2);
            case ALLPAIRS:
                return getDoubleDiffResultAllPairs(testname, features1, features2);
        }
        throw new FRuntimeException("unknown double diff variant: %s", doubleDiffVariant);
    }
    public DoubleDiffResult getDoubleDiffResult(String testname, Collection<String> subMeasurementSet1, Collection<String> subMeasurementSet2,
                                                HashSet<UPair<String>> pairs2skip, Function<String, ErrorEstimationDistribution> errorGetter,
                                                boolean verbose,  DoubleDiffVariant doubleDiffVariant) {

        Vector<FeatureInfo> m1 = map(subMeasurementSet1, (_m) -> new FeatureInfo(_m, diffExpManager.replicateSetFrom, diffExpManager.replicateSetTo));
        Vector<FeatureInfo> m2 = map(subMeasurementSet2, (_m) -> new FeatureInfo(_m, diffExpManager.replicateSetFrom, diffExpManager.replicateSetTo));

        DoubleDiffResult ddr = getDoubleDiffResult(testname, m1, m2, doubleDiffVariant);

        ddr.featureInfos1 = m1;
        ddr.featureInfos2 = m2;

        return ddr;

    }

    int[] tests = new int[2];

    public DoubleDiffResult getDoubleDiffResultFromFeatureInfosByDistributionDifference(String testname, Collection<FeatureInfo> subMeasurementSet1, Collection<FeatureInfo> subMeasurementSet2) {

        log.info("start: %s\n", testname);

        synchronized (tests) {
            tests[0]++;
        }
        Vector<SingleFeatureDiffExp> fromSingles = map(subMeasurementSet1, (_f) -> new SingleFeatureDiffExp(_f.feature,  _f.cond1,  _f.err1,_f.cond2, _f.err2, diffExpManager));
        DiffExpResult fromDiffExp = new DiffExpResult(diffExpManager, testname+".set1", fromSingles, true);
        Vector<SingleFeatureDiffExp> toSingles = map(subMeasurementSet2, (_f) -> new SingleFeatureDiffExp(_f.feature,  _f.cond1,  _f.err1,_f.cond2, _f.err2, diffExpManager));

        DiffExpResult toDiffExp = new DiffExpResult(diffExpManager, testname+".set2", toSingles, true);


        DoubleDiffResult res = new DoubleDiffResult(testname,  map(subMeasurementSet1, (_s) -> _s.feature), diffExpManager.replicateSetFrom, map(subMeasurementSet2, (_s) -> _s.feature), diffExpManager.replicateSetTo, diffExpManager);
        res.estimatedFC = 0.0;
        res.pval = 1.0;

        if(fromDiffExp == null || fromDiffExp.combinedEmpiricalFoldChangeDistrib == null || toDiffExp == null || toDiffExp.combinedEmpiricalFoldChangeDistrib == null) {
            synchronized (tests) {
                tests[1]++;
            }

            System.err.printf("something wrong %d/%d: %s from %s: %.4f to: %s: %.4f",
                    tests[1], tests[0],
                    testname, map(subMeasurementSet1, (_s) -> _s.feature), fromDiffExp.pval,
                    map(subMeasurementSet2, (_s) -> _s.feature), toDiffExp.pval);

            return res;
        }
        log.info("%s:%d after diffexps", testname, tests[0]);
        long t1 = System.currentTimeMillis();
        ErrorEstimationDistribution diffDistribution =  ErrorEstimationDistribution.substract(fromDiffExp.combinedEmpiricalFoldChangeDistrib, toDiffExp.combinedEmpiricalFoldChangeDistrib, false);
        long t2 = System.currentTimeMillis();
        log.info("%s:%d diff distrib took: %.2f sec", testname, tests[0], (t2 - t1) / 1000.0);
        res.estimatedFC = fromDiffExp.estimatedFC - toDiffExp.estimatedFC;
        res.pval = diffDistribution.getCumulativeFrequencyToFoldChange(0.0);
        res.pval = Math.min(res.pval, 1.0 - res.pval);
        System.out.printf("got %s %d,%d %s vs %s :\n%s\t %g\n", testname, subMeasurementSet1.size(), subMeasurementSet2.size(),
                new SparseCumulativeDistribution(fromDiffExp.combinedEmpiricalFoldChangeDistrib).quantiles,
                new SparseCumulativeDistribution(toDiffExp.combinedEmpiricalFoldChangeDistrib).quantiles,
                new SparseCumulativeDistribution(diffDistribution).quantiles,
                res.pval);
        log.info("end test %d : %s\n", tests[0], testname);
        return res;
    }



    static class UsableReplicates {
        Vector<Integer> rep1indeces;
        Vector<Integer> rep2indeces;
        int nc1;
        int nc2;

        UsableReplicates(Vector<Double>  c1data, Vector<Double>  c2data) {
            rep1indeces = filter(rangev(c1data.size()), (_i) -> !Double.isNaN(c1data.get(_i)));
            rep2indeces = filter(rangev(c2data.size()), (_i) -> !Double.isNaN(c2data.get(_i)));
            nc1 = rep1indeces.size();
            nc2 = rep2indeces.size();
        }

    }


    public DoubleDiffResult getDoubleDiffResultAllPairs(String testname, Collection<FeatureInfo> subMeasurementSet1, Collection<FeatureInfo> subMeasurementSet2) {
        DoubleDiffResult result = new DoubleDiffResult(testname, null, fromSet, null, toSet, diffExpManager);

        result.pval = 1.0;
        result.meanFC = 0.0;

        final int nrep1 = fromSet.getNumReplicates();
        final int nrep2 = toSet.getNumReplicates();

        Vector<Integer> rep1range = IteratorUtils.rangev(nrep1);
        Vector<Integer> rep2range = IteratorUtils.rangev(nrep2);

        int numObservations = 0;
        double totalVariance = 0.0;
        double evidenceZTotal = 0.0;

        Vector<UPair<FeatureInfo>> featurePairs = getPairs(toVector(subMeasurementSet1), toVector(subMeasurementSet2));

        HashMap<FeatureInfo, Vector<UPair<FeatureInfo>> > feature2pairs = new HashMap<>();
        HashMap<UPair<FeatureInfo>, Vector<Vector<Integer>>> featurepair2indeces = new HashMap<>();


        HashMap<FeatureInfo, UsableReplicates> feature2usable = new HashMap<>();
        for(int i=0; i<2; i++) {
            for(FeatureInfo f : (i == 0) ? subMeasurementSet1 : subMeasurementSet2) {
                feature2usable.put(f, new UsableReplicates(f.getCond1LogVals(), f.getCond2LogVals()));
            }
        }

        for(UPair<FeatureInfo> featureInfoPair : featurePairs) {

            FeatureInfo f1 = featureInfoPair.getFirst();
            FeatureInfo f2 = featureInfoPair.getSecond();
            UsableReplicates ur1 = feature2usable.get(f1);
            UsableReplicates ur2 = feature2usable.get(f1);

            MapBuilder.updateV(feature2pairs, f1, featureInfoPair);
            MapBuilder.updateV(feature2pairs, f2, featureInfoPair);
            String f1name = f1.feature;
            String f2name = f2.feature;

            ErrorEstimationDistribution bg1 = diffExpManager.getDiffError(f1name);
            ErrorEstimationDistribution bg2 = diffExpManager.getDiffError(f2name);

            double f1c1var = diffExpManager.replicateSetFrom.getError(f1name).getVariance();
            double f1c2var = diffExpManager.replicateSetTo.getError(f1name).getVariance();
            double f2c1var = diffExpManager.replicateSetFrom.getError(f2name).getVariance();
            double f2c2var = diffExpManager.replicateSetTo.getError(f2name).getVariance();

            double corrNorm = (f1c1var + f1c2var + f2c1var + f2c2var);
            double corrSD = Math.sqrt(corrNorm);

            ErrorEstimationDistribution diffBG = getBGDist(bg1, bg2);

            int nObs_pep_pair = (ur1.nc1 * ur1.nc2 * ur2.nc1 * ur2.nc2);
            double pairVariance = nObs_pep_pair * (ur1.nc1 * ur1.nc2) * (ur2.nc2) * (f2c1var) +
                    nObs_pep_pair * (ur1.nc1 * ur1.nc2) * (ur2.nc1) * (f2c2var) +

                    nObs_pep_pair * (ur2.nc1 * ur2.nc2) * (ur1.nc2) * (f1c1var) +
                    nObs_pep_pair * (ur2.nc1 * ur2.nc2) * (ur1.nc1) * (f1c2var);

            // (a, b)
            totalVariance += pairVariance / corrNorm;

            Vector<Vector<Integer>> indeces = toVector(ur1.rep1indeces, ur1.rep2indeces, ur2.rep1indeces, ur2.rep2indeces);

            featurepair2indeces.put(featureInfoPair, indeces);

            for (int rep11idx : ur1.rep1indeces) {
                for (int rep12idx : ur1.rep2indeces) {
                    for (int rep21idx : ur2.rep1indeces) {
                        for (int rep22idx : ur2.rep2indeces) {
                            double diffdiffFC = (f1.getCond1LogVals().get(rep11idx) - f1.getCond2LogVals().get(rep12idx)) - (f2.getCond1LogVals().get(rep21idx) - f2.getCond2LogVals().get(rep22idx));
                            result.meanFC += diffdiffFC;
                            double z0 = diffBG.getConvertedZscore(diffdiffFC);
                            evidenceZTotal += z0;
                            numObservations++;
                        }

                    }
                }
            }

        }


        for(UPair<FeatureInfo> pair : featurePairs) {
            Vector<Vector<Integer>> indeces = featurepair2indeces.get(pair);
            if(indeces == null)
                continue;

            int numObs1 = NumUtils.product(map(indeces, (_v) -> _v.size()));

            ErrorEstimationDistribution pbg1 = diffExpManager.getDiffError(pair.getFirst().feature);
            ErrorEstimationDistribution pbg2 = diffExpManager.getDiffError(pair.getSecond().feature);

            double pair1SD = Math.sqrt(pbg1.getVariance() + pbg2.getVariance());

            for(int overlapIdx = 0; overlapIdx < 2; overlapIdx++) {
                FeatureInfo center = pair.get(overlapIdx == 0);

                double c1var = diffExpManager.replicateSetFrom.getError(center.feature).getVariance();
                double c2var = diffExpManager.replicateSetTo.getError(center.feature).getVariance();


                for(UPair<FeatureInfo> overlap : filter(feature2pairs.get(center), (_pair) -> !_pair.equals(pair))) {

                    ErrorEstimationDistribution obg1 = diffExpManager.getDiffError(overlap.getFirst().feature);
                    ErrorEstimationDistribution obg2 = diffExpManager.getDiffError(overlap.getSecond().feature);

                    double pair2SD = Math.sqrt(obg1.getVariance() + obg2.getVariance());

                    Vector<Vector<Integer>> indeces2 = featurepair2indeces.get(overlap);
                    if(indeces2 == null)
                        continue;

                    int numObs2 = NumUtils.product(map(indeces2, (_v) -> _v.size()));

                    double corrNorm = pair1SD * pair2SD;
                    double p1p2_repOverlap = 0;

                    for(int idxIdx = (overlapIdx == 0) ? 0 : 2, j=0; idxIdx < ((overlapIdx == 0) ? 2 : 4); idxIdx++, j++) {

                        double var = (j == 0) ? c1var : c2var;
                        int ncommon = SetInfo.intersect(indeces.get(idxIdx), indeces2.get(idxIdx)).size();
                        if(ncommon == 0)
                            continue;

                        int noverlaps = (numObs1 / indeces.get(idxIdx).size()) * (numObs2 / indeces2.get(idxIdx).size()) * ncommon;

                        p1p2_repOverlap += (noverlaps * var);
                    }
                    totalVariance += p1p2_repOverlap / corrNorm;

                }

            }

        }

        if(numObservations == 0)
            return  result;

        result.meanFC /= numObservations;

        NormalDistribution pepN = new NormalDistribution(0, Math.sqrt(totalVariance) ); //we will test only the right side (from zero = absolute) -> the variance is the half of it
        result.pval =  2.0 * (1.0 - pepN.cumulativeProbability(Math.abs(evidenceZTotal)));

        return result;
    }

    public DoubleDiffResult getDoubleDiffResultFromFeatureInfosRespectMissingValues(String testname, Collection<FeatureInfo> subMeasurementSet1, Collection<FeatureInfo> subMeasurementSet2) {
        DoubleDiffResult result = new DoubleDiffResult(testname, null, fromSet, null, toSet, diffExpManager);

        result.featureInfos1 = subMeasurementSet1;
        result.featureInfos2 = subMeasurementSet2;

        final int nrep1 = fromSet.getNumReplicates();
        final int nrep2 = toSet.getNumReplicates();

        Vector<Integer> rep1range = IteratorUtils.rangev(nrep1);
        Vector<Integer> rep2range = IteratorUtils.rangev(nrep2);

        double totalVariance = 0.0;
        double evidenceZTotal = 0.0;

        final HashMap<ErrorEstimationDistribution, Integer> EMPTY  = new HashMap<>();
        HashMap<FeatureInfo, HashMap<ErrorEstimationDistribution, Integer>> mes2bg2numUsed = new HashMap<>();

        HashMap<UPair<FeatureInfo>, UPair<Vector<Integer>>> featurepair2usable = new HashMap<>();
        HashMap<FeatureInfo, Vector<UPair<FeatureInfo>>> feature2pairs = new HashMap<>();
        Vector<UPair<FeatureInfo>> fpairs = new Vector<>();

        int numFcsCalced = 0;


        for(FeatureInfo feature1 : subMeasurementSet1) {
            ErrorEstimationDistribution bg1 = diffExpManager.getDiffError(feature1.err1, feature1.err2);

            //log.info("from: %s to: %s feature1: %s normed: %s,%s", fromSet, toSet, feature1, fromSet.getNormed(feature1), toSet.getNormed(feature1));
            Vector<Double> m1_rep1_logvals = feature1.getCond1LogVals();
            Vector<Double> m1_rep2_logvals = feature1.getCond2LogVals();

            for(FeatureInfo feature2 : subMeasurementSet2) {

                ErrorEstimationDistribution bg2 = diffExpManager.getDiffError(feature2.err1, feature2.err2);

                Vector<Double> m2_rep1_logvals = feature2.getCond1LogVals();
                Vector<Double> m2_rep2_logvals = feature2.getCond2LogVals();

                Vector<Integer> rep1_diff_usable = filter(rep1range, (_i) -> !Double.isNaN(m1_rep1_logvals.get(_i)) && !Double.isNaN(m2_rep1_logvals.get(_i)));

                if(rep1_diff_usable.size() == 0)
                    continue;

                Vector<Integer> rep2_diff_usable = filter(rep2range, (_i) -> !Double.isNaN(m1_rep2_logvals.get(_i)) && !Double.isNaN(m2_rep2_logvals.get(_i)));


                if(rep2_diff_usable.size() == 0)
                    continue;


                UPair<FeatureInfo> fpair = UPair.createU(feature1, feature2);
                fpairs.add(fpair);
                featurepair2usable.put(fpair, UPair.createU(rep1_diff_usable, rep2_diff_usable));
                MapBuilder.updateV(feature2pairs, feature1, fpair);
                MapBuilder.updateV(feature2pairs, feature2, fpair);

                double combinedVariance  = bg1.getVariance() + bg2.getVariance();
                double rep1comb = feature1.err1.getVariance() + feature2.err1.getVariance();
                double rep2comb = feature1.err2.getVariance() + feature2.err2.getVariance();

                int n1 = rep1_diff_usable.size();
                int n2 = rep2_diff_usable.size();
                double pairVariance = (n1 * n2) * (combinedVariance);
                pairVariance += (n1 * n2) * (n2 - 1) * (rep1comb) +
                        (n1 * n2) * (n1 - 1) * (rep2comb);

                totalVariance += pairVariance / combinedVariance;


                ErrorEstimationDistribution diffBG = getBGDist(bg2, bg1);

                double combinedSD = Math.sqrt(combinedVariance);

                for(int  rep1Idx : rep1_diff_usable) {
                    for(int rep2Idx : rep2_diff_usable) {

                        double diff1 = m1_rep1_logvals.get(rep1Idx)  - m1_rep2_logvals.get(rep2Idx);
                        double diff2 = m2_rep1_logvals.get(rep1Idx)  - m2_rep2_logvals.get(rep2Idx);

                        double diffFC = (diff2 - diff1);

                        result.meanFC += diffFC;
                        numFcsCalced++;
                        double z0 = diffBG.getConvertedZscore(diffFC);
                        //double z0 = diffFC / combinedSD;

                        evidenceZTotal += z0;
                    }
                }

            }
        }


        result.meanFC /= numFcsCalced;

        for(UPair<FeatureInfo> pair : fpairs) {
            UPair<Vector<Integer>> indeces = featurepair2usable.get(pair);
            if(indeces == null)
                continue;

            ErrorEstimationDistribution bg1 = diffExpManager.getDiffError(pair.getFirst().err1, pair.getFirst().err2);
            double pair1SD = Math.sqrt(bg1.getVariance());

            for(int overlapIdx = 0; overlapIdx < 2; overlapIdx++) {
                FeatureInfo center = pair.get(overlapIdx == 0);

                double c1var = center.err1.getVariance();
                double c2var = center.err2.getVariance();

                for(UPair<FeatureInfo> overlap : filter(feature2pairs.get(center), (_pair) -> !_pair.equals(pair))) {

                    ErrorEstimationDistribution bg2 = diffExpManager.getDiffError(overlap.getFirst().err1, overlap.getFirst().err2);
                    double pair2SD = Math.sqrt(bg2.getVariance());

                    UPair<Vector<Integer>> indeces2 = featurepair2usable.get(overlap);
                    if(indeces2 == null)
                        continue;

                    double corrNorm = pair1SD * pair2SD;
                    int n1 = SetInfo.intersect(indeces.getFirst(), indeces2.getFirst()).size();
                    int n2 = SetInfo.intersect(indeces.getSecond(), indeces2.getSecond()).size();

                    double p1p2_repOverlap = indeces.getSecond().size() * indeces2.getSecond().size() * (n1) * (c1var);
                    p1p2_repOverlap += indeces.getFirst().size() * indeces2.getFirst().size() * (n2) * (c2var);
                    totalVariance += p1p2_repOverlap / corrNorm;

                }

            }

        }

        if(numFcsCalced == 0) {
            result.pval = 1.0;
            return result;
        }

        NormalDistribution pepN = new NormalDistribution(0, Math.sqrt(totalVariance) ); //we will test only the right side (from zero = absolute) -> the variance is the half of it
        result.pval =  2.0 * (1.0 - pepN.cumulativeProbability(Math.abs(evidenceZTotal)));

        return result;

    }


    public DoubleDiffResult getDoubleDiffResultFromFeatureInfos(String testname, Collection<FeatureInfo> subMeasurementSet1, Collection<FeatureInfo> subMeasurementSet2) {

         DoubleDiffResult result = new DoubleDiffResult(testname, null, fromSet, null, toSet, diffExpManager);

         result.featureInfos1 = subMeasurementSet1;
         result.featureInfos2 = subMeasurementSet2;
        /**
         *  solution for auto-check for the degenerate case where one isoform set is differential, the other one is not significantly changing- but not strong
         *  enough to make the double differential change:
         *
         *      calculate the diff exp values ->
         *          1) if it significant for both, everything find
         *          2) if not -> select the highest sum count feature Fmax
         *              and create virtual features for the isoform side not having Fmax
         *              for all features f with count values  f + Fmax
         *              and re-do the test
         */
        if(false) {
            log.info("submes 1: %s\n", subMeasurementSet1);
            log.info("submes 2: %s\n", subMeasurementSet2);
        }



        final int nrep1 = fromSet.getNumReplicates();
        final int nrep2 = toSet.getNumReplicates();

        Vector<Integer> rep1range = IteratorUtils.rangev(nrep1);
        Vector<Integer> rep2range = IteratorUtils.rangev(nrep2);

        double evidenceZTotal = 0.0;

        CorrelationCache correlationCache = new CorrelationCache(subMeasurementSet1, subMeasurementSet2);


        //m1 -> (m2type1 : freq1


        double totalVariance = 0.0;

        final HashMap<ErrorEstimationDistribution, Integer> EMPTY  = new HashMap<>();
        HashMap<FeatureInfo, HashMap<ErrorEstimationDistribution, Integer>> mes2bg2numUsed = new HashMap<>();
        HashSet<UPair<FeatureInfo>> usable = new HashSet<>();


        int numFcsCalced = 0;
        for(FeatureInfo feature1 : subMeasurementSet1) {
            ErrorEstimationDistribution bg1 = diffExpManager.getDiffError(feature1.err1, feature1.err2);

            //log.info("from: %s to: %s feature1: %s normed: %s,%s", fromSet, toSet, feature1, fromSet.getNormed(feature1), toSet.getNormed(feature1));
            Vector<Double> m1_rep1_logvals = feature1.getCond1LogVals();
            Vector<Double> m1_rep2_logvals = feature1.getCond2LogVals();

            for(FeatureInfo feature2 : subMeasurementSet2) {

                ErrorEstimationDistribution bg2 = diffExpManager.getDiffError(feature2.err1, feature2.err2);

                Vector<Double> m2_rep1_logvals = feature2.getCond1LogVals();
                Vector<Double> m2_rep2_logvals = feature2.getCond2LogVals();

                Vector<Integer> rep1_diff_usable = filter(rep1range, (_i) -> !Double.isNaN(m1_rep1_logvals.get(_i)) && !Double.isNaN(m2_rep1_logvals.get(_i)));




                if(rep1_diff_usable.size() == 0)
                    continue;

                Vector<Integer> rep2_diff_usable = filter(rep2range, (_i) -> !Double.isNaN(m1_rep2_logvals.get(_i)) && !Double.isNaN(m2_rep2_logvals.get(_i)));


                if(rep2_diff_usable.size() == 0)
                    continue;



                usable.add(UPair.createU(feature1, feature2));
                totalVariance += 1.0; //the variance of the diagonal element in the covariance matrix (i.e. to (m1,m2) <-> (m1,m2)



                //and update the usages
                MapBuilder.updateMCount(mes2bg2numUsed, feature1, bg2);
                MapBuilder.updateMCount(mes2bg2numUsed, feature2, bg1);

                double combinedVariance  = bg1.getVariance() + bg2.getVariance();
                double rep1comb = feature1.err1.getVariance() + feature2.err1.getVariance();
                double rep2comb = feature1.err2.getVariance() + feature2.err2.getVariance();


                ErrorEstimationDistribution diffBG = getBGDist(bg2, bg1);

                double diffzsum = 0.0;
                for(int  rep1Idx : rep1_diff_usable) {
                    for(int rep2Idx : rep2_diff_usable) {

                        double diff1 = m1_rep1_logvals.get(rep1Idx)  - m1_rep2_logvals.get(rep2Idx);
                        double diff2 = m2_rep1_logvals.get(rep1Idx)  - m2_rep2_logvals.get(rep2Idx);

                        double diffFC = (diff2 - diff1);

                        result.meanFC += diffFC;
                        numFcsCalced++;
                        double z0 = diffBG.getConvertedZscore(diffFC);

                        diffzsum += z0;
                    }
                }

                double ncomb = rep1_diff_usable.size() * rep2_diff_usable.size();
                double totalvar = ncomb
                        + ncomb * (rep1_diff_usable.size() - 1) * (rep2comb / combinedVariance)
                        + ncomb * (rep2_diff_usable.size() - 1) * (rep1comb / combinedVariance);

                double totalSD = Math.sqrt(DistribUtils.getReplicateSummarizationVariance(rep1_diff_usable.size(), rep2_diff_usable.size()));
                double totalSDCheck = Math.sqrt(totalvar);



                if(false && Math.abs(totalSD - totalSDCheck) > 0.01) {
                    System.out.printf("comb: %.3f rep1: %.3f rep2: %.3f total: %.3f\n", combinedVariance, rep1comb, rep2comb, rep1comb + rep2comb);
                    System.out.printf("%d,%d SD: %.2f check: %.2f\n",  rep1_diff_usable.size(), rep2_diff_usable.size(), totalSD, totalSDCheck);
                }
                totalSD = totalSDCheck;
                evidenceZTotal += diffzsum / totalSD;

            }
        }


        result.meanFC /= numFcsCalced;

        for(UPair<FeatureInfo> pair : usable) {
            ErrorEstimationDistribution bg1 = diffExpManager.getDiffError(pair.getFirst().err1, pair.getFirst().err2);
            ErrorEstimationDistribution bg2 = diffExpManager.getDiffError(pair.getSecond().err1, pair.getSecond().err2);


            double selfCorr1 = correlationCache.corr(bg1, bg2, bg2, true);
            double selfCorr2 = correlationCache.corr(bg1, bg2, bg1, false);
            //add row
            totalVariance += NumUtils.sum(map(mes2bg2numUsed.get(pair.getFirst()).entrySet(), (_e) -> _e.getValue() * correlationCache.corr(bg1, bg2, _e.getKey(), true)), 0.0);
            //add column
            totalVariance += NumUtils.sum(map(mes2bg2numUsed.get(pair.getSecond()).entrySet(), (_e) -> _e.getValue() * correlationCache.corr(bg1, bg2, _e.getKey(), false)), 0.0);

            totalVariance -= (selfCorr1  + selfCorr2);


        }


        if(Math.abs(totalVariance) <= 0.01) {
            //no check was possible
            result.pval = 1.0;
            return result;
        }
        NormalDistribution pepN = new NormalDistribution(0, Math.sqrt(totalVariance) ); //we will test only the right side (from zero = absolute) -> the variance is the half of it
        result.pval =  2.0 * (1.0 - pepN.cumulativeProbability(Math.abs(evidenceZTotal)));

        return result;
    }


    static final String FAKE_STRING_PREFIX = "fake.";


    /*static class RawSignalFeature {
        String feature;

        Vector<Double> rep1vals;
        Vector<Double> rep2vals;

        Vector<Double> rep1_logvals;
        Vector<Double> rep2_logvals;
        double sum = 0;

        RawSignalFeature(String feature, NormalizedReplicateSet fromSet, NormalizedReplicateSet toSet) {
            this.feature = feature;

            rep1vals = map(fromSet.getNormed(feature).replicates, (_d) -> toRaw(_d));
            rep2vals = map(toSet.getNormed(feature).replicates, (_d) -> toRaw(_d));
            sum = NumUtils.sum(rep1vals, 0.0) + NumUtils.sum(rep2vals, 0.0);


        }

        RawSignalFeature(String name, RawSignalFeature base, RawSignalFeature toadd) {
            this.feature = name;
            rep1vals = mapIndex(base.rep1vals.size(), (_i) -> base.rep1vals.get(_i) + toadd.rep1vals.get(_i));
            rep2vals = mapIndex(base.rep2vals.size(), (_i) -> base.rep2vals.get(_i) + toadd.rep2vals.get(_i));
        }

        public Vector<Double> getLogVals(boolean onFrom) {
            if(rep1_logvals != null)
                return (onFrom) ? rep1_logvals : rep2_logvals;

            rep1_logvals = map(rep1vals, (_d) -> (_d == 0.0) ? Double.NaN : NumUtils.logN(_d, 2.0));
            rep2_logvals = map(rep2vals, (_d) -> (_d == 0.0) ? Double.NaN : NumUtils.logN(_d, 2.0));
            return (onFrom) ? rep1_logvals : rep2_logvals;
        }

        static double toRaw(double d) {
            return NONMEASURED_RAW_SIGNAL + ((Double.isNaN(d)) ? 0.0: Math.pow(2.0, d));
        }
    }
    */


    class CorrelationCache{
        HashMap<Tuple4<ErrorEstimationDistribution, ErrorEstimationDistribution, ErrorEstimationDistribution, Boolean>, Double> cache = new HashMap<>();

        CorrelationCache(Collection<FeatureInfo> obj1Measurements, Collection<FeatureInfo> obj2Measurements) {
            HashSet<ErrorEstimationDistribution> bgset1 = mapToSet(obj1Measurements, (_m) -> diffExpManager.getDiffError(_m.err1, _m.err2));
            HashSet<ErrorEstimationDistribution> bgset2 = mapToSet(obj2Measurements, (_m) -> diffExpManager.getDiffError(_m.err1, _m.err2));


            /***
             *
             *          f1: (c1f1E  -  c2f1E)   -> var
             *          f2: (c1f2E  -  c2f2E)   -> var
             *
             *
             *
             *          for the same peptide   -> s1-s2 -> diffError(Err(s1) - Err(s2))
             *          first sample-pair  z-value for splicing:  f1: (s1, s2)  f2(s1, s2)   == z-value for (f2-f1)
             *          second sample-pair z-value for splicing:  f1: (s1, s3)  f2(s1, s3)   == z-value for (f2-f1)
             *
             *
             */
            for(ErrorEstimationDistribution bg1 : bgset1) {
                double var1 = bg1.getVariance();


                for(ErrorEstimationDistribution bg2 : bgset2) {
                    double var2 = bg2.getVariance();

                    double divpart1 = Math.sqrt(var1 + var2);

                    for(ErrorEstimationDistribution bg3 : bgset2) {
                        double divpart2 = Math.sqrt(var1 + bg3.getVariance());

                        cache.put(Tuple4.create(bg1, bg2, bg3, true), var1 / (divpart1 * divpart2));
                    }

                    for(ErrorEstimationDistribution bg3 : bgset1) {
                        double divpart2 = Math.sqrt(var2 + bg3.getVariance());

                        cache.put(Tuple4.create(bg1, bg2, bg3, false), var2 / (divpart1 * divpart2));
                    }
                }
            }
            // sqrt(var1 + var2) * sqrt(var1 + var3)
            // (var1 + var2) * (var1 + var3)
            // var1^2 + var1*var2 + var1*var3

        }

        public double corr(ErrorEstimationDistribution bg1, ErrorEstimationDistribution bg2, ErrorEstimationDistribution bg3, boolean on_bg1) {
            return cache.get(Tuple4.create(bg1, bg2, bg3, on_bg1));
        }
    }


    HashMap<UPair<ErrorEstimationDistribution>, ErrorEstimationDistribution> diffBackgrounds = new HashMap<>();


    synchronized ErrorEstimationDistribution getBGDist(ErrorEstimationDistribution bg1, ErrorEstimationDistribution bg2) {

        UPair<ErrorEstimationDistribution> bgIndex = UPair.createU(bg1, bg2);
        ErrorEstimationDistribution bg = diffBackgrounds.get(bgIndex);
        if(bg == null) {

            diffBackgrounds.put(bgIndex, bg = ErrorEstimationDistribution.substract(bg1, bg2));
        }
        return bg;
    }

    public void clear_backgrounds()
    {
        diffBackgrounds = new HashMap<>();
    }

}
