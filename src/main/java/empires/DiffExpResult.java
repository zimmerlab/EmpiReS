package nlEmpiRe;

import lmu.utils.*;
import lmu.utils.fdr.PerformanceResult;
import lmu.utils.fdr.RocInfo;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Vector;
import java.util.function.Function;

import static lmu.utils.ObjectGetter.*;
import org.apache.logging.log4j.Logger;

public class DiffExpResult {

    public static double PERCENT_ROBUST_FEATURE_USED = 0.9;
    Logger log = LogConfig.getLogger();
    DiffExpManager diffExpManager;
    static final NormalDistribution N01 = new NormalDistribution(0, 1);
    static final NormalDistribution NSD2 = new NormalDistribution(0, 2.0);

    public String combinedFeatureName;
    public Vector<String> featureNames = new Vector<>();


    /** the real empirical combination */
    public ErrorEstimationDistribution combinedEmpiricalFoldChangeDistrib;

    /** mean, variance */
    public Vector<UPair<Double>> perFeatureDistribEstimates = new Vector<>();
    //scaled empirical distributions
    public Vector<ErrorEstimationDistribution> perFeatureScaledDistributions = new Vector<>();

    HashMap<Integer, UPair<Double>> confidenceLevel2FoldchangeInterval = null;

    public double estimatedFC = 0.0;

    public double totalZsum = 0.0;
    public double pval = 1.0;

    public double fdr;

    //public double fcDistribPval = 1.0;
    //public double fcDistribFDR = 1.0;

    public double fcEstimatePval = 1.0;
    public double fcEstimateFDR = 1.0;

    public double estimatedVariance = 0.0;
    public double estimatedSD = 0.0;

    public Vector<Double> sdCorrectionScales = new Vector<>();

    public static final double DEFAULT_PVAL_THRESHOLD = 0.05;
    Double nonDiffExpFcThresholdToDefaultThreshold = null;

    HashMap<Double, UPair<Double>> nonDiffExpThreshold2PvalAndFDR = new HashMap<>();

    public Vector<UPair<Double>> perFeaturePvals = new Vector<>();

    public boolean gotCalculatedNonDiffExpFDR(double threshold) {
        return null != nonDiffExpThreshold2PvalAndFDR.get(threshold);
    }


    public double getNonDiffExpLog2FCThresholdToDefaultPvalThreshold() {
        if(nonDiffExpFcThresholdToDefaultThreshold != null)
            return nonDiffExpFcThresholdToDefaultThreshold;

        if(combinedEmpiricalFoldChangeDistrib == null)
            return Double.POSITIVE_INFINITY;

        //look for the log2fc threshold containing 1.0 - DEFAULT_PVAL_THRESHOLD probability mass
        return nonDiffExpFcThresholdToDefaultThreshold = combinedEmpiricalFoldChangeDistrib.getCenteredFCWithTargetProbabilityMass(1.0 - DEFAULT_PVAL_THRESHOLD);

    }

    protected double calcNonDiffExpPval(double threshold) {
        UPair<Double> pre = nonDiffExpThreshold2PvalAndFDR.get(threshold);
        if(pre == null) {

            double pval = 1.0;
            if(combinedEmpiricalFoldChangeDistrib != null) {
                double sumIntegral = combinedEmpiricalFoldChangeDistrib.getCumulativeFrequencyToFoldChange(threshold) - combinedEmpiricalFoldChangeDistrib.getCumulativeFrequencyToFoldChange(-threshold);
                pval = 1.0 - sumIntegral;
            }
            pval = Math.max(Math.pow(10, -8), pval);
            pre = UPair.createU(pval, 1.0);
            nonDiffExpThreshold2PvalAndFDR.put(threshold, pre);
        }
        return pre.getFirst();
    }
    /**
     *
     * @param confidence_level between 0 and 100 (exclusive)
     * @return
     */
    public double getFoldChangeWidthToConfidenceLevel(int confidence_level) {
        UPair<Double> fcInterval = getFoldChangeIntervalToConfidenceLevel(confidence_level);
        return fcInterval.getSecond() - fcInterval.getFirst();

    }

    public boolean isMeasured() {
        return getNumUsedFeatures() > 0;
    }

    static final UPair<Double> UNMEASURED_FC_INTERVAL = UPair.createU(Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);

    /**
     *
     * @param confidence_level between 0 and 100 (exclusive)
     * @return
     */
    public UPair<Double> getFoldChangeIntervalToConfidenceLevel(int confidence_level) {
        if(combinedEmpiricalFoldChangeDistrib == null)
            return UNMEASURED_FC_INTERVAL;

        if(confidence_level < 1 || confidence_level > 99)
            throw new FRuntimeException("unexpected confidence level: %d, it should be 0 < confidence level < 100", confidence_level);
        UPair<Double> fcInterval = (confidenceLevel2FoldchangeInterval == null) ? null : confidenceLevel2FoldchangeInterval.get(confidence_level);
        if(fcInterval != null)
            return fcInterval;

        if(confidenceLevel2FoldchangeInterval == null) {
            confidenceLevel2FoldchangeInterval = new HashMap<>();
        }
        double alpha = 0.5 * ((100 - confidence_level) / 100.0);
        double fc_low = combinedEmpiricalFoldChangeDistrib.getFoldChangeToCumulativeFrequency(alpha);
        double fc_high = combinedEmpiricalFoldChangeDistrib.getFoldChangeToCumulativeFrequency(1.0 - alpha);
        confidenceLevel2FoldchangeInterval.put(confidence_level, fcInterval = UPair.createU(fc_low, fc_high));
        return fcInterval;

    }

    public DiffExpResult(DiffExpManager diffExpManager, String featureName) {
        this(diffExpManager, featureName, toVector(featureName));
    }

    public static Vector<SingleFeatureDiffExp> getSingleDiffExps(DiffExpManager diffExpManager,  Collection<String> featureNames) {
        Vector<SingleFeatureDiffExp> singleFeatureDiffExps = new Vector<>();
        for(String feature : featureNames) {
            SingleFeatureDiffExp singleFeatureDiffExp = new SingleFeatureDiffExp(diffExpManager, feature);
            if(!singleFeatureDiffExp.useable)
                continue;
            singleFeatureDiffExps.add(singleFeatureDiffExp);
        }
        return singleFeatureDiffExps;
    }

    public DiffExpResult(DiffExpManager diffExpManager, String combinedFeatureName, Vector<SingleFeatureDiffExp> tocombine, boolean fromSingleDiffExps) {

        this.diffExpManager = diffExpManager;
        this.combinedFeatureName = combinedFeatureName;

        tocombine = filter(tocombine, (_sfde) -> _sfde != null && _sfde.useable);

        if(tocombine.size() == 0)
            return;

        Vector<ErrorEstimationDistribution> fcDiffErrors = new Vector<>();
        int numFeaturesUsed = tocombine.size();


        if(PERCENT_ROBUST_FEATURE_USED < 1.0) {
            double meanfc = NumUtils.getMedianMean(map(tocombine, (_f) -> _f.featureEstimatedFC), true);
            NumUtils.sort(tocombine, (_t) -> Math.abs(_t.featureEstimatedFC - meanfc), false);

            int selected = (int)Math.ceil(PERCENT_ROBUST_FEATURE_USED * tocombine.size());
            //  System.out.printf("%d -> %d\n", tocombine.size(), selected);
            if(selected > 0) {
                tocombine = VectorUtils.slice(tocombine, 0, selected);
            }

        }

        for(SingleFeatureDiffExp singleFeatureDiffExp : tocombine) {

            this.featureNames.add(singleFeatureDiffExp.feature);

            fcDiffErrors.add(singleFeatureDiffExp.scaledDistribution);
            sdCorrectionScales.add(singleFeatureDiffExp.sdCorrectionScale);
            perFeatureDistribEstimates.add(UPair.createU(singleFeatureDiffExp.getFC(), singleFeatureDiffExp.getSD()));

            estimatedVariance += Math.pow(singleFeatureDiffExp.getSD(), 2.0);

            perFeatureScaledDistributions.add(singleFeatureDiffExp.scaledDistribution);

            if(diffExpManager.empiRe.collect_comparison_pvals) {
                diffExpManager.empiRe.per_comparison_pvals.add(singleFeatureDiffExp.pval);
                diffExpManager.empiRe.per_comparison_fcdist_pvals.add(singleFeatureDiffExp.fcdist_pval);

                perFeaturePvals.add(UPair.createU(singleFeatureDiffExp.pval, singleFeatureDiffExp.fcdist_pval));
            }
            totalZsum += singleFeatureDiffExp.featureZscore;

        }

        /** combine the per feature evidences */
        combinedEmpiricalFoldChangeDistrib = (SingleFeatureDiffExp.USE_NO_COMBINED_DISTRIBS) ? null : ErrorEstimationDistribution.combineAndAverage(perFeatureScaledDistributions);

        estimatedSD = Math.sqrt(estimatedVariance) / numFeaturesUsed;
        estimatedVariance = Math.pow(estimatedSD, 2);
        estimatedFC = (SingleFeatureDiffExp.USE_NO_COMBINED_DISTRIBS) ? NumUtils.getMedianMean(map(tocombine, (_f) -> _f.featureEstimatedFC)) : combinedEmpiricalFoldChangeDistrib.getMostProbableFcWindowCenter();

        if(combinedEmpiricalFoldChangeDistrib != null) {
            fcEstimatePval = combinedEmpiricalFoldChangeDistrib.getCumulativeFrequencyToFoldChange(0.0);
            fcEstimatePval = 2.0 * (Math.min(fcEstimatePval, 1.0 - fcEstimatePval));

        }

        // z-score (a.k.a. Stouffer-like) combination of the evidences -> all independent
        pval = 2.0 * (1.0 - new NormalDistribution(0, Math.sqrt(this.featureNames.size())).cumulativeProbability(Math.abs(totalZsum)));
    }

    public DiffExpResult(DiffExpManager diffExpManager, String combinedFeatureName, Collection<String> featureNames) {
        this(diffExpManager, combinedFeatureName, getSingleDiffExps(diffExpManager, featureNames), true);

    }

    public Vector<ErrorEstimationDistribution> getFeatureErrorDistribs() {
        return map(perFeatureDistribEstimates, (_up) -> ErrorEstimationDistribution.fromNormalDistrib(_up.getFirst(), _up.getSecond()));
    }
    public int getNumUsedFeatures() {
        return featureNames.size();
    }

    public DiffExpManager getDiffExpManager() {
        return diffExpManager;
    }

    public static String getPerformanceString(Vector<DiffExpResult> results, Function<String, Boolean> trueLabeller) {
        PerformanceResult[] result = getPerformances(results, trueLabeller);
        return String.format("%s\n%s\n", result[0], result[1]);
    }

    public static String getPerformanceString(Vector<DiffExpResult> results, HashSet<String> trueLabels) {
        PerformanceResult[] result = getPerformances(results, trueLabels);
        return String.format("%s\n%s\n", result[0], result[1]);
    }
    public static PerformanceResult[] getPerformances(Vector<DiffExpResult> results, HashSet<String> trueLabels) {
        return getPerformances(results, (_n) -> trueLabels.contains(_n), 0.05, 0.0);
    }

    public static PerformanceResult[] getPerformances(Vector<DiffExpResult> results, Function<String, Boolean> trueLabeller) {
        return getPerformances(results, trueLabeller, 0.05, 0.0);
    }

    public static PerformanceResult[] getPerformances(Vector<DiffExpResult> results, HashSet<String> trueLabels, double fdr_treshold, double abs_fc_treshold) {
        return getPerformances(results, (_n) -> trueLabels.contains(_n), fdr_treshold, abs_fc_treshold);
    }

    public static PerformanceResult[] getPerformances(Vector<DiffExpResult> results, Function<String, Boolean> trueLabeller, double fdr_treshold, double abs_fc_treshold) {

        String threshold_suffix = String.format(" fdr<= %.3f abs(log2FC) >= %.2f", fdr_treshold, abs_fc_treshold);
        return new PerformanceResult[] {
                new PerformanceResult("nlEmpiRe.zscore" + threshold_suffix, results, (_d) -> _d.pval, (_d) -> trueLabeller.apply(_d.combinedFeatureName), false,
                        (_d) -> _d.fdr  <= fdr_treshold && Math.abs(_d.estimatedFC) >= abs_fc_treshold, null, false, RocInfo.PerformanceEvaluationStrategy.AVERAGE),

                new PerformanceResult("nlEmpiRe.FCdistrib" + threshold_suffix, results, (_d) -> _d.fcEstimatePval, (_d) -> trueLabeller.apply(_d.combinedFeatureName), false,
                        (_d) -> _d.fcEstimateFDR  <= fdr_treshold && Math.abs(_d.estimatedFC) >= abs_fc_treshold, null, false, RocInfo.PerformanceEvaluationStrategy.OPTIMISTIC)
        };

    }

}
