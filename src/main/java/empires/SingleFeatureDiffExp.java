package empires;

import lmu.utils.LogConfig;
import lmu.utils.Logger;
import lmu.utils.NumUtils;
import lmu.utils.UPair;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.util.Vector;

import static lmu.utils.ObjectGetter.map;

public class SingleFeatureDiffExp {

    public static boolean USE_NO_COMBINED_DISTRIBS = false;
    public static boolean COMBINE_SHIFTED_DISTRIBS_FOR_EACH_REPLICATE_PAIR = false;
    public static boolean CORRECT_OUTLIERS = true;
    Logger log = LogConfig.getLogger();
    public String feature;
    ErrorEstimationDistribution  diffErr;

    public boolean useable = false;
    double sdCorrectionScale = 1.0;
    UPair<Vector<Double>> normedValues;
    double perFeatureSum = 0;
    int numPointsInVariance;

    double featureEstimatedFC;
    double featureEstimatedSD;

    public ErrorEstimationDistribution scaledDistribution;
    double pval;
    double fcdist_pval;
    double featureZscore;


    public SingleFeatureDiffExp(DiffExpManager diffExpManager, String feature) {
        this.feature = feature;
        boolean gotInFrom = diffExpManager.replicateSetFrom.gotFeatureMeasured(feature);
        boolean gotInto = diffExpManager.replicateSetTo.gotFeatureMeasured(feature);

        if(!gotInFrom && !gotInto)
            return;

        if(!diffExpManager.empiRe.compareFeaturesMeasuredOnlyOnOneSide && (!gotInFrom || !gotInto))
            return;

        ReplicatedMeasurement repFrom = diffExpManager.replicateSetFrom.getNormed(feature);
        ReplicatedMeasurement repTo = diffExpManager.replicateSetTo.getNormed(feature);


        int nrep1 = repFrom.nonNanValues.size();
        int nrep2 = repTo.nonNanValues.size();

        if(nrep1 == 0 || nrep2 == 0)
            return;

        init(diffExpManager, repFrom.nonNanValues, repFrom.mean, repTo.nonNanValues, repTo.mean, diffExpManager.replicateSetFrom.getError(feature), diffExpManager.replicateSetTo.getError(feature));
    }

    public SingleFeatureDiffExp(String featureName, DiffExpManager diffExpManager, Vector<Double> reps1, Vector<Double> reps2) {
        this.feature = featureName;
        int nrep1 = reps1.size();
        int nrep2 = reps2.size();

        if(nrep1 == 0 || nrep2 == 0)
            return;

        double repFromMean = NumUtils.getMedianMean(reps1, true); //changed mean to median - effect??? TOCHECK
        double repToMean = NumUtils.getMedianMean(reps2, true);

        ErrorEstimationDistribution err1 = diffExpManager.replicateSetFrom.getErrorByMeanSignalValue(repFromMean);
        ErrorEstimationDistribution err2 = diffExpManager.replicateSetFrom.getErrorByMeanSignalValue(repToMean);

        init(diffExpManager, reps1, repFromMean, reps2, repToMean, err1, err2);

    }

    public SingleFeatureDiffExp(String featureName, Vector<Double> reps1, ErrorEstimationDistribution err1, Vector<Double> reps2, ErrorEstimationDistribution err2, DiffExpManager diffExpManager) {
        this.feature = featureName;
        int nrep1 = reps1.size();
        int nrep2 = reps2.size();

        if(nrep1 == 0 || nrep2 == 0)
            return;

        double repFromMean = NumUtils.getMedianMean(reps1, true); //changed mean to median - effect??? TOCHECK
        double repToMean = NumUtils.getMedianMean(reps2, true);

        init(diffExpManager, reps1, repFromMean, reps2, repToMean, err1, err2);

    }

    private void init(DiffExpManager diffExpManager, Vector<Double> reps1, double repFromMean, Vector<Double> reps2, Double repToMean,
                      ErrorEstimationDistribution rep1error, ErrorEstimationDistribution rep2error) {

        int nrep1 = reps1.size();
        int nrep2 = reps2.size();
        numPointsInVariance = nrep1 * nrep2;

        if(nrep1 == 0 || nrep2 == 0)
            return;

        /** there are valid measurements pairs -> we can use this*/
        useable = true;


        /** the replictes in reps1 and reps2 (normalized log2 signal values) provide differential evidences
         *  for the changes, i.e. all r1 - r2  are evidences for the change between the two replicate sets
         *
         *  as all r1 have a measurement uncertainty (i.e. error background distribution) encoded in the
         *          feature-related ErrorEstimationDistribution in the normalized replicate set "from" in input: "rep1error"
         *   and all r2 have a measurement uncertainty (i.e. error background distribution) encoded in the
         *                    feature-related ErrorEstimationDistribution in the normalized replicate set: "to" in input "rep2error"
         *
         *   -> the values r1 - r2 will have the uncertainty = background distribution
         *    diffErr = rep1error - rep2error
         *
         *      this is provided (and cached) by DiffExpManager (cache as the same difference-distribution might be used
         *      from other features
          */
        diffErr = diffExpManager.getDiffError(rep1error, rep2error);
        double diffErrorVariance = diffExpManager.getDiffError(rep1error, rep2error).getVariance();

        if(CORRECT_OUTLIERS) {
            /** check if this measurement contains outliers -> if so scale up the variance*/
            double calcedFromSD = (reps1.size() == 1) ? rep1error.SD : Math.sqrt(NumUtils.sum(map(reps1, (_d) -> Math.pow(repFromMean - _d, 2.0))) / reps1.size());
            double calcedToSD = (reps2.size() == 1) ? rep2error.SD : Math.sqrt(NumUtils.sum(map(reps2, (_d) -> Math.pow(repToMean - _d, 2.0))) / reps2.size());

            double corrSD1 = Math.max(rep1error.SD,  calcedFromSD);
            double corrSD2 = Math.max(rep2error.SD,  calcedToSD);
            double corrCombinedSD = Math.sqrt(Math.pow(corrSD1, 2.0) + Math.pow(corrSD2, 2.0));
            sdCorrectionScale = Math.max(1.0, corrCombinedSD / diffErr.getSD());

        }

        /**
            independentPart = nrep1 * nrep2 * diffErrorVariance;
            totalVariance = independentPart +
                nrep1 * nrep2 * (nrep2  - 1 ) * diffExpManager.replicateSetFrom.getError(feature).getVariance() +
                nrep1 * nrep2 * (nrep1  - 1 ) * diffExpManager.replicateSetTo.getError(feature).getVariance();

            per distrib back-scale  - divide by (nrep1 * nrep2)
         */

        double perEvidenceVariance = diffErrorVariance
                            + (nrep2  - 1 ) * diffExpManager.replicateSetFrom.getError(feature).getVariance()
                            + (nrep1  - 1 ) * diffExpManager.replicateSetTo.getError(feature).getVariance();


        double totalVariance = perEvidenceVariance * nrep1 * nrep2;

        featureEstimatedSD = Math.sqrt(perEvidenceVariance);

        featureEstimatedSD *= sdCorrectionScale;

        double perEvidenceSDScale = sdCorrectionScale * Math.sqrt(perEvidenceVariance / diffErrorVariance);

        ErrorEstimationDistribution perEvidenceScaledDiffErr = (USE_NO_COMBINED_DISTRIBS) ? null :  (!COMBINE_SHIFTED_DISTRIBS_FOR_EACH_REPLICATE_PAIR) ? null : diffErr.scale(perEvidenceSDScale);
        if(perEvidenceScaledDiffErr != null) {
            perEvidenceScaledDiffErr.getCumulative(); //precalc
        }

        double zSum = 0.0;
        Vector<Double> fromVals = map(reps1, (_d) -> _d + diffExpManager.shifts.get(0));
        Vector<Double> toVals = map(reps2, (_d) -> _d + diffExpManager.shifts.get(1));

        normedValues = UPair.createU(fromVals, toVals);

        Vector<Double> difffcs =  new Vector<>();
        Vector<ErrorEstimationDistribution> single_shifted_errordistribs = new Vector<>();


        for(double logValFrom : fromVals) {
            for(double logValTo : toVals) {
                double diff = logValFrom  - logValTo;
                difffcs.add(diff);
                perFeatureSum += diff;

                zSum += diffErr.getConvertedZscore(diff) / sdCorrectionScale;
                if(!USE_NO_COMBINED_DISTRIBS && COMBINE_SHIFTED_DISTRIBS_FOR_EACH_REPLICATE_PAIR) {
                    single_shifted_errordistribs.add(new ErrorEstimationDistribution(perEvidenceScaledDiffErr, diff));
                }

            }
        }

        featureEstimatedFC = perFeatureSum / (nrep1 * nrep2);


        if(!USE_NO_COMBINED_DISTRIBS) {
            if(!COMBINE_SHIFTED_DISTRIBS_FOR_EACH_REPLICATE_PAIR) {
                //the total variance of the summed evidences is totalVariance -> total SD = Math.sqrt(totalVariance)
                // this value space would be however for nrep1 * nrep2 summed values -> normalize it back to one to have fold-change-space interpretable distribution
                scaledDistribution = diffErr.scale(perEvidenceSDScale / Math.sqrt(nrep1 * nrep2)).shift(featureEstimatedFC);

            } else {
                scaledDistribution = ErrorEstimationDistribution.combineAndAverage(single_shifted_errordistribs);
                featureEstimatedFC = scaledDistribution.getMostProbableFcWindowCenter();
            }
        }


        /* stouffer like calculation */
        double normed_sum_sd = sdCorrectionScale *   Math.sqrt(totalVariance / diffErrorVariance);
        pval = 2.0 * (1.0 - new NormalDistribution(0, normed_sum_sd).cumulativeProbability(Math.abs(zSum)));
        featureZscore = zSum / normed_sum_sd;

        if(scaledDistribution != null) {
            /** direct lookup in fc-distrib */
            fcdist_pval = scaledDistribution.getCumulativeFrequencyToFoldChange(0.0);
            fcdist_pval = 2.0 * (Math.min(fcdist_pval, 1.0 - fcdist_pval));
        }


    }


    public double getFC() {
        return featureEstimatedFC;
    }

    public double getSD() {
        return featureEstimatedSD;
    }
}
