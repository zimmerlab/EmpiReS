package empires;

import lmu.utils.NumUtils;
import lmu.utils.Pair;
import lmu.utils.UPair;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.util.Collections;
import java.util.HashMap;
import java.util.Vector;

import static lmu.utils.ObjectGetter.*;

public class DistribUtils {


    static double NORMALWFACTOR = 4.0;

    static Vector<Pair<Double, Integer>> getNormalAsInt(double step)
    {
        int n = (int)(1.0 / step);
        Vector<Double> LIMITS = mapIndex(n-1, (_i) -> 1.0 - ((_i + 1) / (n + 0.0)));
        NormalDistribution nd = new NormalDistribution(0, 2.0);
        Vector<Double> weights = map(LIMITS, (_d) -> nd.density(nd.inverseCumulativeProbability(_d)));
        double min = NumUtils.min(weights);
        Vector<Integer> iweights = map(weights, (_w) -> (int)Math.round(NORMALWFACTOR * _w / min));

        return mapIndex(LIMITS.size(), (_i) -> Pair.create(LIMITS.get(_i), iweights.get(_i)));
    }

    public static Vector<Pair<Double, Integer>> normalSampling = getNormalAsInt(0.05);
    public static Vector<Double> LIMITS = toVector(0.95, 0.9, 0.80, 0.65, 0.5);


    public static String checkUniform(Vector<Double> pvals, double errAllowed) {
        return checkUniform(pvals, errAllowed, 100);
    }

    /** returns the biggest error and percentage above errAllowed or null
     * if there is none
     * @param pvals
     * @param errAllowed
     * @return
     */
    public static String checkUniform(Vector<Double> pvals, double errAllowed, int maxPercent) {
        Collections.sort(pvals);
        double maxDiff = 0.0;
        double maxDiffPercent = 0.0;
        double diffVal = 0.0;
        for(int i=0; i< maxPercent; i+=2) {
            double perc = i / 100.0;
            double totest = pvals.get((int)(pvals.size() * perc));

            double diff = Math.abs(totest - perc);
            if(diff <= maxDiff)
                continue;

            maxDiff = diff;
            maxDiffPercent = perc;
            diffVal = totest;
        }

        if(Math.abs(maxDiff) < errAllowed)
            return null;

        return String.format("err: %.2f at %.2f%%: %.2f", maxDiff, maxDiffPercent * 100, diffVal);
    }

    static HashMap<UPair<Integer>, NormalDistribution> repsumdistribs = new HashMap<>();

    /**
     *
     * @param nrep1
     * @param nrep2
     * @param diffzsum
     * @return pvalue for unsigned test ( - 4 and 4 will rsult the same pval!)
     */
    public static double getReplicateZSumPvalue(int nrep1, int nrep2, double diffzsum) {
        UPair<Integer> key = UPair.createU(nrep1, nrep2);
        NormalDistribution nd = repsumdistribs.get(key);
        if(nd == null) {
            repsumdistribs.put(key, nd = new NormalDistribution(0.0, Math.sqrt(getReplicateSummarizationVariance(nrep1, nrep2))));
        }

        return (diffzsum < 0 ) ? nd.cumulativeProbability(diffzsum) : 1.0 - nd.cumulativeProbability(diffzsum);

    }

    public static double getReplicateSummarizationVariance(int nrep1, int nrep2) {
        return nrep1 * nrep2 + ((nrep1 * nrep2 == 1) ? 0 : 0.5 * (nrep1 * nrep2 * (nrep1 + nrep2 - 2)));
    }


}
