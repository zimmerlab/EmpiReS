package nlEmpiRe;

import org.apache.logging.log4j.Logger;
import lmu.utils.*;
import lmu.utils.fdr.BenjaminiHochberg;

import java.util.HashMap;
import java.util.Vector;

import static lmu.utils.ObjectGetter.*;

public class DiffExpManager {
    Logger log = LogConfig.getLogger();
    EmpiRe empiRe;
    public Vector<String> totalFeatures;
    public NormalizedReplicateSet replicateSetFrom;
    public NormalizedReplicateSet replicateSetTo;
    Vector<Double> med1;
    Vector<Double> med2;
    public Vector<Double> shifts;
    HashMap<UPair<ErrorEstimationDistribution>, ErrorEstimationDistribution> diffDistributionMap = new HashMap<>();



    public DiffExpManager(NormalizedReplicateSet replicateSetFrom, NormalizedReplicateSet replicateSetTo, Vector<String> totalFeatures, EmpiRe empiRe) {
        this(replicateSetFrom, replicateSetTo, totalFeatures, empiRe, totalFeatures);
    }
    public DiffExpManager(NormalizedReplicateSet replicateSetFrom, NormalizedReplicateSet replicateSetTo, Vector<String> totalFeatures, EmpiRe empiRe, Vector<String> normFeatures) {

        log.info("will norm %s to %s on %d/%d features", replicateSetFrom.getInData().getReplicateSetName(), replicateSetTo.getInData().getReplicateSetName(),
                normFeatures.size(), totalFeatures.size());

        this.totalFeatures = totalFeatures;
        Normalization norm = new Normalization(toVector(med1 = replicateSetFrom.getMedianValues(normFeatures), med2 = replicateSetTo.getMedianValues(normFeatures)));
        shifts = norm.getShifts();

        log.info("shift %s", shifts);

        this.empiRe = empiRe;
        this.replicateSetFrom = replicateSetFrom;
        this.replicateSetTo = replicateSetTo;
    }


    public Vector<Double> getMedianFCDistrib() {
        Vector<Double> normed = new Vector<>();
        double shift = shifts.get(0) - shifts.get(1);

        for(int i=0; i<med1.size(); i++) {
            if(Double.isNaN(med1.get(i)) || Double.isNaN(med2.get(i)))
                continue;

            normed.add(med1.get(i) - med2.get(i) + shift);
        }
        return normed;
    }

    public ErrorEstimationDistribution getDiffError(String feature) {
        return getDiffError(replicateSetFrom.getError(feature), replicateSetTo.getError(feature));
    }

    public synchronized ErrorEstimationDistribution getDiffError(ErrorEstimationDistribution from, ErrorEstimationDistribution to) {
        UPair<ErrorEstimationDistribution> key = UPair.createU(from, to);
        ErrorEstimationDistribution cached = diffDistributionMap.get(key);
        if(cached != null)
            return cached;

        diffDistributionMap.put(key, cached = ErrorEstimationDistribution.substract(from, to));

        return cached;
    }



    /**
     *   correlation is given by
     *
     *   https://en.wikipedia.org/wiki/Pearson_correlation_coefficient
     *
     *                  E[XY]  -    E[X]E[Y]
     *     -----------------------------------------------
     *      sqrt(E[X^2]-[E[X]^2]) * sqrt(E[Y^2]-[E[Y]^2])
     *
     *
     *
     *      in our case X and Y are differences of two empirical distributions with expected value = 0
     *      so where one of the random variables are shared i.e.
     *
     *      either  (i)  X = A - B  and Y = A - C
     *
     *      or      (ii) X = A - C and  Y = B - C
     *
     *
     *      as A,B,C have expected values 0 their differences have also 0 expected values.
     *
     *      so that for (i) the formula above reduces to:
     *
     *                 E[(A-B)(A-C)]                    E(A^2) - E(AC) - E(AB) - E(BC)
     *         ----------------------------    =  ----------------------------------------
     *         sqrt((A-B)^2) * sqrt((A-C)^2)         sqrt(var(A-B)) * sqrt(var(A-C))
     *
     *
     *       due to indepence E(AC) = E(AB) = E(BC) = 0
     *
     *       which leads to:
     *
     *       var(A) / (sqrt(var(A) + var(B)) * sqrt(var(A) + var(C))
     *
     *       as B and C are independent samples from the same distribution i.e. var(B) == var(C)
     *
     *       -> var(A) / (var(A) + var(B))
     *
     *      and
     *
     *      var(C) / (sqrt(var(A) + var(C)) * sqrt(var(B) + var(C))
     *
     *      as A and B are independent samples from the same distribution i.e. var(A) == var(B)
     *
     *      -> var(C) / (var(A) + var(C))
     *
     */


    UPair<Double> getCorrelations(String feature) {
        double var1 = replicateSetFrom.getError(feature).getVariance();
        double var2 = replicateSetTo.getError(feature).getVariance();

        double norm = 1.0 / (var1 + var2);

        /* the first value corresponds on correlations where nrep1 is the same

         */
        return UPair.createU(var1 * norm, var2 * norm);

    }

    public Vector<Pair<DiffExpResult, Double>> calculateNonDifferentialFDRs(Vector<DiffExpResult> diffExpResults, double threshold) {
        if(threshold <= 0.0)
            throw new FRuntimeException("expected some positive number!");

        //check if it is already calculated
        Vector<DiffExpResult> tocalc = filter(diffExpResults, (_r) -> _r.gotCalculatedNonDiffExpFDR(threshold));
        apply(tocalc, (_dr) -> _dr.calcNonDiffExpPval(threshold));

        BenjaminiHochberg.adjust_pvalues(diffExpResults, (_dr) -> _dr.nonDiffExpThreshold2PvalAndFDR.get(threshold).getFirst(),
                (_p) -> _p.getFirst().nonDiffExpThreshold2PvalAndFDR.get(threshold).setSecond(_p.getSecond()));

        return map(diffExpResults, (_dr) -> Pair.create(_dr, _dr.nonDiffExpThreshold2PvalAndFDR.get(threshold).getSecond()));

    }

}
