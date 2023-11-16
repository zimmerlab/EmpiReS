package empires;

import lmu.utils.*;
import lmu.utils.tuple.Tuple3;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.util.*;

import static lmu.utils.IteratorUtils.rangev;
import static lmu.utils.ObjectGetter.*;
import static lmu.utils.ObjectGetter.apply;

public class AutoBackGroundContextProvider implements BackgroundContextFuzzficationStrategyProvider {
    Logger log = LogConfig.getLogger();
    public static double MERGE_BACKGROUNDS_WITH_MAX_SD_DIFF = 0.0001;
    static final double OUTLIER_SIGNIFICANCE_THRESHOLD = 0.0; //0.0000001;
    static final double OUTLIER_CORRECTED_SIGNIFICANCE_TARGET = 0.10;

    public BackGroundContext getContexts(String replicateSetName, Vector<String> sampleNames, Vector<ReplicatedMeasurement> normalized) {

        BackGroundContext backGroundContext = new BackGroundContext();
        log.info("start calc background for replicate set: %s", replicateSetName);
        NumUtils.sort(normalized, (_r) -> _r.mean);



        int[] samplingSizes = new int[normalized.size()];
        int ntotalAvailableEstimates = 0;
        for(int i=0; i<normalized.size(); i++) {
            ntotalAvailableEstimates += normalized.get(i).nonNanValues.size() - 1;
            samplingSizes[i] = ntotalAvailableEstimates >> 1; //we build pair of measurements -> only the half of the points can be used
        }


//
        /** ideally we would have > 200 contexts (overlapping) with each at least 1000 values in it
         * for this we need 101 distinct ones
         */
        int num_contexts = 100;
        int nEstimatesPerContext = Math.max(1000, (int)( ntotalAvailableEstimates / (1.0 + (num_contexts >> 1))));
        int halfcontext_sample_size = nEstimatesPerContext >> 1;
        int[] contextBoundaries = new int[3];

        int hitpos = Arrays.binarySearch(samplingSizes, halfcontext_sample_size);
        contextBoundaries[1] = (hitpos >= 0) ? hitpos : -hitpos - 1;

        log.info("nsample per context: %d total: %d\n", nEstimatesPerContext, ntotalAvailableEstimates);
        hitpos = Arrays.binarySearch(samplingSizes, nEstimatesPerContext);
        contextBoundaries[2] = (hitpos >= 0) ? hitpos : Math.min(samplingSizes.length, -hitpos - 1);

        Vector<Tuple3<Integer, Integer, ErrorEstimationDistribution>> slides = new Vector<>();

        while(contextBoundaries[1] != samplingSizes.length) {

            ErrorEstimationDistribution err = new ErrorEstimationDistribution(contextBoundaries[0], contextBoundaries[2], (_i) -> normalized.get(_i).nonNanValues);
            if(err.min >= 0 || err.max <=0) {
                log.info("invalid boundary %d-%d min/max: %d/%d\n", contextBoundaries[0], contextBoundaries[2], err.min, err.max);
                hitpos = Arrays.binarySearch(samplingSizes, contextBoundaries[1], samplingSizes.length, samplingSizes[contextBoundaries[2]] + nEstimatesPerContext);
                contextBoundaries[2] = (hitpos >= 0) ? hitpos : Math.min(normalized.size(), -hitpos - 1);
                continue;
            }


            slides.add(Tuple3.create(contextBoundaries[0], contextBoundaries[2], err));

            contextBoundaries[0] = contextBoundaries[1];
            contextBoundaries[1] = contextBoundaries[2];
            hitpos = Arrays.binarySearch(samplingSizes, contextBoundaries[1], samplingSizes.length, samplingSizes[contextBoundaries[0]] + nEstimatesPerContext);
            contextBoundaries[2] = (hitpos >= 0) ? hitpos : Math.min(normalized.size(), -hitpos - 1);

        }
        log.info("created %d contextes from %d features, features in slides: %d\n", slides.size(), normalized.size(), nEstimatesPerContext);

        log.info("pre sds: %s", map(slides, (_s) -> String.format("%.3f", _s.get2().getSD())));
        if(slides.size() == 0) {
            throw new FRuntimeException("could not build background distributions for repicate set: %s too few data????", replicateSetName);
        }
        int lastsize = 0;
        int iter = 0;
        while(lastsize != slides.size()) {
            lastsize = slides.size();
            iter++;
            // combine backgrounds with similar SD (or similar overall distrib?)
            //simple solution -> hash them with the resolution and merge all contigs with same values
            Vector<Integer> sdHash = map(slides, (_t) -> (int) (_t.get2().getSD() / MERGE_BACKGROUNDS_WITH_MAX_SD_DIFF));

            Vector<Tuple3<Integer, Integer, ErrorEstimationDistribution>> merged = new Vector<>();
            int lastKey = sdHash.get(0);
            Vector<Tuple3<Integer, Integer, ErrorEstimationDistribution>> tomerge = new Vector<>();


            for (int i = 0; i <= sdHash.size(); i++) {

                if (i == sdHash.size() || sdHash.get(i) != lastKey) {
                    //log.info("%d: merge %d error distribs with key: %d", i, tomerge.size(), lastKey);
                    ErrorEstimationDistribution mergedErr = new ErrorEstimationDistribution(map(tomerge, (_t) -> _t.get2()), map(tomerge, (_t) -> 0.0), false);
                    merged.add(Tuple3.create(tomerge.get(0).get0(), tomerge.get(tomerge.size() - 1).get1(), mergedErr));
                    lastKey = (i == sdHash.size()) ? -1 : sdHash.get(i);
                    tomerge.clear();
                }
                if (i == sdHash.size())
                    break;


                tomerge.add(slides.get(i));
            }

            log.info("iter: %d merged %d contextes to %d contextes max sd: %.3f", iter, slides.size(), merged.size(), MERGE_BACKGROUNDS_WITH_MAX_SD_DIFF);
            slides = merged;
        }

        //log.info("sds: %s", map(slides, (_s) -> String.format("%.3f", _s.get2().getSD())));

        double maxSD = NumUtils.max(slides, (_e) -> _e.get2().getSD());
        backGroundContext.errorBackGrounds = map(slides, (_m) -> _m.get2());
        backGroundContext.meanSignalValue2Error = map(slides, (_m) -> NumUtils.mean(map(rangev(_m.get0(), _m.get1()), (_i) -> normalized.get(_i).mean)));



        HashMap<Integer, Vector<Integer>> outlierSD2outlierFeatures = new HashMap<>();
        HashMap<Integer, Vector<Integer>> feature2reassign = new HashMap<>();
        HashMap<Integer, Integer> outlierSDFreq = new HashMap<>();

        for (int i = 0; i < normalized.size(); i++) {
            Tuple3<Integer, Integer, ErrorEstimationDistribution> q = Tuple3.create(i, 0, null);
            int hitIdx = Collections.binarySearch(slides, q, (_t1, _t2) -> Integer.compare(_t1.get0(), _t2.get0()));
            if (hitIdx < 0) {
                hitIdx = -hitIdx - 1;
            }
            hitIdx = Math.max(0, hitIdx - 1);

            int bestDist = Integer.MAX_VALUE;
            int bestIdx = -1;
            for (int idx = hitIdx; idx < slides.size() && idx < hitIdx + 3; idx++) {
                Tuple3<Integer, Integer, ErrorEstimationDistribution> target = slides.get(idx);
                int dist = Math.max(Math.abs(target.get0() - idx), Math.abs(target.get1() - idx));
                if (dist < bestDist) {
                    bestDist = dist;
                    bestIdx = idx;
                }
            }

            ErrorEstimationDistribution err = backGroundContext.errorBackGrounds.get(bestIdx);

            double median = NumUtils.getMedianMean(normalized.get(i).nonNanValues);

            double sumZ = NumUtils.sum(map(normalized.get(i).nonNanValues, (_v) -> err.getConvertedZscore(_v - median)));
            NormalDistribution nd = new NormalDistribution(0.0, Math.sqrt(normalized.get(i).nonNanValues.size() - 1));

            double pval = 2 * (1.0 - nd.cumulativeProbability(Math.abs(sumZ)));


            if(pval < OUTLIER_SIGNIFICANCE_THRESHOLD) {

                double val = nd.inverseCumulativeProbability(1.0 - ( OUTLIER_CORRECTED_SIGNIFICANCE_TARGET) / 2.0);
                double estimatedSD = Math.abs(sumZ / val);
                double testPval = 2 * (1.0 - new NormalDistribution(0, estimatedSD).cumulativeProbability(Math.abs(sumZ)));

                //System.out.printf("testp: %g target: %g\n", testPval, OUTLIER_CORRECTED_SIGNIFICANCE_TARGET);

                int key = (int)(estimatedSD / MERGE_BACKGROUNDS_WITH_MAX_SD_DIFF);

                //System.out.printf("z: %.2f pval: %.2f zscores: %s\n", z, pval,  map(normalized.get(i).nonNanValues, (_v) -> err.getConvertedZscore(_v - median)));
                //System.out.printf("est SD: %.2f vs %.2f pval of est sd: %.2f\n", estimatedSD, err.getSD(), 2 * (1.0 - nd.cumulativeProbability(Math.abs(sumErr / estimatedSD))));
                if(estimatedSD > maxSD) {
                    MapBuilder.updateV(outlierSD2outlierFeatures, key, i);


                } else {

                    MapBuilder.updateV(feature2reassign, key, i);
                }

            } else {

                backGroundContext.featureIdx2ErrorBackgroundIdx.put(normalized.get(i).featureIdx, bestIdx);


            }




        }



        if(feature2reassign.size() > 0) {
            Vector<UPair<Integer>> sd2erridx = mapIndex(backGroundContext.errorBackGrounds.size(), (_i) -> UPair.createU((int)(backGroundContext.errorBackGrounds.get(_i).getSD() / MERGE_BACKGROUNDS_WITH_MAX_SD_DIFF), _i));
            NumUtils.sort(sd2erridx, (_p) -> _p.getFirst());
            UPair<Integer> q = UPair.createU(-1, -1);
            for(Map.Entry<Integer, Vector<Integer>> sd2outliers  : feature2reassign.entrySet()) {
                q.setFirst(sd2outliers.getKey());
                int hitIdx = Collections.binarySearch(sd2erridx, q, (_p1, _p2) -> Integer.compare(_p1.getFirst(),  _p2.getFirst()));
                if(hitIdx < 0) {
                    hitIdx = Math.max(0, - hitIdx - 2);
                }
                int ERRIDX = hitIdx;
                apply(sd2outliers.getValue(), (_normedIndex) -> backGroundContext.featureIdx2ErrorBackgroundIdx.put(normalized.get(_normedIndex).featureIdx, ERRIDX));
            }
        }

        if(outlierSD2outlierFeatures.size() > 0) {

            for(Map.Entry<Integer, Vector<Integer>> sd2outliers : outlierSD2outlierFeatures.entrySet()) {
                Vector<Integer> outlierFeatures = sd2outliers.getValue();
                int sdKey = sd2outliers.getKey();
                //check if we have enough values
                int numsamples = NumUtils.sum(map(outlierFeatures, (_idx) -> normalized.get(_idx).nonNanValues.size())) >> 1;
                Vector<Double> extra = new Vector<>();
                log.trace("got %d outlier sample values for %d outlier features", numsamples, outlierFeatures.size());
                NormalDistribution sampleSD = new NormalDistribution(0, sdKey * MERGE_BACKGROUNDS_WITH_MAX_SD_DIFF);
                extra.addAll(mapIndex(1000 + 1, (_i) -> sampleSD.sample()));

                ErrorEstimationDistribution outlierErr = new ErrorEstimationDistribution(-1, outlierFeatures.size(),
                        (_i) -> (_i < 0) ? extra : normalized.get(outlierFeatures.get(_i)).nonNanValues);

                int ERR_IDX = backGroundContext.errorBackGrounds.size();
                backGroundContext.errorBackGrounds.add(outlierErr);
                apply(outlierFeatures, (_o) -> backGroundContext.featureIdx2ErrorBackgroundIdx.put(normalized.get(_o).featureIdx, ERR_IDX));


                log.trace("outlier sd: %.2f\n", outlierErr.getSD());
            }



        }
        log.info("got %d/%d outliers , %d reassigned", NumUtils.sum(outlierSD2outlierFeatures.values(), (_v) -> _v.size()), normalized.size(), NumUtils.sum(feature2reassign.values(), (_v) -> _v.size()));

        return backGroundContext;
    }
}
