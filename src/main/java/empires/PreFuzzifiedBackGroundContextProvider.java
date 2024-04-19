package nlEmpiRe;

import lmu.utils.*;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Vector;

import static lmu.utils.IteratorUtils.rangev;
import static lmu.utils.ObjectGetter.*;
import org.apache.logging.log4j.Logger;

public class PreFuzzifiedBackGroundContextProvider implements BackgroundContextFuzzficationStrategyProvider {

    Logger log = LogConfig.getLogger();
    HashMap<String, HashMap<String, String>> feature2sample2fuzzyclass;

    public PreFuzzifiedBackGroundContextProvider(HashMap<String, HashMap<String, String>> feature2sample2fuzzyclass) {
        assert(feature2sample2fuzzyclass.size() > 0);

        this.feature2sample2fuzzyclass = feature2sample2fuzzyclass;

    }

    @Override
    public BackGroundContext getContexts(String replicateSetName, Vector<String> sampleNames, Vector<ReplicatedMeasurement> normalized) {

        BackGroundContext backGroundContext = new BackGroundContext();
        backGroundContext.errorBackGrounds = new Vector<>();
        backGroundContext.meanSignalValue2Error = new Vector<>();

        HashMap<String, HashSet<String>> bg2features = new HashMap<>();
        HashMap<String, String> feature2fuzzyConcept = new HashMap<>();

        for(String feature :  feature2sample2fuzzyclass.keySet()) {
            HashMap<String, Integer> fuzzyConcept2Freq = new HashMap<>();
            HashMap<String, String> mapping = feature2sample2fuzzyclass.get(feature);

            for(String sample : sampleNames) {
                String fuzzyConcept = mapping.get(sample);
                if(fuzzyConcept == null)
                    throw new FRuntimeException("invalid fuzzy concept mapping, no fuzzy concept is mapped to feature: %s and sample: %s", feature, sample);

                MapBuilder.update(fuzzyConcept2Freq, fuzzyConcept);
                MapBuilder.update(bg2features, fuzzyConcept, feature);
            }
            feature2fuzzyConcept.put(feature, Pair.convert_reverse_sorted(fuzzyConcept2Freq, false).get(0).getSecond());
        }

        HashMap<String, ReplicatedMeasurement> normLookup = buildReverseMap(normalized, (_rm) -> _rm.featureName);
        HashMap<String, Integer> fuzzyConcept2BgIdx = new HashMap<>();

        for(String fuzzyConcept : bg2features.keySet()) {
            Vector<ReplicatedMeasurement> features = map_and_filter(bg2features.get(fuzzyConcept), (_fn) -> normLookup.get(_fn), (_r) -> null != _r);
            if(features.size() == 0) {
                log.warn("%s: fuzzy concept: %s has no valid features mapped (%d requested)", replicateSetName, fuzzyConcept, bg2features.get(fuzzyConcept).size());
                continue;
            }

            log.info("%s fuzzyConcept: %s got %d features assigned", replicateSetName, fuzzyConcept, features.size());
            NumUtils.sort(features, (_f) -> _f.mean);
            backGroundContext.meanSignalValue2Error.add(NumUtils.mean(features, (_f) -> _f.mean));
            ErrorEstimationDistribution err = new ErrorEstimationDistribution(0, features.size(), (_i) -> features.get(_i).nonNanValues);
            fuzzyConcept2BgIdx.put(fuzzyConcept, backGroundContext.errorBackGrounds.size());
            backGroundContext.errorBackGrounds.add(err);
        }


        for(ReplicatedMeasurement r : normalized) {
            String fuzzyConcept = feature2fuzzyConcept.get(r.featureName);
            if(fuzzyConcept == null)
                throw new FRuntimeException("invalid fuzzy concept mapping feature: %s has not mapped fuzzy concept in replicate-set: %s", r.featureName, replicateSetName);

            backGroundContext.featureIdx2ErrorBackgroundIdx.put(r.featureIdx, fuzzyConcept2BgIdx.get(fuzzyConcept));
        }

        return backGroundContext;
    }
}

