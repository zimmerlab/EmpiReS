package empires.input;

import lmu.utils.UPair;

import java.util.*;
import java.util.function.BiFunction;

import static lmu.utils.IteratorUtils.rangev;
import static lmu.utils.ObjectGetter.*;

public class ReplicateSetInfo {

    String replicateSetName;

    Vector<String> replicateNames;
    Vector<String> featureNames;

    Vector<Vector<Double>> log2replicate2FeatureData;

    Set<String> featureNamestoTest;
    HashMap<String, Vector<String>> combinedFeatures;

    HashMap<String, Integer> featurename2Idx = null;


    private ReplicateSetInfo(ReplicateSetInfo rsi) {
        this.replicateNames = rsi.replicateNames;
        this.replicateSetName = rsi.replicateSetName;
        this.featureNames = rsi.featureNames;
        this.featureNamestoTest = rsi.featureNamestoTest;
        this.combinedFeatures = rsi.combinedFeatures;
        this.log2replicate2FeatureData = rsi.log2replicate2FeatureData;
    }

    public ReplicateSetInfo(String replicateSetName, Vector<String> replicateNames, Vector<String> featureNames) {
        this.replicateSetName = replicateSetName;
        this.replicateNames = replicateNames;
        this.featureNames = featureNames;

        log2replicate2FeatureData = new Vector<>();
        for(String rep : replicateNames) {
            log2replicate2FeatureData.add(map(featureNames, (_n) -> Double.NaN));
        }
    }


    public ReplicateSetInfo getReplicateSubSet(Vector<Integer> replicates) {
        ReplicateSetInfo rv = new ReplicateSetInfo(this);
        rv.replicateNames = map(replicates, (_i) -> replicateNames.get(_i));
        rv.log2replicate2FeatureData = map(replicates, (_i) -> log2replicate2FeatureData.get(_i));
        return rv;
    }

    public ReplicateSetInfo swapFeatures(Vector<UPair<String>> swaps) {
        ReplicateSetInfo rsi = new ReplicateSetInfo(this);
        rsi.featureNames = cloneVector(featureNames);
        HashMap<String, Integer> lookup = buildIndexMap(rsi.featureNames);
        for(UPair<String> swap : swaps) {
            int idx1 = lookup.get(swap.getFirst());
            int idx2 = lookup.get(swap.getSecond());
            rsi.featureNames.set(idx1, swap.getSecond());
            rsi.featureNames.set(idx2, swap.getFirst());
        }
        rsi.featurename2Idx = null;
        return rsi;
    }

    public int getFeatureIdx(String name) {
        if(featurename2Idx == null) {
            featurename2Idx = buildIndexMap(featureNames);
        }
        return featurename2Idx.get(name);
    }

    public Vector<Double> getReplicateData(String fn) {
        int featureIdx = getFeatureIdx(fn);
        return mapIndex(replicateNames.size(), (_rIdx) -> log2replicate2FeatureData.get(_rIdx).get(featureIdx));
    }
    public HashSet<String> filterFeatureds(BiFunction<String, Vector<Double>, Boolean> excluder) {
        HashSet<String> rv = new HashSet<>();
        Vector<Integer> repRange = rangev(log2replicate2FeatureData.size());
        for(int i=0; i<featureNames.size(); i++) {
            String fn = featureNames.get(i);
            Vector<Double> featureValues = getReplicateData(fn);
            if(excluder.apply(fn, featureValues))
                continue;


            rv.add(fn);
        }
        return rv;
    }

    public ReplicateSetInfo restrictToSubFeatures(Set<String> featureNames) {
        ReplicateSetInfo rsi = new ReplicateSetInfo(this);

        Vector<Integer> featureIndeces = toSortedVector(map(featureNames, (_f) -> getFeatureIdx(_f)), true);
        rsi.featureNames = map(featureIndeces, (_i) -> this.featureNames.get(_i));
        rsi.log2replicate2FeatureData = new Vector<>();
        for(int ri=0; ri<getNumReplicates(); ri++) {

            Vector<Double> v = log2replicate2FeatureData.get(ri);
            rsi.log2replicate2FeatureData.add(map(featureIndeces, (_i) -> v.get(_i)));
        }
        rsi.featurename2Idx = null;


        if(combinedFeatures != null) {
            rsi.combinedFeatures = new HashMap<>();
            for(Map.Entry<String, Vector<String>> e : combinedFeatures.entrySet()) {
                Vector<String> filtered = filter(e.getValue(), (_fn) -> featureNames.contains(_fn));
                if(filtered.size() == 0)
                    continue;

                rsi.combinedFeatures.put(e.getKey(), filtered);
            }
        } else {
            rsi.featureNamestoTest = featureNames;
        }

        return rsi;
    }

    public void setCombinedFeatures(HashMap<String, Vector<String>> combinedFeatures) {
        this.combinedFeatures = combinedFeatures;
    }

    public void setLog2Data(int replicateIdx, Vector<Double> data) {
        log2replicate2FeatureData.set(replicateIdx, data);

    }
    public Vector<String> getReplicateNames() {
        return replicateNames;
    }

    public String getReplicateSetName() {
        return replicateSetName;
    }

    public String getFeatureName(int idx) {
        return featureNames.get(idx);
    }


    public int getNumFeatures() {
        return featureNames.size();
    }

    public int getNumReplicates() {
        return replicateNames.size();
    }

    public Set<String> getFeatureNamesToTest() {
        if(featureNamestoTest != null)
            return featureNamestoTest;

        if(combinedFeatures != null) {
            return featureNamestoTest = combinedFeatures.keySet();
        }
        return featureNamestoTest = toSet(featureNames);
    }

    public Vector<String> getSubFeatures(String combinedName) {
        if(combinedFeatures != null)
            return combinedFeatures.get(combinedName);

        return toVector(combinedName);
    }
    public Vector<String> getFeatureNames() {
        return featureNames;
    }

    public Vector<Vector<Double>> getLog2Data() {
        return log2replicate2FeatureData;
    }



}
