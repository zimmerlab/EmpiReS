package empires.test;


import empires.EmpiRe;
import lmu.utils.UPair;

import java.util.HashMap;
import java.util.Vector;

import static lmu.utils.ObjectGetter.shuffle;
import static lmu.utils.ObjectGetter.toVector;

public class UniformTest
{
    MultiSDSimulation cond1;
    MultiSDSimulation cond2;
    static empires.BackgroundContextFuzzficationStrategyProvider provider = new empires.AutoBackGroundContextProvider();

    empires.NormalizedReplicateSet nrs1;
    empires.NormalizedReplicateSet nrs2;
    empires.EmpiRe empiRe = new EmpiRe();
    private Vector<empires.DiffExpResult> results;
    boolean simulateNans = false;
    Integer numSamples1 = null;
    Integer numSamples2 = null;

    int numFeaturePerCombined = -1;

    int numFeatures = 10_000;
    Double SDFactor = 1.0;

    public HashMap<String, Vector<String>> combinedFeatures = new HashMap<>();

    Vector<empires.DoubleDiffResult> doubleDiffResults;
    Vector<UPair<Vector<String>>> doubleDiffFeatureCombis = null;

    Vector<String> combFeatureVec;
    boolean plotSignals = false;

    public UniformTest() {
        this(null, null, false, 10_000);
    }
    public UniformTest(Integer numSamples1, Integer numSamples2, boolean simulateNans, int numFeatures) {
        this.numFeatures = numFeatures;
        this.simulateNans = simulateNans;
        this.numSamples1 = numSamples1;
        this.numSamples2 = numSamples2;

    }

    public void setNumFeatures(int  numFeatures) {
        this.numFeatures = numFeatures;
    }

    public void setSDFactor(double sdFactor) {
        SDFactor = sdFactor;
    }

    public void setPlotSignals(boolean b) {
        plotSignals = b;
    }

    public void simulate()  {
        empiRe.setCollectComparisonPvals(true);
        cond1 = new MultiSDSimulation(numSamples1, plotSignals, null, simulateNans, numFeatures);
        cond2 = new MultiSDSimulation(numSamples2, plotSignals, cond1.signals, simulateNans, numFeatures);

        cond1.setNumPoints(numFeatures);
        cond2.setNumPoints(numFeatures);

        if(SDFactor != null) {
            cond1.setSDFactor(1.0);
            cond2.setSDFactor(SDFactor);
        }
        cond1.simulate();
        cond2.simulate();
        nrs1 = new empires.NormalizedReplicateSet(cond1.rsi, provider);
        nrs2 = new empires.NormalizedReplicateSet(cond2.rsi, provider);

        Vector<String> fnames = shuffle(nrs1.getInData().getFeatureNames(), true);

        int idx = 0;
        while(idx < fnames.size()) {
            int numSubFeatures = (numFeaturePerCombined > 0) ? numFeaturePerCombined  : 1 + (int)(Math.random() * 10);
            if(idx + numSubFeatures > fnames.size())
                break;

            Vector<String> subvec = new Vector<>();
            for(int j=0; j<numSubFeatures; j++) {
                subvec.add(fnames.get(idx++));
            }
            combinedFeatures.put(String.format("comb%04d", combinedFeatures.size()), subvec);
        }

        combFeatureVec = toVector(combinedFeatures.keySet());
    }
    public void reset() {
        results = null;
    }

    Vector<UPair<Vector<String>>> getDoubleDiffFeatureCombis() {
        if(doubleDiffFeatureCombis != null)
            return doubleDiffFeatureCombis;

        Vector<String> featureNames = nrs1.getInData().getFeatureNames();
        doubleDiffFeatureCombis = new Vector<>();
        int idx = 0;
        while(idx < featureNames.size()) {
            int nfeat1 = 1 + (int)(Math.random() * 10);
            int nfeat2 = 1 + (int)(Math.random() * 10);

            Vector<String> fvec1 = new Vector<>();
            Vector<String> fvec2 = new Vector<>();
            int ntotal = nfeat1 + nfeat2;
            for(int i=0; i<ntotal && idx < featureNames.size(); i++) {
                Vector<String> target = (i < nfeat1) ? fvec1 : fvec2;
                target.add(featureNames.get(idx++));
            }
            if(fvec1.size() == nfeat1 && fvec2.size() == nfeat2) {
                doubleDiffFeatureCombis.add(UPair.createU(fvec1, fvec2));
            }
        }

        return doubleDiffFeatureCombis;
    }

    public Vector<empires.DoubleDiffResult> getDoubleDiffResults() {
        if(doubleDiffResults != null)
            return doubleDiffResults;

        empires.DoubleDiffManager ddm = empiRe.getDoubleDiffManager(nrs1, nrs2);

        doubleDiffResults = new Vector<>();

        int testidx = 0;
        for(UPair<Vector<String>> fcombi : getDoubleDiffFeatureCombis()) {
            testidx++;
            doubleDiffResults.add(ddm.getDoubleDiffResult("test" + testidx, fcombi.getFirst(), fcombi.getSecond()));
        }
        return doubleDiffResults;
    }

    Vector<empires.DiffExpResult> getResults() {
        if(results != null) {
            return results;
        }
        return results = empiRe.getDifferentialResults(nrs1, nrs2, combinedFeatures);
    }





}
