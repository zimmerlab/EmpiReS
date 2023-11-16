package empires;

import lmu.utils.LogConfig;
import lmu.utils.Logger;
import lmu.utils.NumUtils;
import lmu.utils.SetInfo;
import lmu.utils.fdr.BenjaminiHochberg;
import lmu.utils.tuple.Tuple3;

import java.util.*;
import java.util.function.Function;

import static lmu.utils.ObjectGetter.*;

/** improvement on the empiRe model = next level EmpiRe
 *
 *  improvements:
 *      (i) per condition empirical signal fitting
 *      (ii) context "telescoping" (merge contextes with similar error)
 *      (iii) quicker normalization, save signal to context infos if requested for further steps
 *      (iv) outlier detection for single measurments with respect to their variance versus the background context variance
 *          (does it eliminate the need of outlier correction for multiple submeasurements (e.g. multiple peptide, but fold change distorsion )

 */
public class EmpiRe {

    boolean compareFeaturesMeasuredOnlyOnOneSide = false;
    Vector<Double> per_comparison_pvals = new Vector<>();
    Vector<Double> per_comparison_fcdist_pvals = new Vector<>();
    boolean collect_comparison_pvals = false;


    public EmpiRe setCompareFeaturesMeasuredOnlyOneSide(boolean compareFeaturesMeasuredOnlyOnOneSide) {
        this.compareFeaturesMeasuredOnlyOnOneSide = compareFeaturesMeasuredOnlyOnOneSide;
        return this;
    }


    public EmpiRe setCollectComparisonPvals(boolean collect_comparison_pvals) {
        this.collect_comparison_pvals = collect_comparison_pvals;
        return this;
    }


    public Vector<Double> getPerComparisonPvals() {
        return per_comparison_pvals;
    }

    public Vector<Double> getPerComparisonFcDistPvals() {
        return per_comparison_fcdist_pvals;
    }

    public Vector<empires.DiffExpResult> getDifferentialResults(empires.NormalizedReplicateSet replicateSetFrom, empires.NormalizedReplicateSet replicateSetTo, Collection<String> featureNames) {

        return getDifferentialResults(replicateSetFrom, replicateSetTo, buildMap(featureNames, (_n) -> toVector(_n)));
    }

    /** tests all testable features */
    public Vector<empires.DiffExpResult> getDifferentialResults(empires.NormalizedReplicateSet replicateSetFrom, empires.NormalizedReplicateSet replicateSetTo) {
        Set<String> testable = (compareFeaturesMeasuredOnlyOnOneSide) ?
                    SetInfo.union(replicateSetFrom.inData.getFeatureNamesToTest(), replicateSetTo.inData.getFeatureNamesToTest())
                :   SetInfo.intersect(replicateSetFrom.inData.getFeatureNamesToTest(), replicateSetTo.inData.getFeatureNamesToTest());

        if(testable.size() == 0)
            return null;

        return getDifferentialResults(replicateSetFrom, replicateSetTo, testable,
                (_fn) -> {
                        Set<String> s1 = toSet(replicateSetFrom.inData.getSubFeatures(_fn));
                        Set<String> s2 = toSet(replicateSetTo.inData.getSubFeatures(_fn));

                        return (compareFeaturesMeasuredOnlyOnOneSide) ? SetInfo.union(s1, s2) : SetInfo.intersect(s1, s2);
                }
                );
    }

    public Vector<empires.DiffExpResult> getDifferentialResults(empires.NormalizedReplicateSet replicateSetFrom, empires.NormalizedReplicateSet replicateSetTo,
                                                                HashMap<String, Vector<String>> combinedFeature2SubFeatures
                                                               )
    {
        return getDifferentialResults(replicateSetFrom, replicateSetTo, combinedFeature2SubFeatures.keySet(), (_fn) -> combinedFeature2SubFeatures.get(_fn));
    }

    public Vector<empires.DiffExpResult> getDifferentialResults(empires.NormalizedReplicateSet replicateSetFrom, empires.NormalizedReplicateSet replicateSetTo,
                                                                Collection<String> featureNames, Function<String, Collection<String>> feature2subfeature) {

        Logger log = LogConfig.getLogger();

        //calculate shift
        Vector<String> totalFeatures = new Vector<>();
        apply(featureNames, (_fn) -> totalFeatures.addAll(feature2subfeature.apply(_fn)));



        DiffExpManager diffExpManager = new DiffExpManager(replicateSetFrom, replicateSetTo, totalFeatures, this);

        Vector<empires.DiffExpResult> results = new Vector<>();
        Vector<String> gotsubfeatures = filter(featureNames, (_f) -> feature2subfeature.apply(_f).size() > 0);
        log.info("will test diffexp %s - %s on %d/%d features. subfeature distrib: %s", replicateSetFrom.getInData().getReplicateSetName(),
                replicateSetTo.getInData().getReplicateSetName(),
                gotsubfeatures.size(), featureNames.size(),
                NumUtils.getNumInfo(gotsubfeatures, (_f) -> feature2subfeature.apply(_f).size()).getInfoWithQ());

        long t1 = System.currentTimeMillis();
        int ntested = 0;
        for(String combined : gotsubfeatures) {
            Collection<String> subfeatures = feature2subfeature.apply(combined);

            ntested++;
            //log.info("next to test: %s num subfeatures: %d", combined, subfeatures.size());
            empires.DiffExpResult diffExp = new DiffExpResult(diffExpManager, combined, subfeatures);
            if(ntested % 1000 == 0) {
                long t2 = System.currentTimeMillis();
                log.info("diffexp ready: %d/%d took %.4f sec on average", ntested, gotsubfeatures.size(), (t2 - t1) / (1000.0 * ntested));
            }
            if(diffExp.getNumUsedFeatures() == 0)
                continue;
            results.add(diffExp);
        }
        long t2 = System.currentTimeMillis();
        log.info("diffexp ready: %d/%d took %.4f sec on average", ntested, gotsubfeatures.size(), (t2 - t1) / (1000.0 * ntested));

        BenjaminiHochberg.adjust_pvalues(results, (_r) -> _r.pval, (_p) -> _p.getFirst().fdr = _p.getSecond());
        //BenjaminiHochberg.adjust_pvalues(results, (_r) -> _r.fcDistribPval, (_p) -> _p.getFirst().fcDistribFDR = _p.getSecond());
        BenjaminiHochberg.adjust_pvalues(results, (_r) -> _r.fcEstimatePval, (_p) -> _p.getFirst().fcEstimateFDR = _p.getSecond());

        return results;
    }


    public Vector<empires.DoubleDiffResult> getDoubleDiffResults(empires.NormalizedReplicateSet replicateSetFrom, empires.NormalizedReplicateSet replicateSetTo,
                                                                 Vector<Tuple3<String, Collection<String>, Collection<String>>> doubleDiffTests) {
        empires.DoubleDiffManager diffManager = getDoubleDiffManager(replicateSetFrom, replicateSetTo);

        Vector<DoubleDiffResult> results = new Vector<>();
        for(Tuple3<String, Collection<String>, Collection<String>> t : doubleDiffTests) {
            results.add(diffManager.getDoubleDiffResult(t.get0(), t.get1(), t.get2()));
        }
        BenjaminiHochberg.adjust_pvalues(results, (_r) -> _r.pval, (_p) -> _p.getFirst().fdr = _p.getSecond());
        return results;
    }


    public empires.DoubleDiffManager getDoubleDiffManager(empires.NormalizedReplicateSet replicateSetFrom, empires.NormalizedReplicateSet replicateSetTo) {
        HashSet<String> allfeatures = new HashSet<>();
        allfeatures.addAll(replicateSetFrom.getInData().getFeatureNames());
        allfeatures.addAll(replicateSetTo.getInData().getFeatureNames());

        return getDoubleDiffManager(replicateSetFrom, replicateSetTo, toVector(allfeatures));
    }

    public empires.DoubleDiffManager getDoubleDiffManager(empires.NormalizedReplicateSet replicateSetFrom, NormalizedReplicateSet replicateSetTo, Vector<String> allfeatures) {
        return new DoubleDiffManager(replicateSetFrom, replicateSetTo, new DiffExpManager(replicateSetFrom, replicateSetTo, allfeatures, this));
    }
}
