package nlEmpiRe;

import lmu.utils.*;
import nlEmpiRe.input.ReplicateSetInfo;

import java.io.File;
import java.io.PrintWriter;
import java.util.*;

import static lmu.utils.IteratorUtils.rangev;
import static lmu.utils.ObjectGetter.*;
import org.apache.logging.log4j.Logger;

/** describes a set of replicate measurement referring to the same phenotype (a.k.a "a condition")
 *
 *
 */
public class NormalizedReplicateSet {


    Logger log = LogConfig.getLogger();
    ReplicateSetInfo inData;
    Vector<Double> shifts;
    HashMap<Integer, ReplicatedMeasurement> featureIdx2Normalized;
    Vector<ErrorEstimationDistribution> errorBackGrounds = null;
    Vector<Double> meanSignalValue2Error = null;
    HashMap<Integer, Integer> featureIdx2ErrorBackgroundIdx = new HashMap<>();
    HashMap<String, Integer> featurename2idx = null;
    Normalization normalization;

    ErrorEstimationDistribution nonmeasuredDefaultError;
    ReplicatedMeasurement nonMeasuredFeature;


    public NormalizedReplicateSet(ReplicateSetInfo ci) {
        this(ci, new AutoBackGroundContextProvider());
    }
    public NormalizedReplicateSet(ReplicateSetInfo ci, BackgroundContextFuzzficationStrategyProvider contextFuzzficationStrategyProvider) {
        this(ci, contextFuzzficationStrategyProvider, false);
    }
    public NormalizedReplicateSet(ReplicateSetInfo ci, BackgroundContextFuzzficationStrategyProvider contextFuzzficationStrategyProvider, boolean noErrorIfSingleReplicate) {
        inData = ci;
        nonMeasuredFeature = new ReplicatedMeasurement(-1, ci.getNumReplicates());
        try {

            normalization = new Normalization(ci.getLog2Data());

        } catch(Normalization.NonClustarableException nce) {
            BufferPrintWriter bpw = new BufferPrintWriter();
            Vector<String> repnames = ci.getReplicateNames();
            bpw.printf("%s: cannot normalize samples - too many non-overlapping features clustered %d samples into %d groups:\n", ci.getReplicateSetName(),
                    repnames.size(), nce.clusters.size());

            int cidx = 0;
            for(int root : NumUtils.sort(toVector(nce.clusters.keySet()), (_k) -> nce.clusters.get(_k).size(), true)) {
                cidx++;
                bpw.printf("\tcluster: %d size: %d replicates: %s\n", cidx, nce.clusters.get(root).size(),
                        map(nce.clusters.get(root), (_i) -> repnames.get(_i)));
            }
            throw new FRuntimeException(bpw.getBuffer());

        }


        shifts = normalization.getShifts();
        featurename2idx = buildIndexMap(ci.getFeatureNames());

        Vector<Integer> range = rangev(shifts);


        final Vector<Vector<Double>> log2Data = ci.getLog2Data();
        Vector<ReplicatedMeasurement> nanValues = new Vector<>();
        Vector<ReplicatedMeasurement> normalized = new Vector<>();

        for (int i = 0; i < ci.getNumFeatures(); i++) {

            int IDX = i;
            Vector<Double> indata = map(log2Data, (_vd) -> _vd.get(IDX));
            ReplicatedMeasurement r = new ReplicatedMeasurement(ci.getFeatureName(i), i, map(range, (_i) -> indata.get(_i) + shifts.get(_i)));
            if(Double.isNaN(r.mean)) {
                nanValues.add(r);
                continue;
            }
            normalized.add(r);
        }

        if(noErrorIfSingleReplicate && inData.getNumReplicates() < 2) {
            log.warn("unusable for background estimation: %s", inData.getReplicateSetName());
            return;
        }
        BackGroundContext backGroundContext = contextFuzzficationStrategyProvider.getContexts(inData.getReplicateSetName(), inData.getReplicateNames(), normalized);
        featureIdx2ErrorBackgroundIdx = backGroundContext.featureIdx2ErrorBackgroundIdx;
        errorBackGrounds = backGroundContext.errorBackGrounds;
        meanSignalValue2Error = backGroundContext.meanSignalValue2Error;

        if(errorBackGrounds.size() > 0) {
            int maxSDErrIdx = NumUtils.maxObj(rangev(errorBackGrounds), (_idx) -> errorBackGrounds.get(_idx).getSD()).getSecond();

            nonmeasuredDefaultError = errorBackGrounds.get(maxSDErrIdx);

            log.info("reassign %d/%d all-zero measurements to 0.0 max error: %.2f", nanValues.size(), ci.getNumFeatures(), nonmeasuredDefaultError.getSD());

        }

        featureIdx2Normalized = buildReverseMap(normalized, (_n) -> _n.featureIdx);

    }

    public Normalization getNormalization() {
        return normalization;
    }

    public void applyPseudo(double pseudoValue) {
        apply(featureIdx2Normalized.values(), (_n) -> _n.applyPseudo(pseudoValue));
    }
    public void correctBackWards() {
        correctBackWards(0.5);
    }
    public void correctBackWards(double percentile) {
        int selected = (int) (percentile * errorBackGrounds.size());
        ErrorEstimationDistribution ref = errorBackGrounds.get(selected);


        for(int i=selected; i>=0; i--) {
            ErrorEstimationDistribution e = errorBackGrounds.get(i);
            if(e.getSD() > ref.getSD()) {
                ref = e;
            }
            errorBackGrounds.set(i, ref);

            System.out.printf("DEF: %d update: %.2f with %.2f\n", i, e.getSD(), ref.getSD());
        }
    }
    public Vector<Double> getMedianValues(Vector<String> features) {
        Vector<Double> medians = new Vector<>();
        for(String f : features) {
            Integer fidx = featurename2idx.get(f);
            if(fidx == null) {
                medians.add(Double.NaN);
                continue;
            }
            ReplicatedMeasurement rmes = featureIdx2Normalized.get(fidx);
            if(rmes == null) {
                medians.add(Double.NaN);
                continue;
            }
            medians.add(NumUtils.median(rmes.nonNanValues));
        }
        return medians;
    }

    public Vector<Double> getMeanValues(Vector<String> features) {
        Vector<Double> means = new Vector<>();
        for(String f : features) {
            Integer fidx = featurename2idx.get(f);
            if(fidx == null) {
                means.add(Double.NaN);
                continue;
            }
            ReplicatedMeasurement rmes = featureIdx2Normalized.get(fidx);
            if(rmes == null) {
                means.add(Double.NaN);
                continue;
            }
            means.add(rmes.mean);
        }
        return means;
    }
    public ReplicateSetInfo getInData() {
        return inData;
    }

    public boolean gotFeatureMeasured(String feature) {
        Integer fidx = featurename2idx.get(feature);
        return (fidx != null && featureIdx2ErrorBackgroundIdx.get(fidx) != null);
    }

    public ReplicatedMeasurement getNormed(String feature) {
        Integer featureIdx = featurename2idx.get(feature);

        if(featureIdx == null) {
            return nonMeasuredFeature;
        }
        return featureIdx2Normalized.getOrDefault(featureIdx, nonMeasuredFeature);
    }

    public ErrorEstimationDistribution getErrorByMeanSignalValue(double signalValue) {
        int hitidx = Collections.binarySearch(meanSignalValue2Error, signalValue);
        if(hitidx < 0) {
            hitidx = -hitidx - 1;
        }
        hitidx = Math.min(errorBackGrounds.size() - 1, Math.max(hitidx, 0));
        return errorBackGrounds.get(hitidx);
    }

    public ErrorEstimationDistribution getError(String feature) {
        Integer featureIdx = featurename2idx.get(feature);

        if(featureIdx == null) {
            //non measured
            return nonmeasuredDefaultError;
        }
        Integer bgIdx = featureIdx2ErrorBackgroundIdx.get(featureIdx);
        if(bgIdx == null)
            return nonmeasuredDefaultError;

        return errorBackGrounds.get(bgIdx);
    }

    public int getNumReplicates() {
        return inData.getNumReplicates();
    }

    public int getNumRegularBackgrounds(){
        return meanSignalValue2Error.size();
    }

    public ErrorEstimationDistribution getErrorByIdx(int idx) {
        return errorBackGrounds.get(idx);
    }

    public double getMeanSignalErrorForErrorIdx(int idx) {
        return Math.pow(2.0, meanSignalValue2Error.get(idx));
    }


    public void writeErrorDistribInfos(File outdir) {
        NamedFoldChangeDistributionMap nfd = new NamedFoldChangeDistributionMap();
        for(int i=0; i<errorBackGrounds.size(); i++) {
            nfd.addDistrib(""+i, errorBackGrounds.get(i));
        }
        nfd.serialize(new File(outdir, inData.getReplicateSetName()+".errordistribs"));

        PrintWriter pw = FileUtils.getWriter(new File(outdir, inData.getReplicateSetName()+".feature2error.mapping"));
        for(Map.Entry<String, Integer> e : featurename2idx.entrySet()) {
            Integer mapped = featureIdx2ErrorBackgroundIdx.get(e.getValue());
            if(mapped == null)
                continue;
            pw.printf("%s\t%d\n", e.getKey(), mapped);
        }
        pw.close();
    }

}
