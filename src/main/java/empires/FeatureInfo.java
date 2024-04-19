package nlEmpiRe;

import lmu.utils.NumUtils;

import java.util.Vector;

import static lmu.utils.IteratorUtils.rangev;
import static lmu.utils.ObjectGetter.*;

public class FeatureInfo {
    public String feature;
    public Vector<Double> cond1;
    public Vector<Double> cond2;


    Vector<Double> cond1_logvals;
    Vector<Double> cond2_logvals;

    public ErrorEstimationDistribution err1;
    public ErrorEstimationDistribution err2;

    double min1 = 0;
    double min2 = 0;
    double overallmin;
    FeatureInfo(String feature) {
        this.feature = feature;
    }
    public FeatureInfo(String feature, NormalizedReplicateSet rs1, NormalizedReplicateSet rs2) {
        this(feature);
        this.cond1_logvals = rs1.getInData().getReplicateData(feature);
        this.cond2_logvals = rs2.getInData().getReplicateData(feature);
        cond1 = map(cond1_logvals, (_x) -> (Double.isNaN(_x)) ? 0.0 : Math.pow(2.0, _x));
        cond2 = map(cond2_logvals, (_x) -> (Double.isNaN(_x)) ? 0.0 : Math.pow(2.0, _x));

        this.cond1_logvals = rs1.getNormed(feature).replicates;
        this.cond2_logvals = rs2.getNormed(feature).replicates;

        err1 = rs1.getError(feature);
        err2 = rs2.getError(feature);
        setMinInfos();
    }

    static Vector<Double> toLogVals(Vector<Double> vals){
        return map(vals, (_v) -> (_v == 0.0) ? Double.NaN : NumUtils.logN(_v, 2.0));
    }
    public Vector<Double> getCond1LogVals() {
        return (cond1_logvals != null) ? cond1_logvals : (cond1_logvals = toLogVals(cond1));
    }

    public Vector<Double> getCond2LogVals() {
        return (cond2_logvals != null) ? cond2_logvals : (cond2_logvals = toLogVals(cond2));
    }

    public double getMeanLogSignalCond1() {
        return NumUtils.mean(filter(getCond1LogVals(), (_d) -> !Double.isNaN(_d)));
    }

    public double getMeanLogSignalCond2() {
        return NumUtils.mean(filter(getCond2LogVals(), (_d) -> !Double.isNaN(_d)));
    }

    void setMinInfos() {
        min1 = NumUtils.min(cond1);
        min2 = NumUtils.min(cond2);
        overallmin = Math.max(min1, min2);

    }

    void addPseudo(double pseudo) {
        cond1 = map(cond1, (_d) -> (_d <= 0.01) ? 0.0 : _d + pseudo);
        cond2 = map(cond2, (_d) -> (_d <= 0.01) ? 0.0 : _d + pseudo);
        //cond1 = map(cond1, (_d) -> (_d < pseudo) ? 0.0 : _d );
        //cond2 = map(cond2, (_d) -> (_d < pseudo) ? 0.0 : _d );
    }

    public FeatureInfo combine(FeatureInfo featureInfo) {
        FeatureInfo combined = new FeatureInfo(feature+"."+featureInfo.feature);
        combined.cond1 = mapIndex(rangev(cond1.size()), (_i) -> cond1.get(_i) + featureInfo.cond1.get(_i));
        combined.cond2 = mapIndex(rangev(cond2.size()), (_i) -> cond2.get(_i) + featureInfo.cond2.get(_i));
        combined.setMinInfos();
        return combined;
    }


    public static Vector<FeatureInfo> merge(Vector<String> features, NormalizedReplicateSet rs1, NormalizedReplicateSet rs2, double minTreshold, Double pseudo) {
        Vector<FeatureInfo> featureInfos = NumUtils.sort(map(features, (_f) -> new FeatureInfo(_f, rs1, rs2)), (_f) -> _f.overallmin);
        Vector<FeatureInfo> merged = new Vector<>();
        FeatureInfo combined = null;
        Vector<FeatureInfo> combineList = new Vector<>();
        for(FeatureInfo f : featureInfos) {
            if(f.overallmin < minTreshold) {
                combineList.add(f);
                combined = (combined == null) ? f : combined.combine(f);
                continue;
            }
            merged.add(f);
        }
        if(combined != null) {
            if(pseudo == null &&  combined.overallmin  < minTreshold)
                return merged;
            //get  the corresponding error distribs
            if(pseudo != null) {
                combined.addPseudo(pseudo);
            }


            combined.err1 = NumUtils.minObj(combineList, (_f) -> _f.err1.getSD()).getSecond().err1;// cors1.getErrorByMeanSignalValue(combined.getMeanLogSignalCond1());
            //System.out.printf("input sds: %s merged so far: %d  new err1: %.2f\n", map(featureInfos, (_f) -> _f.err1.getSD()), merged.size(), combined.err1.getSD());
            combined.err2 = NumUtils.minObj(combineList, (_f) -> _f.err2.getSD()).getSecond().err2;//
            // combined.err2 = rs1.getErrorByMeanSignalValue(combined.getMeanLogSignalCond2());
            merged.add(combined);
        }
        return merged;
    }


}
