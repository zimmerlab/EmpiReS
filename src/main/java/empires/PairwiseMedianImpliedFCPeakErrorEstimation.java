package nlEmpiRe;

import lmu.utils.MapBuilder;

import java.util.HashMap;
import java.util.Vector;

import static lmu.utils.ObjectGetter.map;

public class PairwiseMedianImpliedFCPeakErrorEstimation {
    public double err = Double.POSITIVE_INFINITY;
    public Double shift = null;
    public ShiftedGroup sg1;
    public ShiftedGroup sg2;

    public static ErrorEstimationDistribution getShiftError(Vector<Double> v1, Vector<Double> v2) {
        assert(v1.size() == v2.size());
        HashMap<Integer, Integer> err2freq = new HashMap<>();
        int ntotal = 0;
        for(int i=0; i< v1.size(); i++)  {
            Double d1 = v1.get(i);
            Double d2 = v2.get(i);
            if(d1 == null || d2 == null || Double.isNaN(d1) || Double.isNaN(d2) || Double.isInfinite(d1) || Double.isInfinite(d2))
                continue;


            int key = (int)((d2 - d1) * ErrorEstimationDistribution.INT_FACTOR);
            MapBuilder.update(err2freq, key);
            ntotal++;

        }

        if(ntotal < 10)
            return null;

        ErrorEstimationDistribution errD = new ErrorEstimationDistribution(err2freq);

        errD.getCumulative(true);
        return errD;

    }

    public PairwiseMedianImpliedFCPeakErrorEstimation(ShiftedGroup sg1, ShiftedGroup sg2) {
        assert(sg1.size() == sg2.size());

        this.sg1 = sg1;
        this.sg2 = sg2;
        ErrorEstimationDistribution errD = getShiftError(sg1.average, sg2.average);
        if(errD == null)
            return;

        shift = -errD.getBestFCShift();
        err = errD.getSD(-shift);

    }

    @Override
    public String toString() {
        return String.format("%s vs %s err: %.2f shift: %.2f", sg1, sg2, err, shift);
    }
}
