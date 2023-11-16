package empires;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.Vector;

import static lmu.utils.ObjectGetter.map;

public class ShiftedGroup {
    int idx;
    HashMap<Integer, Double> idx2shift = new HashMap<>();
    Vector<Double> average;

    public ShiftedGroup(int idx, Vector<Double> data) {
        this.idx = idx;
        idx2shift.put(idx, 0.0);
        average = data;
    }

    public ShiftedGroup(int idx, PairwiseMedianImpliedFCPeakErrorEstimation errorEstimation) {
        this.idx = idx;

        idx2shift.putAll(errorEstimation.sg1.idx2shift);

        System.out.printf("err sg1: %s sg2: %s est shift: %s\n", errorEstimation.sg1.idx2shift, errorEstimation.sg2.idx2shift, errorEstimation.shift);
        for(Map.Entry<Integer, Double> e: errorEstimation.sg2.idx2shift.entrySet()) {
            idx2shift.put(e.getKey(), e.getValue() + errorEstimation.shift);
        }

        final int L = errorEstimation.sg1.average.size();

        average = new Vector<>();
        for(int i=0; i<L; i++) {
            Double d1 = errorEstimation.sg1.average.get(i);
            Double d2 = errorEstimation.sg2.average.get(i);
            if(d1 == null || d2 == null || Double.isNaN(d1) || Double.isNaN(d2) || Double.isInfinite(d1) || Double.isInfinite(d2)) {
                average.add(null);
                continue;
            }


            d2 += errorEstimation.shift;
            double avg = (d1 * errorEstimation.sg1.getNumSamples() + d2 * errorEstimation.sg2.getNumSamples()) / idx2shift.size();
            average.add(avg);
        }
    }

    public Set<Integer> getMembers() {
        return idx2shift.keySet();
    }

    public Vector<Vector<Double>> getShiftedValues(Vector<Vector<Double>> invalues) {
        Vector<Vector<Double>> rv = new Vector<>();
        for(Map.Entry<Integer, Double> e : idx2shift.entrySet()) {
            rv.add(map(invalues.get(e.getKey()), (_d) -> _d + e.getValue()));
        }
        return rv;
    }

    public Vector<Vector<Double>> getRawValues(Vector<Vector<Double>> invalues) {
        Vector<Vector<Double>> rv = new Vector<>();
        for(Map.Entry<Integer, Double> e : idx2shift.entrySet()) {
            rv.add(invalues.get(e.getKey()));
        }
        return rv;
    }
    public int getNumSamples() {
        return idx2shift.size();
    }

    public int size() {
        return average.size();
    }

    @Override
    public String toString() {
        return String.format("sg[%d] %s", idx, map(idx2shift.entrySet(), (_e) -> String.format("%d:%.2f", _e.getKey(), _e.getValue())) );
    }
}
