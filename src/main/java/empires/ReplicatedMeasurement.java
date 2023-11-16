package empires;

import lmu.utils.NumUtils;

import java.util.Vector;

import static lmu.utils.ObjectGetter.*;

public class ReplicatedMeasurement {
    int featureIdx;
    public String featureName;
    public Vector<Double> replicates;
    public Vector<Double> nonNanValues = new Vector<>();
    public Double mean = Double.NaN;

    /** non-measured constructor */
    public ReplicatedMeasurement(int idx, int numReplicates) {
        featureIdx = idx;
        nonNanValues.add(mean = 0.0);
        this.replicates = mapIndex(numReplicates, (_i) -> Double.NaN);

    }

    public ReplicatedMeasurement(String name, int idx, Vector<Double> replicates) {
        this.featureIdx = idx;
        this.featureName = name;
        this.replicates = replicates;
        for(Double d : replicates) {
            if(Double.isNaN(d))
                continue;

            nonNanValues.add(d);

        }


        if(nonNanValues.size() > 1) {
            mean = NumUtils.getMedianMean(nonNanValues, true);
        }
    }


    public void applyPseudo(double pseudo) {

        if(pseudo == 0.0 || nonNanValues.size() == 0)
            return;

        replicates = map(replicates, (_d) -> (Double.isNaN(_d) ? _d : NumUtils.logN(Math.pow(2.0, _d) + pseudo, 2.0)));

        nonNanValues = map(nonNanValues, (_d) -> NumUtils.logN(Math.pow(2.0, _d) + pseudo, 2.0));
        mean = NumUtils.getMedianMean(nonNanValues);
    }

}
