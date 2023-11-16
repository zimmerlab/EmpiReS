package empires;

import java.util.HashMap;
import java.util.Vector;

public class BackGroundContext {
    public Vector<ErrorEstimationDistribution> errorBackGrounds = null;
    public Vector<Double> meanSignalValue2Error = null;
    public HashMap<Integer, Integer> featureIdx2ErrorBackgroundIdx = new HashMap<>();

}
