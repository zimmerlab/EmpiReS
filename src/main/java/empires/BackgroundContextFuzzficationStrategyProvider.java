package nlEmpiRe;

import java.util.Vector;

public interface BackgroundContextFuzzficationStrategyProvider {
    public BackGroundContext getContexts(String replicateSetName, Vector<String> sampleNames, Vector<ReplicatedMeasurement> normalized);
}
