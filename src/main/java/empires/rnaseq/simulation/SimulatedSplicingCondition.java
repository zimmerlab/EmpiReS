package nlEmpiRe.rnaseq.simulation;


import nlEmpiRe.rnaseq.MultiIsoformRegion;

import java.util.HashMap;
import java.util.Vector;

public class SimulatedSplicingCondition {

    public String condition;
    public Vector<HashMap<MultiIsoformRegion, TranscriptSimulation>> replicates = new Vector<>();

    public SimulatedSplicingCondition(String condition) {
        this.condition = condition;
    }
}
