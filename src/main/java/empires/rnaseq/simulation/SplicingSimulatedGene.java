package empires.rnaseq.simulation;

import lmu.utils.Tuple;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Vector;

public class SplicingSimulatedGene {

    public String geneId;
    public HashSet<String> trsWithCounts = new HashSet<>();
    public Vector<Vector<HashMap<Tuple, Double>>> condition2replicate2equivianceClasses = new Vector<>();


}
