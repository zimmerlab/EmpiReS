package empires;

import lmu.utils.UPair;

import java.util.Collection;
import java.util.Vector;

import static lmu.utils.ObjectGetter.map;

public class DoubleDiffSparseDistribution {

    public final static double[] PERCENT_DISTRIBUTION_COVERED_TESTS = {0.01, 0.2, 0.03, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4};
    public Vector<UPair<Double>> jointOverlapInfos = new Vector<>();

    public empires.SparseCumulativeDistribution subSet1FoldChangeDistribution;
    public empires.SparseCumulativeDistribution subSet2FoldChangeDistribution;

    public empires.SparseCumulativeDistribution differenceFoldChangeDistribution;


    public DoubleDiffSparseDistribution(empires.ErrorEstimationDistribution joint1, empires.ErrorEstimationDistribution joint2) {
        empires.ErrorEstimationDistribution diffJoint = empires.ErrorEstimationDistribution.substract(joint1, joint2);

        subSet1FoldChangeDistribution = new empires.SparseCumulativeDistribution(joint1);
        subSet2FoldChangeDistribution = new empires.SparseCumulativeDistribution(joint2);

        differenceFoldChangeDistribution = new SparseCumulativeDistribution(diffJoint);

        calcOverlaps(joint1, joint2);
    }
    public DoubleDiffSparseDistribution(Collection<String> set1, empires.NormalizedReplicateSet replicateSet1,
                                        Collection<String> set2, NormalizedReplicateSet replicateSet2) {

        this(new empires.ErrorEstimationDistribution(map(set1, (_m) -> replicateSet1.getError(_m)), map(set1, (_m) -> replicateSet1.getNormed(_m).mean)),
            new empires.ErrorEstimationDistribution(map(set2, (_m) -> replicateSet2.getError(_m)), map(set1, (_m) -> replicateSet2.getNormed(_m).mean))
        );

    }

    /** get the distribution of percent overlaps of the joint distributions */
    void calcOverlaps(empires.ErrorEstimationDistribution joint1, ErrorEstimationDistribution joint2) {
        /*

                |-----iv1 --------|
                    |------- iv2 -----------|
         */
        for(double percent_of_distribution_covered: PERCENT_DISTRIBUTION_COVERED_TESTS) {
            double fc11 = joint1.getFoldChangeToCumulativeFrequency(percent_of_distribution_covered);
            double fc12 = joint1.getFoldChangeToCumulativeFrequency(1.0 - percent_of_distribution_covered);


            double fc21 = joint2.getFoldChangeToCumulativeFrequency(percent_of_distribution_covered);
            double fc22 = joint2.getFoldChangeToCumulativeFrequency(1.0 - percent_of_distribution_covered);

            double smallerRegion = Math.min(fc22 - fc21, fc12 - fc11);
            double MAXMIN = Math.max(fc11, fc21);
            double MINMAX = Math.min(fc12, fc22);

            double overlap = 0.0;

            //(MINMAX <= MAXMIN) ? 0.0 : (MINMAX - MAXMIN) / smallerRegion;
            if(MINMAX > MAXMIN) {
                double probabilityMass1 = joint1.getCumulativeFrequencyToFoldChange(MINMAX) - joint1.getCumulativeFrequencyToFoldChange(MAXMIN);
                double probabilityMass2 = joint2.getCumulativeFrequencyToFoldChange(MINMAX) - joint2.getCumulativeFrequencyToFoldChange(MAXMIN);
                overlap = Math.min(probabilityMass1, probabilityMass2);
            }
            jointOverlapInfos.add(UPair.createU(percent_of_distribution_covered, overlap));


        }
    }
}

