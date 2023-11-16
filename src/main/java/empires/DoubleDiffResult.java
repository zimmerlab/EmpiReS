package empires;

import lmu.utils.NumUtils;
import lmu.utils.plotting.PlotCreator;

import java.awt.*;
import java.util.Collection;
import java.util.Vector;

import static lmu.utils.ObjectGetter.map;
import static lmu.utils.ObjectGetter.toVector;

public class DoubleDiffResult {

    /**
     * submeasurements along with their applied correction factor and corresponding background distribution
     */

    public String testName = "untested";


    public empires.NormalizedReplicateSet normedReplicateSet1;
    public Collection<String> subFeatures1;

    public empires.NormalizedReplicateSet normedReplicateSet2;
    public Collection<String> subFeatures2;


    public Collection<empires.FeatureInfo> featureInfos1;
    public Collection<empires.FeatureInfo> featureInfos2;

    public DoubleDiffSparseDistribution differenceFoldChangeDistribution;
    public DoubleDiffSparseDistribution correctedDifferenceFoldChangeDistribution;


    public SparseCumulativeDistribution featureSet1DiffDistribution;
    public SparseCumulativeDistribution featureSet2DiffDistribution;
    public SparseCumulativeDistribution unpairedDiffFoldChangeDistribution;

    public DiffExpManager diffExpManager;

    /**
     * estimated by the a fold-change window of size fcWindowForFc
     * -> the one with the highest probability mass is chosen
     */
    public Double estimatedFC = null;

    public double meanFC;
    public double f1meanFC;
    public double f2meanFC;

    public double pval = 1.0;
    public double fdr = 1.0;


    public DoubleDiffResult() {

    }

    DoubleDiffResult(String testName, Collection<String> subFeatures1, empires.NormalizedReplicateSet normedReplicateSet1,
                     Collection<String> subFeatures2, empires.NormalizedReplicateSet normedReplicateSet2, DiffExpManager diffExpManager) {

        this.testName = testName;
        this.subFeatures1 = subFeatures1;
        this.subFeatures2 = subFeatures2;
        this.normedReplicateSet1 = normedReplicateSet1;
        this.normedReplicateSet2 = normedReplicateSet2;
        this.diffExpManager = diffExpManager;

    }



    public String toString() {
        return String.format("%s fdr: %g fc: %.2f %s vs %s", testName, fdr, meanFC, subFeatures1, subFeatures2);
    }

    public void drawDistributions(PlotCreator pc) {
        drawDistributions(pc, false);
    }

    public PlotCreator.BoxPlotBuilder drawDiffFCPlots(PlotCreator pc, boolean corrected) {
        PlotCreator.BoxPlotBuilder bpb = pc.buildBoxPlot();

        pc.setTitle("isoforms1: %d isoforms2: %d fcdiff: %.2f (mean: %.2f)", subFeatures1.size(), subFeatures2.size(), estimatedFC, meanFC);
        Color[] colors = new Color[] {Color.GREEN, Color.BLUE};

        for(int i=0; i<colors.length; i++) {

            Vector<String> src = NumUtils.sort(toVector(((i == 0) ? subFeatures1 : subFeatures2)), (_m) -> normedReplicateSet1.getNormed(_m).mean);

            for(String feature : src) {

                double m1 = normedReplicateSet1.getNormed(feature).mean;
                double m2 = normedReplicateSet2.getNormed(feature).mean;

                empires.ErrorEstimationDistribution err = diffExpManager.getDiffError(feature);

                Vector<Double> quantiles = map(SparseCumulativeDistribution.QUANTILES, (_q) -> m2  - m1 + err.getFoldChangeToCumulativeFrequency(_q));
                bpb.addBox(feature, quantiles, colors[i]);

            }



        }
        bpb.addBox("diff", unpairedDiffFoldChangeDistribution.quantiles);
        pc.setLabels("", "log2 diff fc", null);
        return bpb;


    }

    public PlotCreator.BoxPlotBuilder drawSignalDistributionBoxPlots(PlotCreator pc, boolean corrected) {
        PlotCreator.BoxPlotBuilder bpb = pc.buildBoxPlot();

        pc.setTitle("isoforms1: %d isoforms2: %d fcdiff: %.2f (mean: %.2f)", subFeatures1.size(), subFeatures2.size(), estimatedFC, meanFC);


        Color[] colors = new Color[] {Color.GREEN, Color.BLUE};

        for(int i=0; i<colors.length; i++) {




            Vector<String> src = NumUtils.sort(toVector(((i == 0) ? subFeatures1 : subFeatures2)), (_m) -> normedReplicateSet1.getNormed(_m).mean);

            for(String feature : src) {

                for(int condI = 0; condI < 2; condI++) {

                    empires.NormalizedReplicateSet base = (condI == 0) ? normedReplicateSet1 : normedReplicateSet2;


                    double fc = base.getNormed(feature).mean;
                    ErrorEstimationDistribution bg = base.getError(feature);
                    Vector<Double> quantiles = map(SparseCumulativeDistribution.QUANTILES, (_q) -> fc + bg.getFoldChangeToCumulativeFrequency(_q));
                    Color c = colors[i];
                    if(condI == 0) {
                        c = c.brighter();
                    } else {
                        c = c.darker();
                    }
                    bpb.addBox(feature, quantiles, c);
                }
            }



        }
        pc.setLabels("", "log2 signal", null);
        return bpb;
    }

    public void drawDistributions(PlotCreator pc, boolean corrected) {

        if(featureSet1DiffDistribution != null)
            featureSet1DiffDistribution.drawLine(pc, "fs1");

        if(featureSet2DiffDistribution != null)
            featureSet2DiffDistribution.drawLine(pc, "fs2");

        if(unpairedDiffFoldChangeDistribution != null)
            unpairedDiffFoldChangeDistribution.drawLine(pc, "diff (fs1 -fs2)");

        pc.setTitle("isoforms1: %d isoforms2: %d fcdiff: %.2f (mean: %.2f)", subFeatures1.size(), subFeatures2.size(), estimatedFC, meanFC);

        pc.setLabels("log2(FC)", "frequency", "topright");
    }

    public String getTestName() {
        return testName;
    }

    public empires.NormalizedReplicateSet getNormedReplicateSet1() {
        return normedReplicateSet1;
    }

    public Collection<String> getSubFeatures1() {
        return subFeatures1;
    }

    public NormalizedReplicateSet getNormedReplicateSet2() {
        return normedReplicateSet2;
    }

    public Collection<String> getSubFeatures2() {
        return subFeatures2;
    }

    public Collection<empires.FeatureInfo> getFeatureInfos1() {
        return featureInfos1;
    }

    public Collection<FeatureInfo> getFeatureInfos2() {
        return featureInfos2;
    }

    public DoubleDiffSparseDistribution getDifferenceFoldChangeDistribution() {
        return differenceFoldChangeDistribution;
    }

    public DoubleDiffSparseDistribution getCorrectedDifferenceFoldChangeDistribution() {
        return correctedDifferenceFoldChangeDistribution;
    }

    public SparseCumulativeDistribution getFeatureSet1DiffDistribution() {
        return featureSet1DiffDistribution;
    }

    public SparseCumulativeDistribution getFeatureSet2DiffDistribution() {
        return featureSet2DiffDistribution;
    }

    public SparseCumulativeDistribution getUnpairedDiffFoldChangeDistribution() {
        return unpairedDiffFoldChangeDistribution;
    }

    public DiffExpManager getDiffExpManager() {
        return diffExpManager;
    }

    public Double getEstimatedFC() {
        return estimatedFC;
    }

    public double getMeanFC() {
        return meanFC;
    }

    public double getF1meanFC() {
        return f1meanFC;
    }

    public double getF2meanFC() {
        return f2meanFC;
    }

    public double getPval() {
        return pval;
    }

    public double getFdr() {
        return fdr;
    }
}


