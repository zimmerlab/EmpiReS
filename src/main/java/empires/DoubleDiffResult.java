package nlEmpiRe;

import lmu.utils.NumUtils;
import lmu.utils.plotting.PlotCreator;

import java.awt.*;
import java.util.Collection;
import java.util.Vector;

import static lmu.utils.ObjectGetter.map;
import static lmu.utils.ObjectGetter.toVector;
import static nlEmpiRe.SparseCumulativeDistribution.QUANTILES;

public class DoubleDiffResult {

    /**
     * submeasurements along with their applied correction factor and corresponding background distribution
     */

    public String testName = "untested";


    public NormalizedReplicateSet normedReplicateSet1;
    public Collection<String> subFeatures1;

    public NormalizedReplicateSet normedReplicateSet2;
    public Collection<String> subFeatures2;


    public Collection<FeatureInfo> featureInfos1;
    public Collection<FeatureInfo> featureInfos2;

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

    public double pval = 1.0;
    public double fdr = 1.0;


    public DoubleDiffResult() {

    }

    DoubleDiffResult(String testName, Collection<String> subFeatures1, NormalizedReplicateSet normedReplicateSet1,
                     Collection<String> subFeatures2, NormalizedReplicateSet normedReplicateSet2, DiffExpManager diffExpManager) {

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
        if (!fillFeatures()) {
            return null;
        }
        PlotCreator.BoxPlotBuilder bpb = pc.buildBoxPlot();

        pc.setTitle("isoforms1: %s isoforms2: %s fcdiff: %.2f (mean: %.2f)", subFeatures1, subFeatures2, estimatedFC, meanFC);
        Color[] colors = new Color[] {Color.GREEN, Color.BLUE};

        for(int i=0; i<colors.length; i++) {

            Vector<String> src = NumUtils.sort(toVector(((i == 0) ? subFeatures1 : subFeatures2)), (_m) -> normedReplicateSet1.getNormed(_m).mean);

            for(String feature : src) {

                double m1 = normedReplicateSet1.getNormed(feature).mean;
                double m2 = normedReplicateSet2.getNormed(feature).mean;

                ErrorEstimationDistribution err = diffExpManager.getDiffError(feature);

                Vector<Double> quantiles = map(QUANTILES, (_q) -> m2  - m1 + err.getFoldChangeToCumulativeFrequency(_q));
                bpb.addBox(feature, quantiles, colors[i]);

            }



        }
        if (unpairedDiffFoldChangeDistribution != null) {
            bpb.addBox("diff", unpairedDiffFoldChangeDistribution.quantiles);
        }
        pc.setLabels("", "log2 diff fc", null);
        return bpb;


    }

    public PlotCreator.BoxPlotBuilder drawSignalDistributionBoxPlots(PlotCreator pc, boolean corrected) {
        PlotCreator.BoxPlotBuilder bpb = pc.buildBoxPlot();

        if (!fillFeatures()) {
            return null;
        }
        //pc.setTitle("isoforms1: %d isoforms2: %d fcdiff: %.2f (mean: %.2f)", subFeatures1.size(), subFeatures2.size(), estimatedFC, meanFC);
        pc.setTitle("isoforms1: %s isoforms2: %s fcdiff: %.2f (mean: %.2f)", subFeatures1, subFeatures2, estimatedFC, meanFC);


        Color[] colors = new Color[] {Color.GREEN, Color.BLUE};

        for(int i=0; i<colors.length; i++) {




            Vector<String> src = NumUtils.sort(toVector(((i == 0) ? subFeatures1 : subFeatures2)), (_m) -> normedReplicateSet1.getNormed(_m).mean);

            for(String feature : src) {

                for(int condI = 0; condI < 2; condI++) {

                    NormalizedReplicateSet base = (condI == 0) ? normedReplicateSet1 : normedReplicateSet2;


                    double fc = base.getNormed(feature).mean;
                    ErrorEstimationDistribution bg = base.getError(feature);
                    Vector<Double> quantiles = map(QUANTILES, (_q) -> fc + bg.getFoldChangeToCumulativeFrequency(_q));
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

    boolean fillFeatures() {
        if (subFeatures1 == null && featureInfos1 != null) {
            subFeatures1 = map(featureInfos1, (f) -> f.feature);
        }
        if (subFeatures2 == null && featureInfos2 != null) {
            subFeatures2 = map(featureInfos2, (f) -> f.feature);
        }
        System.out.printf("subfeatures: %s, %s\n", subFeatures1, subFeatures2);
        return subFeatures1 != null && subFeatures2 != null;
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
}


