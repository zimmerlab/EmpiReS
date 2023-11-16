package empires.plotting;

import empires.PairScore;
import lmu.utils.ImageUtils;
import lmu.utils.UPair;
import lmu.utils.plotting.PlotCreator;

import java.awt.image.BufferedImage;
import java.util.Vector;

import static lmu.utils.ObjectGetter.toVector;
import static lmu.utils.ObjectGetter.map;

public class NormalizedReplicateSetPlotting {

    static final Vector<Double> ERRBOXPLOTQUANTILES = toVector(0.05, 0.25, 0.5, 0.75, 0.95);

    empires.NormalizedReplicateSet normalizedReplicateSet;

    public NormalizedReplicateSetPlotting(empires.NormalizedReplicateSet normalizedReplicateSet) {
        this.normalizedReplicateSet = normalizedReplicateSet;
    }

    public BufferedImage plotNormalizationSteps(PlotCreator pc) {
        Vector<BufferedImage> steps = new Vector<>();
        int step = 0;
        for(empires.ShiftedGroup sg : normalizedReplicateSet.getNormalization().getSteps()) {

            step++;
            int idx = 0;
            for(Vector<Double> v : sg.getRawValues(normalizedReplicateSet.getInData().getLog2Data())) {
                pc.cumhist(""+(idx++), v, v.size(), false, true);
            }
            pc.setLabels("log2 signal", "cum. freq", null);
            pc.setTitle("step: %d raw", step);

            BufferedImage rawIm = pc.getImage();

            for(Vector<Double> v : sg.getShiftedValues(normalizedReplicateSet.getInData().getLog2Data())) {
                pc.cumhist(""+(idx++), v, v.size(), false, true);
            }
            pc.setLabels("log2 signal", "cum. freq", null);
            pc.setTitle("step: %d normed", step);
            BufferedImage normedIm = pc.getImage();
            steps.add(ImageUtils.concat(rawIm, normedIm));
        }
        return ImageUtils.vconcat(steps);
    }
    public BufferedImage plotSamplePairHeatMap(PlotCreator pc) {

        PairScore[][] errMatrix = empires.Normalization.getPairWiseErrors(normalizedReplicateSet.getInData().getLog2Data());
        Vector<String> repNames = normalizedReplicateSet.getInData().getReplicateNames();
        PlotCreator.HeatMapBuilder hmp = new PlotCreator.HeatMapBuilder();
        double maxErr = 0.0;
        for(int i=0; i<errMatrix.length; i++) {
            for(int j=i+1; j<errMatrix.length; j++) {
                double err = errMatrix[i][j].err.err;
                maxErr = Math.max(maxErr, err);
                hmp.addValue(repNames.get(i), repNames.get(j), err, true);
            }
        }
        for(int i=0; i<errMatrix.length; i++) {
            hmp.addValue(repNames.get(i), repNames.get(i), maxErr);
        }

        return hmp.getImage();
    }

    public Vector<UPair<Double>> getSignal2Error(double quantile) {
        int numErrs = normalizedReplicateSet.getNumRegularBackgrounds();

        Vector<UPair<Double>> quantiles = new Vector<>();
        for(int i=0; i<numErrs; i++) {
            double signal =  normalizedReplicateSet.getMeanSignalErrorForErrorIdx(i);
            empires.ErrorEstimationDistribution dist = normalizedReplicateSet.getErrorByIdx(i);
            quantiles.add(UPair.createU(signal, dist.getFoldChangeToCumulativeFrequency(quantile)));
        }
        return quantiles;
    }

    public BufferedImage plotBackGroundDistribs(PlotCreator pc) {
        int numErrs = normalizedReplicateSet.getNumRegularBackgrounds();

        int minwidth = numErrs * 20;
        pc.setWidth(Math.max(minwidth, 600));
        pc.reset();
        PlotCreator.BoxPlotBuilder bpb = pc.buildBoxPlot();
        for(int i=0; i<numErrs; i++) {
            String signal = String.format("%.2f", normalizedReplicateSet.getMeanSignalErrorForErrorIdx(i));
            empires.ErrorEstimationDistribution dist = normalizedReplicateSet.getErrorByIdx(i);
            Vector<Double> fcs = map(ERRBOXPLOTQUANTILES, (_p) -> dist.getFoldChangeToCumulativeFrequency(_p));
            bpb.addBox(signal, fcs);
        }
        pc.setTitle(normalizedReplicateSet.getInData().getReplicateSetName());
        pc.setLabels("mean signal", "log2fc error distrib", null);
        return bpb.plot();
    }
}
