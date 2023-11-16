package empires.test;

import empires.input.ReplicateSetInfo;
import lmu.utils.ImageUtils;
import lmu.utils.LogConfig;
import lmu.utils.Logger;
import lmu.utils.NumUtils;
import lmu.utils.plotting.PlotCreator;
import lmu.utils.plotting.CachedPlotCreator;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.awt.image.BufferedImage;
import java.io.File;
import java.util.Vector;
import java.util.function.Function;

import static lmu.utils.IteratorUtils.rangev;
import static lmu.utils.ObjectGetter.*;

public class MultiSDSimulation {
    public empires.input.ReplicateSetInfo rsi;
    public int numSamples;
    public Vector<Double> means;
    public BufferedImage signal2sdPlot;
    public BufferedImage normedSignalsPlot;
    public BufferedImage scatterPlots;
    public Vector<Double> signals;
    public Vector<Double> sds;

    public Double SD_FACTOR = null;
    boolean simulateNANS = true;
    int npoints = 10_000;
    double MAXVAL = 10_000.0;
    boolean withPlots = true;

    public MultiSDSimulation(boolean withPlots) {
        this(null, withPlots, null, false, 10_000);
    }

    public MultiSDSimulation(Integer InNumSamples, boolean withPlots, Vector<Double> inSignals, boolean simulateNANS, int npoints) {

        this.npoints = npoints;
        this.simulateNANS = simulateNANS;
        numSamples = (InNumSamples != null) ? InNumSamples : 2 + (int)(Math.random() * 20);

        GammaDistribution rnd = new GammaDistribution(1.0, 2.0);
        if(inSignals == null) {
            this.signals = new Vector<>();
            while(signals.size() < npoints) {
                double signal = Math.pow(10, 1 + rnd.sample());
                if(signal > MAXVAL)
                    continue;
                signals.add(signal);
            }
        } else {
            signals = inSignals;
        }
        //this.signals = (inSignals != null) ? inSignals : mapIndex(npoints, (_i) -> 10.0 + (Math.random() * MAXVAL));
        //this.signals = (inSignals != null) ? inSignals : mapIndex(npoints, (_i) -> Math.min(MAXVAL, Math.pow(10, 1 + rnd.sample())));
        System.out.printf("signals: %s\n", NumUtils.getNumInfo(this.signals).getInfoWithQ());
        MAXVAL = NumUtils.max(this.signals);
        this.withPlots = withPlots;

    }

    public void setNumPoints(int npoints) {
        this.npoints = npoints;
    }

    public void setSDFactor(double factor) {
        SD_FACTOR = factor;
    }
    double MAXSD;
    Function<Double, Double> signal2sd;


    public double getSDToSignal(double signal) {
        return Math.max(0.2, MAXSD - signal2sd.apply(signal));
        //return 0.3;
    }

    public void writeImage(String prefix, File outdir) {
        ImageUtils.saveImage(signal2sdPlot, new File(outdir, prefix + "_signal2sd.png"));
        ImageUtils.saveImage(normedSignalsPlot, new File(outdir, prefix + "_normed_signals.png"));
        ImageUtils.saveImage(scatterPlots, new File(outdir, prefix + "_scatterplots.png"));

    }
    public void simulate() {


        Logger log = LogConfig.getLogger();
        SD_FACTOR = (SD_FACTOR != null) ? SD_FACTOR : Math.pow(2.0, new NormalDistribution(0, 0.3).sample());

        signal2sd = (_d) -> SD_FACTOR * NumUtils.logN(Math.pow(_d, 1/3.0), 2.0) / 3.0;

        Vector<Vector<Double>> vals = new Vector<>();

        MAXSD = signal2sd.apply(5 * MAXVAL);



        NormalDistribution meanND = new NormalDistribution(0, 5.0);

        Vector<Double> logsignals = map(signals, (_d) -> NumUtils.logN(_d, 2.0));

        sds = map(signals, (_d) -> getSDToSignal(_d));


        PlotCreator pc = (!withPlots) ? null : CachedPlotCreator.getPlotCreator();

        if(pc != null) {
            pc.scatter("", rangev(sds.size()), (_i) -> logsignals.get(_i), (_i) -> sds.get(_i));

            pc.setLabels("logsignal", "sd", null);
            signal2sdPlot = pc.getImage();

        }
        rsi = new ReplicateSetInfo("simulation", mapIndex(numSamples, (_i) -> String.format("sample%02d", _i + 1)),
                mapIndex(npoints, (_i) -> String.format("feat%05d", _i + 1)));


        means = new Vector<>();

        Vector<Integer> idx = rangev(sds.size());
        for(int i=0; i<numSamples; i++) {

            double mean = meanND.sample();
            means.add(mean);

            Vector<Double> samplevals = new Vector<>();
            for(int p = 0; p<npoints; p++) {

                double basesignal = logsignals.get(p);
                double signal = new NormalDistribution( mean + basesignal, sds.get(p)).sample();
                samplevals.add(signal);
            }


            vals.add(samplevals);

            if(simulateNANS) {
                Vector<Integer> nanIndeces = shuffleN(idx, (int)(0.2 * idx.size()));
                log.info("sample: %d got %d/%d NAN indeces", i, nanIndeces.size(), samplevals.size());
                for(int idx2nan : nanIndeces) {
                    samplevals.set(idx2nan, Double.NaN);
                }
            }




            log.trace("%d mean: %.2f got %d/%d nan values\n", i, mean, filteredSize(samplevals, (_d) -> Double.isNaN(_d)), samplevals.size());
            rsi.setLog2Data(i, samplevals);

            Vector<Double> nonnan = filter(samplevals, (_d) -> !Double.isNaN(_d));
            log.trace("sample: %d : %s\n", i, NumUtils.getNumInfo(nonnan).getInfoWithQ());

            if(pc != null) {
                pc.cumhist(String.format("s:%d (%.2f) ", i, mean), nonnan, nonnan.size(), false, true);
            }

        }

        if(pc != null) {
            pc.setLabels("signals", "freq", "bottomright");
            normedSignalsPlot = pc.getImage();

        }


        if(pc != null) {
            Vector<BufferedImage> bims = new Vector<>();
            for (int i = 0; i < numSamples - 1; i++) {
                for (int j = i + 1; j < numSamples; j++) {
                    Vector<Double> v1 = vals.get(i);
                    Vector<Double> v2 = vals.get(j);
                    Vector<Integer> nonnans = filter(rangev(v1.size()), (_i) -> !Double.isNaN(v1.get(_i)) && !Double.isNaN(v2.get(_i)));
                    pc.scatter(i + ":" + j, nonnans, (_i) -> v1.get(_i), (_i) -> v2.get(_i));
                    pc.setLabels("" + i, "" + j, null);
                    ;
                    bims.add(pc.getImage());
                }
            }
            scatterPlots = ImageUtils.vconcat(bims);

        }


    }
}
