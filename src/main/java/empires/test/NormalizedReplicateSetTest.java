package empires.test;

import empires.input.ExpressionSet;
import empires.input.ReplicateSetInfo;
import empires.plotting.NormalizedReplicateSetPlotting;
import lmu.utils.*;
import lmu.utils.plotting.PlotCreator;
import empires.ErrorEstimationDistribution;
import empires.NormalizedReplicateSet;
import lmu.utils.plotting.CachedPlotCreator;
import org.apache.commons.math3.distribution.*;

import java.io.File;
import java.util.HashMap;
import java.util.Vector;
import java.util.function.BiFunction;
import java.util.function.Supplier;

import static lmu.utils.ObjectGetter.toVector;

class NormalizedReplicateSetTest {

    public static void sdCheck(File od) {
        int nreps = 4;
        int nfeatures = 2000;
        NormalDistribution loglevelND = new NormalDistribution(100, 20);
        Vector<Pair<String, BiFunction<Double, Double, Supplier<Double>>>>
                distribs2test = toVector(Pair.create("normal", (_mean, _sd) -> () -> new NormalDistribution(_mean, _sd).sample()),
                Pair.create("chisquared", (_mean, _sd) -> () ->  new ChiSquaredDistribution(Math.pow(_sd, 2.0) / 2.0).sample()),
                Pair.create("exponential", (_mean, _sd) -> () -> new ExponentialDistribution(_sd).sample()),
                Pair.create("poisson", (_mean, _sd) -> () -> (double) (new PoissonDistribution(Math.pow(_sd, 2.0))).sample())
                );

        PlotCreator pc = CachedPlotCreator.getPlotCreator();

        for(int di=0; di<distribs2test.size(); di++) {

            Vector<UPair<Double>> real2sd = new Vector<>();

            String dname = distribs2test.get(di).getFirst();
            for(double sd = 0.5; sd < 6.0; sd += 0.2) {


                Vector<Vector<Double>> logsignals = new Vector<>();
                for(int fi=0; fi<nfeatures; fi++) {
                    double loglevel = 2.0 + Math.abs(loglevelND.sample());

                    Supplier<Double> getter = distribs2test.get(di).getSecond().apply(loglevel, sd);
                    Vector<Double> replicates = new Vector<>();
                    for(int repi = 0; repi < nreps; repi++) {
                        replicates.add(getter.get());
                    }
                    logsignals.add(replicates);
                }
                ErrorEstimationDistribution errD = new ErrorEstimationDistribution(logsignals);
                System.out.printf("%s SD check %.2f = %.2f\n", dname, sd, errD.getSD());

                real2sd.add(UPair.createU(sd, errD.getSD()));
            }
            pc.scatter(dname, real2sd, (_p) -> _p.getFirst(), (_p) -> _p.getSecond());

        }

        pc.setLabels("real SD", "estimated SD", "topleft");
        pc.writeImage(new File(od, "distribestimatecheck.png"));
    }
    public static void writePlots(empires.test.MultiSDSimulation msd, PlotCreator pc, empires.input.ReplicateSetInfo rsi, File od) {
        NormalizedReplicateSet nrs = new NormalizedReplicateSet(rsi);
        empires.plotting.NormalizedReplicateSetPlotting plotting = new NormalizedReplicateSetPlotting(nrs);
        ImageUtils.saveImage(plotting.plotBackGroundDistribs(pc), new File(od, rsi.getReplicateSetName()+"_errdist.png"));

        if(msd == null)
            return;



                /*
        Vector<UPair<Double>> simul = new Vector<>();
        Vector<UPair<Double>> estimated = new Vector<>();
        for(int errIdx = 0 ; errIdx < nrs.getNumRegularBackgrounds(); errIdx++) {
            ErrorEstimationDistribution errD = nrs.getErrorByIdx(errIdx);
            double meanSignal = nrs.getMeanSignalErrorForErrorIdx(errIdx);
            double sd = msd.getSDToSignal(meanSignal);
            simul.add(UPair.createU(meanSignal, sd));
            estimated.add(UPair.createU(meanSignal, Math.sqrt(2) * errD.getSD()));
        }
        pc.line("real", simul, (_p) -> _p.getFirst(), (_p) -> _p.getSecond() );
        pc.line("estimated", estimated, (_p) -> _p.getFirst(), (_p) -> _p.getSecond() );
        pc.setLabels("signal", "SD", "topright");
        pc.writeImage(new File(od, "sd_real_to_est.png"));
        */

        //ImageUtils.saveImage(plotting.plotSamplePairHeatMap(pc), new File(od, rsi.getReplicateSetName()+"_sampleHeatMap.png"));
    }

    public static void main(String[] args) {
        SimpleOptionParser cmd = new SimpleOptionParser("edir", "od", "simul", "sdcheck");
        cmd.setDir("edir");
        cmd.setDir("od");
        cmd.setSwitches("simul", "sdcheck");
        cmd.setOptional("edir");

        if(!OptionParser.parseParams(args, true, false, true, true, cmd))
            return;

        File od = cmd.getFile("od");

        if(cmd.isSet("sdcheck")) {
            sdCheck(od);
            return;
        }

        PlotCreator pc = CachedPlotCreator.getPlotCreator();
        if(cmd.isSet("simul")) {
            empires.test.MultiSDSimulation msd = new MultiSDSimulation(true);
            msd.simulate();
            empires.input.ReplicateSetInfo rsi = msd.rsi;
            writePlots(msd, pc, rsi, od);
            ImageUtils.saveImage(msd.signal2sdPlot, new File(od, "signal2sd.png"));
            ImageUtils.saveImage(msd.normedSignalsPlot, new File(od, "normedsignals.png"));
            ImageUtils.saveImage(msd.scatterPlots, new File(od, "scatters.png"));

            pc.destroy();
            return;

        }

        File edir = cmd.getOptionalFile("edir");

        if(edir != null) {
            HashMap<String, empires.input.ReplicateSetInfo> data =  ExpressionSet.readDataGroupedByCondition(edir);
            for(ReplicateSetInfo rsi : data.values()) {
                writePlots(null, pc, rsi, od);
            }

        }

        pc.destroy();
    }

}