package empires.test;

import empires.input.ExpressionSet;
import empires.input.ReplicateSetInfo;
import empires.plotting.NormalizedReplicateSetPlotting;
import lmu.utils.*;
import lmu.utils.plotting.PlotCreator;
import empires.Normalization;
import empires.NormalizedReplicateSet;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.junit.jupiter.api.RepeatedTest;

import java.awt.image.BufferedImage;
import java.io.File;
import java.util.HashMap;
import java.util.Vector;

import static lmu.utils.ObjectGetter.map;
import static lmu.utils.ObjectGetter.mapIndex;
import static empires.input.ExpressionSet.readDataGroupedByCondition;

import static org.junit.jupiter.api.Assertions.assertTrue;

class NormalizationTest {


    @org.junit.jupiter.api.Test
    @RepeatedTest(value = 20)
    void simpleTest() {
        int npoints = 10_000;
        Vector<Vector<Double>> vals = new Vector<>();
        int nsamples = 2 + (int)(Math.random() * 20);
        NormalDistribution meanND = new NormalDistribution(0, 5.0);
        Vector<Double> means = new Vector<>();
        for(int i=0; i<nsamples; i++) {

            double mean = meanND.sample();
            means.add(mean);
            NormalDistribution nd = new NormalDistribution(mean, 1.0);
            vals.add(mapIndex(npoints, (_i) -> nd.sample()));
        }

        System.out.printf("means: %s\n", means);
        Normalization norm = new Normalization(vals);

        Vector<Double> shifts =  norm.getShifts();


        double tolerance = 0.1;

        int anchorIdx = norm.getAnchorSampleIdx();

        System.out.printf("shifts: %s nachor: %d\n", shifts, anchorIdx);
        double targetMean = means.get(anchorIdx);

        for(int i = 0; i<nsamples; i++) {
            if(anchorIdx == i) {
                assert(Math.abs(shifts.get(i)) < tolerance);
                continue;
            }
            double shift = shifts.get(i);
            Vector<Double> shifted = map(vals.get(i), (_d) -> _d + shift);

            double median = NumUtils.median(shifted, (_d) -> _d);

            System.out.printf("%d: mean was: %.2f target mean: %.2f shift: %.2f median: %.2f (%s)\n", i, means.get(i), targetMean, shift, median, NumUtils.getNumInfo(shifted));
            assertTrue(Math.abs(median - targetMean) < tolerance,
                    String.format("target mean: %.2f shifted %.2f -> %.2f", targetMean, means.get(i), means.get(i) + shift));


        }




    }

    @org.junit.jupiter.api.Test
    @RepeatedTest(value = 10)
    void multiSDTest() {
        multiSDTest(false);
    }


    static void multiSDTest(boolean withplot) {
        MultiSDSimulation simulationInfo = new MultiSDSimulation(withplot);

        Normalization norm = new Normalization(simulationInfo.rsi.getLog2Data());

        Vector<Double> shifts =  norm.getShifts();


        double tolerance = 0.1;

        int anchorIdx = norm.getAnchorSampleIdx();

        System.out.printf("shifts: %s nachor: %d\n", shifts, anchorIdx);
        double targetMean = simulationInfo.means.get(anchorIdx);

        for(int i = 0; i< simulationInfo.numSamples; i++) {
            if(anchorIdx == i) {
                assert(Math.abs(shifts.get(i)) < tolerance);
                continue;
            }
            double shift = shifts.get(i);



            double shiftedmean  = shift + simulationInfo.means.get(i);

            System.out.printf("%d: mean was: %.2f target mean: %.2f shift: %.2f shiftedmean: %.2f\n", i,
                    simulationInfo.means.get(i), targetMean, shift, shiftedmean);
            assertTrue(Math.abs(shiftedmean- targetMean) < tolerance,
                    String.format("target mean: %.2f shifted %.2f -> %.2f (diff: %.2f > %.2f)",
                    targetMean, simulationInfo.means.get(i), simulationInfo.means.get(i) + shift,
                    Math.abs(simulationInfo.means.get(i) + shift - targetMean), tolerance));


        }



    }

    public static void main(String[] args) {
        SimpleOptionParser cmd = new SimpleOptionParser("multisd", "subsample", "targetnumsample");
        cmd.setFile("subsample");
        cmd.setInt("targetnumsample");
        cmd.setSwitches("multisd");
        cmd.setOptional("subsample");
        cmd.setDefault("targetnumsample", "6");

        if(!OptionParser.parseParams(args, true, false, true, true, cmd))
            return;

        if(cmd.isSet("multisd")) {
            multiSDTest(true);
        }

        File subsample = cmd.getOptionalFile("subsample");

        Logger log = LogConfig.getLogger();
        if(subsample != null) {
            int ntarget = cmd.getInt("targetnumsample");
            HashMap<String, empires.input.ReplicateSetInfo> conditions = ExpressionSet.readDataGroupedByCondition(subsample);

            PlotCreator pc = new PlotCreator();

            Vector<BufferedImage> bims = new Vector<>();
            for(ReplicateSetInfo rsi : conditions.values()) {
                if(rsi.getNumReplicates() <= ntarget) {
                    log.info("skip %s: has %d replicates, target is: %d", rsi.getReplicateSetName(), rsi.getNumReplicates(), ntarget);
                    continue;
                }

                NormalizedReplicateSet nrs = new NormalizedReplicateSet(rsi);

                Vector<Integer> replicates = Normalization.getBestSubSample(rsi.getLog2Data(), ntarget);
                NormalizedReplicateSet nrs_sub = new NormalizedReplicateSet(rsi.getReplicateSubSet(replicates));

                BufferedImage bim1 = new empires.plotting.NormalizedReplicateSetPlotting(nrs).plotBackGroundDistribs(pc);
                BufferedImage bim2 = new NormalizedReplicateSetPlotting(nrs_sub).plotBackGroundDistribs(pc);

                bims.add(ImageUtils.concat(bim1, bim2));

            }

            pc.destroy();
            ImageUtils.showImage("subsample test", ImageUtils.vconcat(bims), false);
        }

    }
}