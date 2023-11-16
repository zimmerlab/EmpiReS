package empires.test;

import lmu.utils.*;
import lmu.utils.plotting.PlotCreator;
import empires.DistribUtils;
import empires.DoubleDiffResult;
import org.junit.jupiter.api.RepeatedTest;

import java.awt.image.BufferedImage;
import java.util.Vector;

import static lmu.utils.IteratorUtils.rangev;
import static lmu.utils.ObjectGetter.map;
import static org.junit.jupiter.api.Assertions.assertTrue;

class DoubleDiffManagerTest {
    @RepeatedTest(value = 10)
    public void uniformTest() {

        empires.test.UniformTest uniformTest = new empires.test.UniformTest();
        uniformTest.simulateNans = Math.random() < 0.5;
        uniformTest.simulate();
        Vector<DoubleDiffResult> diffResults = uniformTest.getDoubleDiffResults();
        String err1 = DistribUtils.checkUniform(map(diffResults, (_d) -> _d.pval), 0.02, 10);
        String err2 = DistribUtils.checkUniform(map(diffResults, (_d) -> _d.pval), 0.04, 20);
        assertTrue(err1 == null, String.format("%d x %d: %s", uniformTest.nrs1.getNumReplicates(), uniformTest.nrs2.getNumReplicates(),  err1));
        assertTrue(err2 == null, String.format("%d x %d: %s", uniformTest.nrs1.getNumReplicates(), uniformTest.nrs2.getNumReplicates(),  err2));
    }


    public static void main(String[] args) {
        SimpleOptionParser cmd = new SimpleOptionParser("plotuniform", "nrep1", "nrep2", "simulateNANs", "npoints", "SDfactor");
        cmd.setSwitches("plotuniform");
        cmd.setInt("nrep1", "nrep2");
        cmd.setDouble("SDfactor");
        cmd.setDefault("nrep1", "-1");
        cmd.setDefault("nrep2", "-1");
        cmd.setDefault("npoints", "10000");
        cmd.setOptional("SDfactor");

        cmd.setSwitches("simulateNANs");
        if(!OptionParser.parseParams(args, true, false, true, false, cmd))
            return;


        Logger log = LogConfig.getLogger();

        if(cmd.isSet("plotuniform")) {

            int nrep1 = cmd.getInt("nrep1");
            int nrep2 = cmd.getInt("nrep2");
            boolean withNANs = cmd.isSet("simulateNANs");
            empires.test.UniformTest uniformTest = new UniformTest((nrep1 < 0) ? null : nrep1, (nrep2 < 0) ? null : nrep2, withNANs, 10_000);

            uniformTest.setNumFeatures(cmd.getInt("npoints"));

            if(cmd.isOptionSet("SDfactor")) {
                uniformTest.setSDFactor(cmd.getDouble("SDfactor"));
            }
            uniformTest.simulate();


            Vector<Integer> fidxvec = rangev(uniformTest.cond1.sds.size());


            PlotCreator pc = new PlotCreator();

            pc.scatter("", fidxvec, (_i) -> uniformTest.cond1.sds.get(_i), (_i) -> uniformTest.nrs1.getError(uniformTest.cond1.rsi.getFeatureName(_i)).getSD());
            pc.abline("", null, null, 0.0, 1.0);
            pc.setLabels("simulated SD", "fitted SD", null);
            ImageUtils.showImage("SD scatter cond1", pc.getImage());



            pc.scatter("", fidxvec, (_i) -> uniformTest.cond2.sds.get(_i), (_i) -> uniformTest.nrs2.getError(uniformTest.cond2.rsi.getFeatureName(_i)).getSD());
            pc.abline("", null, null, 0.0, 1.0);
            pc.setLabels("simulated SD", "fitted SD", null);
            ImageUtils.showImage("SD scatter cond2", pc.getImage());


            pc.scatter("", fidxvec, (_i) -> uniformTest.cond1.sds.get(_i), (_i) -> uniformTest.cond2.sds.get(_i));
            pc.abline("", null, null, 0.0, 1.0);
            pc.setLabels("simulated SD 1", "sim SD 2", null);
            ImageUtils.showImage("sim SD s", pc.getImage());






            long t1 = System.currentTimeMillis();
            Vector<DoubleDiffResult> doubleDiffResults = uniformTest.getDoubleDiffResults();
            long t2 = System.currentTimeMillis();
            String err = DistribUtils.checkUniform(map(doubleDiffResults, (_d) -> _d.pval), 0.03);

            log.info("double diff took: %.2f sec (%d tests)", (t2 - t1) / 1000.0, doubleDiffResults.size());
            pc.dens_scatter2("tmp", uniformTest.getDoubleDiffFeatureCombis(),
                    (_p) -> _p.getFirst().size(), (_p) -> _p.getSecond().size());



            pc.setTitle("%d x %d %d\n%s", uniformTest.nrs1.getInData().getReplicateNames().size(), uniformTest.nrs2.getInData().getReplicateNames().size(),
                    uniformTest.nrs1.getInData().getNumFeatures(),
                    err);

            pc.setLabels("isoform set size1", "isoform set size2", null);
            BufferedImage isoformSets = pc.getImage();

            pc.cumhist("pvals", doubleDiffResults, (_d) -> _d.pval, uniformTest.getResults().size(), false, true);
            pc.abline("", null, null, 0.0, 100.0);
            ImageUtils.showImage("test", ImageUtils.concat(isoformSets, pc.getImage()));
            pc.destroy();

        }
    }
}