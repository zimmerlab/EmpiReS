package empires.test;

import empires.EmpiRe;
import empires.input.BackgroundProviderOption;
import empires.input.GeneralOptions;
import lmu.utils.ImageUtils;
import lmu.utils.NumUtils;
import lmu.utils.OptionParser;
import lmu.utils.SimpleOptionParser;
import lmu.utils.plotting.PlotCreator;
import org.junit.jupiter.api.RepeatedTest;

import java.awt.image.BufferedImage;
import java.io.File;
import java.util.Vector;

import static lmu.utils.IteratorUtils.rangev;
import static lmu.utils.ObjectGetter.*;
import static org.junit.jupiter.api.Assertions.assertTrue;

class EmpiReTest {

    @RepeatedTest(value = 10)
    public void uniformTest() {

        empires.test.UniformTest uniformTest = new empires.test.UniformTest();
        uniformTest.simulateNans = Math.random() < 0.5;
        uniformTest.simulate();
        Vector<empires.DiffExpResult> diffResults = uniformTest.getResults();
        String err1 = empires.DistribUtils.checkUniform(map(diffResults, (_d) -> _d.pval), 0.02, 10);
        String err2 = empires.DistribUtils.checkUniform(map(diffResults, (_d) -> _d.pval), 0.04, 20);
        assertTrue(err1 == null, String.format("%d x %d: %s", uniformTest.nrs1.getNumReplicates(), uniformTest.nrs2.getNumReplicates(),  err1));
        assertTrue(err2 == null, String.format("%d x %d: %s", uniformTest.nrs1.getNumReplicates(), uniformTest.nrs2.getNumReplicates(),  err2));
    }


    @RepeatedTest(value = 1)
    public void replicateNumTest() {

        empires.SingleFeatureDiffExp.CORRECT_OUTLIERS = false;
        empires.DiffExpResult.PERCENT_ROBUST_FEATURE_USED = 1.0;

        for(int nrep1 = 2; nrep1 < 20; nrep1++) {
            for(int nrep2 = 2; nrep2 < 20; nrep2++) {
                empires.test.UniformTest uniformTest = new empires.test.UniformTest(nrep1, nrep2, false, 10_000);
                uniformTest.simulate();
                Vector<empires.DiffExpResult> diffResults = uniformTest.getResults();
                String err = empires.DistribUtils.checkUniform(map(diffResults, (_d) -> _d.pval), 0.03, 20);
                assertTrue(err == null, String.format("%d x %d: %s", uniformTest.nrs1.getNumReplicates(), uniformTest.nrs2.getNumReplicates(),  err));
            }
        }

    }

    @RepeatedTest(value = 1)
    public void replicateNumTestSingle() {

        empires.SingleFeatureDiffExp.CORRECT_OUTLIERS = false;
        empires.DiffExpResult.PERCENT_ROBUST_FEATURE_USED = 1.0;

        for(int nrep1 = 2; nrep1 < 20; nrep1++) {
            for(int nrep2 = 3; nrep2 < 20; nrep2++) {
                empires.test.UniformTest uniformTest = new empires.test.UniformTest(nrep1, nrep2, false, 50_000);
                uniformTest.simulate();
                Vector<empires.DiffExpResult> diffResults = uniformTest.getResults();

                String err = empires.DistribUtils.checkUniform(uniformTest.empiRe.getPerComparisonPvals(), 0.03, 20);
                assertTrue(err == null, String.format("%d x %d: %s", uniformTest.nrs1.getNumReplicates(), uniformTest.nrs2.getNumReplicates(),  err));
            }
        }

    }

    public static void main(String[] args) {
        empires.input.GeneralOptions generalOptions = new GeneralOptions();
        empires.input.BackgroundProviderOption backgroundProviderOption = new BackgroundProviderOption();
        SimpleOptionParser cmd = new SimpleOptionParser("plotuniform", "nrep1", "nrep2", "simulateNANs", "npoints", "SDfactor", "featuresToCombine", "plotoutdir");
        cmd.setSwitches("plotuniform");
        cmd.setInt("nrep1", "nrep2", "featuresToCombine");
        cmd.setDouble("SDfactor");
        cmd.setDefault("nrep1", "-1");
        cmd.setDefault("nrep2", "-1");
        cmd.setDefault("npoints", "10000");
        cmd.setDefault("featuresToCombine", "3");
        cmd.setDir("plotoutdir");
        cmd.setOptional("SDfactor", "plotoutdir");

        cmd.setSwitches("simulateNANs");
        if(!OptionParser.parseParams(args, true, false, true, false, cmd, generalOptions, backgroundProviderOption))
            return;

        generalOptions.apply();

        if(cmd.isSet("plotuniform")) {

            int nrep1 = cmd.getInt("nrep1");
            int nrep2 = cmd.getInt("nrep2");
            boolean withNANs = cmd.isSet("simulateNANs");
            empires.test.UniformTest uniformTest = new UniformTest((nrep1 < 0) ? null : nrep1, (nrep2 < 0) ? null : nrep2, withNANs, cmd.getInt("npoints"));

            File plotoutDir = cmd.getOptionalFile("plotoutdir");
            if(plotoutDir != null) {
                uniformTest.setPlotSignals(true);
            }
            uniformTest.provider = backgroundProviderOption.getStrategy();
            uniformTest.numFeaturePerCombined = cmd.getInt("featuresToCombine");
            uniformTest.setNumFeatures(cmd.getInt("npoints"));

            if(cmd.isOptionSet("SDfactor")) {
                uniformTest.setSDFactor(cmd.getDouble("SDfactor"));
            }
            uniformTest.simulate();

            if(plotoutDir != null) {
                uniformTest.cond1.writeImage("cond1", plotoutDir);
                uniformTest.cond2.writeImage("cond2", plotoutDir);
            }
                    /* toremove
        log.info("got %d background distribs\n", errorBackGrounds.size());
        for(int i=0; i<meanSignalValue2Error.size(); i++) {
            ErrorEstimationDistribution e = errorBackGrounds.get(i);
            System.out.printf("next: %d signal: %.2f %s\n", i, meanSignalValue2Error.get(i), new SparseCumulativeDistribution(e).quantiles);
            PlotCreator pc = CachedPlotCreator.getPlotCreator();
            pc.setTitle("%s: %d errbg", ci.getReplicateSetName(), i);
            new SparseCumulativeDistribution(e).drawLine(pc, "1");
            ImageUtils.showImage("t", pc.getImage());
        }

         */


            Vector<Integer> fidxvec = rangev(uniformTest.cond1.sds.size());


            PlotCreator pc = new PlotCreator();

            pc.scatter("", fidxvec, (_i) -> uniformTest.cond1.sds.get(_i), (_i) -> uniformTest.nrs1.getError(uniformTest.cond1.rsi.getFeatureName(_i)).getSD());
            pc.abline("", null, null, 0.0, 1.0);
            pc.setLabels("simulated SD", "fitted SD", null);
            if(plotoutDir != null) {
                ImageUtils.saveImage(pc.getImage(), new File(plotoutDir, "sd_scatter_cond1.png"));

            } else {
                ImageUtils.showImage("SD scatter cond1", pc.getImage());
            }




            pc.scatter("", fidxvec, (_i) -> uniformTest.cond2.sds.get(_i), (_i) -> uniformTest.nrs2.getError(uniformTest.cond2.rsi.getFeatureName(_i)).getSD());
            pc.abline("", null, null, 0.0, 1.0);
            pc.setLabels("simulated SD", "fitted SD", null);
            if(plotoutDir != null) {
                ImageUtils.saveImage(pc.getImage(), new File(plotoutDir, "sd_scatter_cond2.png"));

            } else {

                ImageUtils.showImage("SD scatter cond2", pc.getImage());
            }

            pc.scatter("", fidxvec, (_i) -> uniformTest.cond1.sds.get(_i), (_i) -> uniformTest.cond2.sds.get(_i));
            pc.abline("", null, null, 0.0, 1.0);
            pc.setLabels("simulated SD 1", "sim SD 2", null);
            if(plotoutDir != null) {
                ImageUtils.saveImage(pc.getImage(), new File(plotoutDir, "sim_SDs.png"));
            } else {
                ImageUtils.showImage("sim SD s", pc.getImage());
            }


            uniformTest.empiRe.setCollectComparisonPvals(true);

            String err = empires.DistribUtils.checkUniform(map(uniformTest.getResults(), (_d) -> _d.pval), 0.03);




            pc.setTitle("single feature pvalues");
            pc.cumhist("z-score", uniformTest.empiRe.getPerComparisonPvals(), uniformTest.empiRe.getPerComparisonPvals().size(), false, true);
            pc.cumhist("fcdist", uniformTest.empiRe.getPerComparisonFcDistPvals(), uniformTest.empiRe.getPerComparisonFcDistPvals().size(), false, true);
            pc.abline("", null, null, 0.0, 100.0);
            pc.setLabels("p-value", "cumulative frequency", "bottomright");
            BufferedImage bim1 = pc.getImage();


            pc.setTitle("%d x %d %d/%d features %s", uniformTest.nrs1.getInData().getReplicateNames().size(), uniformTest.nrs2.getInData().getReplicateNames().size(),
                    uniformTest.combFeatureVec.size(), fidxvec.size(),
                    err);
            pc.cumhist("z-score", uniformTest.getResults(), (_d) -> _d.pval, uniformTest.getResults().size(), false, true);
            pc.cumhist("fcdist", uniformTest.getResults(), (_d) -> _d.fcEstimatePval, uniformTest.getResults().size(), false, true);
            pc.abline("", null, null, 0.0, 100.0);
            pc.setLabels("p-value", "cumulative frequency", "bottomright");

            BufferedImage bim2 = pc.getImage();

            if(plotoutDir != null) {
                ImageUtils.saveImage(ImageUtils.concat(bim1, bim2), new File(plotoutDir, "uniformtest.png"));
            } else {
                ImageUtils.showImage("uniform-test", ImageUtils.concat(bim1, bim2));
            }


            empires.SingleFeatureDiffExp.COMBINE_SHIFTED_DISTRIBS_FOR_EACH_REPLICATE_PAIR = true;
            empires.EmpiRe emp2 = new EmpiRe();
            emp2.setCollectComparisonPvals(true);
            Vector<String> allfeatures = new Vector<>();
            apply(uniformTest.combinedFeatures.values(), (_v) -> allfeatures.addAll(_v));
            empires.DiffExpManager diffExpManager = new empires.DiffExpManager(uniformTest.nrs1, uniformTest.nrs2, allfeatures,  emp2);
            Vector<Double> differences = new Vector<>();
            for(empires.DiffExpResult de : uniformTest.getResults()) {

                empires.DiffExpResult diffExp = new empires.DiffExpResult(diffExpManager, de.combinedFeatureName, de.featureNames);
                for(int i=0; i<de.featureNames.size(); i++) {
                    differences.add(de.perFeaturePvals.get(i).getSecond() - diffExp.perFeaturePvals.get(i).getSecond());
                    if(differences.size() % 100 == 0) {
                        System.out.printf("differences: fc-comb vs quick %s\n", NumUtils.getNumInfo(differences).getInfoWithQ());
                    }
                }
            }
            pc.cumhist("z-score", emp2.getPerComparisonPvals(), emp2.getPerComparisonPvals().size(), false, true);
            pc.cumhist("comb-fcdist", emp2.getPerComparisonFcDistPvals(), emp2.getPerComparisonFcDistPvals().size(), false, true);
            pc.abline("", null, null, 0.0, 100.0);
            pc.setLabels("p-value", "cumulative frequency", "bottomright");
            BufferedImage bim3 = pc.getImage();
            ImageUtils.showImage("fcdist-unform test", bim3);



            pc.destroy();
        }
    }
}