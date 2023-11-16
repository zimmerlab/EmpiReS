package empires.input;

import empires.EmpiRe;
import empires.plotting.DiffExpTable;
import empires.plotting.NormalizedReplicateSetPlotting;
import lmu.utils.*;
import lmu.utils.plotting.PlotCreator;
import lmu.utils.tuple.Tuple3;
import lmu.utils.plotting.CachedPlotCreator;

import java.awt.image.BufferedImage;
import java.io.File;
import java.io.PrintWriter;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Vector;
import java.util.function.Function;

import static lmu.utils.ObjectGetter.*;
import static lmu.utils.ObjectGetter.getPairs;

public class ExpressionSet {

    static class ConditionConfiguration{

        String condition;
        Vector<String> replicateNames = new Vector<>();
        HashMap<String, Integer> replicate2InputIndex = new HashMap<>();
        HashMap<String, Integer> replicate2localIndex = new HashMap<>();
        HashMap<Integer, Integer> global2localidx = new HashMap<>();
        empires.input.ReplicateSetInfo data;

        public ConditionConfiguration(String condition) {
            this.condition = condition;
        }

        public void add(String replicateName, int inputIndex) {
            replicate2InputIndex.put(replicateName, inputIndex);
            replicate2localIndex.put(replicateName, replicateNames.size());
            global2localidx.put(inputIndex, replicateNames.size());
            replicateNames.add(replicateName);
        }

        public Collection<Integer> getManagedGlobalIndeces() {
            return global2localidx.keySet();
        }

        public void setFeatureLabels(Vector<String> names) {
            data = new empires.input.ReplicateSetInfo(condition, replicateNames, names);
        }

        public void setFeatureData(int replicateInputIdx, int featureIdx, double val) {
            int local_idx = global2localidx.get(replicateInputIdx);
            data.log2replicate2FeatureData.get(local_idx).set(featureIdx, (val == 0.0) ? Double.NaN : NumUtils.logN(val, 2.0));
        }
    }

    public static HashMap<String, empires.input.ReplicateSetInfo> readDataGroupedByCondition(File exprSetDir) {
        return readDataGroupedByCondition(empires.input.EBROWSER_INPUT.EXPRS.get(exprSetDir), empires.input.EBROWSER_INPUT.FDATA.get(exprSetDir), empires.input.EBROWSER_INPUT.PDATA.get(exprSetDir));
    }

    public static HashMap<String, empires.input.ReplicateSetInfo> readDataGroupedByCondition(File expr, File f, File p_data) {


        HashMap<String, ConditionConfiguration> cond2config = new HashMap<>();


        int cidx = 0;
        for (String[] sp : FileUtils.getFieldSets(p_data, "\t")) {

            cond2config.computeIfAbsent(sp[1], (_c) -> new ConditionConfiguration(_c)).add(sp[0], cidx);
            cidx++;

        }
        ConditionConfiguration[] configs = new ConditionConfiguration[cidx];


        for (ConditionConfiguration cc : cond2config.values()) {
            cc.getManagedGlobalIndeces().stream().forEach((_idx) -> configs[_idx] = cc);
        }


        Vector<String> featureLabels = ObjectGetter.map(FileUtils.readCol(f, "\t", 0), (_s) -> _s.split("\\.")[0]);

        cond2config.values().stream().forEach((_c) -> _c.setFeatureLabels(featureLabels));


        int featureIdx = 0;
        for (String[] sp : FileUtils.getFieldSetsIterable(expr, "\t")) {
            if (featureIdx >= featureLabels.size())
                throw new FRuntimeException("inconsistent data! got %d feature labels, but more data lines", featureLabels.size());
            for (int i = 0; i < sp.length; i++) {
                configs[i].setFeatureData(i, featureIdx, Double.parseDouble(sp[i]));
            }
            featureIdx++;
        }

        return buildMap(cond2config.entrySet(), (_e) -> _e.getKey(), (_e) -> _e.getValue().data);


    }


    public static void main(String[] args) {
        BackgroundProviderOption backgroundProviderOption = new BackgroundProviderOption();
        empires.input.GeneralOptions generalOptions = new GeneralOptions();
        SimpleOptionParser cmd = new SimpleOptionParser("inputdir", "cond1", "cond2", "o", "distribout", "interactive",
                "fcthreshold", "fdrthreshold", "getbestsubN", "normdetails");
        cmd.setDir("inputdir");
        cmd.setOutFile("o");
        cmd.setDouble("fcthreshold", "fdrthreshold");
        cmd.setDefault("fcthreshold", "0.9");
        cmd.setDefault("fdrthreshold", "0.05");


        cmd.setOptional("distribout", "cond1", "cond2");
        cmd.setSwitches("interactive", "normdetails");
        cmd.setInt("getbestsubN");
        cmd.setDefault("getbestsubN", "-1");


        if(!OptionParser.parseParams(args, true, false, true, true, cmd, backgroundProviderOption, generalOptions))
            return;



        generalOptions.apply();

        File inputDir = cmd.getFile("inputdir");
        HashMap<String, empires.input.ReplicateSetInfo> data = readDataGroupedByCondition(inputDir);

        if(cmd.isSet("normdetails")) {
            PlotCreator pc = CachedPlotCreator.getPlotCreator();

            for(String cond : data.keySet())  {
                pc.reset();
                Vector<Vector<Double>> logdata = data.get(cond).getLog2Data();
                Vector<String> repnames = data.get(cond).replicateNames;
                Vector<BufferedImage> images = new Vector<>();
                System.out.printf("start plotting for cond: %s\n",  cond);
                for(int i=0; i<repnames.size(); i++) {
                    for(int j=i+1; j<repnames.size(); j++) {
                        Vector<Double> v1 = logdata.get(i);
                        Vector<Double> v2 = logdata.get(j);

                        System.out.printf("%s %d,%d errd\n", cond, i, j);
                        empires.ErrorEstimationDistribution errD = empires.PairwiseMedianImpliedFCPeakErrorEstimation.getShiftError(v1, v2);

                        errD.getCumulative(true);
                        System.out.printf("%s %d,%d peak\n", cond, i, j);
                        empires.ErrorEstimationDistribution.Peak peak = errD.getBestFCPeak();

                        System.out.printf("%s %d,%d images\n", cond, i, j);
                        Vector<BufferedImage> rowImages = new Vector<>();
                        for(int c = 0; c < 2 ; c++) {

                            System.out.printf("%s %d,%d get line %d\n", cond, i, j, c);
                            Vector<UPair<Double>> line = (c == 0) ? errD.getCumulativeLine() : errD.getNonCumulativeLine();

                            System.out.printf("%s %d,%d line %d: %d\n", cond, i, j, c, line.size());
                            pc.line("", line, (_p) -> _p.getFirst(), (_p) -> _p.getSecond());
                            pc.setTitle("%s vs %s shift: %.2f", repnames.get(i), repnames.get(j), peak.summit * empires.ErrorEstimationDistribution.norm + errD.getMinFC());
                            pc.setLabels("log2fc", "freq", null);

                            pc.abline("", peak.peakStart * empires.ErrorEstimationDistribution.norm + errD.getMinFC(), null, null, null);
                            pc.abline("", peak.summit * empires.ErrorEstimationDistribution.norm + errD.getMinFC(), null, null, null);
                            pc.abline("", peak.peakEnd * empires.ErrorEstimationDistribution.norm + errD.getMinFC(), null, null, null);

                            System.out.printf("will plot line:\n");
                            rowImages.add(pc.getImage());
                            System.out.printf("line plotted!\n");
                        }

                        images.add(ImageUtils.concat(rowImages));
                    }

                }
                System.out.printf("start normalizing for cond: %s\n",  cond);
                empires.NormalizedReplicateSet nrs = new empires.NormalizedReplicateSet(data.get(cond));
                empires.plotting.NormalizedReplicateSetPlotting nrsp = new NormalizedReplicateSetPlotting(nrs);

                ImageUtils.showImage(cond, ImageUtils.vconcat(images));
                ImageUtils.showImage(cond + " norm steps", nrsp.plotNormalizationSteps(CachedPlotCreator.getPlotCreator()));
            }

        }
        int subN = cmd.getInt("getbestsubN");
        if(subN > 0) {
            Tuple3<String, Double, Vector<Integer>> bestSelection = Tuple3.create("", Double.POSITIVE_INFINITY, null);
            for(String cond : data.keySet()) {
                Vector<Vector<Double>> logdata = data.get(cond).getLog2Data();
                if(logdata.size() < subN)
                    continue;
                Vector<Integer> bestSub = empires.Normalization.getBestSubSample(logdata, subN);

                Vector<Vector<Double>> sublogdata = map(bestSub, (_i) -> logdata.get(_i));
                empires.Normalization subNorm = new empires.Normalization(sublogdata);
                int anchor = subNorm.getAnchorSampleIdx();
                System.out.printf("cond: %s\n", cond);
                Vector<Double> sds = new Vector<>();
                for(int i=0; i<sublogdata.size(); i++) {
                    if(i == anchor)
                        continue;

                    empires.ErrorEstimationDistribution err =  empires.PairwiseMedianImpliedFCPeakErrorEstimation.getShiftError(sublogdata.get(i), sublogdata.get(anchor));
                    sds.add(err.getSD());
                    System.out.printf("\terr: %d (anchor: %d) sd: %.2f\n", i, anchor, err.getSD());

                }
                double maxsd = NumUtils.max(sds);

                System.out.printf("%s: %s\n", cond, NumUtils.getNumInfo(sds));
                if(maxsd > bestSelection.get1())
                    continue;

                bestSelection = Tuple3.create(cond, maxsd, bestSub);
            }

            if(bestSelection.get2() == null) {
                System.out.printf("no condition with enough (%d) replicates\n", subN);

                return;
            }

            PrintWriter pw = cmd.getWriter("o");
            int nrep1 = subN >> 1;
            empires.input.ReplicateSetInfo selRepMes = data.get(bestSelection.get0());

            pw.printf("gene");
            Vector<Integer> selIndeces = shuffle(bestSelection.get2(), true);
            int idx = 1;
            for(int i=0; i<selIndeces.size(); i++, idx++) {
                idx = (i == nrep1) ? 1 : idx;
                pw.printf("\tcond%d.rep%d", (i < nrep1) ? 1 : 2, idx);
            }
            pw.println();
            for(String fn : selRepMes.getFeatureNames()) {
                Vector<Double> ldata = selRepMes.getReplicateData(fn);
                if(filteredSize(ldata, (_d) -> !Double.isNaN(_d)) == 0)
                    continue;
                Vector<Integer> intdata = map(map(selIndeces, (_i) -> ldata.get(_i)), (_d) -> Double.isNaN(_d) ? 0 : (int)Math.pow(2.0, _d));
                if(NumUtils.sum(intdata) < subN)
                    continue;

                pw.printf("%s\t%s\n", fn, StringUtils.joinObjects("\t", intdata));
            }
            pw.close();

            return;
        }

        Vector<String> condvec = toVector(data.keySet());
        if(condvec.size() < 2) {
            System.err.printf("too few conditions!\n");
            return;
        }
        String cond1 = cmd.getOptionalValue("cond1", condvec.get(0));
        String cond2 = cmd.getOptionalValue("cond2", condvec.get(1));

        empires.input.ReplicateSetInfo rs1 = data.get(cond1);
        ReplicateSetInfo rs2 = data.get(cond2);

        if(rs1 == null || rs2 == null) {
            Vector<String> unknown = new Vector<>();
            if(rs1 == null) {
                unknown.add(cond1);
            }
            if(rs2 == null) {
                unknown.add(cond2);
            }
            System.err.printf("error, unknown condition(s) requested: %s known conditions: %s\n", unknown, data.keySet());
            return;
        }

        empires.BackgroundContextFuzzficationStrategyProvider backgroundContextFuzzficationStrategyProvider = backgroundProviderOption.getStrategy();
        empires.NormalizedReplicateSet nrs1 = new empires.NormalizedReplicateSet(rs1, backgroundContextFuzzficationStrategyProvider);
        empires.NormalizedReplicateSet nrs2 = new empires.NormalizedReplicateSet(rs2, backgroundContextFuzzficationStrategyProvider);
        Vector<empires.DiffExpResult> diffExpResults = new EmpiRe().getDifferentialResults(nrs1, nrs2);
        empires.plotting.DiffExpTable diffExpTable = new DiffExpTable(nrs1, nrs2, diffExpResults, cmd.getDouble("fcthreshold"), cmd.getDouble("fdrthreshold"));


        diffExpTable.getTable().writeCSV(cmd.getFile("o"));

        File distribout = cmd.getOptionalFile("distribout");
        if(distribout != null) {
            new empires.NamedFoldChangeDistributionMap(diffExpResults).serialize(distribout);
        }

        File trues = EBROWSER_INPUT.TRUES.get(inputDir);
        HashSet<String> truesIds = (!trues.exists()) ? new HashSet<>() : FileUtils.readSet(trues);
        System.out.printf("GOT %d trues!\n", truesIds.size());
        if(truesIds.size() > 0) {
            System.out.println(empires.DiffExpResult.getPerformanceString(diffExpResults, truesIds));
        }

        if(cmd.isSet("interactive")) {
            Vector<Pair<String, Function<empires.DiffExpResult, Object>>> additionalHeaders = new Vector<>();
            if(truesIds.size() > 0) {
                additionalHeaders.add(Pair.create("isTrue", (_d) -> truesIds.contains(_d.combinedFeatureName)));
            }
            diffExpTable.showInteractiveTable(additionalHeaders);
        }
    }
}
