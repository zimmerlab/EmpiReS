package nlEmpiRe.input;

import lmu.utils.*;
import lmu.utils.fdr.PerformanceResult;
import lmu.utils.plotting.PlotCreator;
import lmu.utils.swing.PagedDataTable;
import nlEmpiRe.*;
import lmu.utils.plotting.CachedPlotCreator;
import nlEmpiRe.plotting.DiffExpTable;
import nlEmpiRe.rnaseq.SplicingTest;
import nlEmpiRe.test.rnaseq.BenchmarkGene;

import java.io.File;
import java.util.*;
import java.util.function.Function;

import static lmu.utils.ObjectGetter.*;
import org.apache.logging.log4j.Logger;

public class EQClassInput {

    Logger log = LogConfig.getLogger();
    RNASeqSplicingInfo splicingInfo;
    ExperimentDescriptor experimentDescriptor;
    DSType targetType;

    int limitInput = -1;
    Vector<Double> maxMinCounts = new Vector<>();
    int minCountIntMaxCountCondition;

    public EQClassInput(RNASeqSplicingInfo splicingInfo, ExperimentDescriptor experimentDescriptor,  DSType targetType,
                        int minCountIntMaxCountCondition) {
        this.splicingInfo = splicingInfo;
        this.targetType = targetType;
        this.experimentDescriptor = experimentDescriptor;
        this.minCountIntMaxCountCondition = minCountIntMaxCountCondition;
        splicingInfo.setMinCountIntMaxCountCondition(minCountIntMaxCountCondition);

    }

    static FileUtils.Converter converter = new FileUtils.Converter();
    static final String TOTAL = "TOTAL";
    static class EQInfo {
        String gene;
        Tuple eqClass;
        boolean total;
        HashMap<String, Vector<Double>> cond2values = new HashMap<>();
        DSType targetType;

        EQInfo(Vector<String[]> data, HashMap<Integer, String> idx2cond, DSType targetType) {
            String[] header = data.get(0);
            this.targetType = targetType;
            gene = header[0].substring(1);
            if(! (total = TOTAL.equalsIgnoreCase(header[1]))) {
                eqClass = Tuple.tupleFromCollection(toSortedVector(toVector(header[1].split(",")), true));
            }

            for(int i=1; i<data.size(); i++) {
                converter.set(data.get(i));
                DSType dsType = DSType.get(converter.toString(0));
                if(dsType != targetType)
                    continue;


                for(int j=2, repidx = 0; j < converter.length(); j++, repidx++) {
                    MapBuilder.updateV(cond2values, idx2cond.get(repidx), converter.toDbl(j));
                }

            }
        }

        @Override
        public String toString() {
            return String.format("%s:%s :: %s", gene, eqClass, cond2values);
        }
    }

    public void read(File eqInfo) {
        BlockIterator<String[]>  entryIt = new BlockIterator<>(FileUtils.getFieldSetIterator(eqInfo, "\t"),
                (_block, _next_elem) ->  _block.size() > 0 && _next_elem.length > 0 && _next_elem[0].length() >0  && _next_elem[0].charAt(0) == '>'
        );

        Iterator<EQInfo> eqInfoIterator = filterIterator(getTransformedIterator(entryIt, (_e) -> new EQInfo(_e, experimentDescriptor.getReplicateIndex2Condition(), targetType)), (_e) -> !_e.total);

        BlockIterator<EQInfo> geneIterator = new BlockIterator<>(eqInfoIterator, (_b, _e) -> _b.size() > 0 && !_e.gene.equals(_b.get(0).gene));

        int trueFalseCount[] = new int[2];
        int nadded = 0;



        while(geneIterator.hasNext()) {
            Vector<EQInfo> geneInfo = geneIterator.next();

            HashMap<String, Vector<HashMap<Tuple, Double>>> condition2replicate2EQclassCounts = new HashMap<>();

            EQInfo lead = geneInfo.firstElement();

            for(Map.Entry<String, Vector<Double>> e : lead.cond2values.entrySet()) {
                condition2replicate2EQclassCounts.put(e.getKey(), map(e.getValue(), (_r) -> new HashMap<>()));
            }
            for(EQInfo eq : geneInfo) {
                for(Map.Entry<String, Vector<Double>> e : eq.cond2values.entrySet()) {
                    String cond = e.getKey();
                    Vector<HashMap<Tuple, Double>> target = condition2replicate2EQclassCounts.get(cond);
                    Vector<Double> vals = e.getValue();
                    for(int i=0; i<vals.size(); i++) {
                        target.get(i).put(eq.eqClass, vals.get(i));
                    }
                }
            }

            int trueIdx = experimentDescriptor.isTrue(lead.gene) ? 0 : 1;
            trueFalseCount[trueIdx]++;


            int total = trueFalseCount[0] + trueFalseCount[1];
            if(total % 500 == 0) {
                log.info("read %d trues %d falses limit: %d current: %s so far saved: %d", trueFalseCount[0], trueFalseCount[1], limitInput, lead.gene, nadded);
            }

            if(limitInput < 0 || trueFalseCount[trueIdx] < limitInput) {
                HashMap<String, HashMap<String, Vector<Double>>> reducedCounts = splicingInfo.addGeneInfo(lead.gene, condition2replicate2EQclassCounts);
                for(HashMap<String, Vector<Double>> fc : reducedCounts.values()) {
                    maxMinCounts.add(NumUtils.max(fc.values(), (_v) -> NumUtils.min(_v)));
                }

            }

            if(limitInput > 0 && trueFalseCount[0] >= limitInput && trueFalseCount[1] >= limitInput)
                break;
        }



    }

    public void plotCountDistribution() {
        PlotCreator pc = CachedPlotCreator.getPlotCreator();
        pc.cumhist("", maxMinCounts, maxMinCounts.size(), false, false);
        pc.setLabels("max (per condition) min (in replicates) counts per feature", "num features with <=X counts", null);
        pc.setLog(true, true);
        ImageUtils.showImage("counts", pc.getImage(), false);
    }

    public static void main(String[] args) {
        GeneralOptions generalOptions = new GeneralOptions();
        BackgroundProviderOption backgroundProviderOption = new BackgroundProviderOption();
        SimpleOptionParser cmd = new SimpleOptionParser("eqclasscounts", "samples", "cond2reps", "maxtrnum", "pseudo", "trues", "dstype",
                "showtable", "unmatched","testlimit", "minc", "minctest", "showplots", "cond1", "cond2", "o", "diffexpout", "empireVariant",
                "truediffexp");
        cmd.setFile("eqclasscounts", "samples", "cond2reps", "trues");
        cmd.setInt("maxtrnum", "testlimit", "minc", "minctest");
        cmd.setDefault("maxtrnum", "10");
        cmd.setDefault("testlimit", "-1");
        cmd.setDefault("minc", "5");
        cmd.setDefault("minctest", "5");
        cmd.setOptional("trues", "cond1", "cond2", "o", "truediffexp", "diffexpout", "cond2reps");
        cmd.setDouble("pseudo");

        cmd.setDefault("pseudo", "2.0");
        cmd.setDefault("dstype", DSType.READS.name);
        cmd.setSwitches("showtable", "unmatched", "showplots");
        cmd.setDefault("empireVariant", ""+ DoubleDiffVariant.QUICK_AND_DIRTY.getName());

        if(!OptionParser.parseParams(args, true, false, true, true, cmd, backgroundProviderOption, generalOptions))
            return;


        generalOptions.apply();
        double RAW_PSEUDO = cmd.getDouble("pseudo");
        Logger log = LogConfig.getLogger();
        boolean UNMATCHED = cmd.isSet("unmatched");
        boolean showplots = cmd.isSet("showplots");
        int MAX_TR_NUM = cmd.getInt("maxtrnum");

        DSType targetDSType = DSType.get(cmd.getValue("dstype"));
        ExperimentDescriptor experimentDescriptor = new ExperimentDescriptor(cmd.getOptionalFile("cond2reps"), cmd.getFile("samples"), cmd.getOptionalFile("trues"));

        RNASeqSplicingInfo splicingInfo = new RNASeqSplicingInfo(MAX_TR_NUM, 1.0,  experimentDescriptor.getCond2Reps());

        int minCountIntMaxCountCondition = cmd.getInt("minc");
        EQClassInput eqClassInput = new EQClassInput(splicingInfo, experimentDescriptor, targetDSType, minCountIntMaxCountCondition);


        eqClassInput.limitInput = cmd.getInt("testlimit");
        eqClassInput.read(cmd.getFile("eqclasscounts"));

        if(showplots) {
            eqClassInput.plotCountDistribution();
        }




        Vector<String> conditions = toVector(mapToSet(splicingInfo.getReplicateSetInfos(), (_r) -> _r.getReplicateSetName()));

        String cond1 = cmd.getOptionalValue("cond1", conditions.get(0));
        String cond2 = cmd.getOptionalValue("cond2", conditions.get(1));

        for(String cond : new String[]{cond1, cond2}){
            if(null != filterOne(conditions, (_c) -> _c.equals(cond)))
                continue;

            System.err.printf("invalid condition: >%s< got: %s\n", cond, conditions);
            return;

        }

        if(minCountIntMaxCountCondition < 0) {
            int targetPercentage = 99;
            minCountIntMaxCountCondition = Math.min(splicingInfo.getMinCountLevelForCoveragePercentage(cond1, targetPercentage),
                                                    splicingInfo.getMinCountLevelForCoveragePercentage(cond2, targetPercentage));
            splicingInfo = splicingInfo.restrictToMinCountIntMaxCountCondition(minCountIntMaxCountCondition);
        }

        SplicingTest test = new SplicingTest(splicingInfo, backgroundProviderOption.getStrategy());
        test.doubleDiffVariant = DoubleDiffVariant.get(cmd.getValue("empireVariant"));

        File diffexpOut = cmd.getOptionalFile("diffexpout");

        HashMap<Boolean, HashMap<String, DiffExpResult>> diffexpResults = new HashMap<>();

        HashSet<String> trueDiffexp = (cmd.isOptionSet("truediffexp")) ? FileUtils.readSet(cmd.getFile("truediffexp")) : null;

        Vector<Double> Q2 = toVector(0.01, 0.02, 0.05, 0.1, 0.15, 0.2);
        if(showplots) {
            PlotCreator pc = CachedPlotCreator.getPlotCreator();
            for(String cond : toVector(cond1, cond2)) {
                HashMap<Integer, Vector<Double>> m = test.getEqClassLevel2Count(cond);
                for(int i=0; i<5; i++) {
                    Vector<Double> vals = m.get(i);
                    if(vals == null)
                        continue;
                    Collections.sort(vals);
                    Function<Integer, Double> qF = (_i) -> vals.get((int)((_i / 100.0) * vals.size()));

                    System.out.printf("%s %d at 5%% = %.2f at 10%% = %.2f 20%% = %.2f  25%% = %.2f\n", cond, i+1,
                            qF.apply(5), qF.apply(10), qF.apply(20), qF.apply(25));

                    pc.cumhist(String.format("%s=%d (%d)", cond, i + 1, vals.size()), vals, vals.size(), false, true);
                }

            }
            pc.setLog(true, false);
            pc.setLabels("count", "freq with count <= x", "bottomright");
            ImageUtils.showImage("counts", pc.getImage(), false);

            for(String cond : toVector(cond1, cond2)) {
                Vector<Vector<UPair<Double>>> reps = test.getCountLevel2CoveragePerReplicate(cond);
                for(int i=0; i<reps.size(); i++) {
                    pc.line(String.format("%s:%d", cond, i + 1), reps.get(i), (_p) -> _p.getFirst(), (_p) -> _p.getSecond());
                }

            }

            pc.setLog(true, false);
            pc.setLabels("count", "covered with count >= x", "bottomright");
            ImageUtils.showImage("coverages", pc.getImage(), false);


        }

        if(diffexpOut != null) {
            for(boolean summarize : toVector(false, true)) {
                Vector<DiffExpResult> diff = test.getDiffGenes(cond1, cond2, summarize);
                Vector<DiffExpResult> nosplicdiff = filter(diff, (_d) -> !experimentDescriptor.isTrue(_d.combinedFeatureName));

                if(summarize) {
                    new DiffExpTable(test.getNormalized(cond1), test.getNormalized(cond2), diff).getTable().writeCSV(diffexpOut);
                }


                String prefix = "diffexp summarize:" + summarize+" ";
                System.out.printf("got %d/%d nosplic diffs\n",  nosplicdiff.size(), diff.size());
                diffexpResults.put(summarize, buildReverseMap(nosplicdiff, (_d) -> _d.combinedFeatureName));
                if(trueDiffexp != null) {

                    for(double FC : toVector(0.0, 0.1, 0.25, 0.5)) {
                        PerformanceResult[] results = DiffExpResult.getPerformances(nosplicdiff, trueDiffexp, 0.05, FC);
                        for(PerformanceResult pr: results) {
                            System.out.println(pr.setTitle(prefix + pr.getTitle()));
                        }
                    }


                }
            }

        }

        test.setPseudo(RAW_PSEUDO);
        test.setShowPlots(showplots);

        test.setMinReadCountTresholdInAllReplicatesInAtLeastOneCondition(cmd.getInt("minctest"));




        log.info("calc splicing conditions: %s to test: %s, %s", conditions, cond1, cond2);
        //test.getQuickTest(conditions.get(0), conditions.get(1), trueGenes);
        long t1 = System.currentTimeMillis();
        UPair<Vector<DoubleDiffResult>> splicing = test.getDifferentialAlternativeSplicing(cond1, cond2);
        long t2 = System.currentTimeMillis();
        long DAS_test_time = t2 - t1;
        log.info("DAS test time: %.2f sec.", DAS_test_time / 1000.0);

        PerformanceResult pr = null;
        PerformanceResult prUP = null;

        if(experimentDescriptor.gotTrues()) {

            Function<DoubleDiffResult, Boolean> trueLabeller = (_dr) -> experimentDescriptor.isTrue(_dr.testName.split("\\.")[0]);

            String name = String.format("nlEmpiRe.%s.%d.%d", test.doubleDiffVariant.getName(), (int)cmd.getDouble("pseudo"), cmd.getInt("minc"));
            pr = new PerformanceResult(name, splicing.getFirst(), (_dr) -> _dr.pval, trueLabeller,
                    false, (_dr) -> _dr.fdr <= 0.05, null, false);


            System.out.printf("mincount filter for eq-class input: %d (-minc) settings: %s\n", minCountIntMaxCountCondition, test.getSettings());
            System.out.printf("%s\n", pr);

            /*
            prUP = new PerformanceResult("nlEmpire-up", splicing.getSecond(), (_dr) -> _dr.unpairedPval, trueLabeller,
                    false, (_dr) -> _dr.unpairedFDR <= 0.05, null, false);


            System.out.printf("%s\n", prUP);
            */


            if(showplots) {
                PlotCreator pc = CachedPlotCreator.getPlotCreator();
                pr.drawAUPR(pc, "nlE");
                //prUP.drawAUPR(pc, "nlE-up");
                ImageUtils.showImage("pr", pc.getImage(), false);
            }


        }
        if(cmd.isSet("showtable")) {
            List<Pair<String, Function<BenchmarkGene, Object>>> additionalHeaders = new ArrayList<>();
            if(trueDiffexp != null) {
                additionalHeaders.add(Pair.create("Tdiffexp", (BenchmarkGene _bg) -> trueDiffexp.contains(_bg.geneId)));
            }


            for(Boolean summarized : diffexpResults.keySet()) {
                HashMap<String, DiffExpResult> lookup = diffexpResults.get(summarized);
                String prefix = (summarized) ? "sum.diff" : "comb.diff";

                Function<String, Double> FDR_GETTER = (_s) -> { DiffExpResult der = lookup.get(_s); return (der == null) ? 2.0 : der.fcEstimateFDR;};
                Function<String, Double> FC_GETTER = (_s) -> { DiffExpResult der = lookup.get(_s); return (der == null) ? 0.0 : der.estimatedFC;};
                additionalHeaders.add(Pair.create(prefix+".diff", (BenchmarkGene _bg) -> FDR_GETTER.apply(_bg.geneId)));
                additionalHeaders.add(Pair.create(prefix+".log2FC", (BenchmarkGene _bg) -> FC_GETTER.apply(_bg.geneId)));

                DiffExpManager dm = first(lookup.values()).getDiffExpManager();
                DiffExpTable difft = new DiffExpTable(dm.replicateSetFrom, dm.replicateSetTo, toVector(lookup.values()), 0.3, 0.05, (_s) -> (trueDiffexp == null) ? false : trueDiffexp.contains(_s));

                PagedDataTable.MJFrame mjFrame = difft.showInteractiveTable(Pair.create("diffsplic", (_dr) -> !(experimentDescriptor.gotTrues()) ? false : experimentDescriptor.trues.contains(_dr.combinedFeatureName)));
                mjFrame.setTitle("summ: + " + summarized + " " + mjFrame.getTitle());
            }
            DataTable dt = BenchmarkGene.getTable(splicing.getFirst(), experimentDescriptor.getTrues(), null, splicingInfo, additionalHeaders);
            BenchmarkGene.getWithDetails(dt).setTitle("paired: " + pr);
            //BenchmarkGene.showTable(splicing.getSecond(), trueGenes, null, splicingInfo).setTitle("unpaired: " + prUP);
        }

        File outfile = cmd.getOptionalFile("o");
        if(outfile != null) {
            BenchmarkGene.getTable(splicing.getFirst(), experimentDescriptor.getTrues(), null, splicingInfo).writeCSV(outfile);
        }




    }
}
