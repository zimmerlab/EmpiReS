package empires.release;

import empires.EmpiRe;
import empires.input.DSType;
import lmu.utils.*;
import lmu.utils.fdr.PerformanceResult;
import lmu.utils.fdr.RocInfo;
import lmu.utils.plotting.CachedPlotCreator;
import lmu.utils.plotting.PlotCreator;
import empires.rnaseq.SplicingTest;

import java.io.File;
import java.util.*;

import static lmu.utils.ObjectGetter.*;

public class EQCInput {

    Logger log = LogConfig.getLogger();
    empires.input.RNASeqSplicingInfo splicingInfo;
    empires.input.ExperimentDescriptor experimentDescriptor;
    static empires.input.DSType targetType = empires.input.DSType.READS;

    Vector<Double> maxMinCounts = new Vector<>();
    int minCountIntMaxCountCondition = 5;

    public EQCInput(empires.input.RNASeqSplicingInfo splicingInfo, empires.input.ExperimentDescriptor experimentDescriptor) {
        this.splicingInfo = splicingInfo;
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
        empires.input.DSType targetType;

        EQInfo(Vector<String[]> data, HashMap<Integer, String> idx2cond, empires.input.DSType targetType) {
            String[] header = data.get(0);
            this.targetType = targetType;
            gene = header[0].substring(1);
            if(! (total = TOTAL.equalsIgnoreCase(header[1]))) {
                eqClass = Tuple.tupleFromCollection(toSortedVector(toVector(header[1].split(",")), true));
            } else {
                eqClass = null;
            }

            for(int i=1; i<data.size(); i++) {
                converter.set(data.get(i));
                empires.input.DSType dsType = DSType.get(converter.toString(0));
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

    HashMap<String, empires.NormalizedReplicateSet> cond2normedTotal = new HashMap<>();

    public empires.NormalizedReplicateSet getNormedCondition(String condition) {
        empires.NormalizedReplicateSet nrs = cond2normedTotal.get(condition);
        if(nrs != null)
            return nrs;


        Vector<String> replicates = experimentDescriptor.getCond2Reps().get(condition);
        if(replicates == null)
            throw new FRuntimeException("unknown condition: " + condition);

        Vector<String> genes = new Vector<>();
        Vector<Vector<Double>> logvals = map(replicates, (_r) -> new Vector<>());

        for(String g : gene2totalcountscond2values.keySet()) {
            Vector<Double> vals = map(gene2totalcountscond2values.get(g).get(condition), (_d) -> (_d == 0.0) ? Double.NaN : NumUtils.logN(_d, 2.0));
            if(0 == filteredSize(vals, (_d) -> !Double.isNaN(_d)))
                continue;

            genes.add(g);
            applyIndex(vals.size(), (_i) -> logvals.get(_i).add(vals.get(_i)));
        }

        empires.input.ReplicateSetInfo rsi = new empires.input.ReplicateSetInfo(condition, replicates, genes);
        for(int i=0; i<replicates.size(); i++) {
            rsi.setLog2Data(i, logvals.get(i));
        }

        cond2normedTotal.put(condition, nrs = new empires.NormalizedReplicateSet(rsi));
        return nrs;
    }

    HashMap<String, HashMap<String, Vector<Double>>> gene2totalcountscond2values = new HashMap<>();

    public void read(File eqInfo) {
        BlockIterator<String[]>  entryIt = new BlockIterator<>(FileUtils.getFieldSetIterator(eqInfo, "\t"),
                (_block, _next_elem) ->  _block.size() > 0 && _next_elem.length > 0 && _next_elem[0].length() >0  && _next_elem[0].charAt(0) == '>'
        );

        Iterator<EQInfo> eqInfoIterator = getTransformedIterator(entryIt, (_e) -> new EQInfo(_e, experimentDescriptor.getReplicateIndex2Condition(), targetType));

        BlockIterator<EQInfo> geneIterator = new BlockIterator<>(eqInfoIterator, (_b, _e) -> _b.size() > 0 && !_e.gene.equals(_b.get(0).gene));

        int trueFalseCount[] = new int[2];
        int nadded = 0;



        while(geneIterator.hasNext()) {
            Vector<EQInfo> geneInfo = geneIterator.next();

            EQInfo totalCount = filterOne(geneInfo, (_eq) -> _eq.total);
            if(totalCount == null)
                throw new FRuntimeException("did not found total count for gene: %s", first(geneInfo).gene);

            gene2totalcountscond2values.put(totalCount.gene, totalCount.cond2values);

            HashMap<String, Vector<HashMap<Tuple, Double>>> condition2replicate2EQclassCounts = new HashMap<>();

            Vector<EQInfo> geneInfoWoTotal = filter(geneInfo, (_g) -> !_g.total);

            EQInfo lead = geneInfoWoTotal.firstElement();

            for(Map.Entry<String, Vector<Double>> e : lead.cond2values.entrySet()) {
                condition2replicate2EQclassCounts.put(e.getKey(), map(e.getValue(), (_r) -> new HashMap<>()));
            }
            for(EQInfo eq : geneInfoWoTotal) {
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


            HashMap<String, HashMap<String, Vector<Double>>> reducedCounts = splicingInfo.addGeneInfo(lead.gene, condition2replicate2EQclassCounts);
            for(HashMap<String, Vector<Double>> fc : reducedCounts.values()) {
                maxMinCounts.add(NumUtils.max(fc.values(), (_v) -> NumUtils.min(_v)));
            }
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
        empires.input.GeneralOptions generalOptions = new empires.input.GeneralOptions();
        empires.input.BackgroundProviderOption backgroundProviderOption = new empires.input.BackgroundProviderOption();
        SimpleOptionParser cmd = new SimpleOptionParser("i", "samples", "maxtrnum", "truesplicing", "truediffexp", "cond1", "cond2", "o");
        cmd.setFile("i", "samples", "truesplicing", "truediffexp");
        cmd.setInt("maxtrnum");
        cmd.setDefault("maxtrnum", "10");
        cmd.setOptional("truesplicing", "truediffexp", "cond1", "cond2");
        if(!OptionParser.parseParams(args, true, false, true, true, cmd))
            return;


        generalOptions.apply();
        double RAW_PSEUDO = 2.0;
        Logger log = LogConfig.getLogger();
        int MAX_TR_NUM = cmd.getInt("maxtrnum");

        empires.input.ExperimentDescriptor experimentDescriptor = new empires.input.ExperimentDescriptor(null, cmd.getFile("samples"), cmd.getOptionalFile("trues"));

        empires.input.RNASeqSplicingInfo splicingInfo = new empires.input.RNASeqSplicingInfo(MAX_TR_NUM, 1.0,  experimentDescriptor.getCond2Reps());

        EQCInput eqClassInput = new EQCInput(splicingInfo, experimentDescriptor);

        eqClassInput.read(cmd.getFile("i"));

        Vector<String> conditions = toVector(mapToSet(splicingInfo.getReplicateSetInfos(), (_r) -> _r.getReplicateSetName()));

        String cond1 = cmd.getOptionalValue("cond1", conditions.get(0));
        String cond2 = cmd.getOptionalValue("cond2", conditions.get(1));

        for(String cond : new String[]{cond1, cond2}){
            if(null != filterOne(conditions, (_c) -> _c.equals(cond)))
                continue;

            System.err.printf("invalid condition: >%s< got: %s\n", cond, conditions);
            return;

        }

        SplicingTest test = new SplicingTest(splicingInfo, backgroundProviderOption.getStrategy());
        test.doubleDiffVariant = empires.DoubleDiffVariant.ALLPAIRS;



        HashSet<String> trueDiffexp = FileUtils.readSet(cmd.getOptionalFile("truediffexp"));
        HashSet<String> trueSplic = FileUtils.readSet(cmd.getOptionalFile("truesplicing"));

        log.info("calc diffexp conditions: %s to test: %s, %s", conditions, cond1, cond2);

        HashMap<String, empires.DiffExpResult> gene2diffexp = buildReverseMap(new EmpiRe().getDifferentialResults(eqClassInput.getNormedCondition(cond1), eqClassInput.getNormedCondition(cond2)), (_e) -> _e.combinedFeatureName);

        if(trueDiffexp != null) {
            System.out.printf("%s\n", new PerformanceResult("EmpiRe-diffexp", filter(gene2diffexp.values(), (_d) -> (trueSplic == null) ? true : !trueSplic.contains(_d.combinedFeatureName)), (_d) -> _d.pval, (_d) -> trueDiffexp.contains(_d.combinedFeatureName), false, (_d) -> _d.fdr <= 0.05,
                    null, false, RocInfo.PerformanceEvaluationStrategy.OPTIMISTIC));
        }
        log.info("calc splicing conditions: %s to test: %s, %s", conditions, cond1, cond2);
        //test.getQuickTest(conditions.get(0), conditions.get(1), trueGenes);
        long t1 = System.currentTimeMillis();
//        UPair<Vector<DoubleDiffResult>> allSplicingResults = test.getDifferentialAlternativeSplicing(cond1, cond2);

//        HashMap<String, DoubleDiffResult> splicingFull = buildReverseMap(allSplicingResults.getFirst(), (_ddr) -> _ddr.testName);

        HashMap<String, empires.DoubleDiffResult> splicing = new HashMap<>();
//        splicing = buildReverseMap(allSplicingResults.getFirst(), (_ddr) -> _ddr.testName.split("\\.")[0]);


        long t2 = System.currentTimeMillis();
        long DAS_test_time = t2 - t1;
        log.info("DAS test time: %.2f sec.", DAS_test_time / 1000.0);


//        if(trueSplic != null) {
//            System.out.printf("%s\n", new PerformanceResult("EmpiRe-diffsplic", toVector(splicing.values()), (_d) -> _d.pval, (_d) -> trueSplic.contains(_d.testName.split("\\.")[0]), false, (_d) -> _d.fdr <= 0.05,
//                        null, false, RocInfo.PerformanceEvaluationStrategy.OPTIMISTIC));
//
//
//        }
        empires.DoubleDiffResult MISSING_DDR = new empires.DoubleDiffResult();
        Vector<empires.DoubleDiffResult> MISSING_VDDR = new Vector<>();
        MISSING_VDDR.add(MISSING_DDR);
        Vector<String> genelist = NumUtils.sort(toVector(gene2diffexp.keySet()), (_g) -> Math.min(gene2diffexp.get(_g).fdr, splicing.getOrDefault(_g, MISSING_DDR).fdr), false);



        DataTable.HeaderGetterManager<String> hgm = new DataTable.HeaderGetterManager<>("gene", (_g) -> _g);
        if(trueDiffexp != null) {
            hgm.add("true.diffexp", (_g) -> trueDiffexp.contains(_g));
        }
        if(trueSplic != null) {
            hgm.add("true.splicing", (_g) -> trueSplic.contains(_g));
        }
        hgm.add("diffexp.fdr", (_g) -> gene2diffexp.get(_g).fdr)
                .add("diffexp.log2fc", (_g) -> gene2diffexp.get(_g).estimatedFC)
                .add("diffsplic.most.signif.test", (_g) -> splicing.getOrDefault(_g, MISSING_DDR).testName)
                .add("diffsplic.fdr", (_g) -> splicing.getOrDefault(_g, MISSING_DDR).fdr)
                .add("diffsplic.difflog2fc", (_g) -> splicing.getOrDefault(_g, MISSING_DDR).meanFC);


        DataTable.buildTable(genelist, hgm).writeCSV(cmd.getFile("o"));


        HashMap<String, String> typeMap = new HashMap<>();
        UPair<Vector<empires.DoubleDiffResult>> diffTransResults = test.getDifferentialTranscriptExpression(cond1, cond2);
        HashMap<String, empires.DoubleDiffResult> transcript2diffexp = buildReverseMap(diffTransResults.getFirst(), (_ddr) -> _ddr.testName.split("\\.")[0]);
//        HashMap<String, DoubleDiffResult> transcript2diffexp = new HashMap<>();
        genelist.forEach(_g -> transcript2diffexp.put(_g, new empires.DoubleDiffResult()));

        final double CUTOFF = 0.051; //fixme this should be param
        genelist.forEach(_g -> {
            if (gene2diffexp.get(_g).fdr <= CUTOFF && splicing.getOrDefault(_g, MISSING_DDR).testName.equals("untested")) {
                typeMap.put(_g, "A1:single_transcript_DE");
            } else if (gene2diffexp.get(_g).fdr <= CUTOFF && splicing.getOrDefault(_g, MISSING_DDR).fdr > CUTOFF) {
                typeMap.put(_g, "A2:similar_for_all_trs");
            } else if (transcript2diffexp.get(_g).fdr <= CUTOFF) {  //here look if 1 transcript is DE
                if(gene2diffexp.get(_g).fdr <= CUTOFF) {
                    typeMap.put(_g, "B1:big_change_for_tr");
                } else if (splicing.getOrDefault(_g, MISSING_DDR).fdr <= CUTOFF) {
                    typeMap.put(_g, "C1:minor_change_in_tr");
                } else {
                    typeMap.put(_g, "UNKNOWN:not_DAS_not_DE_but_difftrans");
                }
            } else if (splicing.getOrDefault(_g, MISSING_DDR).fdr <= CUTOFF) {
                if (gene2diffexp.get(_g).fdr <= CUTOFF) {
                    typeMap.put(_g, "B2:change_in_composition_and_overall");
                } else {
                    typeMap.put(_g, "C2:minor_change_in_composition");
                }
            }
        });



        DataTable.HeaderGetterManager<String> hgm2 = new DataTable.HeaderGetterManager<>("gene", (_g) -> _g);
        if(trueDiffexp != null) {
            hgm2.add("true.diffexp", (_g) -> trueDiffexp.contains(_g));
        }
        if(trueSplic != null) {
            hgm2.add("true.splicing", (_g) -> trueSplic.contains(_g));
        }
        hgm2.add("diffType", typeMap::get)
                .add("diffexp.fdr", (_g) -> gene2diffexp.get(_g).fdr)
                .add("diffexp.log2fc", (_g) -> gene2diffexp.get(_g).estimatedFC)
                .add("diffsplic.most.signif.test", (_g) -> splicing.getOrDefault(_g, MISSING_DDR).testName)
                .add("diffsplic.fdr", (_g) -> splicing.getOrDefault(_g, MISSING_DDR).fdr)
                .add("diffsplic.difflog2fc", (_g) -> splicing.getOrDefault(_g, MISSING_DDR).meanFC)
                .add("diffsplic.f1.difflog2fc("+cond1+"-"+cond2+")", (_g) -> splicing.getOrDefault(_g, MISSING_DDR).f1meanFC)
                .add("diffsplic.f2.difflog2fc("+cond1+"-"+cond2+")", (_g) -> splicing.getOrDefault(_g, MISSING_DDR).f2meanFC);


        DataTable.buildTable(genelist, hgm2).writeCSV(new File(cmd.getFile("o").getAbsolutePath() + ".complete.info"));


//        DataTable.HeaderGetterManager<String> hgm3 = new DataTable.HeaderGetterManager<>("gene", (_g) -> _g.split("\\.")[0]);
//        hgm3.add("diffsplic.most.signif.test", (_g) -> splicingFull.getOrDefault(_g, MISSING_DDR).testName)
//                .add("diffsplic.fdr", (_g) -> splicingFull.getOrDefault(_g, MISSING_DDR).fdr)
//                .add("diffsplic.difflog2fc", (_g) -> splicingFull.getOrDefault(_g, MISSING_DDR).meanFC)
//                .add("diffsplic.f1.difflog2fc("+cond1+"-"+cond2+")", (_g) -> splicingFull.getOrDefault(_g, MISSING_DDR).f1meanFC)
//                .add("diffsplic.f2.difflog2fc("+cond1+"-"+cond2+")", (_g) -> splicingFull.getOrDefault(_g, MISSING_DDR).f2meanFC);
//
//        DataTable.buildTable(splicingFull.keySet(), hgm3).writeCSV(new File(cmd.getFile("o").getAbsolutePath() + ".splicing.full.info"));
    }
}

