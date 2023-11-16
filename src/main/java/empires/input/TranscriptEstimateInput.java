package empires.input;

import empires.rnaseq.GFFBasedIsoformRegionGetter;
import empires.rnaseq.SplicingTest;
import lmu.utils.*;
import lmu.utils.fdr.PerformanceResult;
import lmu.utils.plotting.PlotCreator;
import empires.DoubleDiffResult;
import lmu.utils.plotting.CachedPlotCreator;
import empires.test.rnaseq.BenchmarkGene;

import java.io.File;
import java.util.*;
import java.util.function.Function;

import static lmu.utils.ObjectGetter.*;

public class TranscriptEstimateInput {

    Logger log = LogConfig.getLogger();
    public static final String ABUNDANCE_FILE_KALLISTO = "abundance.tsv";
    public static final String ABUNDANCE_FILE_SALMON = "quant.sf";
    public static final String ABUNDANCE_FILE_STRINGTIE = "t_data.ctab";

    static String[] ABUNDANCE_FILES = new String[]{ABUNDANCE_FILE_KALLISTO, ABUNDANCE_FILE_SALMON, ABUNDANCE_FILE_STRINGTIE};

    empires.input.ExperimentDescriptor experimentDescriptor;
    HashMap<String, String> transcript2gene;
    empires.input.RNASeqSplicingInfo rsi;

    public TranscriptEstimateInput(File root, empires.input.ExperimentDescriptor experimentDescriptor, HashMap<String, String> transcript2gene) {
        this.experimentDescriptor = experimentDescriptor;
        rsi = new empires.input.RNASeqSplicingInfo(1000, 0, experimentDescriptor.getCond2Reps());
        HashMap<String, HashSet<String>> gene2trs = new HashMap<>();

        HashMap<String, Vector<HashMap<String, Double>>> cond2rep2tr2counts = new HashMap<>();
        for(Map.Entry<String, Vector<String>> e: experimentDescriptor.getCond2Reps().entrySet()) {
            //ReplicateSetInfo rsi = new ReplicateSetInfo("")

            Vector<HashMap<String, Double>> rep2tr2counts = new Vector<>();
            cond2rep2tr2counts.put(e.getKey(), rep2tr2counts);
            String selectedAbundanceFile = null;

            HashMap<String, Vector<String>> errors = new HashMap<>();

            for(String abundanceFile : ABUNDANCE_FILES) {
                Vector<File> files = map(e.getValue(), (_r) -> new File(new File(root, _r), abundanceFile));

                if(filteredSize(files, (_f) -> !_f.exists()) == 0) {
                    selectedAbundanceFile = abundanceFile;
                    break;
                }
                errors.put(abundanceFile, map(files, (_f) -> _f.getAbsolutePath()));

            }
            if(selectedAbundanceFile == null) {
                throw new FRuntimeException("don't find estimate files for condition: %s tested: %s", e.getKey(), errors);
            }
            for(String replicate : e.getValue()) {
                    File repFile = new File(new File(root, replicate), selectedAbundanceFile);
                    rep2tr2counts.add(readTrEstimates(repFile));

            }

            HashSet<String> trs = new HashSet<>();
            apply(rep2tr2counts, (_m) -> trs.addAll(_m.keySet()));
            apply(trs, (_t) -> MapBuilder.update(gene2trs, transcript2gene.get(_t), _t));
        }

        log.info("tr2gene distirb: %s", NumUtils.getNumInfo(gene2trs.values(), (_v) -> _v.size()).getInfoWithQ());


        for(String g : gene2trs.keySet()) {
            HashMap<String, Vector<HashMap<Tuple, Double>>> condition2replicate2EQclassCounts = new HashMap<>();
            HashSet<String> measuredTrs = gene2trs.get(g);
            for(String cond : cond2rep2tr2counts.keySet()) {
                Vector<HashMap<String, Double>> condvals = cond2rep2tr2counts.get(cond);
                Vector<HashMap<Tuple, Double>> eq2counts = map(condvals, (_c) -> new HashMap<>());
                for(String tr : measuredTrs) {
                    Tuple trEQ = new Tuple(tr);
                    applyIndex(condvals.size(), (_idx) -> eq2counts.get(_idx).put(trEQ, condvals.get(_idx).getOrDefault(tr, 0.0)));

                }

                condition2replicate2EQclassCounts.put(cond, eq2counts);
            }
            rsi.addGeneInfo(g, condition2replicate2EQclassCounts);
        }
    }

    static final String[] KALLISTO_EST_HEADERS = new String[]{"target_id", "est_counts"};
    static final String[] SALMON_EST_HEADERS = new String[]{"Name", "NumReads"};
    static final String[] STRINGTIE_EST_HEADERS = new String[]{"t_name", "FPKM"};
    static final String[][] TEST_HEADER_COMBIS = new String[][]{KALLISTO_EST_HEADERS, SALMON_EST_HEADERS, STRINGTIE_EST_HEADERS};

    static HashMap<String, Double> readTrEstimates(File f) {
        HashSet<String> headers = toSet(FileUtils.getHeaders(f));
        String[] touse = null;
        for(String[] totest : TEST_HEADER_COMBIS) {
            if(filteredSize(toVector(totest), (_s) -> headers.contains(_s)) != totest.length)
                continue;

            touse = totest;
            break;
        }
        if(touse == null)
            throw new FRuntimeException("could not guess headers from %s headers: %s checked: %s", f.getAbsolutePath(), headers,
                        map(TEST_HEADER_COMBIS, (_t) -> Arrays.toString(_t)));


        return buildMap(FileUtils.getHrIterator(f, FileUtils.getHeaderedReader("\t").add(touse[0], "id").add(touse[1], "est")),
                (_hr) -> _hr.getString("id"), (_hr) -> _hr.getDbl("est"));
    }

    static HashMap<String, String> getTranscript2Gene(File gtf) {
        empires.rnaseq.GFFBasedIsoformRegionGetter gff = new GFFBasedIsoformRegionGetter(gtf, null, null);
        HashMap<String, String> transcript2gene = new HashMap<>();
        apply(gff.getRegions(), (_m) ->
            apply(_m.isoforms.keySet(), (_tr) -> transcript2gene.put(_tr, _m.id))
        );
        return transcript2gene;
    }
    public static void main(String[] args) {

        empires.input.BackgroundProviderOption backgroundProviderOption = new BackgroundProviderOption();
        SimpleOptionParser cmd = new SimpleOptionParser("cond2reps", "gtf", "trestimateroot", "pseudo", "trues",
                "showtable",  "minc", "minctest", "showplots", "cond1", "cond2", "o");
        cmd.setFile("cond2reps", "trues");
        cmd.setDir("trestimateroot");
        cmd.setInt( "minc", "minctest");

        cmd.setDefault("minc", "0");
        cmd.setDefault("minctest", "0");
        cmd.setOptional("trues", "cond1", "cond2");
        cmd.setDouble("pseudo");
        cmd.setDefault("pseudo", "2.0");
        cmd.setSwitches("showtable", "showplots");

        if(!OptionParser.parseParams(args, true, false, true, true, cmd, backgroundProviderOption))
            return;

        double RAW_PSEUDO = cmd.getDouble("pseudo");
        Logger log = LogConfig.getLogger();
        boolean showplots = cmd.isSet("showplots");

        HashMap<String, String> tr2genes = TranscriptEstimateInput.getTranscript2Gene(cmd.getFile("gtf"));
        empires.input.ExperimentDescriptor experimentDescriptor = new ExperimentDescriptor(cmd.getFile("cond2reps"), null, cmd.getOptionalFile("trues"));

        TranscriptEstimateInput transcriptEstimateInput = new TranscriptEstimateInput(cmd.getFile("trestimateroot"), experimentDescriptor, tr2genes);

        RNASeqSplicingInfo splicingInfo = transcriptEstimateInput.rsi;

        empires.rnaseq.SplicingTest test = new SplicingTest(splicingInfo, backgroundProviderOption.getStrategy());
        test.setPseudo(RAW_PSEUDO);
        test.setShowPlots(showplots);

        test.setMinReadCountTresholdInAllReplicatesInAtLeastOneCondition(cmd.getInt("minctest"));
        Vector<String> conditions = test.getConditions();

        String cond1 = cmd.getOptionalValue("cond1", conditions.get(0));
        String cond2 = cmd.getOptionalValue("cond2", conditions.get(1));

        for(String cond : new String[]{cond1, cond2}){
            if(null != filterOne(conditions, (_c) -> _c.equals(cond)))
                continue;

            System.err.printf("invalid condition: >%s< got: %s\n", cond, conditions);
            return;

        }


        log.info("calc splicing conditions: %s to test: %s, %s", conditions, cond1, cond2);
        //test.getQuickTest(conditions.get(0), conditions.get(1), trueGenes);
        long t1 = System.currentTimeMillis();
        UPair<Vector<DoubleDiffResult>> splicing = test.getDifferentialAlternativeSplicing(cond1, cond2);
        long t2 = System.currentTimeMillis();
        long DAS_test_time = t2 - t1;
        log.info("DAS test time: %.2f sec.", DAS_test_time / 1000.0);

        PerformanceResult pr = null;

        if(experimentDescriptor.gotTrues()) {

            Function<DoubleDiffResult, Boolean> trueLabeller = (_dr) -> experimentDescriptor.isTrue(_dr.testName.split("\\.")[0]);

            pr = new PerformanceResult("nlEmpire", splicing.getFirst(), (_dr) -> _dr.pval, trueLabeller,
                    false, (_dr) -> _dr.fdr <= 0.05, null, false);


            System.out.printf("settings: %s\n", test.getSettings());
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
        DataTable dataTable = BenchmarkGene.getTable(splicing.getFirst(), experimentDescriptor.getTrues(), null, splicingInfo);
        dataTable.writeCSV(cmd.getFile("o"));
        if(cmd.isSet("showtable")) {
            BenchmarkGene.showTable(splicing.getFirst(), experimentDescriptor.getTrues(), null, splicingInfo).setTitle("paired: " + pr);
            //BenchmarkGene.showTable(splicing.getSecond(), trueGenes, null, splicingInfo).setTitle("unpaired: " + prUP);
        }




    }

}
