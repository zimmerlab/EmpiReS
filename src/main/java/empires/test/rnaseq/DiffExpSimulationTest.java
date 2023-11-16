package empires.test.rnaseq;

import empires.input.EBROWSER_INPUT;
import empires.input.ExpressionSet;
import empires.input.ReplicateSetInfo;
import empires.plotting.DiffExpTable;
import empires.rnaseq.DiffExpSimulation;
import lmu.utils.*;
import lmu.utils.fdr.PerformanceResult;
import lmu.utils.plotting.PlotCreator;
import lmu.utils.swing.PagedDataTable;
import lmu.utils.tuple.Tuple3;
import empires.DiffExpResult;
import empires.EmpiRe;
import empires.Normalization;
import empires.NormalizedReplicateSet;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;
import java.util.function.Function;

import static lmu.utils.ObjectGetter.*;
import static empires.input.ExpressionSet.readDataGroupedByCondition;

public class DiffExpSimulationTest {

    static final int MAX_NUM_INVALIDS = 800;

    static empires.rnaseq.DiffExpSimulation simulateDiffExp(HashMap<String, empires.input.ReplicateSetInfo> conditions, int num_replicates, double meanFc, double minFC) {

        empires.input.ReplicateSetInfo maxRepl = NumUtils.maxObj(conditions.values(), (_rsi) -> _rsi.getNumReplicates()).getSecond();

        Vector<Integer> replicates = Normalization.getBestSubSample(maxRepl.getLog2Data(), num_replicates);

        NormalDistribution diffExp = new NormalDistribution(meanFc, 0.3 * meanFc);

        return new empires.rnaseq.DiffExpSimulation(maxRepl.getReplicateSubSet(replicates), diffExp, minFC, 0.2);

    }

    static void writeEBFile(empires.rnaseq.DiffExpSimulation de, File outdir, empires.input.EBROWSER_INPUT type) {
        File of = type.get(outdir);

        PrintWriter pw = FileUtils.getWriter(of);
        UPair<empires.input.ReplicateSetInfo> diffs = de.getSimulation();
        Vector<String> features = filter(diffs.getFirst().getFeatureNames(), (_s) -> de.getUsedFeatures().contains(_s));

        switch (type) {
            case EXPRS:
                for(String f : features) {
                    Vector<Double> v = diffs.get(true).getReplicateData(f);
                    v.addAll(diffs.get(false).getReplicateData(f));

                    pw.println(StringUtils.join("\t", v, (_d) -> String.format("%.0f", Math.pow(2.0, _d))));
                }
                break;
            case FDATA:
                apply(features, (_f) -> pw.printf("%s\t%s\n", _f, _f));
                break;
            case PDATA:
                for(int i=0; i<2; i++) {
                    empires.input.ReplicateSetInfo rsi = diffs.get(i == 0);
                    String n = (i == 0) ? "control" : "diff";
                    for(int j=0; j<rsi.getNumReplicates(); j++) {
                        pw.printf("%s.rep%d\t%d\n", n, j + 1, i);
                    }

                }
                break;
        }
        pw.close();
    }

    static void writeEBFiles(empires.rnaseq.DiffExpSimulation de, File outdir) {
        for(empires.input.EBROWSER_INPUT eb : empires.input.EBROWSER_INPUT.values()) {
            writeEBFile(de, outdir, eb);
        }

        PrintWriter pw = FileUtils.getWriter(outdir, "trues.txt");
        for(int i=0; i<2; i++) {
            apply(de.getSimulatedFeatures(i == 0), (_s) -> pw.println(_s));
        }
        pw.close();

    }

    static empires.rnaseq.DiffExpSimulation getSimulatedDiffExp(HashMap<String, empires.input.ReplicateSetInfo> conditions, double meanFc, double minFC, int num_replicates) {
        Vector<empires.input.ReplicateSetInfo> v = filter(conditions.values(), (_r) -> _r.getNumReplicates() >= 6);
        empires.input.ReplicateSetInfo base = v.get((int)(v.size() * Math.random()));
        Vector<Integer> replicates = Normalization.getBestSubSample(base.getLog2Data(), num_replicates);
        NormalDistribution diffExp = new NormalDistribution(meanFc, 0.3 * meanFc);

        empires.rnaseq.DiffExpSimulation diffExpSimul = new empires.rnaseq.DiffExpSimulation(base.getReplicateSubSet(replicates), diffExp, minFC, 0.1);

        while(diffExpSimul.getNumInvalids() > MAX_NUM_INVALIDS) {
            System.out.printf("resimulate got %d invalids\n", diffExpSimul.getNumInvalids());
            diffExpSimul = getSimulatedDiffExp(conditions, meanFc, minFC, num_replicates);
        }
        return diffExpSimul;
    }


    static class EmpiReResult {
        Vector<Double> shifts;
        empires.rnaseq.DiffExpSimulation diffExpSimulation;
        NormalizedReplicateSet control;
        NormalizedReplicateSet diff;
        Vector<DiffExpResult> diffExpResults;
        PerformanceResult pr;
        Function<String, Boolean> isTrue;
        double minFC;

        HashMap<String, PerformanceResult> methodPerf = new HashMap<>();


        public void showTable() {
            empires.plotting.DiffExpTable table = new DiffExpTable(control, diff, diffExpResults, minFC * 0.7, 0.05, (_id) -> isTrue.apply(_id));
            PagedDataTable.MJFrame mj = table.showInteractiveTable(Pair.create("usedFeature", (_de) -> diffExpSimulation.getUsedFeatures().contains(_de.combinedFeatureName)),
                    Pair.create("simulatedFC", (_de) -> diffExpSimulation.getSimulatedFC(_de.combinedFeatureName)),
                    Pair.create("testInfo", (_de) -> diffExpSimulation.getNormTestInfo(_de.combinedFeatureName))

            );
            mj.setTitle("inv" + diffExpSimulation.getNumInvalids() +":" + pr);
            mj.setVisible(true);

        }
    }

    static EmpiReResult runEmpiRe(empires.rnaseq.DiffExpSimulation diffExpSimulation, double minFC) {


        UPair<empires.input.ReplicateSetInfo> simul = diffExpSimulation.getSimulation();

        EmpiReResult empiReResult = new EmpiReResult();
        empiReResult.diffExpSimulation = diffExpSimulation;
        empiReResult.minFC = minFC;
        empiReResult.control = new NormalizedReplicateSet(simul.getFirst());
        empiReResult.diff = new NormalizedReplicateSet(simul.getSecond());

        EmpiRe emp = new EmpiRe();

        HashSet<String> allfeatures = new HashSet<>();
        allfeatures.addAll(empiReResult.control.getInData().getFeatureNames());
        allfeatures.addAll(empiReResult.diff.getInData().getFeatureNames());

        empiReResult.diffExpResults = emp.getDifferentialResults(empiReResult.control, empiReResult.diff, allfeatures);
        empiReResult.shifts = first(empiReResult.diffExpResults).getDiffExpManager().shifts;

        empiReResult.isTrue = (_s) ->  diffExpSimulation.getSimulatedFeatures(true).contains(_s)
                || diffExpSimulation.getSimulatedFeatures(false).contains(_s);

        empiReResult.pr = new PerformanceResult("nlEmpire", empiReResult.diffExpResults, (_d) -> _d.pval,
                (_d) -> empiReResult.isTrue.apply(_d.combinedFeatureName)
                , false, (_d) -> _d.fdr <= 0.05, null);

        System.out.printf("fc: %s\n",  new PerformanceResult("nlEmpire-fc", empiReResult.diffExpResults, (_d) -> _d.fcEstimatePval,
                (_d) -> empiReResult.isTrue.apply(_d.combinedFeatureName)
                , false, (_d) -> _d.fcEstimateFDR <= 0.05, null)
        );

        empiReResult.methodPerf.put("empiRe", empiReResult.pr);
        return empiReResult;
    }

    static EmpiReResult testdiffExp(HashMap<String, empires.input.ReplicateSetInfo> conditions, double meanFc, double minFC, int num_replicates, File od,
                                    File diffexp_script) {
        empires.rnaseq.DiffExpSimulation diffExpSimulation = getSimulatedDiffExp(conditions, meanFc, minFC, num_replicates);

        EmpiReResult empiReResult = runEmpiRe(diffExpSimulation, minFC);

        if(diffexp_script != null) {
            apply(runDEMethods(diffExpSimulation, diffexp_script, od).entrySet(),
                    (_e) -> empiReResult.methodPerf.put(_e.getKey().n, _e.getValue())
            );
        }
        return empiReResult;
    }




    enum DEMethod {
        DESEQ("DESeq"),
        EDGER("edgeR"),
        limma("limma")
        ;
        String n;

        DEMethod(String n) {
            this.n = n;
        }

    }


    static HashMap<DEMethod, PerformanceResult> runDEMethods(empires.rnaseq.DiffExpSimulation diffExpSimulation, File script, File outdir) {
        writeEBFiles(diffExpSimulation, outdir);

        HashMap<DEMethod, PerformanceResult> rv = new HashMap<>();

        if(script == null) {
            return rv;
        }



        for(DEMethod dm : DEMethod.values()) {
            File outf = new File(outdir, dm.n+".out");

            String invoke = String.format("%s %s %s %s %s %s", script.getAbsolutePath(),
                    empires.input.EBROWSER_INPUT.EXPRS.get(outdir).getAbsolutePath(),
                    empires.input.EBROWSER_INPUT.PDATA.get(outdir).getAbsolutePath(),
                    EBROWSER_INPUT.FDATA.get(outdir).getAbsolutePath(),
                    dm.n,
                    outf.getAbsolutePath()
            );

            try
            {
                Process proc = Runtime.getRuntime().exec(invoke);
                proc.waitFor();

            }
            catch (IOException ie) {
                throw new FRuntimeException("i/o error: %s while running: %s!", ie, ie.getMessage(), invoke);
            }
            catch(InterruptedException iex) {
                throw new FRuntimeException("interrupted!");
            }


            Vector<Tuple3<String, Double, Double>> res = toVector(FileUtils.getConvertedHrIterator(outf, FileUtils.getHeaderedReader("\t").add("GENE.ID", "gene")
                            .add("log2FC", "fc").add("ADJ.PVAL", "fdr"),
                    false, (_hr) -> Tuple3.create(_hr.getString("gene"), _hr.getDbl("fc"), _hr.getDbl("fdr"))));

            HashSet<String> diffRes = mapToSet(res, (_r) -> _r.get0());
            Set<String> missing = SetInfo.minus(diffExpSimulation.getUsedFeatures(), diffRes);
            apply(missing, (_m) -> res.add(Tuple3.create(_m, 0.0, 1.0)));

            System.out.printf("%s: add %d missing\n", dm.n, missing.size());

            PerformanceResult pr = new PerformanceResult(dm.n, res, (_t) -> _t.get2(), (_t) -> diffExpSimulation.isSimulated(_t.get0()), false,
                    (_t) -> _t.get2() <= 0.05, null, false);

            System.out.printf("%s\n", pr);
            rv.put(dm, pr);
        }

        return rv;
    }


    public static void main(String[] args) {

        SimpleOptionParser cmd = new SimpleOptionParser("ebrowser", "numreps", "meanfc", "minfc", "ntests", "od", "descript", "perfout");
        cmd.setFile("ebrowser", "descript");
        cmd.setDir("od");
        cmd.setDefault("numreps", "6");
        cmd.setInt("ntests");
        cmd.setDefault("ntests", "-1");
        cmd.setDouble("meanfc");
        cmd.setDouble("minfc");
        cmd.setDefault("meanfc", "1.0");
        cmd.setDefault("minfc", "0.3");
        cmd.setOptional("od", "descript", "perfout");

        if(!OptionParser.parseParams(args, true, false, true, true, cmd))
            return;


        HashMap<String, empires.input.ReplicateSetInfo> conditions = ExpressionSet.readDataGroupedByCondition(cmd.getFile("ebrowser"));
        double meanFC = cmd.getDouble("meanfc");
        double minFC = cmd.getDouble("minfc");
        int numReps = cmd.getInt("numreps");

        int ntests = cmd.getInt("ntests");

        File od = cmd.getOptionalFile("od");
        File script = cmd.getOptionalFile("descript");
        PrintWriter perfout = cmd.getOptionalWriter("perfout");

        if(ntests > 0) {

            if(perfout != null) {
                perfout.printf("testIdx\tmethod\tprec\trec\tauroc\taupr\n");
            }

            for(int i=0; i<ntests; i++) {

                EmpiReResult result = testdiffExp(conditions, meanFC, minFC, numReps, od, script);

                if(perfout != null) {


                    for(Map.Entry<String, PerformanceResult> e : result.methodPerf.entrySet()) {

                        PerformanceResult pr = e.getValue();
                        perfout.printf("%d\t%s\t%.4f\t%.4f\t%.4f\n", i + 1, e.getKey(), pr.PREC, pr.RECALL, pr.AUROC, pr.AUPRC);
                    }
                    perfout.flush();
                }

                System.out.printf("%s\n", result.pr);
                System.err.printf("%s\n", result.pr);
                //assertTrue(pr.PREC >= 0.93, "perf too bad: " + pr);

                if(result.pr.PREC < 0.8) {

                    result.showTable();
                    break;
                }

            }
            if(perfout != null) {
                perfout.close();
            }
            return;
        }

        PlotCreator pc = new PlotCreator();


        DiffExpSimulation diffExpSimulation = getSimulatedDiffExp(conditions, meanFC, minFC, numReps);



        if(od != null) {

            runDEMethods(diffExpSimulation, script, od);

        }

        EmpiReResult result = runEmpiRe(diffExpSimulation, minFC);


        UPair<ReplicateSetInfo> simul = diffExpSimulation.getSimulation();

        NormalizedReplicateSet control = result.control;
        NormalizedReplicateSet diff = result.diff;

        Vector<Double> shifts = result.shifts;

        HashMap<UPair<String>, Vector<Double>> replicatePair2Simulations = new HashMap<>();

        Vector<Double> foldchanges_up = new Vector<>();
        Vector<Double> foldchanges_dn = new Vector<>();

        for(int i=0; i<2; i++) {
            Vector<Double> target = (i == 0) ? foldchanges_up : foldchanges_dn;
            for(String feature : diffExpSimulation.getSimulatedFeatures(i == 0)) {

                double m1 = control.getNormed(feature).mean + shifts.get(0);
                double m2 = diff.getNormed(feature).mean + shifts.get(1);

                target.add(m1 - m2);

            }
        }



        pc.hist(String.format("up %d/%d", foldchanges_up.size(), result.diffExpResults.size()), foldchanges_up, 100, false, false);
        pc.hist(String.format("dn %d/%d", foldchanges_dn.size(), result.diffExpResults.size()), foldchanges_dn, 100, false, false);
        pc.setLimits(- meanFC -2.0, meanFC + 2.0, null, null);
        pc.setLabels("foldchange", "freq", "topright");
        ImageUtils.showImage("test", pc.getImage());


        System.out.printf("perf: %s\n", result.pr);

        result.showTable();
    }

}
