package empires.test.rnaseq;

import empires.input.EBROWSER_INPUT;
import lmu.utils.*;
import lmu.utils.fdr.PerformanceResult;
import lmu.utils.tuple.Tuple3;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;
import java.util.Vector;

import static lmu.utils.ObjectGetter.*;

public class TestStandardDEMethods {

    public static void main(String[] args) {
        SimpleOptionParser cmd = new SimpleOptionParser("dir", "script", "outfile");
        cmd.setDir("dir");
        cmd.setFile("script", "outfile");
        cmd.setOptional("script", "outfile");
        if(!OptionParser.parseParams(args, true, false, true, true, cmd))
            return;

        File script = cmd.getOptionalFile("script");
        File ebdir = cmd.getFile("dir");

        HashSet<String> trues = FileUtils.readSet(empires.input.EBROWSER_INPUT.TRUES.get(ebdir));

        File testFile = cmd.getOptionalFile("outfile");

        if(testFile != null) {
            Vector<Tuple3<String, Double, Double>> res = toVector(FileUtils.getConvertedHrIterator(testFile, FileUtils.getHeaderedReader("\t").add("GENE.ID", "gene")
                            .add("log2FC", "fc").add("ADJ.PVAL", "fdr"),
                    false, (_hr) -> Tuple3.create(_hr.getString("gene"), _hr.getDbl("fc"), _hr.getDbl("fdr"))));

            PerformanceResult pr = new PerformanceResult(testFile.getName(), res, (_t) -> _t.get2(), (_t) -> trues.contains(_t.get0()), false, (_t) -> _t.get2() <= 0.05, null, false);
            System.out.printf("%s\n", pr);
            return;
        }
        if(script == null) {
            System.out.printf("need a scriptfile to perform tests!\n");
            return;
        }
        for(empires.test.rnaseq.DiffExpSimulationTest.DEMethod dm : DiffExpSimulationTest.DEMethod.values()) {
            File outf = new File(ebdir, dm.n + ".out");

            String invoke = String.format("%s %s %s %s %s %s", script.getAbsolutePath(),
                    empires.input.EBROWSER_INPUT.EXPRS.get(ebdir).getAbsolutePath(),
                    empires.input.EBROWSER_INPUT.PDATA.get(ebdir).getAbsolutePath(),
                    EBROWSER_INPUT.FDATA.get(ebdir).getAbsolutePath(),
                    dm.n,
                    outf.getAbsolutePath()
            );

            try {
                Process proc = Runtime.getRuntime().exec(invoke);
                proc.waitFor();

            } catch (IOException ie) {
                throw new FRuntimeException("i/o error: %s while running: %s!", ie, ie.getMessage(), invoke);
            } catch (InterruptedException iex) {
                throw new FRuntimeException("interrupted!");
            }


            Vector<Tuple3<String, Double, Double>> res = toVector(FileUtils.getConvertedHrIterator(outf, FileUtils.getHeaderedReader("\t").add("GENE.ID", "gene")
                            .add("log2FC", "fc").add("ADJ.PVAL", "fdr"),
                    false, (_hr) -> Tuple3.create(_hr.getString("gene"), _hr.getDbl("fc"), _hr.getDbl("fdr"))));

            HashSet<String> diffRes = mapToSet(res, (_r) -> _r.get0());
            Set<String> missing = SetInfo.minus(trues, diffRes);
            apply(missing, (_m) -> res.add(Tuple3.create(_m, 0.0, 1.0)));

            System.out.printf("%s: add %d missing\n", dm.n, missing.size());

            PerformanceResult pr = new PerformanceResult(dm.n, res, (_t) -> _t.get2(), (_t) -> trues.contains(_t.get0()), false, (_t) -> _t.get2() <= 0.05, null, false);

            System.out.printf("%s\n", pr);


        }
    }
}
