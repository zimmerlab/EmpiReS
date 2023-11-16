package empires.evaluation;

import lmu.utils.FileUtils;
import lmu.utils.OptionParser;
import lmu.utils.Pair;
import lmu.utils.SimpleOptionParser;
import empires.DiffExpResult;
import empires.EmpiRe;
import empires.NormalizedReplicateSet;
import empires.input.ReplicateSetInfo;
import empires.plotting.DiffExpTable;

import java.io.File;
import java.util.HashMap;
import java.util.Vector;
import java.util.function.Function;

import static lmu.utils.ObjectGetter.buildReverseMap;
import static empires.input.ExpressionSet.readDataGroupedByCondition;

public class CompareToDiffExpTools {


    static class DiffResult {
        String id = "";
        double fdr = 1.1;
        double fc = 0.0;

        DiffResult() {

        }
        DiffResult(FileUtils.HeaderedReader hr) {
            this.id = hr.getString("gene");
            this.fdr = hr.getDbl("fdr");
            this.fc = hr.getDbl("fc");
        }
        public static HashMap<String, DiffResult> read(File diffout) {
            return buildReverseMap(
                    FileUtils.getConvertedHrIterator(diffout,
                            FileUtils.getHeaderedReader("\t").add("GENE.ID", "gene").add("log2FC", "fc").add("ADJ.PVAL", "fdr"),
                            false,
                            (_hr) -> new DiffResult(_hr)
                        ), (_d) -> _d.id);
        }
    }


    enum DiffExpMethod{
        DESEQ("DESeq.out"),
        EDGER("edgeR.out"),
        LIMMA("limma.out")
        ;

        String filename;

        DiffExpMethod(String fn) {
            this.filename = fn;
        }

        public HashMap<String, DiffResult> read(File dir) {
            File dfile = new File(dir, filename);
            if(!dfile.exists())
                return null;

            return DiffResult.read(dfile);
        }
    }

    HashMap<String, ReplicateSetInfo> dataGroupedByCondition;

    NormalizedReplicateSet cond1;
    NormalizedReplicateSet cond2;
    Vector<DiffExpResult> diffResults;

    HashMap<DiffExpMethod, HashMap<String, DiffResult>> method2diff = new HashMap<>();

    CompareToDiffExpTools(File EBrowserDir) {
        dataGroupedByCondition = readDataGroupedByCondition(EBrowserDir);

        for(DiffExpMethod dmethod : DiffExpMethod.values()) {
            HashMap<String, DiffResult> m = dmethod.read(EBrowserDir);
            if(m == null)
                continue;

            method2diff.put(dmethod, m);
        }

        cond1 = new NormalizedReplicateSet(dataGroupedByCondition.get("0"));
        cond2 = new NormalizedReplicateSet(dataGroupedByCondition.get("1"));
        EmpiRe empiRe = new EmpiRe();
        diffResults = empiRe.getDifferentialResults(cond1, cond2);

        //Pair<String, Function<DiffExpResult, Object>> ... additionalHeaders


    }

    public void showTable() {
        DiffExpTable diffExpTable = new DiffExpTable(cond1, cond2, diffResults);

        DiffResult def = new DiffResult();

        Vector<Pair<String, Function<DiffExpResult, Object>>> headers = new Vector<>();

        for(DiffExpMethod method : method2diff.keySet()) {

            HashMap<String, DiffResult> m2d = method2diff.get(method);

            headers.add(Pair.create(method+".fdr", (_de) -> m2d.getOrDefault(_de.combinedFeatureName, def).fdr));
            headers.add(Pair.create(method+".log2fc", (_de) -> m2d.getOrDefault(_de.combinedFeatureName, def).fc));
        }


        diffExpTable.showInteractiveTable(headers);
    }

    public static void main(String[] args) {
        SimpleOptionParser cmd = new SimpleOptionParser("dir");
        cmd.setDir("dir");
        if(!OptionParser.parseParams(args, true, false, true, true, cmd))
            return;

        CompareToDiffExpTools comp = new CompareToDiffExpTools(cmd.getFile("dir"));
        comp.showTable();
    }
}
