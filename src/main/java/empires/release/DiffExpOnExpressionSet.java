package empires.release;

import lmu.utils.DataTable;
import lmu.utils.NumUtils;
import lmu.utils.OptionParser;
import lmu.utils.SimpleOptionParser;
import empires.BackgroundContextFuzzficationStrategyProvider;
import empires.DiffExpResult;
import empires.EmpiRe;
import empires.NormalizedReplicateSet;
import empires.input.BackgroundProviderOption;
import empires.input.GeneralOptions;
import empires.input.ReplicateSetInfo;

import java.io.File;
import java.util.HashMap;
import java.util.Vector;

import static lmu.utils.ObjectGetter.toVector;
import static empires.input.ExpressionSet.readDataGroupedByCondition;

public class DiffExpOnExpressionSet {

    public static void main(String[] args) {
        BackgroundProviderOption backgroundProviderOption = new BackgroundProviderOption();
        GeneralOptions generalOptions = new GeneralOptions();
        SimpleOptionParser cmd = new SimpleOptionParser("inputdir", "cond1", "cond2", "o");
        cmd.setDir("inputdir");
        cmd.setOutFile("o");
        cmd.setOptional("cond1", "cond2");

        if (!OptionParser.parseParams(args, true, false, true, true, cmd))
            return;


        generalOptions.apply();
        File inputDir = cmd.getFile("inputdir");
        HashMap<String, ReplicateSetInfo> data = readDataGroupedByCondition(inputDir);
        Vector<String> condvec = toVector(data.keySet());
        if(condvec.size() < 2) {
            System.err.printf("too few conditions!\n");
            return;
        }
        String cond1 = cmd.getOptionalValue("cond1", condvec.get(0));
        String cond2 = cmd.getOptionalValue("cond2", condvec.get(1));

        ReplicateSetInfo rs1 = data.get(cond1);
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

        BackgroundContextFuzzficationStrategyProvider backgroundContextFuzzficationStrategyProvider = backgroundProviderOption.getStrategy();
        NormalizedReplicateSet nrs1 = new NormalizedReplicateSet(rs1, backgroundContextFuzzficationStrategyProvider);
        NormalizedReplicateSet nrs2 = new NormalizedReplicateSet(rs2, backgroundContextFuzzficationStrategyProvider);
        Vector<DiffExpResult> diffExpResults = new EmpiRe().getDifferentialResults(nrs1, nrs2);

        NumUtils.sort(diffExpResults, (_d) -> _d.pval, false);
        DataTable.buildTable(diffExpResults, DataTable.buildHeader("gene", (DiffExpResult _d) -> _d.combinedFeatureName)
                .add("log2FC", (_d) -> _d.estimatedFC)
                .add("fdr", (_d) -> _d.fdr)
                //.add("log2signal."+cond1+".mean", (_d) -> nrs1.getNormed(_d.combinedFeatureName).mean)
                //.add("log2signal."+cond2+".mean", (_d) -> nrs2.getNormed(_d.combinedFeatureName).mean)
        ).writeCSV(cmd.getFile("o"));

    }
}
