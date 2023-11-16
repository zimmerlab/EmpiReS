package empires.input;

import empires.plotting.DiffExpTable;
import lmu.utils.FRuntimeException;
import lmu.utils.OptionParser;
import lmu.utils.SimpleOptionParser;
import empires.DiffExpResult;
import empires.EmpiRe;
import empires.NormalizedReplicateSet;

import java.util.HashMap;
import java.util.Vector;

import static empires.input.ExpressionSet.readDataGroupedByCondition;

public class RNASeqCountInfo {

    public static void main(String[] args) {

        SimpleOptionParser cmd = new SimpleOptionParser("ebrowser", "cond1", "cond2");

        if(!OptionParser.parseParams(args, true, false, true, true, cmd))
            return;


        HashMap<String, ReplicateSetInfo> conditions = ExpressionSet.readDataGroupedByCondition(cmd.getFile("ebrowser"));

        String cond1 = cmd.getValue("cond1");
        String cond2 = cmd.getValue("cond2");

        if(conditions.get(cond1) == null)
            throw new FRuntimeException("unknown condition: %s known: %s", cond1, conditions.keySet());

        if(conditions.get(cond2) == null)
            throw new FRuntimeException("unknown condition: %s known: %s", cond2, conditions.keySet());


        NormalizedReplicateSet nrs1 = new NormalizedReplicateSet(conditions.get(cond1));
        NormalizedReplicateSet nrs2 = new NormalizedReplicateSet(conditions.get(cond2));

        EmpiRe empiRe = new EmpiRe();

        Vector<DiffExpResult> diffExpResults = empiRe.getDifferentialResults(nrs1, nrs2);

        empires.plotting.DiffExpTable table = new DiffExpTable(nrs1, nrs2, diffExpResults);

        table.showInteractiveTable();


    }


}
