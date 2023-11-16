package empires.input;

import empires.EmpiRe;
import empires.plotting.DiffExpTable;
import lmu.utils.*;

import java.io.File;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Vector;


import static lmu.utils.ObjectGetter.apply;

public class MicroArrayInput {
    HashMap<String, Vector<String>> condition2replicatenames = new HashMap<>();

    HashMap<String, ReplicateSetInfo> replicateSets = new HashMap<>();
    HashMap<String, empires.NormalizedReplicateSet> normedReplicateSets = new HashMap<>();
    HashMap<String, Vector<String>> gene2oligos = new HashMap<>();
    Vector<String> featureNames = new Vector<>();
    empires.BackgroundContextFuzzficationStrategyProvider backgroundContextFuzzficationStrategyProvider;

    public MicroArrayInput(File rawdata, File labels, File mapping, empires.BackgroundContextFuzzficationStrategyProvider backgroundContextFuzzficationStrategyProvider) {

        this.backgroundContextFuzzficationStrategyProvider = backgroundContextFuzzficationStrategyProvider;

        apply(FileUtils.getHrIterator(labels, FileUtils.getHeaderedReader("\t").add("id").add("condition")),
                (_hr) -> MapBuilder.updateV(condition2replicatenames, _hr.getString("condition"), _hr.getString("id")));


        apply(FileUtils.getFieldSetIterator(mapping, "\t"), (_sp) -> MapBuilder.updateV(gene2oligos, _sp[1], _sp[0]));


        HashMap<String, Vector<Double>> replicate2Data = new HashMap<>();
        Iterator<String[]> dataIterator = FileUtils.getFieldSetIterator(rawdata, "\t");
        String[] header = dataIterator.next();

        FileUtils.Converter c = new FileUtils.Converter();
        while(dataIterator.hasNext()) {
            c.set(dataIterator.next());
            featureNames.add(c.toString(0));
            for(int i=1; i<header.length; i++) {
                double val = c.toDbl(i);
                MapBuilder.updateV(replicate2Data, header[i], (val == 0.0) ? Double.NaN : NumUtils.logN(val, 2.0));
            }
        }

        for(String cond : condition2replicatenames.keySet()) {
            Vector<String> replicates =  condition2replicatenames.get(cond);
            ReplicateSetInfo rsi = new ReplicateSetInfo(cond,replicates, featureNames);
            for(int i=0; i<replicates.size(); i++) {
                rsi.setLog2Data(i, replicate2Data.get(replicates.get(i)));
            }
            rsi.setCombinedFeatures(gene2oligos);
            replicateSets.put(cond, rsi);
        }
    }

    public empires.NormalizedReplicateSet getNormalized(String condition) {
        empires.NormalizedReplicateSet normed = normedReplicateSets.get(condition);
        if(normed == null) {
            normedReplicateSets.put(condition, normed = new empires.NormalizedReplicateSet(replicateSets.get(condition), backgroundContextFuzzficationStrategyProvider));
        }
        return normed;
    }

    public Vector<empires.DiffExpResult> getDiffResults(String cond1, String cond2) {
        empires.NormalizedReplicateSet nrs1 = getNormalized(cond1);
        empires.NormalizedReplicateSet nrs2 = getNormalized(cond2);
        return (new EmpiRe()).getDifferentialResults(nrs1, nrs2);
    }

    public static void main(String[] args) {
        BackgroundProviderOption backgroundProviderOption = new BackgroundProviderOption();
        SimpleOptionParser cmd = new SimpleOptionParser("raw", "labels", "mapping", "cond1", "cond2", "o", "name2distribout", "interactive");
        cmd.setFile("raw", "labels", "mapping");
        cmd.setOptional("o", "name2distribout");
        cmd.setSwitches("interactive");
        if(!OptionParser.parseParams(args, true, false, true, true, cmd, backgroundProviderOption))
            return;


        MicroArrayInput microArrayInput = new MicroArrayInput(cmd.getFile("raw"), cmd.getFile("labels"), cmd.getFile("mapping"), backgroundProviderOption.getStrategy());

        String cond1 = cmd.getValue("cond1");
        String cond2 = cmd.getValue("cond2");
        Vector<empires.DiffExpResult> diffExpResults = microArrayInput.getDiffResults(cond1, cond2);
        empires.plotting.DiffExpTable diffExpTable = new DiffExpTable(microArrayInput.getNormalized(cond1), microArrayInput.getNormalized(cond2),
                diffExpResults);

        if(cmd.isOptionSet("o")) {
            diffExpTable.getTable().writeCSV(cmd.getFile("o"));
        }
        if(cmd.isOptionSet("name2distribout")) {
            new empires.NamedFoldChangeDistributionMap(diffExpResults).serialize(cmd.getFile("name2distribout"));
        }
        if(cmd.isSet("interactive")) {
            diffExpTable.showInteractiveTable();
        }

    }
}
