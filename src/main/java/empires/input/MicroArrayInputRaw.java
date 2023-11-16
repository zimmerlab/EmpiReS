package empires.input;

import empires.EmpiRe;
import empires.plotting.DiffExpTable;
import lmu.utils.*;

import java.io.File;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Vector;

import static lmu.utils.ObjectGetter.*;

public class MicroArrayInputRaw {
    static final String CELLSUFFIX = ".CEL.gz";

    Logger log = LogConfig.getLogger();
    HashMap<String, Vector<String>> condition2replicatenames = new HashMap<>();

    HashMap<String, empires.input.ReplicateSetInfo> replicateSets = new HashMap<>();
    HashMap<String, empires.NormalizedReplicateSet> normedReplicateSets = new HashMap<>();
    HashMap<String, String> probeset2gene = new HashMap<>();
    HashMap<String, Vector<String>> gene2oligos = new HashMap<>();
    Vector<String> featureNames = new Vector<>();
    empires.BackgroundContextFuzzficationStrategyProvider backgroundContextFuzzficationStrategyProvider;

    public MicroArrayInputRaw(File rawdata, File labels, File mapping, boolean nolog, empires.BackgroundContextFuzzficationStrategyProvider backgroundContextFuzzficationStrategyProvider) {

        apply(FileUtils.getHrIterator(labels, FileUtils.getHeaderedReader("\t").add("id").add("condition")),
                (_hr) -> MapBuilder.updateV(condition2replicatenames, _hr.getString("condition"), _hr.getString("id")));


        this.backgroundContextFuzzficationStrategyProvider = backgroundContextFuzzficationStrategyProvider;
        probeset2gene  = FileUtils.readMap(mapping, "\t", 0, 1, "string", "string");

        HashMap<String, Vector<Double>> replicate2Data = new HashMap<>();
        Iterator<String[]> dataIterator = FileUtils.getFieldSetIterator(rawdata, " ");
        Vector<String> headers = map(dataIterator.next(), (_s) -> (!_s.endsWith(CELLSUFFIX)) ? _s : _s.split(CELLSUFFIX)[0]);

        FileUtils.Converter c = new FileUtils.Converter();
        int numUnmapped = 0;
        int nAll = 0;
        while(dataIterator.hasNext()) {
            nAll++;
            c.set(dataIterator.next());
            String oligoName = c.toString(0);
            String probeSetName = (probeset2gene.containsKey(oligoName)) ? oligoName : null;
            for(int i=oligoName.length() - 1; probeSetName == null && i>=0; i--){
                if(!Character.isDigit(oligoName.charAt(i)))
                    break;

                String testProbeSetName = oligoName.substring(0, i);
                if(!probeset2gene.containsKey(testProbeSetName))
                    continue;

                probeSetName = testProbeSetName;

            }
            if(probeSetName == null) {
                //no gene mapping
                numUnmapped++;
                continue;
            }


            MapBuilder.updateV(gene2oligos, probeset2gene.get(probeSetName), oligoName);
            featureNames.add(oligoName);

            for(int i=0; i<headers.size(); i++) {
                double val = c.toDbl(i + 1);
                //System.out.printf("%s->%s %s:%.2f\n", oligoName, probeset2gene.get(probeSetName), headers.get(i), val);
                MapBuilder.updateV(replicate2Data, headers.get(i), (nolog) ? val : (val == 0.0) ? Double.NaN : NumUtils.logN(val, 2.0));
            }
        }

        log.info("got %d/%d unmapped oligos", numUnmapped, nAll);
        for(String cond : condition2replicatenames.keySet()) {
            Vector<String> replicates =  condition2replicatenames.get(cond);
            empires.input.ReplicateSetInfo rsi = new ReplicateSetInfo(cond,replicates, featureNames);
            for(int i=0; i<replicates.size(); i++) {
                rsi.setLog2Data(i, replicate2Data.get(replicates.get(i)));
            }
            rsi.setCombinedFeatures(gene2oligos);
            replicateSets.put(cond, rsi);
        }
    }


    public HashMap<String, Vector<String>> getGene2Oligos() {
        return gene2oligos;
    }

    public HashMap<Integer, Vector<String>> getNumOligos2Genes() {
        HashMap<Integer, Vector<String>> m = new HashMap<>();
        apply(gene2oligos.entrySet(), (_e) -> MapBuilder.updateV(m, _e.getValue().size(), _e.getKey()));
        return m;
    }

    public HashMap<String, Integer> getCondition2NumReplicates() {
        return buildMap(replicateSets.values(), (_e) -> _e.getReplicateSetName(), (_e) -> _e.getNumReplicates());
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
        SimpleOptionParser cmd = new SimpleOptionParser("raw", "labels", "mapping", "cond1", "cond2", "o", "name2distribout", "interactive", "od", "minreps", "nolog");
        cmd.setFile("raw", "labels", "mapping");
        cmd.setDir("od");
        cmd.setInt("minreps");
        cmd.setDefault("minreps", "3");
        cmd.setOptional("o", "name2distribout", "cond1", "cond2", "od", "o");
        cmd.setSwitches("interactive", "nolog");
        if(!OptionParser.parseParams(args, true, false, true, true, cmd, backgroundProviderOption))
            return;

        boolean nolog = cmd.isSet("nolog");


        MicroArrayInputRaw microArrayInput = new MicroArrayInputRaw(cmd.getFile("raw"), cmd.getFile("labels"), cmd.getFile("mapping"), nolog, backgroundProviderOption.getStrategy());

        String cond1 = cmd.getOptionalValue("cond1", null);
        String cond2 = cmd.getOptionalValue("cond2", null);
        File outDir = cmd.getOptionalFile("od");


        if(outDir == null && (cond1 == null || cond2 == null )) {
            System.err.printf("provide output dir (-od) for all cond comparisons\n");
            return;
        }
        File outF = cmd.getOptionalFile("o");

        if(cond1 != null && cond2 != null && (outF == null && !cmd.isSet("interactive"))) {
            System.err.printf("provide output file (-o) \n");
            return;
        }
        if(cond1 != null && cond2 != null && ((outF != null) || cmd.isSet("interactive"))) {
            Vector<empires.DiffExpResult> diffExpResults = microArrayInput.getDiffResults(cond1, cond2);
            empires.plotting.DiffExpTable diffExpTable = new empires.plotting.DiffExpTable(microArrayInput.getNormalized(cond1), microArrayInput.getNormalized(cond2),
                    diffExpResults);

            if(outF != null) {
                diffExpTable.getTable().writeCSV(outF);
            }
            if(cmd.isOptionSet("name2distribout")) {
                new empires.NamedFoldChangeDistributionMap(diffExpResults).serialize(cmd.getFile("name2distribout"));
            }
            if(cmd.isSet("interactive")) {
                diffExpTable.showInteractiveTable();
            }
        }
        if(outDir != null) {
            int minReps = cmd.getInt("minreps");
            Vector<String> conds = filter(microArrayInput.condition2replicatenames.keySet(), (_c) -> microArrayInput.condition2replicatenames.get(_c).size() >= minReps);
            if(conds.size() < 2) {
                System.err.printf("there are no condition pair with enough replicates (-minreps = %d) : %s", minReps,
                        map(microArrayInput.condition2replicatenames.keySet(), (_c) -> String.format("%s: %d",_c, microArrayInput.condition2replicatenames.get(_c).size())));
                return;
            }
            for(UPair<String> cp : getPairs(conds, true)) {
                File cpdir = new File(outDir, String.format("%s_VS_%s", cp.getFirst(), cp.getSecond()));
                cpdir.mkdirs();
                Vector<empires.DiffExpResult> diffExpResults = microArrayInput.getDiffResults(cp.getFirst(),cp.getSecond());
                empires.plotting.DiffExpTable diffExpTable = new DiffExpTable(microArrayInput.getNormalized(cp.getFirst()), microArrayInput.getNormalized(cp.getSecond()),
                        diffExpResults);

                diffExpTable.getTable().writeCSV(new File(cpdir, "diffexp.table"));
                new empires.NamedFoldChangeDistributionMap(diffExpResults).serialize(new File(cpdir, "diffexp.distribout"));
            }
        }

    }
}
