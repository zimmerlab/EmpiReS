package empires.test;

import empires.EmpiRe;
import empires.input.BackgroundProviderOption;
import empires.input.GeneralOptions;
import empires.input.ReplicateSetInfo;
import empires.test.rnaseq.BenchmarkGene;
import lmu.utils.*;

import java.io.File;
import java.util.*;
import java.util.function.BiFunction;
import java.util.function.Function;
import static lmu.utils.ObjectGetter.*;
import static lmu.utils.ObjectGetter.filteredSize;

public class CPTACSplicing {
    public static void main(String[] args) {
        empires.input.GeneralOptions generalOptions = new GeneralOptions();
        empires.input.BackgroundProviderOption backgroundProviderOption = new BackgroundProviderOption();
        SimpleOptionParser cmd = new SimpleOptionParser("matrix", "labels", "gpeptidome", "doublediffvariant", "o");
        cmd.setFile("matrix", "labels", "gpeptidome");
        cmd.setDefault("doublediffvariant", empires.DoubleDiffVariant.ALLPAIRS.getName());


        if (!OptionParser.parseParams(args, true, false, true, true, cmd, backgroundProviderOption, generalOptions))
            return;

        generalOptions.apply();

        File gpeptidome = cmd.getFile("gpeptidome");

        HashMap<UPair<String>, Boolean> peptide2junction = new HashMap<>();
        apply(FileUtils.getHrIterator(gpeptidome, FileUtils.getHeaderedReader("\t").add("gene", "peptide", "junction")),
                (_hr) -> peptide2junction.put(UPair.createU(_hr.getString("gene"), _hr.getString("peptide")), _hr.getBool("junction")));

        empires.DoubleDiffVariant doubleDiffVariant = empires.DoubleDiffVariant.get(cmd.getValue("doublediffvariant"));
        HashMap<String, String> replicate2condition = FileUtils.readMap(cmd.getFile("labels"), "\t", 0, 1);
        HashMap<String, Vector<String>> condition2replicates = new HashMap<>();
        apply(replicate2condition.entrySet(), (_e) -> MapBuilder.updateV(condition2replicates, _e.getValue(), _e.getKey()));

        boolean noplots = cmd.isSet("noplots");
        File matrix = cmd.getFile("matrix");

        Iterator<String[]> it = FileUtils.getFieldSetIterator(matrix, "\t");
        Vector<String> header = map(it.next(), (_s) -> _s.trim());
        String idfield = "pep";

        HashSet<String> headerSet = toSet(header);
        if (filteredSize(replicate2condition.keySet(), (_n) -> !headerSet.contains(_n)) > 0) {
            throw new FRuntimeException("some replicates missing in matrix file: %s\nheaders: %s\n", SetInfo.minus(replicate2condition.keySet(), headerSet), headerSet);
        }
        Set<String> notyetmapped = SetInfo.minus(headerSet, replicate2condition.keySet());
        if (!notyetmapped.contains(idfield)) {
            System.err.printf("could not find idfield: %s in not replicated headers: %s! please provide a valid feature id field with -idfield!", idfield, notyetmapped);
            return;
        }

        /** read-in matrix data - skip idfield column */
        int idIdx = header.indexOf(idfield);
        Vector<String> features = new Vector<>();
        Vector<Vector<Double>> label2data = map(header, (_i) -> new Vector<>());

        FileUtils.Converter converter = new FileUtils.Converter();
        while (it.hasNext()) {
            converter.set(it.next());
            features.add(converter.toString(idIdx));
            if (converter.length() != header.size())
                throw new FRuntimeException("invalid line: %s expected %d fields, got %d!", toVector(converter.getRaw()), converter.length(), header.size() + 1);

            for (int i = 0; i < header.size(); i++) {
                if (idIdx == i)
                    continue;
                double v = converter.toDbl(i);
                label2data.get(i).add((v == 0.0) ? Double.NaN : NumUtils.logN(v, 2.0));
            }
        }

        HashMap<String, empires.NormalizedReplicateSet> conditions = new HashMap<>();

        Logger log = LogConfig.getLogger();
        for (String cond : condition2replicates.keySet()) {

            Vector<String> replicates = condition2replicates.get(cond);
            empires.input.ReplicateSetInfo rsi = new ReplicateSetInfo(cond, replicates, features);
            for (int i = 0; i < replicates.size(); i++) {
                String rep = replicates.get(i);
                int idx = first(filterIndex(header, (_s) -> rep.equals(_s)));
                //System.out.printf("add %s:%s -> %d\n", cond, rep, idx);
                rsi.setLog2Data(i, label2data.get(idx));
            }
            conditions.put(cond, new empires.NormalizedReplicateSet(rsi));
        }

        Vector<String> condvec = toSortedVector(conditions.keySet(), true);

        String cond1 = condvec.get(0);
        String cond2 = condvec.get(1);

        empires.EmpiRe emp = new EmpiRe();
        empires.NormalizedReplicateSet n1 = conditions.get(cond1);
        empires.NormalizedReplicateSet n2 = conditions.get(cond2);


        int[] numtests = {0, 0};

        empires.FeatureBasedSplicingTest.SplicingTestFilter testFilter =
                (String gene, Vector<String> isoform_ids_eqclass1, Vector<String> isoform_ids_eqclass2, Vector<String> peptide_ids1, Vector<String> peptide_ids2)
                        ->
                {
                    numtests[0]++;
                    if(numtests[0] % 100 == 0) {
                        log.info("at splicing test %d ...", numtests[0]);
                    }
                    Set<String> uniqPeptides1 = mapToSet(peptide_ids1, (_p) -> _p.split("_")[0]);
                    Set<String> uniqPeptides2 = mapToSet(peptide_ids2, (_p) -> _p.split("_")[0]);

                    int numJunctions1 = filteredSize(uniqPeptides1, (_p) -> peptide2junction.get(UPair.createU(gene, _p)));
                    int numJunctions2 = filteredSize(uniqPeptides2, (_p) -> peptide2junction.get(UPair.createU(gene, _p)));

                    if ((numJunctions1 == 0 && uniqPeptides1.size() < 2) || (numJunctions2 == 0 && uniqPeptides2.size() < 2)) {

                        numtests[1]++;
                        log.info("skip test %d/%d: %s %s -> %s (%d junctions) %s -> %s (%d junctions)", numtests[1], numtests[0],
                                gene, peptide_ids1, uniqPeptides1, numJunctions1, peptide_ids2, uniqPeptides2, numJunctions2);
                        return false;
                    }

                    return true;

                };

        Vector<empires.DoubleDiffResult> results = new empires.FeatureBasedSplicingTest(n1, n2, gpeptidome, doubleDiffVariant, 2, testFilter).getSplicingResults();

        /** benchmarkgenes wraps the doublediff results -> this wrapping is here irrelevant
            f g.ddr is the DoubleDiffResult it has pval, fdr and meanfc
            the used feature infos are in g.ddr.featureInfos1 and g.ddr.featureInfos2
            of type FeatureInfo
            in FeatureInfo the field: feature is the name/id of the feature
            and via the methods getCond1LogVals(), getCond2LogVals() one can retrieve the normalized log signal values
            in the analyzed conditions
         */
        Function<Collection<empires.FeatureInfo>, Set<String>> getUniqPeps = (_v) -> mapToSet(_v, (_f) ->  _f.feature.split("_")[0]);
        BiFunction<String, Set<String>, Integer> getNumJunctions = (_g, _s) -> (_s == null) ? 0 : filteredSize(_s, (_f) -> peptide2junction.getOrDefault(UPair.createU(_g, _s), false));

        DataTable table = empires.test.rnaseq.BenchmarkGene.getTable(results, null, null, null,
                Pair.create("num.uniqpeps1", (empires.test.rnaseq.BenchmarkGene g) -> getUniqPeps.apply(g.ddr.featureInfos1).size()),
                Pair.create("num.uniqpeps2", (empires.test.rnaseq.BenchmarkGene g) -> getUniqPeps.apply(g.ddr.featureInfos2).size()),
                Pair.create("num.junctions1", (empires.test.rnaseq.BenchmarkGene g) -> getNumJunctions.apply(g.geneId, getUniqPeps.apply(g.ddr.featureInfos1))),
                Pair.create("num.junctions2", (BenchmarkGene g) -> getNumJunctions.apply(g.geneId, getUniqPeps.apply(g.ddr.featureInfos2)))
        );
        table.writeCSV(cmd.getFile("o"));
    }
}
