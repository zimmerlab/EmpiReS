package empires.input;

import empires.EmpiRe;
import empires.plotting.DiffExpTable;
import lmu.utils.*;

import java.io.File;
import java.util.*;
import java.util.function.BiFunction;

import static lmu.utils.ObjectGetter.*;
import static lmu.utils.ObjectGetter.filter;

public class LFQInput {
    public static void main(String[] args) {
        Locale.setDefault(Locale.UK);
        empires.input.GeneralOptions generalOptions = new GeneralOptions();
        BackgroundProviderOption backgroundProviderOption = new BackgroundProviderOption();
        SimpleOptionParser cmd = new SimpleOptionParser("trues", "falses", "pgroup", "labels", "cond1", "cond2", "peptides",
                "useonlypeps",  "minfc", "o", "interactive", "minpeps", "minnummeasuredpercondition", "useonlyproteins",
                "gpeptidome", "splicingout", "maxisoforms", "minpepsforsplicingtest");

        cmd.setFile("trues", "pgroup", "labels", "peptides", "useonlypeps", "gpeptidome");

        cmd.setDouble("minfc");
        cmd.setDefault("minfc", "0.2");
        cmd.setInt("minpeps", "minnummeasuredpercondition", "maxisoforms", "minpepsforsplicingtest");
        cmd.setDefault("minpeps", "1");
        cmd.setDefault("maxisoforms", "10");
        cmd.setDefault("minpepsforsplicingtest", "1");


        cmd.setDefault("minnummeasuredpercondition", "2");
        cmd.setOptional("peptides",  "useonlypeps", "useonlyproteins", "trues", "falses", "splicingout", "gpeptidome");

        cmd.setSwitches("interactive");

        if (!OptionParser.parseParams(args, true, false, true, true, cmd, backgroundProviderOption, generalOptions))
            return;

        Logger log = LogConfig.getLogger();

        int minnummeasuredpercondition = cmd.getInt("minnummeasuredpercondition");
        int MINPEPS = cmd.getInt("minpeps");
        int MINPEPS_SPLICING = cmd.getInt("minpepsforsplicingtest");

        generalOptions.apply();

        File gpeptidome = cmd.getOptionalFile("gpeptidome");
        File splicingout = cmd.getOptionalFile("splicingout");

        if((gpeptidome == null) != (splicingout == null)) {
            System.err.printf("the parameter -gpeptidome and -splicingout make only sense together!\n");
            return;
        }
        Set<String> NEEDED_GPEPTIDOME_HEADERS = toSet("gene", "id", "peptide");
        if(gpeptidome != null && splicingout != null) {
            HashSet<String> gpheaders = toSet(FileUtils.getHeaders(gpeptidome));

            if(SetInfo.intersect(NEEDED_GPEPTIDOME_HEADERS, gpheaders).size() != NEEDED_GPEPTIDOME_HEADERS.size()) {
                System.err.printf("some needed headers from gpeptidome file (%s) are missing: %s!", gpeptidome.getAbsolutePath(), SetInfo.minus(NEEDED_GPEPTIDOME_HEADERS, gpheaders));
                return;
            }

        }

        double MINFC = cmd.getDouble("minfc");
        HashSet<String> trues = (!cmd.isOptionSet("trues")) ? null : FileUtils.readSet(cmd.getFile("trues"));
        HashSet<String> falses = (!cmd.isOptionSet("falses")) ? null : FileUtils.readSet(cmd.getFile("falses"));
        File pgroup = cmd.getOptionalFile("pgroup");

        HashMap<String, String> labelmapping = FileUtils.readMap(cmd.getFile("labels"), "\t", 0, 1);



        empires.input.ProteinGroupInfo pgi = new ProteinGroupInfo(pgroup, labelmapping);

        log.info("COND2LABELS: %s\n", pgi.getCondition2ReplicateMap());
        log.info("samplelist: %s\n", pgi.getSampleList());
        File peps = cmd.getOptionalFile("peptides");
        if (peps != null) {
            pgi.readPeptides(peps, false);
        }

        if (cmd.isOptionSet("useonlypeps")) {
            pgi.restrictToPeptides(FileUtils.readSet(cmd.getFile("useonlypeps")));
        }

        if (trues != null && falses != null) {
            apply(pgi.id2info.values(), (_h) -> _h.checkMultiOrg(trues, falses));
            pgi.setFalses(falses);
            pgi.setTrues(trues);
            pgi.setTrueFalseLabels();
        }


        BiFunction<UPair<Vector<Double>>, Integer, UPair<Vector<String>>> toNV = (_t, minidx) ->
        {
            double min2 = _t.getSecond().get(minidx);
            //Math.min(NumUtils.min(filter(_t.getFirst(), (_d) -> _d > 0.0)), NumUtils.min(filter(_t.getSecond(), (_d) -> _d > 0.0)));

            return UPair.createU(map(_t.getFirst(), (_d) -> String.format("%.2f", _d / min2))
                    , map(_t.getSecond(), (_d) -> String.format("%.2f", _d / min2)));

        };

        /* OPT CHECK* */

        log.info("pgi id2info sizes: %s", NumUtils.getNumInfo(pgi.id2info.values(), (_p) -> _p.getNumPepIds()).getInfoWithQ());


        Set<String> ids = filterToSet(pgi.id2info.keySet(), (_i) -> !pgi.id2info.get(_i).toFilter());
        log.info("got %d peps filtered: %d", pgi.id2info.size(), ids.size());

        String cond1 = cmd.getValue("cond1");
        String cond2 = cmd.getValue("cond2");

        for(String cond : toVector(cond1, cond2)) {
            if (pgi.getCondition2ReplicateMap().get(cond) != null)
                continue;

            System.err.printf("unknown condition: %s! known conditions: %s\n", cond, pgi.getCondition2ReplicateMap().keySet());
            return;
        }

        ReplicateSetInfo cond1repset = pgi.getProteinToPeps(cond1, false);
        ReplicateSetInfo cond2repset = pgi.getProteinToPeps(cond2, false);

        Set<String> c1peps = cond1repset.filterFeatureds((_id, _values) -> filteredSize(_values, (_d) -> !Double.isNaN(_d)) < minnummeasuredpercondition);
        Set<String> c2peps = cond2repset.filterFeatureds((_id, _values) -> filteredSize(_values, (_d) -> !Double.isNaN(_d)) < minnummeasuredpercondition);

        Set<String> oneSideOkPeps = filterToSet(cond1repset.getFeatureNames(),
                (_n) ->
                {
                    int numMeasured1 = filteredSize(cond1repset.getReplicateData(_n), (_d) -> !Double.isNaN(_d));
                    int numMeasured2 = filteredSize(cond2repset.getReplicateData(_n), (_d) -> !Double.isNaN(_d));

                    return (numMeasured1 > 0 && numMeasured2 > 0 && Math.max(numMeasured1, numMeasured2) > 1);
                });

        Set<String> commonPeps = oneSideOkPeps; //PepsSetInfo.intersect(c1peps, c2peps);

        log.info("c1peps: %d/%d c2 peps: %d/%d common: %d\n", c1peps.size(), cond1repset.getFeatureNames().size(), c2peps.size(), cond2repset.getFeatureNames().size(),
                commonPeps.size());

        ReplicateSetInfo restricted = cond1repset.restrictToSubFeatures(commonPeps);


        Set<String> selids = restricted.getFeatureNamesToTest();

        if(cmd.isOptionSet("useonlyproteins")) {
            selids = FileUtils.readSet(cmd.getFile("useonlyproteins"));
        }

        log.info("sel ids: %d/%d\n", selids.size(), cond1repset.getFeatureNamesToTest().size());

        log.info("first ids to test: %s\n", VectorUtils.slice(toVector(restricted.getFeatureNamesToTest()), 0, 10));
        empires.EmpiRe empiRe = new EmpiRe();
        empires.BackgroundContextFuzzficationStrategyProvider backgroundContextFuzzficationStrategyProvider = backgroundProviderOption.getStrategy();
        empires.NormalizedReplicateSet normc1 = new empires.NormalizedReplicateSet(cond1repset, backgroundContextFuzzficationStrategyProvider);
        empires.NormalizedReplicateSet normc2 = new empires.NormalizedReplicateSet(cond2repset, backgroundContextFuzzficationStrategyProvider);



        if(gpeptidome != null && splicingout != null) {

            new empires.FeatureBasedSplicingTest(normc1, normc2, gpeptidome, empires.DoubleDiffVariant.ALLPAIRS, MINPEPS_SPLICING, null).getSplicingResultTable().writeCSV(splicingout);
            /*
            Set<String> allpeps = SetInfo.intersect(toSet(cond1repset.featureNames), toSet(cond2repset.featureNames));

            HashMap<String, Tuple> peptide2eqclass = new HashMap<>();
            HashMap<String, HashSet<String>> gene2peptides = new HashMap<>();

            HashSet<String> nonUniqPeptides = new HashSet<>();

            apply(FileUtils.getHrIterator(gpeptidome, FileUtils.getHeaderedReader("\t").add("gene", "id", "peptide")),
                    (_hr) -> {
                            String pep = _hr.getString("peptide");
                            if(!allpeps.contains(pep))
                                return;

                            if(peptide2eqclass.containsKey(pep)) {
                                nonUniqPeptides.add(pep);
                                return;
                            }
                            Tuple eqC = Tuple.tupleFromCollection(toSortedVector(toVector(_hr.getString("id").split(",")), true));
                            peptide2eqclass.put(pep, eqC);
                            MapBuilder.update(gene2peptides, _hr.getString("gene"), pep);
                    });

            if(nonUniqPeptides.size() > 0) {
                log.warn("found %d/%d non-gene uniqly mapped peptides (%d genes)! (ignore them)", nonUniqPeptides.size(), peptide2eqclass.size(), gene2peptides.size());
            }

            RNASeqSplicingInfo splicingInfo = new RNASeqSplicingInfo(cmd.getInt("maxisoforms"), -1, pgi.cond2labels);


            DoubleDiffManager diffManager = new EmpiRe().getDoubleDiffManager(normc1, normc2);


            Vector<Integer> eqlassSizes = new Vector<>();

            Vector<DoubleDiffResult> splicingResults = new Vector<>();
            for(String gene : gene2peptides.keySet()) {

                Vector<String> uniq_peps  = filter(gene2peptides.get(gene), (_p) -> !nonUniqPeptides.contains(_p));
                System.out.printf("%s got %d/%d uniq peps\n", gene, uniq_peps.size(), gene2peptides.get(gene).size());
                if(uniq_peps.size() == 0)
                    continue;

                HashMap<Tuple, Vector<String>> eq2peps = new HashMap<>();
                for(String up : uniq_peps) {
                    MapBuilder.updateV(eq2peps, peptide2eqclass.get(up), up);
                }

                eqlassSizes.add(eq2peps.size());
                if(eq2peps.size() < 2)
                    continue;


                System.out.printf("%s: %d eqclasses: %s\n", gene, eq2peps.size(), eq2peps);

                Vector<Tuple> tkeys = toVector(eq2peps.keySet());

                DoubleDiffResult bestHit = null;
                for(UPair<Tuple> t : getPairs(tkeys, true)) {
                    DoubleDiffResult ddr = diffManager.getDoubleDiffResult(gene+"."+t, eq2peps.get(t.getFirst()), eq2peps.get(t.getSecond()));
                    System.out.printf("%s vs %s: %g\n", eq2peps.get(t.getFirst()), eq2peps.get(t.getSecond()), ddr.pval);
                    if(bestHit == null || bestHit.pval > ddr.pval) {
                        bestHit = ddr;
                    }
                }
                splicingResults.add(bestHit);
            }


            System.out.printf("eqclass distrib: %s\n", NumUtils.getNumInfo(eqlassSizes).getInfoWithQ());
            BenchmarkGene.getTable(splicingResults, null, null, splicingInfo).writeCSV(splicingout);

             */

        }
        Set<String> proteinsToTest = SetInfo.intersect(ids, selids);
        Vector<empires.DiffExpResult> res = empiRe.getDifferentialResults(normc1, normc2, proteinsToTest, (_c) -> restricted.getSubFeatures(_c));

        res = filter(res, (_de) -> _de.featureNames.size() >= MINPEPS);


        if(cmd.isOptionSet("trues")) {

            Vector<empires.DiffExpResult> trueResults = filter(res, (_r) -> pgi.id2info.get(_r.combinedFeatureName).is_true);
            System.out.println(empires.DiffExpResult.getPerformanceString(res, (_n) -> pgi.id2info.get(_n).is_true));


        }

        empires.plotting.DiffExpTable table = new DiffExpTable(normc1, normc2, res, MINFC, 0.05, (!cmd.isOptionSet("trues")) ? (_id) -> false : (_id) -> pgi.id2info.get(_id).is_true);

        table.getTable().writeCSV(cmd.getFile("o"));

        if(cmd.isSet("interactive")) {
            table.showInteractiveTable().setVisible(true);
        }
    }
}

