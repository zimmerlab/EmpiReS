package empires.test;

import empires.input.BackgroundProviderOption;
import empires.input.ProteinGroupInfo;
import empires.input.ReplicateSetInfo;
import empires.plotting.DiffExpTable;
import lmu.utils.*;
import lmu.utils.fdr.PerformanceResult;
import empires.BackgroundContextFuzzficationStrategyProvider;
import empires.DiffExpResult;
import empires.EmpiRe;
import empires.NormalizedReplicateSet;

import java.io.File;
import java.util.*;
import java.util.function.BiFunction;

import static lmu.utils.ObjectGetter.*;

public class MaxLFQTest {

    public static void main(String[] args) {
        Locale.setDefault(Locale.UK);
        empires.input.BackgroundProviderOption backgroundProviderOption = new BackgroundProviderOption();
        SimpleOptionParser cmd = new SimpleOptionParser("trues", "falses", "pgroup", "labels", "condpair", "peptides", "od", "export",
                "useonlypeps", "onlyreport", "minfc", "tableout", "precheck", "msqrobids", "precalc", "ids", "showtable", "minpeps");

        cmd.setFile("trues", "pgroup", "labels", "peptides", "useonlypeps", "msqrobids");

        cmd.setDir("od", "export", "precalc");

        cmd.setDouble("minfc");
        cmd.setDefault("minfc", "0.2");
        cmd.setInt("minpeps");
        cmd.setDefault("minpeps", "1");
        cmd.setOptional("peptides",  "od", "export", "useonlypeps", "tableout", "msqrobids", "precalc", "ids");

        cmd.setSwitches("onlyreport", "precheck", "showtable");
        if (!OptionParser.parseParams(args, true, false, true, true, cmd, backgroundProviderOption))
            return;

        int MINPEPS = cmd.getInt("minpeps");

        File msqrobids = cmd.getOptionalFile("msqrobids");

        HashSet<String> msqrobid = FileUtils.readSet(msqrobids);

        File tout = cmd.getOptionalFile("tableout");

        double MINFC = cmd.getDouble("minfc");
        HashSet<String> trues = (!cmd.isOptionSet("trues")) ? null : FileUtils.readSet(cmd.getFile("trues"));
        HashSet<String> falses = (!cmd.isOptionSet("falses")) ? null : FileUtils.readSet(cmd.getFile("falses"));
        File pgroup = cmd.getOptionalFile("pgroup");

        HashMap<String, String> labelmapping = FileUtils.readMap(cmd.getFile("labels"), "\t", 0, 1);
        String[] cp = cmd.getValue("condpair").split(":");
        UPair<String> CP = UPair.createU(cp[0], cp[1]);
        empires.input.ProteinGroupInfo pgi = new ProteinGroupInfo(pgroup, labelmapping);

        System.out.printf("COND2LABELS: %s\n", pgi.getCondition2ReplicateMap());
        System.out.printf("samplelist: %s\n", pgi.getSampleList());
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

        System.out.printf("pgi id2info sizes: %s\n", NumUtils.getNumInfo(pgi.id2info.values(), (_p) -> _p.getNumPepIds()).getInfoWithQ());


        Set<String> ids = filterToSet(pgi.id2info.keySet(), (_i) -> !pgi.id2info.get(_i).toFilter());
        System.out.printf("got %d peps filterd: %d\n", pgi.id2info.size(), ids.size());

        empires.input.ReplicateSetInfo cond1repset = pgi.getProteinToPeps(CP.getFirst(), false);
        empires.input.ReplicateSetInfo cond2repset = pgi.getProteinToPeps(CP.getSecond(), false);

        Set<String> c1peps = cond1repset.filterFeatureds((_id, _values) -> filteredSize(_values, (_d) -> !Double.isNaN(_d)) < 2);
        Set<String> c2peps = cond2repset.filterFeatureds((_id, _values) -> filteredSize(_values, (_d) -> !Double.isNaN(_d)) < 2);

        Set<String> oneSideOkPeps = filterToSet(cond1repset.getFeatureNames(),
                (_n) ->
                {
                    int numMeasured1 = filteredSize(cond1repset.getReplicateData(_n), (_d) -> !Double.isNaN(_d));
                    int numMeasured2 = filteredSize(cond2repset.getReplicateData(_n), (_d) -> !Double.isNaN(_d));

                    return (numMeasured1 > 0 && numMeasured2 > 0 && Math.max(numMeasured1, numMeasured2) > 1);
                });

        Set<String> commonPeps = oneSideOkPeps; //PepsSetInfo.intersect(c1peps, c2peps);

        System.out.printf("c1peps: %d/%d c2 peps: %d/%d common: %d\n", c1peps.size(), cond1repset.getFeatureNames().size(), c2peps.size(), cond2repset.getFeatureNames().size(),
                commonPeps.size());

        ReplicateSetInfo restricted = cond1repset.restrictToSubFeatures(commonPeps);


        Set<String> selids = restricted.getFeatureNamesToTest();

        if(cmd.isOptionSet("ids")) {
            selids = FileUtils.readSet(cmd.getFile("ids"));
        }

        System.out.printf("sel ids: %d/%d\n", selids.size(), cond1repset.getFeatureNamesToTest().size());

        System.out.printf("first ids to test: %s\n", VectorUtils.slice(toVector(restricted.getFeatureNamesToTest()), 0, 10));
        EmpiRe empiRe = new EmpiRe();
        BackgroundContextFuzzficationStrategyProvider backgroundContextFuzzficationStrategyProvider = backgroundProviderOption.getStrategy();
        NormalizedReplicateSet normc1 = new NormalizedReplicateSet(cond1repset, backgroundContextFuzzficationStrategyProvider);
        NormalizedReplicateSet normc2 = new NormalizedReplicateSet(cond2repset, backgroundContextFuzzficationStrategyProvider);

        Set<String> proteinsToTest = SetInfo.intersect(ids, selids);
        Vector<DiffExpResult> res = empiRe.getDifferentialResults(normc1, normc2, proteinsToTest, (_c) -> restricted.getSubFeatures(_c));

        res = filter(res, (_de) -> _de.featureNames.size() >= MINPEPS);



        Vector<DiffExpResult> trueResults = filter(res, (_r) -> pgi.id2info.get(_r.combinedFeatureName).is_true);

        PerformanceResult pr = new PerformanceResult("nlEmpire", res,
                (_r) -> _r.pval, (_r) -> pgi.id2info.get(_r.combinedFeatureName).is_true,
                false, (_r) -> _r.fdr <= 0.05 && Math.abs(_r.estimatedFC) >= MINFC , null);

        PerformanceResult prFC = new PerformanceResult("nlEmpire.FCDist", res,
                (_r) -> _r.fcEstimatePval, (_r) -> pgi.id2info.get(_r.combinedFeatureName).is_true,
                false, (_r) -> _r.fcEstimateFDR <= 0.05 && Math.abs(_r.estimatedFC) >= MINFC , null);

        System.out.printf("res: %s\n", pr);
        System.out.printf("resFC: %s\n", prFC);


        if(cmd.isSet("showtable")) {
            empires.plotting.DiffExpTable table = new DiffExpTable(normc1, normc2, res, MINFC, 0.05, (_id) -> pgi.id2info.get(_id).is_true);
            table.showInteractiveTable().setVisible(true);
        }
    }
}
