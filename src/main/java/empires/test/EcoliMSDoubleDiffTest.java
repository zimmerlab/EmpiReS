package empires.test;

import empires.EmpiRe;
import empires.input.BackgroundProviderOption;
import empires.input.GeneralOptions;
import empires.input.ProteinGroupInfo;
import empires.input.ReplicateSetInfo;
import lmu.utils.*;
import lmu.utils.fdr.BenjaminiHochberg;

import java.io.File;
import java.io.PrintWriter;
import java.util.*;
import java.util.function.BiFunction;

import static lmu.utils.ObjectGetter.*;
import static lmu.utils.ObjectGetter.filter;

public class EcoliMSDoubleDiffTest {

    public static void main(String[] args) {
        Locale.setDefault(Locale.UK);
        empires.input.GeneralOptions generalOptions = new GeneralOptions();
        empires.input.BackgroundProviderOption backgroundProviderOption = new BackgroundProviderOption();
        SimpleOptionParser cmd = new SimpleOptionParser("pgroup", "labels", "peptides", "o", "variant", "minpeps");

        cmd.setFile("pgroup", "labels", "peptides");
        cmd.setInt("minpeps");
        cmd.setDefault("minpeps", "8");
        cmd.setDefault("variant", empires.DoubleDiffVariant.MISSINGVAL.getName());

        if (!OptionParser.parseParams(args, true, false, true, true, cmd, backgroundProviderOption, generalOptions))
            return;

        int MINPEPS = cmd.getInt("minpeps");
        Logger log = LogConfig.getLogger();
        empires.DoubleDiffVariant doubleDiffVariant = empires.DoubleDiffVariant.get(cmd.getValue("variant"));

        generalOptions.apply();
        File pgroup = cmd.getOptionalFile("pgroup");

        HashMap<String, String> labelmapping = FileUtils.readMap(cmd.getFile("labels"), "\t", 0, 1);
        empires.input.ProteinGroupInfo pgi = new ProteinGroupInfo(pgroup, labelmapping);

        log.info("COND2LABELS: %s\n", pgi.getCondition2ReplicateMap());
        log.info("samplelist: %s\n", pgi.getSampleList());
        File peps = cmd.getOptionalFile("peptides");
        if (peps != null) {
            pgi.readPeptides(peps, false);
        }

        Vector<empires.input.ReplicateSetInfo> replicateSetInfos = map(pgi.getCondition2ReplicateMap().keySet(), (_c) -> pgi.getProteinToPeps(_c, false));


        for (ReplicateSetInfo rsi : replicateSetInfos) {
            System.out.printf("%s %d proteins pep distrib: %s\n", rsi.getReplicateSetName(), rsi.getFeatureNamesToTest().size(),
                    NumUtils.getNumInfo(rsi.getFeatureNamesToTest(), (_s) -> rsi.getSubFeatures(_s).size()).getInfoWithQ());
        }

        Vector<empires.NormalizedReplicateSet> normalizedReplicateSets = map(replicateSetInfos, (_rsi) -> new empires.NormalizedReplicateSet(_rsi));


        int NUM_TESTS = 10_000;

        PrintWriter pw = cmd.getWriter("o");

        Vector<Pair<empires.DoubleDiffManager, String>> diffsplicprots = new Vector<>();
        for(UPair<empires.NormalizedReplicateSet> up : getPairs(normalizedReplicateSets, true)) {
            Set<String> commonProts = SetInfo.intersect(up.getFirst().getInData().getFeatureNamesToTest(), up.getSecond().getInData().getFeatureNamesToTest());
            Vector<String> protsToSelect = filter(commonProts, (_f) -> up.getFirst().getInData().getSubFeatures(_f).size() >= MINPEPS && up.getSecond().getInData().getSubFeatures(_f).size() >= MINPEPS);

            empires.EmpiRe emp = new EmpiRe();

            empires.DoubleDiffManager ddm = emp.getDoubleDiffManager(up.getFirst(), up.getSecond());
            diffsplicprots.addAll(map(protsToSelect, (_p) -> Pair.create(ddm, _p)));
        }

        log.info("got %d diffsplic selections!", diffsplicprots.size());


        for(int iter=0; iter < 50; iter++) {
            Vector<Double> tests = new Vector<>();
            Vector<Pair<empires.DoubleDiffManager, String>> shuffled = shuffleN(diffsplicprots, NUM_TESTS);
            for(Pair<empires.DoubleDiffManager, String> up : shuffled) {
                Set<String> pepset = new HashSet<>();
                String prot = up.getSecond();
                pepset.addAll(up.getFirst().getDiffExpManager().replicateSetFrom.getInData().getSubFeatures(prot));
                pepset.addAll(up.getFirst().getDiffExpManager().replicateSetTo.getInData().getSubFeatures(prot));

                Vector<String> pepvec = shuffle(toVector(pepset), false);
                int split_pos = (int)(pepvec.size() * (0.3 + (Math.random() * 0.4)));
                Vector<String> peps1 = VectorUtils.slice(pepvec, 0, split_pos);
                Vector<String> peps2 = VectorUtils.slice(pepvec, split_pos, pepvec.size());

                tests.add(up.getFirst().getDoubleDiffResult("test"+ tests.size(), peps1, peps2, doubleDiffVariant).pval);
            }

            Vector<Double> fdrs = new Vector<>();
            BenjaminiHochberg.adjust_pvalues(tests, (_pval) -> _pval, (_p) -> fdrs.add(_p.getSecond()));

            int nsignif = filteredSize(fdrs, (_fdr) -> _fdr <= 0.05);
            int nsignif01 = filteredSize(fdrs, (_fdr) -> _fdr <= 0.01);
            pw.printf("iter %d %d/%d = %.2f%% signifs at 0.05 %d/%d = %.2f%% at 0.01\n", iter+1, nsignif, tests.size(), (100.0 * nsignif) / tests.size(),
                    nsignif01, tests.size(), (100.0 * nsignif01) / tests.size());

            pw.flush();
            log.info("iter:%d %d/%d = %.2f%% signifs\n", iter + 1, nsignif, tests.size(), (100.0 * nsignif) / tests.size());
        }

        Vector<Double> tests = new Vector<>();
        for(Pair<empires.DoubleDiffManager, String> up : diffsplicprots) {
            Set<String> pepset = new HashSet<>();
            String prot = up.getSecond();
            pepset.addAll(up.getFirst().getDiffExpManager().replicateSetFrom.getInData().getSubFeatures(prot));
            pepset.addAll(up.getFirst().getDiffExpManager().replicateSetTo.getInData().getSubFeatures(prot));

            Vector<String> pepvec = shuffle(toVector(pepset), false);
            int split_pos = (int)(pepvec.size() * (0.3 + (Math.random() * 0.4)));
            Vector<String> peps1 = VectorUtils.slice(pepvec, 0, split_pos);
            Vector<String> peps2 = VectorUtils.slice(pepvec, split_pos, pepvec.size());

            tests.add(up.getFirst().getDoubleDiffResult("test"+ tests.size(), peps1, peps2, doubleDiffVariant).pval);
        }
        Vector<Double> fdrs = new Vector<>();
        BenjaminiHochberg.adjust_pvalues(tests, (_pval) -> _pval, (_p) -> fdrs.add(_p.getSecond()));

        int nsignif = filteredSize(fdrs, (_fdr) -> _fdr <= 0.05);
        int nsignif01 = filteredSize(fdrs, (_fdr) -> _fdr <= 0.01);
        pw.printf("merged %d/%d = %.2f%% signifs at 0.05 %d/%d = %.2f%% at 0.01\n",  nsignif, tests.size(), (100.0 * nsignif) / tests.size(),
                nsignif01, tests.size(), (100.0 * nsignif01) / tests.size());

        pw.flush();
        log.info("merged %d/%d = %.2f%% signifs\n", nsignif, tests.size(), (100.0 * nsignif) / tests.size());

        pw.close();
    }
}


