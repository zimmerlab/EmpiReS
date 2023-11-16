package empires.test;

import empires.SparseCumulativeDistribution;
import lmu.utils.*;
import lmu.utils.plotting.PlotCreator;
import lmu.utils.tuple.Tuple4;
import empires.ErrorEstimationDistribution;
import lmu.utils.plotting.CachedPlotCreator;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.awt.image.BufferedImage;
import java.util.HashMap;
import java.util.Vector;

import static lmu.utils.IteratorUtils.rangev;
import static lmu.utils.ObjectGetter.*;

public class PepVarTest {

    static double MISSING_VALUE_PROB = 0.2;
    static NormalDistribution SDND = new NormalDistribution(0.4, 0.4);
    static class PepSource {
        NormalDistribution ndC1;
        NormalDistribution ndC2;

        Vector<Boolean> missing_rep1;
        Vector<Boolean> missing_rep2;

        double c1var;
        double c2var;
        double var;

        public PepSource(double sd1, double sd2, int nrep1, int nrep2) {
            missing_rep1 = mapIndex(nrep1, (_r) -> Math.random() < MISSING_VALUE_PROB);
            missing_rep2 = mapIndex(nrep2, (_r) -> Math.random() < MISSING_VALUE_PROB);
            ndC1 = new NormalDistribution(0.0, sd1);
            ndC2 = new NormalDistribution(0.0, sd2);
            c1var = Math.pow(sd1, 2.0);
            c2var = Math.pow(sd2, 2.0);
            var= c1var + c2var;
        }

        public PepSource(int nrep1, int nrep2) {
            this(0.2 + Math.abs(SDND.sample()), 0.2 + Math.abs(SDND.sample()), nrep1, nrep2);
        }


        UPair<Vector<Double>> sample() {
            return UPair.createU(map(missing_rep1, (_b) -> (_b) ? Double.NaN : ndC1.sample()), map(missing_rep2, (_b) -> (_b) ? Double.NaN : ndC2.sample()));
        }

        public String getShortInfo() {
            //return String.format("<m:%d,%d,sd:%.2f,%.2f>", filteredSize(missing_rep1, (_r) -> _r), filteredSize(missing_rep2, (_r) -> _r), Math.sqrt(c1var), Math.sqrt(c2var));
            return String.format("<%d,%d>", filteredSize(missing_rep1, (_r) -> _r), filteredSize(missing_rep2, (_r) -> _r));
        }

        public String toString() {
            return String.format("p: %.2f %s + %.2f %s", ndC1.getStandardDeviation(), missing_rep1,  ndC2.getStandardDeviation(), missing_rep2);
        }

    }

    static Vector<Integer> getCommonIndeces(Vector<Double> v1, Vector<Double> v2) {
        Vector<Integer> rv = new Vector<>();
        for(int i=0; i<v1.size(); i++) {
            if(Double.isNaN(v1.get(i)) || Double.isNaN(v2.get(i)))
                continue;

            rv.add(i);
        }

        return rv;
    }

    static void corrTestAllPairs(int numPeps1, int numPeps2, int nreps1, int nreps2, boolean use_correlation, PlotCreator pc) {
        Vector<PepSource> iso1 = mapIndex(numPeps1, (_i) -> new PepSource(nreps1, nreps2));
        Vector<PepSource> iso2 = mapIndex(numPeps2, (_i) -> new PepSource(nreps1, nreps2));

        System.out.printf("%s\n%s\n", iso1, iso2);

        int NTESTS = 30000;


        Vector<UPair<PepSource>> pepPairs = getPairs(iso1, iso2);
        Vector<Double> pvals = new Vector<>();

        for(int testi=0; testi<NTESTS; testi++) {
            Vector<UPair<Vector<Double>>> iso1_samples = map(iso1, (_i) -> _i.sample());
            Vector<UPair<Vector<Double>>> iso2_samples = map(iso2, (_i) -> _i.sample());


            double sumEvidences = 0;
            int numObservations = 0;

            double totalVariance = 0;

            HashMap<PepSource, Vector<UPair<PepSource>> > pep2pairs = new HashMap<>();
            HashMap<UPair<PepSource>, Vector<Vector<Integer>>> peppair2indeces = new HashMap<>();

            class UsableReplicates {
                Vector<Integer> rep1indeces;
                Vector<Integer> rep2indeces;
                int nc1;
                int nc2;

                UsableReplicates(UPair<Vector<Double>> data) {
                    rep1indeces = filter(rangev(data.getFirst().size()), (_i) -> !Double.isNaN(data.getFirst().get(_i)));
                    rep2indeces = filter(rangev(data.getSecond().size()), (_i) -> !Double.isNaN(data.getSecond().get(_i)));
                    nc1 = rep1indeces.size();
                    nc2 = rep2indeces.size();
                }

            }
            HashMap<Pair<Boolean, Integer>, UsableReplicates> pep2usable = new HashMap<>();
            for(int pep1idx = 0; pep1idx < iso1.size(); pep1idx++) {
                pep2usable.put(Pair.create(true, pep1idx), new UsableReplicates(iso1_samples.get(pep1idx)));
            }
            for(int pep2idx = 0; pep2idx < iso2.size(); pep2idx++) {
                pep2usable.put(Pair.create(false, pep2idx), new UsableReplicates(iso2_samples.get(pep2idx)));
            }


            HashMap<UPair<PepSource>, Vector<Vector<Integer>>> peppair2used = new HashMap<>();


            for(int pep1idx = 0; pep1idx < iso1.size(); pep1idx++) {
                UPair<Vector<Double>> pep1 = iso1_samples.get(pep1idx);
                PepSource pepSource1 = iso1.get(pep1idx);

                UsableReplicates ur1 = pep2usable.get(Pair.create(true, pep1idx));

                for(int pep2idx = 0; pep2idx < iso2.size(); pep2idx++) {
                    PepSource pepSource2 = iso2.get(pep2idx);
                    UPair<Vector<Double>> pep2 = iso2_samples.get(pep2idx);
                    UsableReplicates ur2 = pep2usable.get(Pair.create(false, pep2idx));

                    UPair<PepSource> pairKey = UPair.createU(pepSource1, pepSource2);

                    MapBuilder.updateV(pep2pairs, pairKey.getFirst(), pairKey);
                    MapBuilder.updateV(pep2pairs, pairKey.getSecond(), pairKey);

                    double corrNorm = (pepSource1.var + pepSource2.var);
                    double corrSD = Math.sqrt(corrNorm);

                    int nObs_pep_pair = (ur1.nc1 * ur1.nc2 * ur2.nc1 * ur2.nc2);
                    double pairVariance = nObs_pep_pair * (ur1.nc1 * ur1.nc2) * (ur2.nc2 ) * (pepSource2.c1var) +
                            nObs_pep_pair * (ur1.nc1 * ur1.nc2) * (ur2.nc1 ) * (pepSource2.c2var) +

                            nObs_pep_pair * (ur2.nc1 * ur2.nc2) * (ur1.nc2 ) * (pepSource1.c1var) +
                            nObs_pep_pair * ( ur2.nc1 * ur2.nc2) * (ur1.nc1 ) * (pepSource1.c2var);

                    // (a, b)
                    totalVariance += (use_correlation) ? pairVariance / corrNorm : pairVariance;

                    Vector<Vector<Integer>>  indeces = toVector(ur1.rep1indeces, ur1.rep2indeces, ur2.rep1indeces, ur2.rep2indeces);

                    peppair2indeces.put(pairKey, indeces);

                    for (int rep11idx : ur1.rep1indeces) {
                        for (int rep12idx : ur1.rep2indeces) {
                            for (int rep21idx : ur2.rep1indeces) {
                                for (int rep22idx : ur2.rep2indeces) {
                                    double diffdiffFC = (pep1.getFirst().get(rep11idx) - pep1.getSecond().get(rep12idx)) - (pep2.getFirst().get(rep21idx) - pep2.getSecond().get(rep22idx));
                                    sumEvidences += (use_correlation) ? diffdiffFC / corrSD : diffdiffFC;
                                    numObservations++;
                                }

                            }
                        }
                    }
                }
            }



            for(UPair<PepSource> pair : pepPairs) {
                Vector<Vector<Integer>> indeces = peppair2indeces.get(pair);
                if(indeces == null)
                    continue;

                int numObs1 = NumUtils.product(map(indeces, (_v) -> _v.size()));

                double pair1SD = Math.sqrt(pair.getFirst().var + pair.getSecond().var);

                for(int overlapIdx = 0; overlapIdx < 2; overlapIdx++) {
                    PepSource center = pair.get(overlapIdx == 0);

                    double c1var = center.c1var;
                    double c2var = center.c2var;


                    for(UPair<PepSource> overlap : filter(pep2pairs.get(center), (_pair) -> !_pair.equals(pair))) {



                        double pair2SD = Math.sqrt(overlap.getFirst().var + overlap.getSecond().var);

                        Vector<Vector<Integer>> indeces2 = peppair2indeces.get(overlap);
                        if(indeces2 == null)
                            continue;

                        int numObs2 = NumUtils.product(map(indeces2, (_v) -> _v.size()));

                        double corrNorm = pair1SD * pair2SD;
                        double p1p2_repOverlap = 0;

                        for(int idxIdx = (overlapIdx == 0) ? 0 : 2, j=0; idxIdx < ((overlapIdx == 0) ? 2 : 4); idxIdx++, j++) {

                            double var = (j == 0) ? c1var : c2var;
                            int ncommon = SetInfo.intersect(indeces.get(idxIdx), indeces2.get(idxIdx)).size();
                            if(ncommon == 0)
                                continue;

                            int noverlaps = (numObs1 / indeces.get(idxIdx).size()) * (numObs2 / indeces2.get(idxIdx).size()) * ncommon;

                            p1p2_repOverlap += (noverlaps * var);
                        }
                        totalVariance += (use_correlation) ? p1p2_repOverlap / corrNorm : p1p2_repOverlap;

                    }

                }

            }



            if(numObservations == 0)
                continue;

            //System.out.printf("peppairs:  %s overlaps: %.3f\n", pepPairs, peptidePairOverlaps);
            //System.out.printf("nobs: %d var: %.2f + %.2f + %.2f\n", numObservations, diagonalVariance,  peptidePairOverlaps, replicatePairOverlaps);
            //System.out.printf("observ: %d total var: %.2f sumE: %.2f\n", numObservations, totalVariance, sumEvidences);
            NormalDistribution pND = new NormalDistribution(0, Math.sqrt(totalVariance));
            //pvals.add(totalND.cumulativeProbability(sumEvidences));
            pvals.add(pND.cumulativeProbability(sumEvidences));


        }

        if(pvals.size() == 0)
            return;
        System.out.printf("pvals: %s\n", NumUtils.getNumInfo(pvals).getInfoWithQ());
        //pc.cumhist(String.format("p: %dx%d r: %dx%d", numPeps1, numPeps2, nreps1, nreps2), pvals, pvals.size(), false, true);
        pc.cumhist(String.format("%s, %s"
                , StringUtils.joinObjects(",", iso1, (_i) -> _i.getShortInfo())
                , StringUtils.joinObjects(",", iso2, (_i) -> _i.getShortInfo())),
                pvals, pvals.size(), false, true);



    }


    static void corrTest(int numPeps1, int numPeps2, int nreps1, int nreps2, boolean use_correlation, PlotCreator pc) {
        Vector<PepSource> iso1 = mapIndex(numPeps1, (_i) -> new PepSource(nreps1, nreps2));
        Vector<PepSource> iso2 = mapIndex(numPeps2, (_i) -> new PepSource(nreps1, nreps2));

        System.out.printf("%s\n%s\n", iso1, iso2);

        int NTESTS = 30000;


        Vector<UPair<PepSource>> pepPairs = getPairs(iso1, iso2);
        Vector<Double> pvals = new Vector<>();

        for(int testi=0; testi<NTESTS; testi++) {
            Vector<UPair<Vector<Double>>> iso1_samples = map(iso1, (_i) -> _i.sample());
            Vector<UPair<Vector<Double>>> iso2_samples = map(iso2, (_i) -> _i.sample());


            double sumEvidences = 0;
            int numObservations = 0;

            double totalVariance = 0;

            HashMap<UPair<PepSource>, UPair<Vector<Integer>>> peppair2commonindeces = new HashMap<>();

            HashMap<PepSource, Vector<UPair<PepSource>>> pep2pairs = new HashMap<>();


            for(int pep1idx = 0; pep1idx < iso1.size(); pep1idx++) {
                UPair<Vector<Double>> pep1 = iso1_samples.get(pep1idx);
                PepSource pepSource1 = iso1.get(pep1idx);

                for(int pep2idx = 0; pep2idx < iso2.size(); pep2idx++) {
                    PepSource pepSource2 = iso2.get(pep2idx);
                    UPair<Vector<Double>> pep2 = iso2_samples.get(pep2idx);

                    Vector<Integer> cond1_indeces = getCommonIndeces(pep1.getFirst(), pep2.getFirst());
                    Vector<Integer> cond2_indeces = getCommonIndeces(pep1.getSecond(), pep2.getSecond());

                    int n1 = cond1_indeces.size();
                    int n2 = cond2_indeces.size();

                    if(n1 == 0 || n2 == 0)
                        continue;

                    UPair<PepSource> pairKey = UPair.createU(pepSource1, pepSource2);
                    peppair2commonindeces.put(pairKey, UPair.createU(cond1_indeces, cond2_indeces));

                    MapBuilder.updateV(pep2pairs, pairKey.getFirst(), pairKey);
                    MapBuilder.updateV(pep2pairs, pairKey.getSecond(), pairKey);

                    double corrNorm = (pepSource1.var + pepSource2.var);
                    double corrSD = Math.sqrt(corrNorm);
                    double pairVariance = (n1 * n2) * (pepSource1.var + pepSource2.var);
                    pairVariance += (n1 * n2) * (n2 - 1) * (pepSource1.c1var + pepSource2.c1var) +
                                              (n1 * n2) * (n1 - 1) * (pepSource1.c2var + pepSource2.c2var);

                    totalVariance += (use_correlation) ? pairVariance / corrNorm : pairVariance;

                    for(int rep1idx : cond1_indeces) {
                        for(int rep2idx : cond2_indeces) {
                            double diffdiffFC = (pep1.getFirst().get(rep1idx) - pep1.getSecond().get(rep2idx)) - (pep2.getFirst().get(rep1idx) - pep2.getSecond().get(rep2idx));
                            sumEvidences += (use_correlation) ? diffdiffFC / corrSD : diffdiffFC;
                            numObservations++;
                        }
                    }
                }
            }

            for(UPair<PepSource> pair : pepPairs) {
                UPair<Vector<Integer>> indeces = peppair2commonindeces.get(pair);
                if(indeces == null)
                    continue;

                double pair1SD = Math.sqrt(pair.getFirst().var + pair.getSecond().var);

                for(int overlapIdx = 0; overlapIdx < 2; overlapIdx++) {
                    PepSource center = pair.get(overlapIdx == 0);

                    double c1var = center.c1var;
                    double c2var = center.c2var;

                    for(UPair<PepSource> overlap : filter(pep2pairs.get(center), (_pair) -> !_pair.equals(pair))) {

                        double pair2SD = Math.sqrt(overlap.getFirst().var + overlap.getSecond().var);

                        UPair<Vector<Integer>> indeces2 = peppair2commonindeces.get(overlap);
                        if(indeces2 == null)
                            continue;

                        double corrNorm = pair1SD * pair2SD;
                        int n1 = SetInfo.intersect(indeces.getFirst(), indeces2.getFirst()).size();
                        int n2 = SetInfo.intersect(indeces.getSecond(), indeces2.getSecond()).size();

                        double p1p2_repOverlap = indeces.getSecond().size() * indeces2.getSecond().size() * (n1) * (c1var);
                        p1p2_repOverlap += indeces.getFirst().size() * indeces2.getFirst().size() * (n2) * (c2var);
                        totalVariance += (use_correlation) ? p1p2_repOverlap / corrNorm : p1p2_repOverlap;

                    }

                }

            }

            if(numObservations == 0)
                continue;

            //System.out.printf("peppairs:  %s overlaps: %.3f\n", pepPairs, peptidePairOverlaps);
            //System.out.printf("nobs: %d var: %.2f + %.2f + %.2f\n", numObservations, diagonalVariance,  peptidePairOverlaps, replicatePairOverlaps);
            NormalDistribution pND = new NormalDistribution(0, Math.sqrt(totalVariance));
            //pvals.add(totalND.cumulativeProbability(sumEvidences));
            pvals.add(pND.cumulativeProbability(sumEvidences));


        }

        //pc.cumhist(String.format("p: %dx%d r: %dx%d", numPeps1, numPeps2, nreps1, nreps2), pvals, pvals.size(), false, true);
        pc.cumhist(String.format("%s, %s"
                , StringUtils.joinObjects(",", iso1, (_i) -> _i.getShortInfo())
                , StringUtils.joinObjects(",", iso2, (_i) -> _i.getShortInfo())),
                pvals, pvals.size(), false, true);



    }

    static double sample(Vector<NormalDistribution> d) {
        return NumUtils.mean(map(d, (_nd) -> _nd.sample()));
    }

    static double sampleE(Vector<Double> vals) {
        return vals.get((int)(Math.random() * vals.size()));
    }

    public static void checkDependenceDistrib() {
        HashMap<Integer, Long> obs = new HashMap<>();
        ErrorEstimationDistribution.setFcWithBin(0.0005);
        Vector<NormalDistribution> tests = toVector(new NormalDistribution(0, 1.2), new NormalDistribution(0, 1.2),
                new NormalDistribution(-1.0, 0.5), new NormalDistribution(1.0, 0.5),
                new NormalDistribution(-3.0, 0.6), new NormalDistribution(3.0, 0.6));

        Vector<Double> tvals = new Vector<>();
        for(NormalDistribution nd : tests) {
            for(int i=0; i<10_000; i++) {
                double sample = nd.sample();
                tvals.add(sample);
                MapBuilder.update(obs, (int)(sample * ErrorEstimationDistribution.INT_FACTOR), 1l);
            }
        }
        ErrorEstimationDistribution errd = new ErrorEstimationDistribution(obs, true);
        errd.getCumulative(true);
        ErrorEstimationDistribution diffErr = ErrorEstimationDistribution.substract(errd, errd);
        diffErr.getCumulative(true);
        PlotCreator pc = CachedPlotCreator.getPlotCreator();
        new empires.SparseCumulativeDistribution(errd).drawLine(pc, "errd");
        new empires.SparseCumulativeDistribution(diffErr).drawLine(pc, "differr");

        int nreps1 = 3;
        int nreps2 = 4;

        double perEvidenceVariance = diffErr.getVariance()
                + (nreps2  - 1 ) * errd.getVariance()
                + (nreps1  - 1 ) * errd.getVariance();


        double totalVariance = perEvidenceVariance * nreps1 * nreps2;

        ErrorEstimationDistribution scaled = diffErr.scale(Math.sqrt(totalVariance / diffErr.getVariance()));

        scaled.getCumulative(true);

        ImageUtils.showImage("test + " + errd.getVariance() , pc.getImage());

        Vector<Double> vals = new Vector<>();
        for(int i=0; i<10_000; i++) {
            Vector<Double> rep1 = mapIndex(nreps1, (_d) -> sampleE(tvals));
            Vector<Double> rep2 = mapIndex(nreps2, (_d) -> sampleE(tvals));

            vals.add(NumUtils.sum(getPairs(rep1, rep2), (_p) -> _p.getFirst() - _p.getSecond()));
        }

        System.out.printf("vals info: %s\n", NumUtils.getNumInfo(vals));
        System.out.printf("sd of scaled: %.2f\n", scaled.getSD());
        pc.cumhist("sampled", vals, vals.size(), false, true);
        new empires.SparseCumulativeDistribution(diffErr).drawCumulative(pc, "diffErr", 100);
        new SparseCumulativeDistribution(scaled).drawCumulative(pc, "diffErr-scaled", 100);
        pc.setLabels("sum of evidences", "cumfreq", "bottomright");
        ImageUtils.showImage("check  "  , pc.getImage());
    }

    public static void main(String[] args) {

        SimpleOptionParser cmd = new SimpleOptionParser("usecorr", "checkdep", "allpairs", "singlepeps", "runspertest");
        cmd.setSwitches("usecorr", "checkdep", "allpairs", "singlepeps");
        cmd.setInt("runspertest");
        cmd.setDefault("runspertest", "10");

        if(!OptionParser.parseParams(args, true, false, true, true, cmd))
            return;

        if(cmd.isSet("checkdep")) {
            checkDependenceDistrib();
            return;
        }

        int RUNSPERTEST = cmd.getInt("runspertest");
        PlotCreator pc = CachedPlotCreator.getPlotCreator();

        Vector<Tuple4<Integer, Integer, Integer, Integer>> single_peppair_tests =  toVector(
                Tuple4.create(1, 1, 40, 3),
                Tuple4.create(1, 1, 5, 8),
                Tuple4.create(1, 1, 1, 2),
                Tuple4.create(1, 1, 1, 5),
                Tuple4.create(1, 1, 3, 3),
                Tuple4.create(1, 1, 7, 3),
                Tuple4.create(1, 1, 4, 4)
        );

        Vector<Tuple4<Integer, Integer, Integer, Integer>> tests = toVector(
                Tuple4.create(4, 2, 1, 2),
                Tuple4.create(4, 2, 1, 5)

                ,Tuple4.create(3, 3, 7, 3)
                ,Tuple4.create(3, 3, 3, 3)
                ,Tuple4.create(3, 13, 3, 13),
                Tuple4.create(3, 10, 40, 3)

        );

        tests.addAll(single_peppair_tests);
        Vector<BufferedImage> bims = new Vector<>();

        Vector<Tuple4<Integer, Integer, Integer, Integer>> TODO = (cmd.isSet("singlepeps")) ? single_peppair_tests : tests;
        boolean use_correlation = cmd.isSet("usecorr");
        //for(boolean use_correlation : toVector(true)) {
        for(Tuple4<Integer, Integer, Integer, Integer> t : TODO) {

            for(int i=0; i<RUNSPERTEST; i++) {
                if(cmd.isSet("allpairs")) {
                    corrTestAllPairs(t.get0(), t.get1(), t.get2(), t.get3(), use_correlation, pc);
                } else {
                    corrTest(t.get0(), t.get1(), t.get2(), t.get3(), use_correlation, pc);
                }

            }
            pc.setTitle("peps: %d x %d reps: %d % d use correlation: %s\nmissing prob: %.2f", t.get0(), t.get1(), t.get2(), t.get3(), use_correlation, MISSING_VALUE_PROB);
            pc.abline("", null, null, 0.0, 100.0);
            //pc.setLabels("pval", "cum.freq", "bottomright");
            pc.setLabels("pval", "cum.freq", null);

            bims.add(pc.getImage());
        }

        ImageUtils.showImage("test", ImageUtils.concat(bims));
    }
}
