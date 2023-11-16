package empires.test;


import lmu.utils.DataTable;
import lmu.utils.ImageUtils;
import lmu.utils.NumUtils;
import lmu.utils.StringUtils;
import lmu.utils.plotting.PlotCreator;
import lmu.utils.swing.PagedDataTable;
import empires.DistribUtils;
import empires.ErrorEstimationDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.junit.jupiter.api.Test;

import java.util.TreeMap;
import java.util.Vector;
import java.util.function.Function;

import static lmu.utils.ObjectGetter.map;
import static lmu.utils.ObjectGetter.mapIndex;
import static org.junit.jupiter.api.Assertions.assertTrue;

public class DistSumTest {

    @Test
    public void TreeMapTest() {

        TreeMap<Long, Integer> m = new TreeMap<>();

        for (int i = 0; i < 1_000_000; i++) {
            long val = (long) (Math.random() * Long.MAX_VALUE);
            m.put(val, 0);
        }

        System.out.printf("got %d vals\n", m.size());
        long t1 = System.currentTimeMillis();
        double navg_hits = 0.0;
        int NTESTS = 1_0;
        for (int i = 0; i < NTESTS; i++) {

            long v1 = (long) (Math.random() * Long.MAX_VALUE);
            long v2 = (long) (Math.random() * Long.MAX_VALUE);
            long minv = Math.min(v1, v2);
            long maxv = Math.max(v1, v2);

            //int numHits = m.subMap(minv, true, maxv, true).size();
            int numHits = m.subMap(minv, true, maxv, true).keySet().size();
            navg_hits += (numHits / (0.0 + NTESTS));

        }
        long t2 = System.currentTimeMillis();
        System.out.printf("%d tests avg: %.2f took: %.2f sec\n", NTESTS, navg_hits, (t2 - t1) / 1_000.0);

    }

    @Test
    public void testSum() {


        NormalDistribution meanND = new NormalDistribution(13.0, 1.0);
        NormalDistribution sdND = new NormalDistribution(2.0, 1.7);
        Vector<NormalDistribution> normals = mapIndex(154, (_i) -> new NormalDistribution(meanND.sample(), Math.abs(sdND.sample())));
        Vector<ErrorEstimationDistribution> errs = map(normals, (_n) -> ErrorEstimationDistribution.fromNormalDistrib(_n));

        NormalDistribution sumND = new NormalDistribution(NumUtils.sum(normals, (_n) -> _n.getMean()), Math.sqrt(NumUtils.sum(normals, (_n) -> Math.pow(_n.getStandardDeviation(), 2.0))));
        NormalDistribution normed = new NormalDistribution(sumND.getMean() / normals.size(), sumND.getStandardDeviation() / normals.size());

        ErrorEstimationDistribution summedUp = ErrorEstimationDistribution.combineAndAverage(errs);

        double alpha = 0.01;
        double minval = normed.inverseCumulativeProbability(alpha);
        double maxval = normed.inverseCumulativeProbability(1.0 - alpha);
        double step = (maxval - minval) / 20.0;

        double sumDiff = 0.0;
        int ntests = 0;
        for (double v = minval; v < maxval; v += step) {

            double p1 = normed.cumulativeProbability(v);
            double p2 = summedUp.getCumulativeFrequencyToFoldChange(v);
            sumDiff += Math.abs(p1 - p2);
            System.out.printf("%.2f sumnd: %g  summedup: %g\n", v, p1, p2);
            ntests++;
        }
        System.out.printf("avg err: %g\n", sumDiff / ntests);
        double avgErr = sumDiff / ntests;
        assert (avgErr < 0.02);

    }

    static class TestResult{
        int idx;
        Vector<Double> reps1;
        Vector<Double> reps2;

        double npval;
        double epval;

        double nlog10_npval;
        double nlog10_epval;
        double avg;
        double sum;

        TestResult(int idx, Vector<Double> reps1, Vector<Double> reps2, ErrorEstimationDistribution scaled, double totalVariance ) {
            this.idx = idx; this.reps1 = reps1; this.reps2 = reps2;



            for (double d1 : reps1) {
                for (double d2 : reps2) {
                    double diff = d1 - d2;
                    sum += diff;


                }
            }


            avg = sum / (reps1.size() * reps2.size());
            ErrorEstimationDistribution shifted = new ErrorEstimationDistribution(scaled, avg);
            shifted.getSD(avg);
            double N1 = reps1.size() * reps1.size();
            NormalDistribution combined = new NormalDistribution(sum / N1, Math.sqrt(totalVariance) / N1);

            npval = combined.cumulativeProbability(0.0);
            npval = 2.0 * (Math.min(npval, 1.0 - npval));

            epval = shifted.getCumulativeFrequencyToFoldChange(0.0);
            epval = 2.0 * (Math.min(epval, 1.0 - epval));


            nlog10_epval = - NumUtils.logN(epval, 10.0);
            nlog10_npval = - NumUtils.logN(npval, 10.0);

        }

        public static DataTable toTable(Vector<TestResult> results) {
            Function<Vector<Double>, String> vfmt = (_v) -> StringUtils.join(",", _v,  (_d) -> String.format("%.2f", _d));
            return DataTable.buildTable(results,
                    DataTable.buildHeader("idx", (TestResult _r) -> _r.idx)
                    .add("reps1", (_r) -> vfmt.apply(_r.reps1))
                    .add("reps2", (_r) -> vfmt.apply(_r.reps2))
                    .add("avg", (_r) -> _r.avg)
                    .add("npval", (_r) -> _r.npval)
                    .add("epval", (_r) -> _r.epval)

                    .add("nlog10npval", (_r) -> _r.nlog10_npval)
                    .add("nlog10epval", (_r) -> _r.nlog10_epval)
                    .add("diff.nlog10", (_r) -> _r.nlog10_epval - _r.nlog10_npval)

            );
        }
    }

    @Test
    public void testDependency() {
        NormalDistribution nd1 = new NormalDistribution(7.0, 1.7);
        NormalDistribution nd2 = new NormalDistribution(7.0, 0.3);

        ErrorEstimationDistribution err1 = ErrorEstimationDistribution.fromNormalDistrib(nd1);
        ErrorEstimationDistribution err2 = ErrorEstimationDistribution.fromNormalDistrib(nd2);


        int nrep1 = 4;
        int nrep2 = 7;

        double var1 = Math.pow(nd1.getStandardDeviation(), 2.0);
        double var2 = Math.pow(nd2.getStandardDeviation(), 2.0);
        double sumVar = var1 + var2;


        double independentVariance = nrep1 * nrep2 * sumVar;

        double totalVariance = independentVariance;

        totalVariance += nrep1 * nrep2 * (nrep2 - 1) * var1;

        //every rep from the secondside has nrep1 -1 partner
        totalVariance += nrep1 * nrep2 * (nrep1 - 1) * var2;


        Vector<TestResult> results = new Vector<>();


        /**
         *  nrep1 * nrep2 * sumVar * ( 1 +  (nrep2 - 1) * var1)  + ((nrep1 - 1)) * var2))
         *
         *
         */

        double totalSD = Math.sqrt(totalVariance);
        double indepSD = Math.sqrt(independentVariance);

        //double dependentFactor = totalSD / indepSD;

        double scaleFactor = Math.sqrt(1.0 + ((nrep2 - 1) * err1.getVariance()) + ((nrep1 - 1) * err2.getVariance()));


        double renormed = totalSD / (nrep1 * nrep2);

        double current = Math.sqrt(err1.getVariance() + err2.getVariance());

        double dependentFactor = renormed / current;

        System.out.printf("err1: %.2f err2: %.2f current: %.2f renormed: %.2f factor: %.2f\n", err1.getSD(), err2.getSD(), current, renormed, dependentFactor);


        ErrorEstimationDistribution diffDistrib = ErrorEstimationDistribution.substract(err1, err2);
        ErrorEstimationDistribution scaleddiffDistrib = diffDistrib.scale(dependentFactor);



        for (int i = 0; i < 100_000; i++) {
            Vector<Double> v1 = mapIndex(nrep1, (_i) -> nd1.sample());
            Vector<Double> v2 = mapIndex(nrep2, (_i) -> nd2.sample());

            results.add(new TestResult(i, v1, v2, scaleddiffDistrib, totalVariance));
        }

        String err = DistribUtils.checkUniform(map(results, (_r) -> _r.npval), 0.01);

        String erremp = DistribUtils.checkUniform(map(results, (_r) -> _r.epval), 0.02);

        PlotCreator pc = new PlotCreator();
        pc.cumhist("normal", results, (_r) -> _r.npval,  results.size(), false, true);
        pc.cumhist("emp",  results, (_r) -> _r.epval, results.size(), false, true);
        pc.abline("", null, null, 0.0,100.0);
        ImageUtils.showImage("uniform?", pc.getImage(), false);
        pc.scatter("", results, (_r) -> _r.nlog10_npval, (_r) -> _r.nlog10_epval);
        pc.setLabels("normal", "emprical", null);
        ImageUtils.showImage("test", pc.getImage());

        DataTable dt = TestResult.toTable(results);
        PagedDataTable.getDetailViewFrame(dt, null, null);
        assertTrue(err == null, String.format("normal: %s", err));
        assertTrue(erremp == null, String.format("eno: %s", erremp));


    }

    public static void main(String[] args) {
        DistSumTest dst = new DistSumTest();
        dst.testDependency();
    }
}
