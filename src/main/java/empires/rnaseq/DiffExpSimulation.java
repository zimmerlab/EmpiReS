package empires.rnaseq;

import lmu.utils.*;
import lmu.utils.tuple.Tuple3;
import lmu.utils.tuple.Tuple4;
import empires.Normalization;
import empires.NormalizedReplicateSet;
import empires.input.ReplicateSetInfo;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.util.*;
import java.util.function.BiFunction;

import static lmu.utils.IteratorUtils.rangev;
import static lmu.utils.ObjectGetter.*;

public class DiffExpSimulation {

    ReplicateSetInfo input;
    int nrep1;
    int nrep2;
    NormalDistribution foldchangeND;
    double minfc;
    double percentDifferential;

    Vector<UPair<String>> swaps = new Vector<>();

    HashSet<String> upRegulated = new HashSet<>();
    HashSet<String> dnRegulated = new HashSet<>();
    HashMap<String, Double> feature2fc = new HashMap<>();
    HashMap<String, String> feature2testInfo = new HashMap<>();
    Set<String> usedFeatures = new HashSet<>();

    ReplicateSetInfo control;
    ReplicateSetInfo diff;

    static final double MINCOUNT_FILTER = 10.0;

    Set<String> invalids;
    public DiffExpSimulation(ReplicateSetInfo replicates, NormalDistribution foldchangeND, double minfc, double percentDifferential) {
        this(replicates, replicates.getNumReplicates() >> 1, foldchangeND, minfc, percentDifferential);

    }

    public DiffExpSimulation(ReplicateSetInfo replicates, int simul_reps1, NormalDistribution foldchangeND, double minfc, double percentDifferential) {
        this.input = replicates;
        if(replicates.getNumReplicates() - simul_reps1 < 2)
            throw new FRuntimeException("won't/cant simulate replicates %d,%d (input: %d)", replicates.getNumReplicates(), simul_reps1, replicates.getNumReplicates() - simul_reps1);


        Logger log = LogConfig.getLogger();
        this.foldchangeND = foldchangeND;
        this.minfc = minfc;
        this.percentDifferential = percentDifferential;

        Vector<Integer> repidx = shuffle(rangev(replicates.getNumReplicates()), false);
        Vector<Integer> reps1 = VectorUtils.slice(repidx, 0, simul_reps1);
        Vector<Integer> reps2 = VectorUtils.slice(repidx, simul_reps1, repidx.size());



        control = replicates.getReplicateSubSet(reps1);
        ReplicateSetInfo subset = replicates.getReplicateSubSet(reps2);

        BiFunction<ReplicateSetInfo, String, Boolean> stableMeasured = (_rsi, _fn) -> filteredSize(_rsi.getReplicateData(_fn), (_d) -> !Double.isNaN(_d) && _d > 2.0) >= Math.min(3, _rsi.getNumReplicates());

        ReplicateSetInfo C = control;
        ReplicateSetInfo D = subset;
        HashSet<String> measuredFeatures1 = filterToSet(control.getFeatureNames(), (_c) -> stableMeasured.apply(C, _c));
        HashSet<String> measuredFeatures2 = filterToSet(subset.getFeatureNames(), (_c) -> stableMeasured.apply(D, _c));



        for(String f : SetInfo.intersect(measuredFeatures1, measuredFeatures2)) {
            Vector<Double> cdata = map(control.getReplicateData(f), (_d) -> Double.isNaN(_d) ? 0.0 : Math.pow(2.0, _d));
            Vector<Double> ddata = map(subset.getReplicateData(f), (_d) -> Double.isNaN(_d) ? 0.0 : Math.pow(2.0, _d));
            double cmin = NumUtils.min(cdata);
            double dmin = NumUtils.min(ddata);


            double maxmin = Math.max(cmin, dmin);

            usedFeatures.add(f);

            if(cmin < MINCOUNT_FILTER || dmin < MINCOUNT_FILTER)
                continue;



        }

        control = control.restrictToSubFeatures(usedFeatures);
        subset = subset.restrictToSubFeatures(usedFeatures);

        NormalizedReplicateSet nr1 = new NormalizedReplicateSet(control);
        NormalizedReplicateSet nr2 = new NormalizedReplicateSet(subset);

        Vector<String> totalFeatures = toVector(usedFeatures);
        Normalization norm = new Normalization(toVector(nr1.getMedianValues(totalFeatures), nr2.getMedianValues(totalFeatures)));
        Vector<Double> shifts = norm.getShifts();
        double controlShift = shifts.get(0);
        double diffShift = shifts.get(1);

        double SHIFT = controlShift - diffShift;

        double MAXDIFF = minfc * 0.4;

        int numRepDiff = Math.abs(control.getNumReplicates() - subset.getNumReplicates());
        int minNumRep = Math.min(control.getNumReplicates() , subset.getNumReplicates());


        Vector<Tuple4<String, Double, Double, Boolean>> pre_feature2signalInfo = NumUtils.sort(
                map(usedFeatures,
                (_s) ->
                {
                    //calculate num inversions
                    Vector<Pair<Double, Integer>> v = map(nr1.getNormed(_s).nonNanValues, (_d) -> Pair.create(_d + controlShift, -1));
                    v.addAll(map(nr2.getNormed(_s).nonNanValues, (_d) -> Pair.create(_d + diffShift, 1)));
                    NumUtils.sort(v, (_p) -> _p.getFirst());

                    int maxDiff = 0;
                    int sum = 0;
                    for(Pair<Double, Integer> p : v) {
                        sum += p.getSecond();
                        maxDiff = Math.max(maxDiff, Math.abs(sum));
                    }

                    feature2testInfo.put(_s, String.format("%s maxdiff: %d repdiff: %s", v, maxDiff, numRepDiff));
                    //System.out.printf("%s maxdiff: %d repdiff: %s\n", v, maxDiff, numRepDiff);
                    //should be around rep diff

                    double cmedian = NumUtils.getMedianMean(nr1.getNormed(_s).nonNanValues);
                    double dmedian = NumUtils.getMedianMean(nr2.getNormed(_s).nonNanValues);

                    return Tuple4.create(_s,  cmedian, Math.abs(cmedian - dmedian + SHIFT), maxDiff < minNumRep);
                }), (_t) -> _t.get1());



        Vector<Tuple4<String, Double, Double, Boolean>> feature2signalInfo = filter(pre_feature2signalInfo, (_t) ->  _t.get2() < MAXDIFF );

        invalids = mapToSet(filter(pre_feature2signalInfo, (_p) -> !_p.get3() && _p.get2() >= minfc), (_t) -> _t.get0());

        //usedFeatures = mapToSet(feature2signalInfo, (_f) -> _f.get0());
        //control = control.restrictToSubFeatures(usedFeatures);
        //subset = subset.restrictToSubFeatures(usedFeatures);
        boolean[] masked = new boolean[feature2signalInfo.size()];

        int entities_to_simulate = feature2signalInfo.size();

        int num_up_simul = Math.max(1, (int)(0.5 * percentDifferential * entities_to_simulate));


        if(entities_to_simulate < num_up_simul * 2)
            throw new FRuntimeException("something is wrong... should simulate %d entities differential, but got only %d", num_up_simul * 2, entities_to_simulate);



        int up_step = (int)(entities_to_simulate / (0.0 + (num_up_simul)));
        //mark the entitites to simulate first
        Vector<Integer> up_selected = new Vector<>();
        for(int i=0; i< feature2signalInfo.size(); i+=up_step) {
            up_selected.add(i);
            masked[i] = true;
        }

        for(int idx : up_selected) {

            String fn = feature2signalInfo.get(idx).get0();
            double cmedian = NumUtils.getMedianMean(nr1.getNormed(fn).nonNanValues);
            double dmedian = NumUtils.getMedianMean(nr2.getNormed(fn).nonNanValues);

            double fc = 0.0;
            while(fc < minfc) {
                fc = Math.abs(foldchangeND.sample());
            }

            fc += Math.max(0.0, dmedian - cmedian);

            double baseSignal = feature2signalInfo.get(idx).get1();
            Tuple3<String, Double, Double> q = Tuple3.create("", baseSignal + fc, 0.0);
            int hitIdx = Collections.binarySearch(feature2signalInfo, q, (_p1, _p2) -> Double.compare(_p1.get1(), _p2.get1()));
            if(hitIdx < 0) {
                hitIdx = - hitIdx - 1;
            }
            while(hitIdx < feature2signalInfo.size() && masked[hitIdx]) {
                hitIdx++;
            }

            if(hitIdx >= feature2signalInfo.size())
                continue;

            double diffSignal = feature2signalInfo.get(hitIdx).get1();
            String fn1 = feature2signalInfo.get(idx).getFirst();
            String fn2 = feature2signalInfo.get(hitIdx).getFirst();

            feature2fc.put(fn1, diffSignal - baseSignal);
            feature2fc.put(fn2, baseSignal - diffSignal);
            masked[hitIdx] = true;
            swaps.add(UPair.createU(fn1, fn2));

        }



        diff = subset.swapFeatures(swaps);

        apply(swaps, (_s) -> upRegulated.add(_s.getFirst()));
        apply(swaps, (_s) -> dnRegulated.add(_s.getSecond()));



    }

    public Set<String> getUsedFeatures() {
        return usedFeatures;
    }

    public double getSimulatedFC(String feature) {
        return feature2fc.getOrDefault(feature, 0.0);
    }

    public HashSet<String> getSimulatedFeatures(boolean up) {
        return (up) ? upRegulated : dnRegulated;
    }

    public UPair<ReplicateSetInfo> getSimulation() {
        return UPair.createU(control, diff);
    }

    public String getNormTestInfo(String feature) {
        return feature2testInfo.getOrDefault(feature, "");
    }

    public int getNumInvalids() {
        return invalids.size();
    }

    public boolean isSimulated(String feature) {
        return upRegulated.contains(feature) || dnRegulated.contains(feature);
    }
}
