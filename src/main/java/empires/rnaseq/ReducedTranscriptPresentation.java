package nlEmpiRe.rnaseq;

import lmu.utils.*;
import lmu.utils.tuple.Tuple4;
import nlEmpiRe.input.EQClassInput;

import java.util.*;

import static lmu.utils.IteratorUtils.rangev;
import static lmu.utils.ObjectGetter.*;

import org.apache.logging.log4j.Logger;
public class ReducedTranscriptPresentation {

    Logger log = LogConfig.getLogger();
    static final int DEFAULT_MINTRANSCRIPTS_TO_KEEP_IN_RESTRICT = 10;
    HashMap<Tuple, Double> inputEqClassToCounts;

    HashMap<String, Vector<Tuple>> tr2eqs = new HashMap<>();
    HashMap<String, Double> tr2Count;

    /* remove the ones with too few reads to distinct from other transcripts */
    HashMap<String, Double> restrictedTr2Count;
    HashMap<Tuple, Double> restrictedEQCounts;

    Vector<String> leadTranscriptsToTest = new Vector<>();
    HashMap<String, String>  tr2ToTestTr = new HashMap<>();
    HashMap<String, HashSet<String>> tr2cluster = new HashMap<>();
    HashSet<String> removedTranscripts = new HashSet<>();

    Vector<Tuple4<String, String, Double, Double>> removes = new Vector<>();
    Vector<Tuple> orderedEqClasses = new Vector<>();
    HashMap<Tuple, Integer> eqClass2Idx;

    public interface DetailProcesser {
        public void init(ReducedTranscriptPresentation reducer, HashMap<String, Double> tr2Count, HashMap<String, Vector<Tuple>> tr2eqs, HashMap<Tuple, Double> EqClassToCounts);
        public void addIter(int iter, HashMap<String, Double> tr2maxInfo, String trToRemove, HashMap<String, PriorityQueue<Pair<String, Double>>> tr2mostSim);
        public void addRemove(String transcript, HashMap<Tuple, Double> restrictedEQCounts, HashMap<String, Double> restrictedTr2Count);
        public void addRemoves(int iter ,Vector<Tuple4<String, String, Double, Double>> removes);
    }

    static double getInfoScore(double value) {
        //return (value == 0.0) ? 1.0 : NumUtils.logN(value, 10.0);
        return NumUtils.logN(2 + value, 10.0);
    }

    void updateSimScore(HashMap<String, PriorityQueue<Pair<String, Double>>> tr2mostSim, String tr1, String tr2, double score) {
        PriorityQueue<Pair<String, Double>> old = tr2mostSim.get(tr1);
        if(old == null) {

            tr2mostSim.put(tr1, old = new PriorityQueue<>((_p1, _p2) -> Double.compare(_p1.getSecond(), _p2.getSecond())));
        }

        old.add(Pair.create(tr2, score));

    }

    public ReducedTranscriptPresentation(HashMap<Tuple, Double> eqClasses) {
        this(eqClasses, DEFAULT_MINTRANSCRIPTS_TO_KEEP_IN_RESTRICT);
    }
    public ReducedTranscriptPresentation(HashMap<Tuple, Double> EqClassToCounts, int minTranscriptsToClean) {
        this(EqClassToCounts, minTranscriptsToClean, null);
    }

    public ReducedTranscriptPresentation(HashMap<Tuple, Double> EqClassToCounts, int minTranscriptsToClean, DetailProcesser detailProcesser) {
        this.inputEqClassToCounts = EqClassToCounts;

        tr2eqs = new HashMap<>();
        tr2Count = getPartialCounts(EqClassToCounts.keySet(), EqClassToCounts, tr2eqs);

        if(detailProcesser != null) {
            detailProcesser.init(this, tr2Count, tr2eqs, EqClassToCounts);
        }
        removedTranscripts.addAll(tr2Count.keySet());

        restrictedTr2Count = tr2Count;
        restrictedEQCounts = EqClassToCounts;

        HashMap<String, Double> tr2maxInfo = new HashMap<>();
        HashMap<String, PriorityQueue<Pair<String, Double>>> tr2mostSim = new HashMap<>();

        Vector<String> trnames = toSortedVector(tr2eqs.keySet(), true);
        HashMap<String, Integer> tr2idx = buildIndexMap(trnames);

        UnionFind uf = new UnionFind(trnames.size());

        HashMap<String, Vector<Tuple>> restrictedTR2eqs = new HashMap<>();
        restrictedTR2eqs.putAll(tr2eqs);


        int iter = 0;
        HashMap<UPair<String>, UPair<Double>> infos = new HashMap<>();
        while(uf.getNumRoots() > minTranscriptsToClean) {

            iter++;
            tr2maxInfo.clear();
            tr2mostSim.clear();
            for(Map.Entry<String, Double> e : restrictedTr2Count.entrySet()) {
                String tr1 = e.getKey();
                HashMap<String, Double> conditional = getPartialCounts(restrictedTR2eqs.get(tr1), restrictedEQCounts, null);
                for(String tr2: conditional.keySet()) {

                    if(tr1.equals(tr2))
                        continue;

                    double tr1uniq = e.getValue() - conditional.get(tr2);
                    double tr2uniq = restrictedTr2Count.get(tr2) - conditional.get(tr2);

                    double simScore = getInfoScore(tr1uniq) * getInfoScore(tr2uniq);

                    updateSimScore(tr2mostSim, tr1, tr2, simScore);
                    updateSimScore(tr2mostSim, tr2, tr1, simScore);

                    infos.put(UPair.createUSorted(tr1, tr2), UPair.createU(tr1uniq, tr2uniq));

                    tr2maxInfo.put(tr1, Math.max(tr1uniq, tr2maxInfo.getOrDefault(tr1, 0.0)));
                    tr2maxInfo.put(tr2, Math.max(tr2uniq, tr2maxInfo.getOrDefault(tr2, 0.0)));
                }
            }


            //remove transcript with the minimal partial score
            if(tr2maxInfo.size() == 0) {
                log.warn("empty tr2max at iter: %d roots: %d", iter, uf.getNumRoots());
                break;
            }
            String transcriptToRemove = NumUtils.minObj(tr2maxInfo.entrySet(), (_e) -> _e.getValue()).getSecond().getKey();

            if(detailProcesser != null) {
                detailProcesser.addIter(iter, tr2maxInfo, transcriptToRemove, tr2mostSim);
            }


            HashSet<String> preTranscripts = new HashSet<>();

            preTranscripts.addAll(restrictedTr2Count.keySet());


            restrictedEQCounts = delete(restrictedEQCounts, toSet(transcriptToRemove));
            restrictedTR2eqs.clear();
            restrictedTr2Count = getPartialCounts(restrictedEQCounts.keySet(), restrictedEQCounts, restrictedTR2eqs);

            Vector<Tuple4<String, String, Double, Double>> iterRemoves = new Vector<>();

            for(String preTr : preTranscripts) {
                if(restrictedTr2Count.containsKey(preTr))
                    continue;

                //connect to the best existing
                PriorityQueue<Pair<String, Double>> bestSims = tr2mostSim.get(preTr);


                while(bestSims.size() != 0) {
                    Pair<String, Double> p = bestSims.poll();
                    if(!restrictedTr2Count.containsKey(p.getFirst()))
                        continue;

                    UPair<Double> c = infos.get(UPair.createUSorted(preTr, p.getFirst()));
                    Tuple4<String, String, Double, Double> remove = Tuple4.create(preTr, p.getFirst(), c.getFirst(), c.getSecond());
                    iterRemoves.add(remove);
                    removes.add(remove);

                    uf.connect(tr2idx.get(preTr), tr2idx.get(p.getFirst()));
                    break;
                }
            }

            if(detailProcesser != null) {
                detailProcesser.addRemoves(iter, iterRemoves);
                detailProcesser.addRemove(transcriptToRemove, restrictedEQCounts, restrictedTr2Count);
            }

        }

        //System.out.printf("got %d clusters (target: %d) roots: %d\n", uf.getClusters().size(),  minTranscriptsToClean, uf.getNumRoots());

        for(Map.Entry<Integer, HashSet<Integer>> e : uf.getClusters().entrySet()) {
            HashSet<String> cluster = mapToSet(e.getValue(), (_idx) -> trnames.get(_idx));
            String parent = trnames.get(e.getKey());
            leadTranscriptsToTest.add(parent);
            apply(cluster, (_s) -> tr2cluster.put(_s, cluster));
            for(String s : cluster) {
                tr2ToTestTr.put(s, parent);
            }
        }

        removedTranscripts.removeAll(leadTranscriptsToTest);
        orderedEqClasses.addAll(restrictedEQCounts.keySet());

        eqClass2Idx = buildIndexMap(orderedEqClasses);

    }


    public String getLeadTr(String tr) {
        return tr2ToTestTr.get(tr);
    }
    public int getEqClassIdx(Tuple t) {
        return eqClass2Idx.get(t);
    }


    public Tuple getEqClassByIdx(int idx) {
        return orderedEqClasses.get(idx);
    }

    
    public HashSet<String> getTrCluster(String tr) {
        return tr2cluster.get(tr);
    }

    public Vector<String> getTranscriptsToTest() {
        return leadTranscriptsToTest;
    }

    public HashSet<String> getRemovedTranscripts() {
        return removedTranscripts;
    }


    public Tuple reduceTuple(Tuple inTuple) {
        return Tuple.tupleFromCollection(toSortedVector(mapToSet(rangev(inTuple.cardinality()), (_i) -> getLeadTr(inTuple.getAsString(_i))), true));
    }
    public HashMap<Tuple, Double> reduce(HashMap<Tuple, Double> eqClasses) {
        return delete(eqClasses, removedTranscripts);
    }

    public Vector<Tuple4<String, String, Double, Double>> getRemoves() {
        return removes;
    }

    public boolean isTestable(String tr1, String tr2) {
        return tr2cluster.get(tr1) != tr2cluster.get(tr2);
    }


    public UPair<Double> getPartialInfos(String tr1, String tr2) {
        double countTr1 = tr2Count.get(tr1);
        HashMap<String, Double> conditional = getPartialCounts(tr2eqs.get(tr1), inputEqClassToCounts, null);

        double tr1uniq = countTr1 - conditional.getOrDefault(tr2, 0.0);
        double tr2uniq = tr2Count.get(tr2) - conditional.getOrDefault(tr2, 0.0);

        return UPair.createU(tr1uniq, tr2uniq);
    }

    public static HashMap<String, Double> getPartialCounts(Collection<Tuple> tuples, HashMap<Tuple, Double> counts, HashMap<String, Vector<Tuple>> tr2eqs) {
        HashMap<String, Double> tr2count = new HashMap<>();

        for(Tuple t : tuples) {
            double c = counts.get(t);

            for(int ti=0; ti<t.cardinality(); ti++) {
                MapBuilder.update(tr2count, t.getAsString(ti), c);
                if(tr2eqs != null) {
                    MapBuilder.updateV(tr2eqs, t.getAsString(ti), t);
                }
            }
        }
        return tr2count;

    }

    public Collection<Tuple> getRestrictedEQClasses() {
        return orderedEqClasses;
    }

    public int getNumRestrictedEQClasses() {
        return orderedEqClasses.size();
    }

    public int getNumClusters() {
        return leadTranscriptsToTest.size();
    }

    public static HashMap<Tuple, Double> delete(HashMap<Tuple, Double> counts, HashSet<String> tr2remove) {
        HashMap<Tuple, Double> rv = new HashMap<>();
        for(Map.Entry<Tuple, Double>  e : counts.entrySet()) {
            Tuple t = e.getKey();
            Vector<String> newKey = filter(mapIndex(t.cardinality(), (_idx) -> t.getAsString(_idx)), (_trId) -> !tr2remove.contains(_trId));
            if(newKey.size() == 0)
                continue; //lost them...
            Tuple reduced = Tuple.tupleFromCollection(newKey);
            MapBuilder.update(rv, reduced, e.getValue());
        }
        return rv;

    }

    public static HashMap<Tuple, Integer> reduce(HashMap<Tuple, Integer> counts, HashMap<String, String> tr2leadTr) {
        HashMap<Tuple, Integer> rv = new HashMap<>();
        for(Map.Entry<Tuple, Integer>  e : counts.entrySet()) {
            Tuple t = e.getKey();
            Vector<String> newKey = toSortedVector(toSet(mapIndex(t.cardinality(), (_idx) -> tr2leadTr.get(t.getAsString(_idx)))), true);
            if(newKey.size() == 0)
                continue; //lost them...
            Tuple reduced = Tuple.tupleFromCollection(newKey);
            MapBuilder.update(rv, reduced, e.getValue());
        }
        return rv;
    }


}
