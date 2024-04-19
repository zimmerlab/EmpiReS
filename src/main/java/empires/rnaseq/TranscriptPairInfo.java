package nlEmpiRe.rnaseq;

import lmu.utils.Pair;
import lmu.utils.Tuple;
import lmu.utils.tuple.Tuple3;

import java.util.Collection;
import java.util.HashSet;
import java.util.Vector;

import static lmu.utils.ObjectGetter.*;

public class TranscriptPairInfo {
    public String tr1;
    public String tr2;

    public Vector<Tuple> eqClassesInBoth = new Vector<>();
    public Vector<Tuple> exclTr1 = new Vector<>();
    public Vector<Tuple> exclTr2 = new Vector<>();


    public Vector<Tuple3<String, Vector<Tuple>, Vector<Tuple>>> testableCombis = new Vector<>();


    protected TranscriptPairInfo(String tr1, String tr2) {
        this.tr1 = tr1; this.tr2 = tr2;
    }

    int calcTestable() {
        Vector<String> names = toVector("merged."+tr1+"."+tr2, "excl."+tr1, "excl."+tr2);
        Vector<Vector<Tuple>> eqs = toVector(eqClassesInBoth, exclTr1, exclTr2);

        for(int i=0; i<eqs.size(); i++) {
            for(int j=i+1; j<eqs.size(); j++) {
                if(eqs.get(i).size() == 0 || eqs.get(j).size() == 0)
                    continue;

                testableCombis.add(Tuple3.create(names.get(i)+"_VS_"+names.get(j), eqs.get(i), eqs.get(j)));
            }
        }
        return testableCombis.size();
    }


    public static Vector<TranscriptPairInfo> getTrPairInfos(Vector<Pair<Tuple, String>> usedFeatures) {

        Vector<Tuple> tuples = map(usedFeatures, (_p) -> _p.getFirst());
        Vector<HashSet<String>> lookups = map(tuples, (_t) -> toSet(mapIndex(_t.cardinality(), (_i) -> _t.getAsString(_i))));
        HashSet<String> trSet = new HashSet<>();
        apply(lookups, (_h) -> trSet.addAll(_h));

        Vector<String> trVec = toSortedVector(trSet, true);
        Vector<TranscriptPairInfo> pairs = new Vector<>();

        for(int i=0; i<trVec.size(); i++) {
            for(int j=i+1; j<trVec.size(); j++) {
                TranscriptPairInfo tp = new TranscriptPairInfo(trVec.get(i), trVec.get(j));

                for(int ti = 0; ti < tuples.size(); ti++) {
                    boolean t1 = lookups.get(ti).contains(tp.tr1);
                    boolean t2 = lookups.get(ti).contains(tp.tr2);

                    if(!t1 && !t2)
                        continue;

                    Tuple eqClass = tuples.get(ti);
                    if(t1 && t2) {
                        tp.eqClassesInBoth.add(eqClass);
                        continue;
                    }
                    if(t1) {
                        tp.exclTr1.add(eqClass);
                    } else {
                        tp.exclTr2.add(eqClass);
                    }
                }

                if(tp.calcTestable() == 0)
                    continue;

                pairs.add(tp);
            }
        }

        return pairs;
    }
}
