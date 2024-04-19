package nlEmpiRe;

import lmu.utils.*;

import java.util.*;

import static lmu.utils.IteratorUtils.rangev;
import static lmu.utils.ObjectGetter.*;
import org.apache.logging.log4j.Logger;

public class Normalization {


    ShiftedGroup merged;
    Vector<Double> shifts;
    int anchorIdx = -1;

    Vector<ShiftedGroup> steps = new Vector<>();

    static class NonClustarableException extends RuntimeException{
        public HashMap<Integer, HashSet<Integer>> clusters;
        public NonClustarableException(HashMap<Integer, HashSet<Integer>> clusters) {
            this.clusters = clusters;
        }
    }


    /**
     *
     * @param replicate2measurementlogvalues
     *      first dimension: num of replicates
     *      second dimension: num of measurements. null-s NaN-s are allowed
     *
     *      assumption: for every measured object (second dimension) every replicate measures the same.
     *          the precision of the measurement is however is influenced by the a signal level dependent error (normal value)
     *
     *          goal is to find a set of shifts (factor sizes) so that
     *          the sum(SD) is minimal
     */
    public Normalization(Vector<Vector<Double>> replicate2measurementlogvalues) {

        Logger log = LogConfig.getLogger();
        int N = replicate2measurementlogvalues.size();

        PairScore[][] errMatrix = getPairWiseErrors(replicate2measurementlogvalues);
        HashMap<Integer, ShiftedGroup> root2group = buildMap(rangev(N), (_i) -> errMatrix[_i][_i].sg1);

        // calc all pairs
        // then take most similar pair
        // merge them -> update the scores to the other ones (priority queue with invalidation)
        PriorityQueue<PairScore> pq = new PriorityQueue<>((_ps1, _ps2) -> Double.compare(_ps1.err.err, _ps2.err.err));
        for(int i=0; i < N - 1; i++) {
            for(int j=i+1; j<N; j++) {
                if(errMatrix[i][j].err == null || errMatrix[i][j].err.shift == null)
                    continue;
                pq.add(errMatrix[i][j]);
            }
        }

        UnionFind uf = new UnionFind(N);

        while(!pq.isEmpty() && uf.getNumRoots() > 1) {
            PairScore ps = pq.poll();

            if(ps.sg1 != root2group.get(uf.root(ps.sg1.idx)) || ps.sg2 != root2group.get(uf.root(ps.sg2.idx)))
                continue; //invalidated


            Integer root = uf.connect(ps.sg1.idx, ps.sg2.idx);
            if(root == null) //already merged
                continue;



            if(ps.err == null)
                throw new FRuntimeException("err is null for: %d %d\n", ps.sg1.idx, ps.sg2.idx);
            ShiftedGroup sg = new ShiftedGroup(root, ps.err);


            steps.add(sg);
            log.debug("next merge: %s", sg);
            root2group.put(root, sg);

            //add scores to the remaining roots
            for(int otherRoot : uf.getRootSet()) {
                if(otherRoot - root == 0)
                    continue;

                PairScore nps = new PairScore(sg, root2group.get(otherRoot));
                if(nps.err == null || nps.err.shift == null)
                    continue;
                pq.add(nps);
            }


        }
        if(uf.getNumRoots() > 1) {
            throw new NonClustarableException(uf.getClusters());
        }
        int finalRoot = uf.getRootSet().iterator().next();
        merged = root2group.get(finalRoot);

        anchorIdx = 0;
        for(Map.Entry<Integer, Double> e : merged.idx2shift.entrySet()) {
            if(Math.abs(merged.idx2shift.get(anchorIdx)) <= Math.abs(e.getValue()))
                continue;

            anchorIdx = e.getKey();
        }
        shifts = mapIndex(N, (_i) -> merged.idx2shift.get(_i));
    }


    public Vector<ShiftedGroup> getSteps() {
        return steps;
    }

    public static Vector<Integer> getBestSubSample(Vector<Vector<Double>> replicate2measurementlogvalues, int num_required_sammples) {
        Logger log = LogConfig.getLogger();
        int N = replicate2measurementlogvalues.size();
        if(num_required_sammples >= replicate2measurementlogvalues.size())
            return rangev(replicate2measurementlogvalues.size());

        PairScore[][] errMatrix = getPairWiseErrors(replicate2measurementlogvalues);
        HashMap<Integer, ShiftedGroup> root2group = buildMap(rangev(N), (_i) -> errMatrix[_i][_i].sg1);

        // calc all pairs
        // then take most similar pair
        // merge them -> update the scores to the other ones (priority queue with invalidation)
        PriorityQueue<PairScore> pq = new PriorityQueue<>((_ps1, _ps2) -> Double.compare(_ps1.err.err, _ps2.err.err));
        for(int i=0; i < N - 1; i++) {
            for(int j=i+1; j<N; j++) {
                pq.add(errMatrix[i][j]);
            }
        }

        log.info("will get best %d subsamples from %d samples", replicate2measurementlogvalues.size(), num_required_sammples);

        UnionFind uf = new UnionFind(N);
        int[] roots = null;
        while(!pq.isEmpty() && uf.getNumRoots() > 1) {
            PairScore ps = pq.poll();

            if(ps.sg1 != root2group.get(uf.root(ps.sg1.idx)) || ps.sg2 != root2group.get(uf.root(ps.sg2.idx)))
                continue; //invalidated

            int r1 = uf.root(ps.sg1.idx);
            int r2 = uf.root(ps.sg2.idx);

            log.info("%d -> %d, %d ->%d", ps.sg1.idx, r1, ps.sg2.idx, r2);
            if(r1 == r2)
                continue;

            int s1 = uf.getClusterSize(r1);
            int s2 = uf.getClusterSize(r2);

            log.info("sizes: %d, %d", s1, s2);
            if( (s1 + s2) >= num_required_sammples) {
                log.info("ready, break");
                roots = new int[]{r1, r2};
                break;
            }
            Integer root = uf.connect(ps.sg1.idx, ps.sg2.idx);
            if(root == null) //already merged
                continue;

            log.debug("merged clusters of sizes: %d (root: %d), %d(root: %d) to new root: %d", s1, r1, s2, r2, root);

            if(ps.err == null)
                throw new FRuntimeException("err is null for: %d %d\n", ps.sg1.idx, ps.sg2.idx);
            ShiftedGroup sg = new ShiftedGroup(root, ps.err);

            //log.debug("next merge: %s err-SD: %.3f", sg, ps.err);
            root2group.put(root, sg);

            //add scores to the remaining roots
            for(int otherRoot : uf.getRootSet()) {
                if(otherRoot - root == 0)
                    continue;

                pq.add(new PairScore(sg, root2group.get(otherRoot)));
            }
        }



        HashMap<Integer, HashSet<Integer>> clusters = uf.getClusters();
        HashSet<Integer> c1 = clusters.get(roots[0]);
        HashSet<Integer> c2 = clusters.get(roots[1]);

        log.info("found subclusters: %d, %d", c1.size(), c2.size());

        HashSet<Integer> bigger = (c1.size() >= c2.size()) ? c1 : c2;
        HashSet<Integer> smaller = (c1 == bigger) ? c2 : c1;

        Vector<Integer> selected = new Vector<>();
        selected.addAll(bigger);
        if(bigger.size() + smaller.size() == num_required_sammples) {
            selected.addAll(smaller);
            return selected;
        }


        pq.clear();
        for(int idx1 : bigger) {
            for(int idx2 : smaller) {
                pq.add(errMatrix[idx1][idx2]);
            }
        }
        //add the ones from the smaller with the smallest distance to the current bigger
        while(selected.size() < num_required_sammples) {


            while(!pq.isEmpty()) {
                PairScore ps = pq.poll();
                if(!bigger.contains(ps.sg1.idx) || !smaller.contains(ps.sg2.idx))
                    continue;

                int smallerIdx = ps.sg2.idx;
                smaller.remove(smallerIdx);
                bigger.add(smallerIdx);
                selected.add(smallerIdx);
                for(int idx2 : smaller) {
                    pq.add(errMatrix[smallerIdx][idx2]);
                }

            }
        }

        return  selected;

    }
    public static PairScore[][] getPairWiseErrors(Vector<Vector<Double>> replicate2measurementlogvalues) {
        int N = replicate2measurementlogvalues.size();

        Vector<ShiftedGroup> groups = map(rangev(N), (_i) -> new ShiftedGroup(_i, replicate2measurementlogvalues.get(_i)));

        PairScore[][] scoreMatrix = new PairScore[N][N];



        for(int i=0; i < N; i++) {
            scoreMatrix[i][i] = new PairScore(groups.get(i), groups.get(i));
            for(int j=i+1; j<N; j++) {
                PairScore ps = new PairScore(groups.get(i), groups.get(j));
                scoreMatrix[i][j] = ps;
                scoreMatrix[j][i] = ps.flip();
            }
        }
        return scoreMatrix;
    }


    public int getAnchorSampleIdx() {
        return anchorIdx;
    }

    public Vector<Double> getShifts() {
        return shifts;
    }
}
