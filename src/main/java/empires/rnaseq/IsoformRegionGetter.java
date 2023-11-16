package empires.rnaseq;

import lmu.utils.*;
import lmu.utils.tuple.Tuple3;

import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Vector;

public abstract class IsoformRegionGetter
{
    Logger log = LogConfig.getLogger();
    public KeyedIntervalForest<Tuple3<String, Boolean, Region1D<String>>, String> region_tree = null;

    HashMap<String, UPair<RegionSearcher>> chr2searcher = new HashMap<>();

    private HashMap<String, Integer> chr2end = new HashMap<>();

    static class RegionSearcher {
        boolean strand;
        public Vector<Region1D<String>> start_sorted = new Vector<>();
        public Vector<Region1D<String>> end_sorted = new Vector<>();

        public RegionSearcher(boolean strand, Vector<Region1D<String>> infos) {
            this.strand = strand;
            if (infos == null)
                return;

            this.start_sorted.addAll(infos);
            this.end_sorted.addAll(infos);
            NumUtils.sort(start_sorted, (gli) -> gli.getX1());
            NumUtils.sort(end_sorted, (gli) -> gli.getX2());
        }

        public Pair<Integer, Region1D<String>> getClosestRegion(int pos, boolean end_based) {
            int mindist = Integer.MAX_VALUE;
            Region1D<String> rv = null;

            Vector<Region1D<String>> v = (!end_based) ? start_sorted : end_sorted;
            lmu.utils.VectorUtils.Comparator2<Region1D<String>, Integer> comparer = (!end_based) ? (_gli, _p) -> Integer.compare(_gli.getX1(), _p) : (_gli, _p) -> Integer.compare(_gli.getX2(), _p);

            Pair<Boolean, Integer> r = VectorUtils.binarySearch(v, pos, comparer);


            for (int j = 0; j < 2; j++) {
                if ((j == 0 && r.getSecond() == 0) || (j == 1 && r.getSecond() >= v.size()))
                    continue;

                int idx = r.getSecond() - ((j == 0) ? 1 : 0);
                Region1D<String> gli = v.get(idx);

                if (gli.getX1() <= pos && gli.getX2() >= pos)
                    return Pair.create(0, gli);

                int dist = Math.min(Math.abs(gli.getX1() - pos), Math.abs(gli.getX2() - pos));
                if (dist >= mindist)
                    continue;

                mindist = dist;
                rv = gli;
            }

            return Pair.create(mindist,rv);

        }

        public Pair<Integer, Region1D<String>> getClosestRegion(int pos)
        {
            Pair<Integer, Region1D<String>> pre_start = getClosestRegion(pos, false);
            Pair<Integer, Region1D<String>> pre_end = getClosestRegion(pos, true);

            return (pre_start.getFirst() < pre_end.getFirst()) ? pre_start : pre_end;
            /*
            int mindist = Integer.MAX_VALUE;
            Region1D<String> rv = null;

            for (int i=0; i<2; i++)
            {
                Vector<Region1D<String>> v = (i==0) ? start_sorted : end_sorted;
                lmu.utils.VectorUtils.Comparator2<Region1D<String>, Integer> comparer = (i==0) ? (_gli, _p) -> Integer.compare(_gli.getX1(), _p) : (_gli, _p) -> Integer.compare(_gli.getX2(), _p);

                Pair<Boolean, Integer> r = VectorUtils.binarySearch(v, pos, comparer);

                for (int j = 0; j < 2; j++)
                {
                    if ((j == 0 && r.getSecond() == 0) || (j == 1 && r.getSecond() >= v.size()))
                        continue;

                    int idx = r.getSecond() - ((j == 0) ? 1 : 0);
                    Region1D<String> gli = v.get(idx);

                    if (gli.getX1() <= pos && gli.getX2() >=pos)
                        return Pair.create(0, gli);

                    int dist = Math.min(Math.abs(gli.getX1() - pos), Math.abs(gli.getX2() - pos));
                    if (dist >= mindist)
                        continue;

                    mindist = dist;
                    rv = gli;
                }
            }
            return Pair.create(mindist, rv);*/
        }
    }

    public Integer getChrLength(String chr)
    {
        init();
        return chr2end.get(chr);
    }
    synchronized void init()
    {
        if (region_tree != null)
            return;


        region_tree = new KeyedIntervalForest<>((t) -> t.get0(), (t) -> t.get2().getX1(), (t) -> t.get2().getX2());

        HashMap<String, Vector<Integer>> chr2ends = new HashMap<>();

        HashMap<String, Vector<Region1D<String>>> pos2regs = new HashMap<>();
        HashMap<String, Vector<Region1D<String>>> neg2regs = new HashMap<>();


        int nregs = 0;
        for (Tuple3<String, Boolean, Region1D<String>> t : getRegionInfos())
        {
            nregs++;
            region_tree.add(t);
            if(t.get1() == null)
                throw new FRuntimeException("null for strand: %s", t);
            HashMap<String, Vector<Region1D<String>>> target = (t.get1()) ? pos2regs : neg2regs;
            MapBuilder.updateV(target, t.get0(), t.get2());
            MapBuilder.updateV(chr2ends, t.get0(), t.get2().getX2());
        }
        ObjectGetter.apply(chr2ends.entrySet(), (e) -> chr2end.put(e.getKey(), NumUtils.max(e.getValue()) + 1000));
        ObjectGetter.apply(pos2regs.values(), (c) -> Collections.sort(c));
        ObjectGetter.apply(chr2ends.keySet(), (c) -> chr2searcher.put(c, UPair.createU(new RegionSearcher(true, pos2regs.get(c)), new RegionSearcher(false, neg2regs.get(c)))));

        log.info("inited with %d regions chr2ends: %s", nregs, chr2end);
    }

    /** chr, strand, region<id> */
    public abstract Vector<Tuple3<String, Boolean, Region1D<String>>> getRegionInfos();

    public abstract empires.rnaseq.MultiIsoformRegion getRegionById(String id);

    static Vector<Tuple3<String, Boolean, Region1D<String>>> EMPTY  = new Vector<>();

    public Vector<Tuple3<String, Boolean, Region1D<String>>> getRegionInfos(String chr, Integer start , Integer end)
    {
        init();
        if (chr == null)
            return getRegionInfos();


        Integer chrend = null;
        if (null == (chrend = chr2end.get(chr)))
        {
            if (chr2end.get("chr"+chr) != null)
            {
                chr = "chr"+chr;
                chrend = chr2end.get(chr);
            }
            else
            {
                return EMPTY;
            }
        }


        start = (start != null) ? start : 0;
        end = (end != null) ? end : chrend;


        return region_tree.getIntersecting(chr, start, end);
    }

    public Iterator<empires.rnaseq.MultiIsoformRegion> getRegionsAround(String chr, Integer start , Integer end)
    {
        return ObjectGetter.getTransformedIterator(region_tree.getAround(chr, start, end).iterator(), (t) -> getRegionById(t.get2().getObject()));
    }


    public Iterator<empires.rnaseq.MultiIsoformRegion> getRegions(String chr, Integer start , Integer end)
    {
        return ObjectGetter.getTransformedIterator(getRegionInfos(chr, start, end).iterator(), (t) -> getRegionById(t.get2().getObject()));
    }

    public void clearCache()
    {

    }

    public Pair<Integer, Region1D<String>> getPreRegion(Boolean strand, String chr, int pos)
    {
        init();
        UPair<RegionSearcher> searcher = chr2searcher.get(chr);

        if (searcher == null)
            return null;

        return searcher.get(strand).getClosestRegion(pos, false);
    }

    public Pair<Integer, Region1D<String>> getClosestRegion(Boolean strand, String chr, int pos)
    {
        init();
        UPair<RegionSearcher> searcher = chr2searcher.get(chr);

        if (searcher == null)
            return null;

        if (strand != null)
        {
            RegionSearcher rs = (strand) ? searcher.getFirst() : searcher.getSecond();
            return rs.getClosestRegion(pos);
        }
        Pair<Integer, Region1D<String>> r1 = searcher.getFirst().getClosestRegion(pos);
        Pair<Integer, Region1D<String>> r2 = searcher.getSecond().getClosestRegion(pos);

        return (r1.getFirst() < r2.getFirst()) ? r1 : r2;
    }

    public Pair<Integer, Region1D<String>> getClosestRegion(String chr, Boolean strand,  int x1, int x2)
    {
        init();
        Pair<Integer, Region1D<String>> info1 = getClosestRegion(strand, chr, x1);
        Pair<Integer, Region1D<String>> info2 = getClosestRegion(strand, chr, x2);
        if (info1 == null || info2 == null)
        {
            log.warn("requested %s%c %d,%d info1 null: %s info2 null: %s", chr, GenomicUtils.getStrand(strand), x1, x2, info1 == null, info2 == null);

            return null;
        }
        return (info1.getFirst() > info2.getFirst()) ? info2 : info1;
    }



    public Integer getClosestRegionDistance(String chr, Boolean strand, RegionVector merged)
    {
        Pair<Integer, Region1D<String>> p = getClosestRegion(chr, strand, merged.getX1(), merged.getX2());

        return (p == null) ? null : p.getFirst();
    }

    public Iterable<empires.rnaseq.MultiIsoformRegion> getRegionsIteratble() {
        return () -> getRegions();
    }
    public abstract Iterator<empires.rnaseq.MultiIsoformRegion> getRegions();

    public abstract Pair<empires.rnaseq.MultiIsoformRegion, RegionVector> getIsoform(String id);
    public abstract Pair<MultiIsoformRegion, RegionVector> getCoding(String id);


}


