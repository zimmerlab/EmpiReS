package empires.rnaseq;


import lmu.utils.*;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Vector;

import static lmu.utils.ObjectGetter.filter;

public class MultiIsoformRegion//<A>
{
    //A obj;

    public IsoformRegionGetter factory = null;
    public String src = null;
    public String src_version = null;

    public String id;
    public String name;
    public String chr;
    public String biotype;
    public Boolean strand;

    public int start = -1;
    public int end = -1;

    public HashMap<String, RegionVector> isoforms = new HashMap<>();

    public HashSet<String> erroneous_isoforms = new HashSet<>();
    public HashMap<String, String> iso2biotype = new HashMap<>();
    public HashMap<String, String> iso2coding = new HashMap<>();
    public HashMap<String, String> coding2iso = new HashMap<>();
    public HashMap<String, RegionVector> coding = new HashMap<>();

    HashSet<Region1D> exons = new HashSet<>();

    public MultiIsoformRegion()
    {

    }


    public String getDisplayName()
    {
        return (name != null) ? name : id;
    }

    public String getLocationInfo()
    {
        return String.format("%s (%s):%s %s%c:%d-%d", name, id, biotype, chr, GenomicUtils.getStrand(strand), start, end);
    }

    public int getLength()
    {
        return end - start;
    }

    public GenomicRegionVector getTrGRV(String id)
    {
        GenomicRegionVector grv = new GenomicRegionVector(id, chr, strand, isoforms.get(id));
        grv.gene = this.id;
        grv.type = biotype;
        return grv;
    }

    public GenomicRegionVector getCodingGRV(String id)
    {
        return new GenomicRegionVector(this.id + ":" + id, chr, strand, coding.get(id));
    }

    public String getIsoBioType(String id)
    {
        String bt = iso2biotype.get(id);
        return (bt == null) ? "unknown" : bt;
    }


    public UPair<RegionVector> getUTR(String trid)
    {
        RegionVector tr = isoforms.get(trid);

        if(tr == null)
            throw new FRuntimeException("unknown tr id: %s requested for %s", trid, id);

        String codingid = iso2coding.get(tr);

        if(codingid == null)
            throw new FRuntimeException("tr id: %s has no coding part in %s", trid, id);

        RegionVector coding_rv = coding.get(codingid);



        RegionVector utr = tr.substract(coding_rv);
        Vector<Region1D> pre = (strand) ? filter(utr.getRegions(), (_r) -> _r.getX1() < coding_rv.getX1()) : filter(utr.getRegions(), (_r) -> _r.getX2() > coding_rv.getX2());
        Vector<Region1D> post = (strand) ? filter(utr.getRegions(), (_r) -> _r.getX2() > coding_rv.getX2()) : filter(utr.getRegions(), (_r) -> _r.getX1() < coding_rv.getX1());

        return UPair.createU(new RegionVector(pre, true), new RegionVector(post, true));
    }


    public void addRegionVector(String id, String biotype, RegionVector rv)
    {
        start = (start < 0 ) ? rv.getX1() : Math.min(start, rv.getX1());
        end = (end < 0 ) ? rv.getX2() : Math.max(end, rv.getX2());
        rv.setObject(id);
        isoforms.put(id, rv);
        iso2biotype.put(id, biotype);
        ObjectGetter.apply(rv.getRegions(), (_r) -> exons.add(_r));
    }
    RegionVector merged_tr;

    public synchronized RegionVector getMergedTranscript()
    {
        if (merged_tr != null)
            return merged_tr;

        return merged_tr = RegionVector.merge(isoforms.values());
    }

    HashMap<UPair<Integer>, HashSet<String>> introns = null;

    public synchronized HashMap<UPair<Integer>, HashSet<String>> getIntrons()
    {
        if (introns != null)
            return introns;

        HashMap<UPair<Integer>, HashSet<String>> rv = new HashMap<>();

        for (Map.Entry<String, RegionVector> e : isoforms.entrySet())
        {
            for(Region1D r : e.getValue().getInverseRegionVector().getRegions())
            {
                MapBuilder.update(rv, UPair.createU(r.getX1(), r.getX2()), e.getKey());
            }
        }
        return introns = rv;

    }

    public String toString()
    {
        return String.format("%s:%s (%s) %s(%c) %d-%d (%d trs)", id, biotype, name, chr, GenomicUtils.getStrand(strand), start, end, isoforms.size());
    }


    HashMap<RegionVector, TranscriptMapper> id2mapper = null;

    public RegionVector getGeneRegionVector(RegionVector rv, int start, int end)
    {
        return new RegionVector(getGeneRegions(rv, start, end), true);
    }

    public synchronized  Vector<Region1D> getGeneRegions(RegionVector rv, int start, int end)
    {
        if(id2mapper == null)
        {
            id2mapper = new HashMap<>();
        }
        TranscriptMapper tm = id2mapper.get(rv);
        if(tm == null)
        {

            id2mapper.put(rv, tm = new TranscriptMapper(rv, this.start, this.end, strand));
        }
        return tm.getGeneRegions(start, end);
    }

    static class TranscriptMapper
    {
        RegionVector rv;
        boolean strand;
        int gene_start;
        int gene_end;

        Vector<Pair<Integer, Region1D>> start2exon = null;

        TranscriptMapper(RegionVector rv, int gene_start, int gene_end, boolean strand)
        {
            this.rv = rv;
            this.strand = strand;
            this.gene_start = gene_start;
            this.gene_end = gene_end;
        }

        synchronized Vector<Pair<Integer, Region1D>> getStart2Exon()
        {
            if (start2exon != null)
                return start2exon;

            start2exon = new Vector<>();
            int start = 0;
            int x = (strand) ? 0 : rv.getNumRegions() - 1;
            int step = (strand) ? 1 : -1;

            for(int i=x; i>=0 && i<rv.getNumRegions(); i+= step)
            {
                Region1D src = rv.getRegion(i);
                Region1D r = (strand) ? src.translate(-gene_start) : new Region1D(gene_end- src.getX2(), gene_end - src.getX1());
                start2exon.add(Pair.create(start, r));
                start += rv.getRegion(i).getLength();
            }
            //System.out.printf("strand: %s :: %s -> %s\n", strand, rv, start2exon);
            return start2exon;
        }


        Pair<Integer,Integer> getExonIdxAndGenePosition(int pos, boolean end)
        {
            Pair<Boolean, Integer> p = VectorUtils.binarySearch(getStart2Exon(), pos, (_p, _i) -> Integer.compare(_p.getFirst(), _i));
            int exidx = (p.getFirst() && !end) ? p.getSecond() : p.getSecond() -1;
            if (exidx < 0 )
                throw new FRuntimeException("requested pos %d in transcript: %s, start2exon: %s binary search: %s\n",
                        pos, rv, getStart2Exon(), p);

            Region1D ex = start2exon.get(exidx).getSecond();
            int rel_pos = pos - start2exon.get(exidx).getFirst();


            //System.out.printf("pos: %d exons: %s p: %s  ex: %s REL POS: %d\n",pos,  getStart2Exon(), p, ex, rel_pos);
            return Pair.create(exidx, ex.getX1() + rel_pos);

        }

        public Vector<Region1D> getGeneRegions(int start, int end)
        {
            Pair<Integer, Integer> ex1 = getExonIdxAndGenePosition(start, false);
            Pair<Integer, Integer> ex2 = getExonIdxAndGenePosition(end, true);

            if (ex2.getSecond() == start2exon.get(ex2.getFirst()).getSecond().getX1())
            {
                ex2 = Pair.create(ex2.getFirst() -1, rv.getRegion(ex2.getFirst()-1).getX2());
            }

            if (ex1.getFirst() - ex2.getFirst() == 0)
                return ObjectGetter.toVector(new Region1D(ex1.getSecond() , ex2.getSecond(), null));

            Vector<Region1D> rv = ObjectGetter.toVector(new Region1D(ex1.getSecond(), start2exon.get(ex1.getFirst()).getSecond().getX2()));

            for (int i = ex1.getFirst() + 1; i < ex2.getFirst(); i++)
            {
                rv.add(start2exon.get(i).getSecond());

            }
            rv.add(new Region1D(start2exon.get(ex2.getFirst()).getSecond().getX1(), ex2.getSecond(), null));

            return rv;

        }
    }
}

