package empires.rnaseq;


import lmu.utils.ObjectGetter;
import lmu.utils.Region1D;
import lmu.utils.RegionVector;

import java.util.Collection;

public class GenomicRegionVector extends RegionVector
{
    public String id = null;
    public String symbol = null;
    public String src = null;
    public String type = null;
    public String chr;
    public Boolean strand;
    public String gene = null;

    public GenomicRegionVector(String id, String chr, Boolean strand, Collection<Region1D> regions)
    {
        super();
        this.id = id;
        this.chr =chr;
        this.strand =strand;
        ObjectGetter.apply(regions, (_r) -> addRegion(_r));
    }


    public GenomicRegionVector(String id, String chr, boolean strand, RegionVector rv)
    {
        this(id, chr, strand, rv.getRegions());
    }

    public void setType(String type)
    {
        this.type = type;
    }

    public String getType()
    {
        return type;
    }

    public void setSrc(String src)
    {
        this.src = src;
    }

    public String getSrc()
    {
        return src;
    }

    public void setId(String id)
    {
        this.id = id;
    }

    public String getGene()
    {
        return gene;
    }

    public String getId()
    {
        return id;
    }
    public String getChr()
    {
        return chr;
    }

    public Boolean getStrand()
    {
        return strand;
    }

    public String toString()
    {
        return String.format("%s on %s(%c) %s", id, chr, GenomicUtils.getStrand(strand), super.toString());
    }







}

