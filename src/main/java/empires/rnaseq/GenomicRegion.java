package nlEmpiRe.rnaseq;

import lmu.utils.*;
import lmu.utils.tuple.*;


public class GenomicRegion extends Region1D
{
    String chr;
    Boolean strand;
    public GenomicRegion(String chr, Boolean strand, int start, int end)
    {
        super(start, end, null);
        this.chr =chr;
        this.strand =strand;
    }

    public GenomicRegion(Tuple4<String, String, Boolean, Region1D> t)
    {
        this(t.get1(), t.get2(), t.get3().getX1(), t.get3().getX2());
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
        return String.format("%s%c:%d-%d", chr, GenomicUtils.getStrand(strand), getX1(), getX2());
    }
}

