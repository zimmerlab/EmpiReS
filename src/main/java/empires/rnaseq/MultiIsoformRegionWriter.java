package empires.rnaseq;

import lmu.utils.NumUtils;
import lmu.utils.ObjectGetter;
import lmu.utils.RegionVector;
import lmu.utils.StringUtils;

import java.io.PrintWriter;
import java.util.Map;
import java.util.Vector;

public class MultiIsoformRegionWriter
{
    static String DEFAULT_SRC =  "ENSEMBL";
    PrintWriter pw;
    final public static String PROTEIN_RV_TYPE = "*CODING*";

    boolean write_merged = false;
    public MultiIsoformRegionWriter(PrintWriter pw, boolean merged)
    {
        this.pw = pw;
        write_merged = merged;
        pw.printf("regionid\tname\tregiontype\tsrc\tchr\tstrand\tstart\tend\trvid\trvtype\tregvec\n");
    }

    public void write(empires.rnaseq.MultiIsoformRegion g, String src)
    {
        String prefix = String.format("%s\t%s\t%s\t%s\t%s\t%c\t%d\t%d", g.id, g.name, g.biotype, src, g.chr, empires.rnaseq.GenomicUtils.getStrand(g.strand), g.start, g.end);
        if (write_merged)
        {
            pw.printf("%s\t%s\t%s\t%s\n", prefix, g.id, g.biotype, g.getMergedTranscript().getSimpleRepresentation());
            return;
        }
        for (Map.Entry<String, RegionVector> e : g.isoforms.entrySet())
        {
            pw.printf("%s\t%s\t%s\t%s\n", prefix, e.getKey(), g.getIsoBioType(e.getKey()), e.getValue().getSimpleRepresentation());
        }
        for(Map.Entry<String, RegionVector> e : g.coding.entrySet())
        {
            pw.printf("%s\t%s\t%s\t%s\n", prefix, e.getKey(), PROTEIN_RV_TYPE + "::" + e.getKey()+ "::" + g.coding2iso.get(e.getKey()), e.getValue().getSimpleRepresentation());
        }
    }

    public void write(Vector<empires.rnaseq.GenomicRegionVector> grv)
    {
        empires.rnaseq.GenomicRegionVector g = grv.get(0);
        int start = NumUtils.min(grv, (rv) -> rv.getX1());
        int end = NumUtils.max(grv, (rv) -> rv.getX2());


        String prefix = String.format("%s\t%s\t%s\t%s\t%s\t%c\t%d\t%d", g.getGene(), g.symbol, g.getType(), g.getSrc(), g.getChr(), empires.rnaseq.GenomicUtils.getStrand(g.getStrand()), start, end);

        if (write_merged)
        {
            RegionVector rv = RegionVector.merge(ObjectGetter.convert(grv));
            pw.printf("%s\t%s\t%s\t%s\n", prefix, g.getGene(), g.getType(), rv.getSimpleRepresentation());
            return;
        }

        for (GenomicRegionVector t: grv)
        {
            pw.printf("%s\t%s\t%s\t%s\n", prefix, t.getId(), t.getType(), t.getSimpleRepresentation());
        }
    }

    public static void writeBED(MultiIsoformRegion mir, PrintWriter pw)
    {
        for(Map.Entry<String, RegionVector> e: mir.isoforms.entrySet())
        {
            String sizes = StringUtils.joinObjects(",", e.getValue().getRegions(), (_r) -> ""+ _r.getLength());
            String starts = StringUtils.joinObjects(",", e.getValue().getRegions(), (_r) -> ""+ _r.getX1());

            pw.printf("chr%s\t%d\t%d\t%s\t0\t%c\t%d\t%d\t0,0,0\t%d\t%s\t%s\n",
                    mir.chr, e.getValue().getX1(), e.getValue().getX2(),
                    mir.id+"."+e.getKey(),
                    GenomicUtils.getStrand(mir.strand),
                    e.getValue().getX1(), e.getValue().getX2(),
                    e.getValue().getNumRegions(),
                    sizes,
                    starts


            );
        }
    }

    public void close()
    {
        pw.close();
    }



}

