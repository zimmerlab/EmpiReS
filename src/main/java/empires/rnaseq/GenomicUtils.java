package empires.rnaseq;



import lmu.utils.*;
import lmu.utils.tuple.Tuple3;

import java.util.Collections;
import java.util.Vector;
import java.util.regex.Pattern;

public class  GenomicUtils
{

    static Pattern chrP = Pattern.compile(":-\\(\\)\\s+");
    static FileUtils.StringSplitter chrInfoSplitter = FileUtils.getTokenizer(":- ()/,;");
    static FileUtils.StringSplitter chrInfoSplitterChrom = FileUtils.getTokenizer(": ()/,;");

    static Tuple3<String, Integer,Integer> parseRegion(String chr, String r1, String r2)
    {
        try
        {
            return Tuple3.create(chr, (r1.equals(".")) ? 0 : Integer.parseInt(r1), (r2.equals(".")) ? Integer.MAX_VALUE : Integer.parseInt(r2));
        }
        catch(NumberFormatException nfe)
        {

        }
        return null;
    }

    static public GenomicRegion parseStrandedChrRegion(String line)
    {
        Tuple3<String, Integer,Integer> t = parseChrRegion(line);

        Boolean strand = (t.get0().length() == 0) ? null : GenomicUtils.getStrand(t.get0().substring(t.get0().length() -1));
        String chr = (strand == null) ? t.get0() : StringUtils.substring(t.get0(), -1);
        return new GenomicRegion(chr, strand, t.get1(), t.get2());
    }

    static public Tuple3<String, Integer,Integer> parseChrRegion(String line)
    {
        String chr = chrInfoSplitterChrom.split(line)[0];
        String patched = "X:" + line.substring(chr.length());
        String sp[] = chrInfoSplitter.split(patched);

        if (sp.length != 3)
            return null;

        try
        {
            return parseRegion(chr, sp[1], sp[2]);
        }
        catch(NumberFormatException nfe)
        {

        }
        return null;
    }


    public static char getStrand(Boolean b)
    {
        return (b == null) ? ' ' : (b) ? '+' : '-';
    }


    public static Boolean getStrand(String s)
    {
        if (s == null  || s.length() == 0)
            return null;


        if (s.charAt(0) == '+' ||  s.charAt(0) == '-')
            return s.charAt(0) == '+';

        try
        {
            return Boolean.parseBoolean(s);
        }
        catch (Exception e)
        {
            //silent error
        }
        System.err.printf("strand: >%s< converted to null", s);
        return null;
    }

    public static final String CANONICAL_SPLICE_MOTIF = "GTAG";


    static public Vector<Region1D> convertTranscriptRegions2GenomicRegions(Vector<Region1D> genomic_coding_regions, Vector<Region1D> toconvert)
    {
        Vector<Region1D> rv = new Vector<Region1D>();

        int start_exon=0;
        Region1D ce = genomic_coding_regions.get(0);

        int tr_start = 0;
        int tr_end = ce.getLength();

        for(Region1D r : toconvert)
        {
            int r_trstart = r.getX1()*3;
            int r_trend = r.getX2()*3;
            while (tr_end < r_trstart)
            {
                start_exon++;
                ce = genomic_coding_regions.get(start_exon);
                tr_start = tr_end;
                tr_end = tr_start+ ce.getLength();

            }


            while (r_trend > tr_end)
            {
                int offset = Math.max(0, r_trstart-tr_start);
                rv.add(new Region1D(ce.getX1()+offset, ce.getX2(), null));

                start_exon++;
                if (start_exon == genomic_coding_regions.size())
                {
                    return ObjectGetter.filter(rv, (_r) -> _r.getLength() > 0);
                }
                ce = genomic_coding_regions.get(start_exon);
                tr_start = tr_end;
                tr_end = tr_start+ce.getLength();

            }

            int offset = Math.max(0, r_trstart-tr_start);
            rv.add(new Region1D(ce.getX1()+offset, ce.getX1()+(r_trend-tr_start), null));
        }

        return ObjectGetter.filter(rv, (r) -> r.getLength() > 0);
    }


    static void transcribed2RegionVector(Pair<Integer,Integer> tr_region, RegionVector global_regions, RegionVector addto)
    {
        int prelength=0;
        int trlength=0;
        boolean in_region=false;
        for (Region1D r : global_regions.getRegions())
        {
            prelength = trlength;
            trlength += r.getLength();

            if (tr_region.getSecond()<prelength)
                return;

            if (tr_region.getFirst()>=trlength)
                continue;


            int in_exstart = (tr_region.getFirst()<=prelength)?r.getX1():r.getX1()+tr_region.getFirst()-prelength;

            int in_exend = (trlength < tr_region.getSecond())?r.getX2(): r.getX1()+tr_region.getSecond()-prelength;

            addto.addRegion(new Region1D(in_exstart, in_exend, r.getObject()));

        }
    }

    static public RegionVector transcribed2RegionVector(Pair<Integer,Integer> tr_region, RegionVector global_regions)
    {
        RegionVector rv = new RegionVector();
        transcribed2RegionVector(tr_region,global_regions, rv);
        return rv;
    }
    static public RegionVector transcribed2RegionVector(Vector<Pair<Integer,Integer>> tr_regions, RegionVector global_regions)
    {
        RegionVector rv = new RegionVector();
        for (Pair<Integer,Integer> p: tr_regions)
        {
            transcribed2RegionVector(p,global_regions, rv);
        }

        return rv;
    }

    public static int toGlobal(int local, int g_start, int g_end,  boolean strand)
    {
        return (strand) ? g_start+local : g_end+1-local;
    }

    public static Pair<Integer,Integer> toGlobal(Pair<Integer,Integer> local, int g_start, int g_end,  boolean strand)
    {
        if (strand)
        {
            return new Pair<Integer,Integer>(g_start+local.getFirst(), g_start+local.getSecond());
        }
        return new Pair<Integer,Integer>(g_end+1-local.getSecond(), g_end+1-local.getFirst());
    }

    public static Region1D toGlobal(Region1D local, int g_start, int g_end,  boolean strand)
    {
        if (strand)
        {
            return new Region1D(g_start+local.getX1(), g_start+local.getX2(), null);
        }
        return new Region1D(g_end+1-local.getX2(), g_end+1-local.getX1(), null);
    }

    public static RegionVector toGlobal(RegionVector local, int g_start, int g_end,  boolean strand)
    {
        RegionVector global = new RegionVector();


        if (strand)
        {
            for (Region1D r : local.getRegions())
            {
                if (local.getLength()==0)
                {
                    continue;
                }
                global.addRegion(r.getX1()+g_start, r.getX2()+g_start);
            }

            return global;
        }
        //reverse
        Vector<Region1D> v = local.getRegions();
        final int E = g_end + 1;
        for (int i=v.size()-1; i>=0; i--)
        {
            Region1D r = v.get(i);

            if (r.getLength()==0)
            {
                continue;
            }

            global.addRegion(E - r.getX2(), E - r.getX1());
        }

        return global;

    }

    public static int toLocal(int pos, int g_start, int g_end, boolean strand)
    {
        return (strand) ? pos - g_start : g_end + 1 - pos;
    }

    public static Pair<Integer,Integer> toLocal(Pair<Integer,Integer> global, int g_start, int g_end, boolean strand)
    {
        if(strand)
        {
            return new Pair<Integer,Integer>(global.getFirst()-g_start, global.getSecond()-g_start);
        }
        return new Pair<Integer,Integer>(g_end+1-global.getSecond(), g_end+1-global.getFirst());
    }

    public static Region1D toLocal(Region1D global, int g_start, int g_end, boolean strand)
    {
        if(strand)
        {
            return new Region1D(global.getX1()-g_start, global.getX2()-g_start, global.getObject());
        }
        return new Region1D(g_end+1-global.getX2(), g_end+1-global.getX1(), global.getObject());
    }

    public static RegionVector toLocal(RegionVector global, int g_start, int g_end, boolean strand)
    {
        RegionVector rv = new RegionVector();
        rv.setObject(global.getObject());
        if (strand)
        {
            for (Region1D r : global.getRegions())
            {
                rv.addRegion(toLocal(r,g_start,g_end,strand));
            }
        }
        else
        {
            Vector<Region1D> v = global.getRegions();
            for (int i=v.size()-1; i>=0; i--)
            {
                Region1D r = v.get(i);
                rv.addRegion(toLocal(r,g_start,g_end,strand));
            }


        }
        return rv;
    }

    /*


    public static class ProteinSeqCheck
    {
        Logger log = LogConfig.getLogger();
        Protein p;
        public Vector<Integer> shifts = new Vector<Integer>();
        public Vector<Pair<Integer,Integer>> cuts =new Vector<Pair<Integer,Integer>> ();

        public Vector<String> translated = new Vector<String> ();

        public boolean multi_src;
        public boolean single_src;

        public int getNumSolutions()
        {
            return shifts.size();
        }

        public String getTranslated()
        {
            if (!single_src)
            {
                throw new RuntimeException("no single valid solution!");
            }

            return getTranslated(0);
        }
        public String getTranslated(int sol)
        {
            return translated.get(sol);
        }
        public Vector<Pair<Integer,Integer>> getFixedProteinRegionsInGene()
        {
            if (!single_src)
            {
                throw new RuntimeException("no single valid solution!");
            }
            return getFixedProteinRegionsInGene(0);
        }
        /*
        public Vector<Pair<Integer,Integer>> getFixedProteinRegionsInGene(int sol)
        {
            Vector<Pair<Integer,Integer>> rv = new Vector<Pair<Integer,Integer>>();
            int shift=shifts.get(sol);
            Pair<Integer,Integer> cut = cuts.get(sol);
            String ps = translated.get(sol);




            int n2skip = shift+cut.getFirst()*3;

            int coding_length=0;
            for (CodingExon ce: p.coding_exons)
            {
                coding_length+=ce.length;
            }
            int endcut=coding_length-n2skip-3*ps.length();
            int lastpos=coding_length-endcut;

            int celength=0;
            for (CodingExon ce: p.coding_exons)
            {
                int start=ce.start_in_gene;
                int end = ce.end_in_gene+1;

                if (n2skip>0)
                {
                    start=Math.min(start+n2skip, end);
                    n2skip-=(start-n2skip);
                }
                celength+=(end-start);
                boolean stophere=false;
                if (celength>lastpos)
                {
                    end-=(celength-lastpos);
                    stophere=true;
                }


                rv.add(new Pair<Integer,Integer>(start,end));

                if (stophere)
                {
                    break;
                }


            }
            return rv;
        }
        public ProteinSeqCheck(Protein p, int max_no_aa_accepted, GenomeSequenceExtractor extractor, String gene_seq)
        {
            this.p=p;
            String prot_dna_seq = null;
            String coding_prot_dna_seq = null;

            if (extractor!=null)
            {
                prot_dna_seq=GenomicUtils.getProteinSeq(p, extractor, true);
            }
            else
            {
                RegionVector v = p.getCodingRegionVector(false);
                prot_dna_seq = GenomicUtils.getSplicedSeq(gene_seq, v, true);
            }



            for (int i=0; i<3; i++)
            {
                String ps=GenomicUtils.translate(prot_dna_seq, 0, i);

                int nx=0;
                Pair<Integer,Integer> cutPos = protein.DigeUtils.getXCut(ps);

                if (cutPos.getSecond()<0)
                {
                    continue;
                }

                for (int j=cutPos.getFirst(); j<cutPos.getSecond()+1; j++)
                {
                    if (ps.charAt(j)!='X')
                        continue;

                    nx++;
                }
                if (nx>max_no_aa_accepted)
                {
                    //too much x , skip it
                    continue;
                }

                translated.add(ps.substring(cutPos.getFirst(),cutPos.getSecond()));
                shifts.add(i);
                cuts.add(cutPos);

            }

            multi_src=translated.size()>1;
            single_src=translated.size()==1;

        }



    }
    */
    public static char REVERSE_DNA[]=new char[Character.MAX_VALUE];
    public static char DNA2RNA[]=new char[Character.MAX_VALUE];

    static
    {
        for (char c=0; c<Character.MAX_VALUE; c++)
        {
            REVERSE_DNA[c]='N';
            DNA2RNA[c]='N';
        }
        REVERSE_DNA['T']='A';
        REVERSE_DNA['A']='T';
        REVERSE_DNA['G']='C';
        REVERSE_DNA['C']='G';

        REVERSE_DNA['t']='A';
        REVERSE_DNA['a']='T';
        REVERSE_DNA['g']='C';
        REVERSE_DNA['c']='G';


        DNA2RNA['T']='A';
        DNA2RNA['A']='U';
        DNA2RNA['G']='C';
        DNA2RNA['C']='G';

        DNA2RNA['t']='A';
        DNA2RNA['a']='U';
        DNA2RNA['g']='C';
        DNA2RNA['c']='G';
    }

    public static String antisense(String s)
    {
        return antisense(s,false);
    }
    public static String antisense(String s, boolean get_rna)
    {
        final char[] conv = (get_rna)? DNA2RNA : REVERSE_DNA;

        StringBuffer sb = new StringBuffer();
        for (int i=0; i<s.length(); i++)
        {
            sb.append(conv[s.charAt(i)]);

        }
        return sb.toString();
    }
    public static String reverse_complement(String s)
    {
        StringBuffer sb = new StringBuffer();
        for (int i=s.length()-1; i>=0; i--)
        {
            sb.append(REVERSE_DNA[s.charAt(i)]);

        }
        return sb.toString();
    }


    public static String getSplicedSeq(final String chr, boolean strand, RegionVector regvec, GenomeSequenceExtractor fromdna, boolean dna)
    {
        return getSplicedSeq(chr,strand,regvec.getRegionPairs(),fromdna,dna);
    }

    public static String getSplicedSeq(final String chr, boolean strand, Vector<Pair<Integer,Integer>> regions, GenomeSequenceExtractor fromdna, boolean dna)
    {
        return getSplicedSeq(strand,chr,regions,fromdna,null,dna);
    }

    public static String getSplicedSeq(String seq, RegionVector regionvec, boolean dna)
    {
        return getSplicedSeq(seq, regionvec.getRegionPairs(), dna);
    }

    public static String getSplicedSeq(String seq, Vector<Pair<Integer,Integer>> regions, boolean dna)
    {
        return getSplicedSeq(true,null,regions, null, seq,dna);
    }

    public static String getSplicedSeq(boolean strand, final String chr, Vector<Pair<Integer,Integer>> regions, GenomeSequenceExtractor fromdna,
                                       String relative, boolean dna)
    {

        Collections.sort(regions);
        StringBuffer sb = new StringBuffer();
        int real_length=0;
        for(int i=0; i<regions.size(); i++)
        {
            Pair<Integer,Integer> region = regions.get(i);
            if (region.getSecond()<=region.getFirst())
            {
                continue; //nothing to extract
            }
            String genomeseq = null;
            if (relative!=null)
            {
                //checks
                if(region.getFirst()<0 || region.getFirst()>relative.length())
                {
                    throw new RuntimeException("Invalid region: "+region+" (outof bounds in relative of length: "+relative.length()+") regions: "+regions);

                }
                if (region.getSecond()>relative.length())
                {
                    throw new RuntimeException("Invalid region: "+region+" (outof bounds in relative of length: "+relative.length()+") regions: "+regions);

                }

                if (region.getFirst()<0 || region.getSecond()<0 || region.getSecond()<=region.getFirst() || region.getSecond()>relative.length())
                {
                    throw new RuntimeException("invalid region: "+region+" in regions: "+
                            regions+" relative length: "+relative.length());
                }
                genomeseq = relative.substring(region.getFirst(), region.getSecond());
            }
            else
            {
                genomeseq = fromdna.getGenomeSeq(chr, region.getFirst(), region.getSecond());
            }

            if (genomeseq==null)
            {
                throw new RuntimeException(String.format("error by extracting %s for %s\n", ""+region, chr));

            }
            //System.out.println("APPEND:\n"+genomeseq.toString());
            sb.append(genomeseq);
        }

        if (!strand)
        {
            String revsb = reverse_complement(sb.toString());
            sb.setLength(0);
            sb.append(revsb);
        }

        String protdna = sb.toString();
        if (dna)
        {

            return protdna;
        }
        throw new FRuntimeException("translate is not yet migrated to nlEmpiRe...");

    }




}
