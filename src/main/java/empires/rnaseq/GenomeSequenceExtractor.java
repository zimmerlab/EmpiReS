package nlEmpiRe.rnaseq;

import lmu.utils.*;

import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.apache.logging.log4j.Logger;



public class GenomeSequenceExtractor
{

    final String MAPPER_INDEX = ".chr.mapper.index";
    Logger log = LogConfig.getLogger();
    /** allows to read regions from genome
     * from a single fasta file with entry per chromosome
     * or a directory containing such
     */


    public interface ChromosomeNameMapper
    {
        //returns null if not a chromosome file, otherwise the name
        public String getName(String filename);
    }
    static Pattern std_chrom_name_pat = Pattern.compile(".+.dna.chromosome.([^\\.]+).fa");

    static class DefaultChromosomeNameMapper implements ChromosomeNameMapper
    {
        public String getName(String filename)
        {
            Matcher m = std_chrom_name_pat.matcher(filename);
            if (!m.matches())
                return null;

            return "chr"+m.group(1);

        }
    }
    static DefaultChromosomeNameMapper defaultMapper = new DefaultChromosomeNameMapper();
    static class ChrLookup
    {
        String chrname;
        File f;
        BufferedRandomAccessFile raf;
        long pos;

        int line_width = 60;
        int line_width_c = line_width+1;

        public ChrLookup(File f, BufferedRandomAccessFile raf, String name, long p)
        {
            this(f, raf, name, p, 60, 61);
        }

        public ChrLookup(File f, BufferedRandomAccessFile raf, String name, long p, int lw, int lwc)
        {

            this.raf=raf; this.chrname=name; pos=p; this.f = f;
            line_width = lw; line_width_c = lwc;
        }
    }



    HashMap<String, ChrLookup> chr2lookup = null;


    synchronized void initChrLookup(File root, File index, ChromosomeNameMapper mapper)
    {
        if(root == null || !root.exists())
            throw new FRuntimeException("invalid dir/toplevel file: %s for genome sequence", root.getAbsolutePath());

        if (index==null && (root != null && root.exists() && root.isDirectory()))
        {
            index = new File(root, MAPPER_INDEX);
        }
        //check if there is an index File

        if (index != null && index.exists())
        {
            chr2lookup = new HashMap<String, ChrLookup>();

            BufferedRandomAccessFile single_raf = null;
            //try to read
            HashMap<String, BufferedRandomAccessFile> path2buff = new HashMap<String, BufferedRandomAccessFile>();
            for(String[] sp : FileUtils.getFieldSets(index, "\t"))
            {
                if (sp.length < 3) {
                    log.info("skip line: %s", StringUtils.join("\t", sp));
                    continue;
                }
                String chr = sp[0];
                //check format -> if second is a file then own variant, otherwise fi
                File fastafile = new File(root, sp[1]);
                if(fastafile.exists())
                {
                    long fp = Long.parseLong(sp[2]);
                    log.trace("%s:%s:%d\n", sp[0], fastafile.getAbsolutePath(), fp);
                    BufferedRandomAccessFile raf = path2buff.get(fastafile.getAbsolutePath());
                    if (raf == null)
                    {
                        try {
                            raf = new BufferedRandomAccessFile(fastafile, "r", 30000);
                        }
                        catch(IOException ie)
                        {
                            throw new FRuntimeException("i/o error at opening RAF on %s", ie, ie.getMessage(), fastafile.getAbsolutePath());
                        }
                        path2buff.put(fastafile.getAbsolutePath(), raf);
                    }

                    //System.out.printf("%s\n", line);
                    chr2lookup.put(chr, new ChrLookup(fastafile, raf, chr, fp));
                    continue;
                }

                //.fi file : (chr, start, length, line length
                //but if so -> the  root must be a file
                if(!root.isFile())
                    throw new FRuntimeException("invalid index format %s for %s", index.getAbsolutePath(), root.getAbsolutePath());

                if(sp.length != 5)
                    throw new FRuntimeException("invalid fai index format expected 5 fields, but got: %d", sp.length);
                try {
                    single_raf = (single_raf != null) ? single_raf : new BufferedRandomAccessFile(root, "r", 30000);
                }
                catch(IOException ie)
                {
                    throw new FRuntimeException("i/o error at opening RAF on %s", ie, ie.getMessage(), root.getAbsolutePath());
                }
                long fp = Long.parseLong(sp[2]);
                int line_length = Integer.parseInt(sp[3]);

                int line_length_c = Integer.parseInt(sp[4]);
                chr2lookup.put(chr, new ChrLookup(root, single_raf, chr, fp, line_length, line_length_c));

            }
            log.info("read %d chromosomes from index file", chr2lookup.size(), index.getAbsolutePath());
        }
        if (chr2lookup!=null)
            return;

        //load 
        chr2lookup = new HashMap<String, ChrLookup>();


        File toplevel = (root.isDirectory()) ? null : root;
        File fna = null;

        if(toplevel == null) {


            for (String s : root.list()) {
                if (s.contains("dna.toplevel") && (s.endsWith(".fa") || s.endsWith(".fasta"))) {
                    toplevel = new File(root, s);
                    continue;
                }
                if (s.endsWith(".fna")) {
                    fna = new File(root, s);
                }
            }
            if (toplevel == null && fna == null) {
                //try to read per chromosome fasta files
                for (String s : root.list()) {
                    String chrname = mapper.getName(s);
                    if (chrname == null) {
                        log.debug("skip file: " + s);
                        continue;
                    }
                    File f = new File(root, s);
                    try {
                        BufferedRandomAccessFile raf = new BufferedRandomAccessFile(f, "r", 10000);
                        String header = raf.readLine();
                        chr2lookup.put(chrname, new ChrLookup(f, raf, chrname, header.length() + 1));
                        log.debug("add chr: " + chrname + ": " + s);
                    } catch (IOException e) {
                        throw new RuntimeException("error while initializing raf from " + f.getAbsolutePath() + " for fasta seek", e);
                    }
                }

            }
        }
        File target = (toplevel!=null)?toplevel:fna;

        if (toplevel==null)
        {
            if (chr2lookup.size()==0)
            {
                throw new RuntimeException("could not initalize genome from "+root.getAbsolutePath()+" neither chromosme fastas, nor toplevel fasta, nor .fna file found!");
            }

        }
        if (target!=null)
        {
            try
            {

                BufferedRandomAccessFile raf = new BufferedRandomAccessFile(target, "r", 30000);
                String line=null;
                int nadded=0;
                while (null!=(line=raf.readLine()))
                {
                    if (line.length()==0 || line.charAt(0)!='>')
                        continue;

                    String chr = null;
                    if (fna == target)
                    {
                        String sp[] = line.split("\\|");
                        if (sp.length<4)
                        {
                            log.debug("unexcpected .fna entry: "+line+" could not extract chr ...");
                            continue;
                        }
                        chr = sp[3];
                    }
                    else
                    {
                        chr = line.split("\\s+")[0].substring(1);
                    }


                    chr2lookup.put(chr, new ChrLookup(target, raf, chr, raf.getFilePointer()));

                    log.info("next chr from toplevel: %s: %s\n%s\n%s\n", chr, target.getAbsolutePath(), line, raf.readLine());
                    nadded++;
                    //System.out.println("line: "+line+" fp: "+(raf.getFilePointer()-(long)line.length()));
                }

                log.info("added "+ nadded+" chrs from "+target.getAbsolutePath());
            }
            catch(IOException e)
            {
                throw new FRuntimeException("i/o error: %s while reading genome from %s", e, e.getMessage(), target.getAbsolutePath());
            }
        }


        //cannot write index file 
        try
        {
            if (!index.exists() && !index.createNewFile())
            {
                log.warn("cannot create index file to "+index.getAbsolutePath());
                return;
            }
        }
        catch(IOException e)
        {
            log.warn("i/o error while creating index file: "+index.getAbsolutePath());
            return;
        }


        if (!index.canWrite())
        {
            log.warn("cannot write index file to "+index.getAbsolutePath());
            return;
        }

        //write index file
        try
        {
            PrintWriter pw = new PrintWriter(new FileWriter(index));
            for (ChrLookup cl: chr2lookup.values())
            {
                pw.printf("%s\t%s\t%d\n", cl.chrname, cl.f.getName(), cl.pos);
            }
            pw.close();
        }
        catch(IOException ie)
        {
            log.warn("i/o error while writing index file: "+index.getAbsolutePath());
        }


    }

    public HashMap<String, Long> getChr2Length()
    {
        HashMap<String, Long> rv = new HashMap<>();
        HashMap<File, Vector<ChrLookup>> f2chr = new HashMap<>();

        ObjectGetter.apply(chr2lookup.values(), (_c) -> MapBuilder.updateV(f2chr, _c.f, _c));

        for (Vector<ChrLookup> v : f2chr.values())
        {
            NumUtils.sort(v, (c) -> c.pos, false);

            long last = v.get(0).f.length();
            for (int i=0; i<v.size(); i++)
            {
                String chr = v.get(i).chrname;
                long p = v.get(i).pos;
                long next = (i == v.size() -1) ? last : v.get(i+1).pos;
                rv.put(chr, next - p);
            }
        }

        return rv;

    }
    File dna_dir;

    /** if matcher does not match the file is not considered,
     * if matches, group(1) should return the chromosome name 
     *  if not specified the standard pattern is used
     */
    public GenomeSequenceExtractor(File fastadir)
    {
        this(fastadir, null);
    }

    public GenomeSequenceExtractor(File fastadir, File indexfile)
    {
        this(fastadir, indexfile, (ChromosomeNameMapper)null);
    }
    public GenomeSequenceExtractor(File fastadir, File indexfile, ChromosomeNameMapper mapper)
    {

        dna_dir = fastadir;
        if (mapper==null)
        {
            mapper=defaultMapper;
        }
        /*if (!fastadir.isDirectory())
        {
            throw new RuntimeException("GenomeSequenceExtractor constructor expected a directory as input, got "+fastadir.getAbsolutePath());
        }
        */

        initChrLookup(fastadir, indexfile, mapper);


    }

    public boolean gotChr(String chr)
    {
        return null!=getChr(chr,false);
    }

    public Collection<String> getChromosomes()
    {
        return chr2lookup.keySet();
    }

    ChrLookup getChr(String chr, boolean throw_error)
    {
        ChrLookup chrl = chr2lookup.get(chr);

        if (chrl==null)
        {
            chrl=chr2lookup.get("chr"+chr);
        }

        if (chrl==null && throw_error)
        {

            throw new RuntimeException("unknown chr:"+chr+" got followings: "+chr2lookup.keySet());
        }

        return chrl;
    }

    StringBuffer intern_buffer = new StringBuffer();

    boolean verbose = false;

    void setVerbose(boolean b)
    {
        verbose = b;
    }

    public synchronized String getGenomeSeq(String chr, int from, int to)
    {
        if (to<=from)
        {
            return null;
        }

        ChrLookup chrl = getChr(chr,true);
        int jump = (int)(from/(double)chrl.line_width);
        int inline_jump = from % chrl.line_width;

        int relative = chrl.line_width_c * jump+inline_jump;
        if (inline_jump==0)
        {
            relative--;
        }
        if (verbose)
        {
            System.out.printf("from: %d jump: %d inline nump: %d relative: %d\n", from, jump, inline_jump, relative);
        }
        try
        {

            intern_buffer.setLength(0);
            //StringBuffer rv =new StringBuffer();
            long p = chrl.pos-1+relative;
            p = Math.max(p, 0);
            chrl.raf.seek(p);
            int rlength=to-from;
            String line;
            int nread=0;
            while (intern_buffer.length()<rlength)
            {
                line = chrl.raf.readLine();
                if (line==null)
                {
                    log.error("invalid region requested: "+chr+" "+from+"-"+to+" after reading "+nread+" chars end of file dna dir: ", dna_dir.getAbsolutePath());
                    break;
                }
                intern_buffer.append(line);
                nread+=line.length();
            }
            while (intern_buffer.length()<rlength)
            {
                intern_buffer.append('N'); //fill up missing length
            }

            intern_buffer.setLength(rlength);
            return intern_buffer.toString();
        }
        catch(IOException e)
        {
            throw new RuntimeException("error while reading region "+chr+" "+from+"-"+to+" from "+chrl.f.getAbsolutePath()+":"+e.getMessage(), e);
        }
    }

    public String getGenomeSeq(String chr, int from, int to, boolean strand)
    {
        String seq = getGenomeSeq(chr,from,to);
        if (strand)
        {
            return seq;
        }
        return GenomicUtils.reverse_complement(seq);
    }
    /*public String getGenomeSeq(String chr, int from, int to, boolean strand)
    {
        return getGenomeSeq(chr,from,to,strand,false);
        
    }
    */
    /*public String getGenomeSeq(String chr, int from, int to)
    {
        StringBuffer sb = getGenomeSeqBuff(chr,from,to);
        if (sb==null)
        {
            return null;
        }
        return sb.toString();
        
    }
    */


    /*public Vector<GenomicRegion> readRegions(Iterator<GenomicRegion> regIt)
{
    Vector<GenomicRegion> v = new Vector<GenomicRegion>();
    while (regIt.hasNext())
    {
        v.add(regIt.next());
    }
    Collections.sort(v);


    for (GenomicRegion gr : v)
    {
        gr.sequence = getGenomeSeq(gr.chr, gr.region.getX1(), gr.region.getX2());
    }
    return v;
}
*/
    public static void main (String args[]) throws Exception
    {

        SimpleOptionParser cmd = new SimpleOptionParser("extractchrs", "test", "i", "index", "org","onlyindex");
        cmd.setDir( "extractchrs");
        cmd.setFile("i");
        //cmd.setFile("index");
        cmd.setSwitches("test", "onlyindex");
        cmd.setOptional("extractchrs", "index", "i", "org");

        if (!OptionParser.parseParams(args, true, false, true, cmd))
            return;


        //new GenomeSequenceExtractor(new File(args[0]), (File)null);
        GenomeSequenceExtractor gse = new GenomeSequenceExtractor(cmd.getFile("i"), cmd.getOptionalFile("index"));



        if (cmd.isSet("onlyindex"))
            return;


        if (cmd.isOptionSet("extractchrs"))
        {
            File od = cmd.getFile("extractchrs");

            HashMap<String, Long> ch2l = new HashMap<>();

            for(String chr : gse.getChromosomes())
            {
                File of = new File(od, chr+".fasta");
                long chrlength = 0;
                System.out.printf("write %s -> %s\n", chr, of.getAbsolutePath());
                try
                {
                    ChrLookup cl = gse.getChr(chr, false);
                    cl.raf.seek(cl.pos);


                    PrintWriter pw = new PrintWriter(of);
                    pw.printf(">%s\n", chr);
                    String line = null;
                    while (null != ( line = cl.raf.readLine()))
                    {
                        if (line.length()>0 && line.charAt(0) == '>')
                            break;

                        chrlength += line.length();
                        pw.println(line);
                    }
                    pw.close();
                }
                catch(IOException ie)
                {
                    System.err.printf("error: %s while writing chr file for %s (%s)\n", ie.getMessage(), chr, of.getAbsolutePath());
                }

                ch2l.put(chr, chrlength);
            }
            FileUtils.writeMap(new File(od, "chr.lengths"), "\t", ch2l);
            return;
        }

        System.out.printf(">");
        /*
        for(String chr : gse.getChromosomes())
        {
            System.out.println("check first 120 chars of chr: "+chr);
            //System.out.println("chck:\n"+gse.fl.get(chr)+"\n");
            System.out.println(gse.getGenomeSeq(chr, 0, 120)+"\n");
            
        }
        */
        BufferedReader bfr = new BufferedReader(new InputStreamReader(System.in));

        String last = null;
        String chr = "";
        int start = -1;
        int end = -1;
        String gseq = null;
        String line = null;
        while (null!=(line=bfr.readLine()))
        {
            if (line.equals("rev"))
            {
                if (gseq == null)
                {
                    System.out.println("no sequence to reverse complement...");
                    System.out.printf(">");
                    continue;

                }

                String rc = GenomicUtils.reverse_complement(gseq);
                System.out.printf("%s %d-%d reverse complement:\n%s\n", chr, start, end, rc);
                System.out.printf(">");
                continue;
            }
            String[] sp = line.split("\\s+|-|:|,");

            try
            {
                chr = sp[0];
                start = Integer.parseInt(sp[1]);
                end = Integer.parseInt(sp[2]);
                gseq = gse.getGenomeSeq(chr,start, end);
                System.out.printf("%s %d-%d:\n%s\n", chr, start, end, gseq);
            }
            catch(Exception e)
            {
                System.err.printf("invalid format: %s should be: chr start end\n", e.getMessage());
            }
            System.out.printf(">");

        }

    }




}
