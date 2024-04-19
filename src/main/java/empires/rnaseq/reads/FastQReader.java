package nlEmpiRe.rnaseq.reads;

import lmu.utils.*;
import lmu.utils.index.Indexer;
import lmu.utils.index.StringBasedDoubleRAFIndex;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.BitSet;
import java.util.Locale;
import java.util.Vector;


public class FastQReader
{

    public int lineno=0;
    public int read_num=0;

    File f;
    //BufferedRandomAccessFile bfr = null;
    BufferedReader bfr = null;
    boolean gotnext=false;
    FastQRecord record = new FastQRecord();

    int num_skipped_lines = 0;

    public FastQReader(File f)
    {
        this.f=f;
        bfr = FileUtils.getReader(f);
        /*
        try
        {
            bfr = new BufferedRandomAccessFile(f, "r", 1000000);
        }
        catch(IOException e)
        {
            throw new RuntimeException("cannot open fastq file: "+f.getAbsolutePath());
        }
        */
    }

    public long getFP()
    {
        return 0;
        /*
        try
        {
            return bfr.getFilePointer();
        }
        catch(IOException ie)
        {
            throw new FRuntimeException("i/o error: %s while getting fp", ie, ie.getMessage());
        }
        */

    }
    public FastQRecord current()
    {
        return record;
    }



    //public static FastQRecord readRecord(BufferedRandomAccessFile bfr)
    public static FastQRecord readRecord(BufferedReader bfr)
    {
        return next(bfr, new FastQRecord(), null);
    }

    public static FastQRecord parse(String header, String read, String qual)
    {
        return parse(new FastQRecord() ,header, read, qual);
    }

    public static FastQRecord parse(FastQRecord record, String header, String read, String qual)
    {
        record.setHeader(header);
        record.setSeq(read);
        record.setQual(qual);
        record.valid=true;
        return record;
    }

    static FastQRecord nextWithHeaderSearch(BufferedReader bfr, FastQRecord record, FastQReader fqr, String header_prefix) {
        record.valid = false;
        if (bfr == null) {
            return null;
        }
        try {
            //long fp = bfr.getFilePointer();
            long fp = 0;
            String header = null;
            while (true) {

                header = bfr.readLine();
                if (header == null) {
                    //end of file
                    fqr = null;
                    bfr.close();
                    bfr = null;
                    return null;

                }
                if (header.length() == 0 || !header.startsWith(header_prefix)) {
                    header = null;
                    fqr.num_skipped_lines++;
                    continue;
                }


                break;
            }

            if(header == null) {
                fqr = null;
                bfr.close();
                bfr = null;
                return null;
            }


            return next(bfr, record, fqr, header);


        } catch (IOException e) {
            if (fqr != null) {
                throw new RuntimeException("i/o error:" + e.getMessage() + " after read: " + fqr.read_num + " line: " + fqr.lineno + " while reading fastq file:" + fqr.f.getAbsolutePath());
            }
            throw new RuntimeException("i/o error:" + e.getMessage() + " while reading next record from " + bfr);

        }
    }

    static FastQRecord next(BufferedReader bfr, FastQRecord record, FastQReader fqr)
    {
        record.valid=false;
        if (bfr==null)
        {
            return null;
        }
        try {
            //long fp = bfr.getFilePointer();
            long fp = 0;
            String header = bfr.readLine();
            return next(bfr, record, fqr, header);
        }
        catch(IOException e)
        {
            if (fqr != null)
            {
                throw new RuntimeException("i/o error:"+e.getMessage()+" after read: "+ fqr.read_num+" line: "+fqr.lineno+" while reading fastq file:" + fqr.f.getAbsolutePath());
            }
            throw new RuntimeException("i/o error:"+e.getMessage()+" while reading next record from " + bfr);

        }

    }
    //static FastQRecord next(BufferedRandomAccessFile bfr, FastQRecord record, FastQReader fqr)
    static FastQRecord next(BufferedReader bfr, FastQRecord record, FastQReader fqr, String header)
    {
        record.valid=false;
        if (bfr==null)
        {
            return null;
        }

        try
        {
            //long fp = bfr.getFilePointer();
            long fp = 0;
            //String header = bfr.readLine();
            if (header==null)
            {
                //end of file
                bfr.close();
                bfr=null;
                return null;

            }

            String read = bfr.readLine();
            String dummy = bfr.readLine();
            String qual = bfr.readLine();

            if (qual==null)
            {
                if (fqr != null)
                {
                    throw new RuntimeException("invalid fastq entry, could not parse 4 lines after read: "+ fqr.read_num+" line: "+fqr.lineno+" while reading fastq file:" + fqr.f.getAbsolutePath());
                }
                throw new FRuntimeException("invalid fastq entry, could not parse 4 lines!\nheader:%s\nread: %s\nline3: %s\nqual: %s\n", header, read, dummy, qual);

            }


            if (fqr!=null)
            {
                fqr.lineno+=4;
                fqr.read_num++;
            }

            parse(record,header,read,qual);
            record.fp=fp;

        }
        catch(IOException e)
        {
            if (fqr != null)
            {
                throw new RuntimeException("i/o error:"+e.getMessage()+" after read: "+ fqr.read_num+" line: "+fqr.lineno+" while reading fastq file:" + fqr.f.getAbsolutePath());
            }
            throw new RuntimeException("i/o error:"+e.getMessage()+" while reading next record from " + bfr);

        }

        return record;

    }

    public FastQRecord next(FastQRecord r)
    {
        return next(bfr,r, this);
    }

    public FastQRecord next(String check_header_prefix) {
        return (check_header_prefix == null) ? next() : nextWithHeaderSearch(bfr, new FastQRecord(), this, check_header_prefix);
    }
    public FastQRecord next()
    {
        return next(bfr,record, this);
    }

    public static Indexer<FastQRecord> getIndexer(File f, boolean gzipped)
    {
        return (f == null ) ? null : new FastQReader.FastQIndex(f, gzipped);
    }

    public static class FastQIndex implements Indexer<FastQRecord>
    {
        StringBasedDoubleRAFIndex sindex = null;
        FileUtils.RecordIterator rit = new FileUtils.RecordIterator('\n');
        public FastQIndex(File f, boolean gzipped)
        {
            sindex = new StringBasedDoubleRAFIndex(f,gzipped);
        }


        public FastQRecord read(int idx)
        {
            String[] sp = rit.split(sindex.read(idx));
            return FastQReader.parse(sp[0], sp[1], sp[3]);
        }

        public void write(FastQRecord fq, int idx)
        {
            sindex.write(String.format("%s\n%s\n%s\n%s", fq.header.toString(), fq.readseq.toString(), fq.header.toString(), fq.qualstring.toString()), idx);
        }

        public int size()
        {
            return sindex.size();
        }
    }


}



