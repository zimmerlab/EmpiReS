package nlEmpiRe.rnaseq.reads;

import lmu.utils.StringUtils;

import java.io.PrintWriter;
import lmu.utils.*;
import java.util.*;
import java.io.*;

public class FastQRecord
{
    static final int NORMQUAL = 33;

    public long fp;
    boolean valid=false;
    public int readid=-1;
    public StringBuffer header = new StringBuffer();
    public StringBuffer readseq = new StringBuffer();
    public StringBuffer qualstring = new StringBuffer();
    int[] quality = null;
    boolean qual_calced=false;

    public int readlength=0;

    public void setHeader(String h)
    {
        header.setLength(0);
        header.append(h);
    }

    public int getCombinedQual(int pos)
    {
        if (readseq.charAt(pos)=='N')
            return 0;

        if (!qual_calced)
        {
            getQuality();
        }

        return quality[pos];
    }
    void realloc(int length)
    {
        quality=new int[length];
    }

    public StringBuffer setQualInfo(int level, char below_mark, StringBuffer in)
    {
        int[] q = getQuality();
        for (int i=0; i<q.length; i++)
        {
            if (q[i] >= level)
                continue;

            in.setCharAt(i, below_mark);
        }
        return in;
    }

    public StringBuffer getQualInfo(int level, char below_mark)
    {
        return setQualInfo(level, below_mark, StringUtils.initBuffer(getQuality().length));
    }




    public int[] getQuality()
    {
        if (quality==null || quality.length<qualstring.length())
        {
            int nlength=(int)(100*(1+(qualstring.length()/100.0)));
            realloc(nlength);
            qual_calced=false;
        }

        if (qual_calced)
        {
            return quality;
        }

        final int ql=qualstring.length();
        for (int i=0; i<ql; i++)
        {
            quality[i]=qualstring.charAt(i)-NORMQUAL;
        }
        qual_calced=true;
        return quality;
    }
    public void setQual(String qseq)
    {
        qual_calced=false;
        qualstring.setLength(0);
        qualstring.append(qseq);
    }
    public void setSeq(String seq)
    {
        readseq.setLength(0);
        readseq.append(seq);
        readlength=readseq.length();
    }

    public void write(PrintWriter pw, String nid)
    {
        nid = (nid == null) ? header.substring(1) : nid;
        pw.printf("@%s\n%s\n+%s\n%s\n", nid, readseq.toString(), nid, qualstring.toString());
    }

    public void writeTrimmed(PrintWriter pw, String nid, int left, int right)
    {
        if (right <= 0 )
        {
            right = readlength + right;
        }

        String nrs = readseq.substring(left, right);
        String nqual = qualstring.substring(left, right);
        if (nid==null)
        {
            nid=header.toString();
        }
        pw.printf("@%s\n%s\n+%s\n%s\n", nid, nrs, nid, nqual);
    }

    public String getFirstHeaderPart()
    {
        int wsp = header.indexOf(" ");
        wsp = (wsp < 0) ? header.length() : wsp;
        return header.substring(1, wsp);
    }
    public void trim(int left, int right)
    {
        setSeq(readseq.substring(left, right));
        setQual(qualstring.substring(left, right));
    }

    public String toString()
    {
        return String.format("%s\n%s\n%s\n%s\n", header.toString(), readseq.toString(), header.toString(), qualstring.toString());
    }
}



