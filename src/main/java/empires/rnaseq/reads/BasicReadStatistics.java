package empires.rnaseq.reads;

import lmu.utils.*;
import lmu.utils.index.BitSetIndex;

import javax.xml.stream.XMLStreamConstants;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamReader;
import javax.xml.stream.XMLStreamWriter;
import java.io.File;
import java.util.BitSet;
import java.util.Locale;
import java.util.Vector;

public class BasicReadStatistics
{
    int n=0;
    int lowqual_treshold = 20;
    static final int INIT_BITSETSIZE = 20_000_000;
    static final int BITSETSIZE_INCREMENT = 1_000_000;

    AutoResizeVector<AutoResizeBitSet> ns = new AutoResizeVector<>();

    AutoResizeVector<AutoResizeBitSet> lq = new AutoResizeVector<>();

    Vector<AutoResizeBitSet.IndexedBitSetContainer> ns_info = new Vector<>();
    Vector<AutoResizeBitSet.IndexedBitSetContainer> lq_info = new Vector<>();


    public BasicReadStatistics()
    {

    }

    public BasicReadStatistics(File xml, File bs_indexfile)
    {
        BitSetIndex bsi = (bs_indexfile == null) ? null : new BitSetIndex(bs_indexfile);

        try
        {
            XMLStreamReader xsr = XMLUtils.getReader(FileUtils.getReader(xml));
            while (xsr.hasNext())
            {
                int event = xsr.next();

                if (event != XMLStreamConstants.START_ELEMENT)
                    continue;

                if (xsr.getLocalName().equals("infos"))
                {
                    n = XMLUtils.getIntAttrib(xsr, "n", null);
                    continue;
                }
                if (xsr.getLocalName().equals("n"))
                {
                    int num_n = XMLUtils.getIntAttrib(xsr, "value", null);
                    Integer count = XMLUtils.getIntAttrib(xsr, "count", null);
                    Integer bs_index = XMLUtils.getIntAttrib(xsr, "bs_index", null);

                    if (bsi != null)
                    {
                        System.out.println("count: "+ count +" bs index: "+bs_index);
                        ns_info.add(new AutoResizeBitSet.IndexedBitSetContainer(count, bsi, bs_index));
                        continue;
                    }
                    try
                    {
                        BitSet bs = XMLUtils.readEncodedCData(xsr);
                        ns.set(num_n, new AutoResizeBitSet(bs, 1_000_000));
                    }
                    catch(ClassCastException ce)
                    {
                        throw new FRuntimeException("error in %s: n:%d cannot decode bitset from %d long text", xml, num_n, xsr.getText().length());
                    }

                    continue;
                }

                if (xsr.getLocalName().equals("lq"))
                {
                    int num_lq = XMLUtils.getIntAttrib(xsr, "value", null);
                    Integer count = XMLUtils.getIntAttrib(xsr, "count", null);
                    Integer bs_index = XMLUtils.getIntAttrib(xsr, "bs_index", null);

                    if (bsi != null)
                    {
                        lq_info.add(new AutoResizeBitSet.IndexedBitSetContainer(count, bsi, bs_index));
                    }

                    BitSet bs = XMLUtils.readEncodedCData(xsr);
                    lq.set(num_lq, new AutoResizeBitSet(bs, 1_000_000));
                    continue;
                }

            }
        }
        catch(XMLStreamException xse)
        {
            throw new FRuntimeException("xml error: %s while reading from %s", xse, xse.getMessage(), xml);
        }

    }

    public int getN()
    {
        return n;
    }


    Vector<BitSet> getCumulativeBitSets(AutoResizeVector<AutoResizeBitSet> arv)
    {
        Vector<BitSet> v = new Vector<BitSet>();
        v.add(null);
        BitSet last = null;
        for (int i=1; i<arv.size(); i++)
        {
            BitSet current = (arv.get(i) == null) ? null : arv.get(i).getBitSet();

            BitSet bs = (current == null) ? last : (last == null) ? current : SetInfo.union(last, arv.get(i).getBitSet());

            v.add(bs);
            last = bs;

        }
        BitSet all = new BitSet(n);



        for(int i=0; i<n; all.set(i++));

        if (last != null)
        {
            all = SetInfo.minus(all, last);
        }
        v.set(0, all);
        return v;
    }


    public Vector<BitSet> getNBitSets()
    {
        return getCumulativeBitSets(ns);
    }

    public Vector<BitSet> getLQBitSets()
    {
        return getCumulativeBitSets(lq);
    }

    public void update(int rid, FastQRecord rec)
    {
        n = Math.max(n, rid);
        int qual[] = rec.getQuality();
        final int L = rec.readlength;
        int numN = 0;
        int numLQ = 0;
        for (int i=0; i<L; i++)
        {
            if(rec.readseq.charAt(i) == 'N')
            {
                numN++;
            }

            if(qual[i] < lowqual_treshold)
            {
                numLQ++;
            }
        }

        if (numN != 0)
        {
            AutoResizeBitSet bs = ns.get(numN);
            if (bs == null)
            {
                ns.set(numN, bs = new AutoResizeBitSet(INIT_BITSETSIZE, BITSETSIZE_INCREMENT));
            }
            bs.set(rid);
        }

        if (numLQ != 0)
        {
            AutoResizeBitSet bs = lq.get(numLQ);
            if (bs == null)
            {
                lq.set(numLQ, bs = new AutoResizeBitSet(INIT_BITSETSIZE, BITSETSIZE_INCREMENT));
            }
            bs.set(rid);
        }

    }

    void printInfo()
    {
        System.out.printf("reads: %d\n", n);
        System.out.printf("num n-s\n");
        Vector<AutoResizeBitSet.IndexedBitSetContainer> info = getNInfo();


        double norm = 100.0 / n;

        for (int i=0; i<info.size(); i++)
        {
            int c = info.get(i).getCardinality();
            System.out.printf("\t%d: %d (%2.2f%%)\n", i, c, c * norm);
        }

        System.out.printf("num lq-s (treshold: %d)\n", lowqual_treshold);

        info = getLQInfo();


        for (int i=0; i<info.size(); i++)
        {
            int c = info.get(i).getCardinality();
            System.out.printf("\t%d: %d (%2.2f%%)\n", i, c, c * norm);
        }


    }


    public double getPercentNUnder(int p)
    {
        Vector<AutoResizeBitSet.IndexedBitSetContainer> v = getNInfo();
        if (p >= v.size())
            return 100.0;

        return v.get(p).getCardinality() / (double)n;
    }


    public double getPercentLQUnder(int p)
    {
        Vector<AutoResizeBitSet.IndexedBitSetContainer> v = getLQInfo();
        if (p >= v.size())
            return 100.0;

        return v.get(p).getCardinality() / (double)n;
    }

    public Vector<AutoResizeBitSet.IndexedBitSetContainer> getNInfo()
    {
        if(ns_info.size() != 0)
            return ns_info;

        Vector<BitSet> nbits = AutoResizeBitSet.getRelevantBitSet(n, AutoResizeBitSet.toBitSets(ns), 0.99, true);
        for (int i=0; i<nbits.size(); i++)
        {
            ns_info.add(new AutoResizeBitSet.IndexedBitSetContainer(nbits.get(i)));
        }
        return ns_info;
    }

    public Vector<AutoResizeBitSet.IndexedBitSetContainer> getLQInfo()
    {
        if(lq_info.size() != 0)
            return lq_info;

        Vector<BitSet> lqbits = AutoResizeBitSet.getRelevantBitSet(n, AutoResizeBitSet.toBitSets(lq), 0.99, true);
        for (int i=0; i<lqbits.size(); i++)
        {
            lq_info.add(new AutoResizeBitSet.IndexedBitSetContainer(lqbits.get(i)));
        }
        return lq_info;
    }


    public void writeXML(File f, File bitsetfile)
    {
        int bsi_index = 0;
        try
        {
            BitSetIndex bsi = null;
            if (bitsetfile !=null)
            {
                bsi = new BitSetIndex(bitsetfile);
            }
            //BitSetIndex bsi = new BitSetIndex(bitsetfile);

            XMLStreamWriter xsw = XMLUtils.getWriter(FileUtils.getWriter(f));
            xsw.writeStartDocument();
            xsw.writeStartElement("infos");
            xsw.writeAttribute("n", ""+n);
            xsw.writeAttribute("bs_indexed", ""+(bsi!=null));

            Vector<BitSet> nbits = ObjectGetter.map(getNInfo(), (AutoResizeBitSet.IndexedBitSetContainer ibc) -> ibc.getBitSet());

            //Vector<Integer> cumulative_up = AutoResizeBitSet.getAutoCumulativeCardinality(ns, true);
            //Vector<Integer> cumulative_down = AutoResizeBitSet.getAutoCumulativeCardinality(ns, false);


            /*BitSet full = AutoResizeBitSet.getFullBS(n);
            Vector<BitSet> v = AutoResizeBitSet.getCumulativeBitSets(AutoResizeBitSet.toBitSets(ns), true);

            v.insertElementAt(SetInfo.minus(full, v.get(v.size()-1)), 0);
            Vector<BitSet> v2 = v;

            ObjectGetter.apply(IteratorUtils.getIntIterator(1, v.size()), (Integer idx) -> v2.set(idx, SetInfo.union(v2.get(0), v2.get(idx))));

            int maxlevel = (int)(0.9999 * n);
            v = ObjectGetter.filter(v, (BitSet bs) -> bs.cardinality() <= maxlevel);
            */
            bsi_index = AutoResizeBitSet.writeBitSetData(xsw, "n", nbits, bsi,bsi_index, AutoResizeBitSet.DIRECT_BS_FLAG);

            Vector<BitSet> lqbits = ObjectGetter.map(getLQInfo(), (AutoResizeBitSet.IndexedBitSetContainer ibc) -> ibc.getBitSet());
            //bsi_index = AutoResizeBitSet.writeBitSetData(xsw, "lq", AutoResizeBitSet.toBitSets(lq),bsi,bsi_index, AutoResizeBitSet.CUM_UP_BS_FLAG);
            bsi_index = AutoResizeBitSet.writeBitSetData(xsw, "lq", lqbits, bsi, bsi_index, AutoResizeBitSet.DIRECT_BS_FLAG);

            xsw.writeEndElement();
            xsw.writeEndDocument();
            xsw.close();
        }
        catch(XMLStreamException xse)
        {
            throw new FRuntimeException("xml error: %s while writing %s", xse, xse.getMessage(), f);
        }
    }


    public static void main (String args[])
    {
        SimpleOptionParser cmd = new SimpleOptionParser("i", "convert");
        cmd.setFile("i");
        cmd.setOptional("convert");
        Locale.setDefault(Locale.UK);

        if (!OptionParser.parseParams(args, true, false, true, cmd))
            return;


        File xml = cmd.getFile("i");
        File bs_index = new File(xml.getAbsolutePath()+".bs_index");

        BasicReadStatistics brs = new BasicReadStatistics(xml, (bs_index.exists()) ? bs_index : null);

        if (cmd.isOptionSet("convert"))
        {
            File f = cmd.getFile("convert");
            File f_i = new File(f.getAbsolutePath()+".bs_index");
            brs.writeXML(f, f_i);
        }


        brs.printInfo();
    }
}


