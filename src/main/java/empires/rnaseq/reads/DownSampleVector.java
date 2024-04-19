package nlEmpiRe.rnaseq.reads;

import lmu.utils.ObjectGetter;

import java.util.Vector;

public class DownSampleVector
{
    public double[] downsampled;
    public double[] reads;
    public double[] frags;
    public double[] rv_dn;

    public DownSampleVector(int SN)
    {
        downsampled = new double[SN];
        reads = new double[SN];
        frags = new double[SN];
        rv_dn = new double[SN];
    }
    void reset()
    {
        for(int i=0; i<downsampled.length; downsampled[i++]=0);
        for(int i=0; i<reads.length; reads[i++]=0);
        for(int i=0; i<frags.length; frags[i++]=0);
    }

    public void update(DownSampleVector dsv) {
        for(int i=0; i<downsampled.length; i++) {
            downsampled[i] += dsv.downsampled[i];
        }
        for(int i=0; i<reads.length; i++) {
            reads[i] += dsv.reads[i];
        }
        for(int i=0; i<frags.length; i++) {
            frags[i] += dsv.frags[i];
        }
    }


    public static DownSampleVector update(DownSampleVector dsv, Vector<Integer> values, ObjectGetter.MapGetter<Integer, Double> downsampler)
    {

        double sum = 0;
        double maxval = 0.0;
        final int NS = values.size();
        for (int i=0; i< NS; i++)
        {

            int count = values.get(i);
            dsv.frags[i] += (count > 0 ) ? 1 : 0;
            dsv.reads[i] += count;


            sum += (dsv.rv_dn[i] = downsampler.get(count));
            maxval = Math.max(maxval, dsv.rv_dn[i]);
        }

        if (sum < 0.0001)
            return dsv;

        for (int i = 0; i < NS; dsv.downsampled[i] += (dsv.rv_dn[i++] / maxval)) ;

        return dsv;
    }
}


