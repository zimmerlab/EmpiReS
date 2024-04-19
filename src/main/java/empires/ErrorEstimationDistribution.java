package nlEmpiRe;

import lmu.utils.*;
import lmu.utils.tuple.Tuple3;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.util.*;
import java.util.function.Function;

import static lmu.utils.IteratorUtils.rangev;
import static lmu.utils.ObjectGetter.*;
import org.apache.logging.log4j.Logger;

public class ErrorEstimationDistribution {

    Logger log = LogConfig.getLogger();
    public static double INT_FACTOR = 100;
    static double MINFC = 1.0 / INT_FACTOR;
    public static double norm = 1.0 / (INT_FACTOR + 0.0);
    static NormalDistribution N01 = new NormalDistribution(0.0, 1.0);

    public static void setFcWithBin(double fcWithBin) {
        MINFC = norm = fcWithBin;
        INT_FACTOR = (int)(1.0 / fcWithBin);
    }
    int min = 0;
    int max = 0;
    HashMap<Integer, Long> vals = null;
    long[] cumulative = null;
    double[] zscores = null;

    Double SD = null;
    Double variance  = null;

    static Double NORM_FACTOR =  (double)(1 << 30); // long max val 2^63  if we add up 5 distributions  with values 2^10 -> max value is 2^50 -> how to deal with it
    //if we have 2^30



    long[] normed;
    double max_z = 0.0;
    long TOTAL = 0;

    public static ErrorEstimationDistribution fromNormalDistrib(double mean, double sd) {
        return fromNormalDistrib(new NormalDistribution(mean, sd));
    }

    static Tuple3<Integer, Integer, long[]> cutToPercent(Tuple3<Integer, Integer, long[]> t, double percent) {
        long sumL = NumUtils.sum(t.get2());
        long throwaway = (long)(percent * sumL * 0.5);
        int start = 0;
        long pre = 0;

        while(pre + t.get2()[start] < throwaway) {
            pre += t.get2()[start];
            start++;
        }

        long post = 0;
        int end = t.get2().length - 1;
        //System.out.printf("cut length: %d throw away %d/%d start: %d\n", t.get2().length, throwaway, sumL, start);
        while(post + t.get2()[end] < throwaway) {
            post += t.get2()[end];
            end--;
        }
        long[] shrinked = Arrays.copyOfRange (t.get2(), start, end + 1);
        shrinked[0] += pre;
        shrinked[shrinked.length - 1] += post;
        //System.out.printf("sum before: %d after: %d\n", NumUtils.sum(t.get2()), NumUtils.sum(shrinked));
        /*System.out.printf("cut %d-%d to: %d-%d sum: %d throw away: %d pre: %d: %d post: %d: %d\n", t.get0(), t.get1(),
                t.get0() + start, t.get0() + end,
                sumL, throwaway, start, pre, end, post
                );

         */
        return Tuple3.create(t.get0() + start, t.get0() + start + shrinked.length, shrinked);
    }

    static Tuple3<Integer, Integer, long[]> sumUpTmp(Tuple3<Integer, Integer, long[]> d1, Tuple3<Integer, Integer, long[]> d2) {
        int minsum = Math.min(d1.get0(), d2.get0());
        int maxsum = Math.max(d1.get1(), d2.get1());

        Logger log = LogConfig.getLogger();

        long[] summed = new long[maxsum - minsum + 1];

        for(Tuple3<Integer, Integer, long[]> d : toVector(d1, d2)) {
            long[] raw = d.get2();
            for(int i=0; i< raw.length; i++) {
                summed[i + d.get0() - minsum] += raw[i];
            }
        }
        return Tuple3.create(minsum, maxsum, summed);
    }

    static Tuple3<Integer, Integer, long[]> sumUp(Tuple3<Integer, Integer, long[]> d1, Tuple3<Integer, Integer, long[]> d2) {
        return sumUp(d1, d2, 1);
    }

    static ErrorEstimationDistribution fromTuple(Tuple3<Integer, Integer, long[]> t) {
        return new ErrorEstimationDistribution(t.get0(), t.get1(), toCumulative(t.get2()), true);
    }

    static String toQuantiles(Tuple3<Integer, Integer, long[]> t) {
        return "" + new SparseCumulativeDistribution(fromTuple(t)).quantiles;
    }

    static Tuple3<Integer, Integer, long[]> sumUp(Tuple3<Integer, Integer, long[]> d1, Tuple3<Integer, Integer, long[]> d2, int factorD1) {

        Logger log = LogConfig.getLogger();
        long[] n1 = d1.get2();
        long[] n2 = d2.get2();
        int minsum = d1.get0() + d2.get0();
        int maxsum = d1.get1() + d2.get1();


        if(factorD1 > 1) {
            double n2factor = 1.0 / (factorD1);
            for(int i=0; i<n2.length; i++) {
                n2[i] *= n2factor;
            }
        }

        long[] summed = new long[maxsum - minsum + 1];

        for(int i=0; i < n1.length; i++) {
            for(int j=0; j < n2.length; j++) {
                int sumIdx = (d1.get0()) + d2.get0() + i + j - minsum;
                summed[sumIdx] += n1[i] * n2[j];
            }
        }
        long[] normed = getNormedFreqs(toCumulative(summed));
        Tuple3<Integer, Integer, long[]> combined = Tuple3.create(minsum, minsum  + normed.length, normed);
        return (combined.get2().length < 1000) ? combined : cutToPercent(combined, 0.001);
    }

    public static class ErrorEstimationCombiner{
        Logger log = LogConfig.getLogger();
        int numDistribs = 0;
        Tuple3<Integer,Integer, long[]> last;
        ErrorEstimationDistribution rv = null;
        public boolean closed = false;
        ErrorEstimationDistribution normed = null;
        Vector<ErrorEstimationDistribution> inputs = new Vector<>();

        ErrorEstimationCombiner() {

        }


        public void add(ErrorEstimationDistribution errorEstimationDistribution) {


            if(closed)
                throw new FRuntimeException("invalid state: closed. You cannot add further distributions after invoking get()");

            inputs.add(errorEstimationDistribution);
            Tuple3<Integer,Integer, long[]> normed = Tuple3.create(errorEstimationDistribution.min, errorEstimationDistribution.max, getNormedFreqs(errorEstimationDistribution.getCumulative()));
            //log.info("normed of %s max: %s", errorEstimationDistribution, NumUtils.max(normed.get2()));
            if(numDistribs == 0) {
                numDistribs = 1;
                last = normed;
                rv = errorEstimationDistribution;
                return;
            }

            numDistribs++;
            rv = null;
            last = sumUp(last, normed, numDistribs  - 1);
        }

        static Tuple3<Integer,Integer, long[]> toTuple(ErrorEstimationDistribution d) {
            return Tuple3.create(d.min, d.max, getNormedFreqs(d.getCumulative()));
        }
        public ErrorEstimationDistribution getNormed() {
            if(normed != null)
                return normed;

            if(numDistribs < 2)
                return inputs.get(0);


            int min = last.get0();
            long[] normed = last.get2();
            HashMap<Integer, Integer> scaled = new HashMap<>();
            double norm = 1.0 / (0.0 + numDistribs);

            for(int i=0; i<normed.length; i++) {
                int fc = (int)((min + i) * norm);

                MapBuilder.update(scaled, fc, (int)normed[i]);

            }

            ErrorEstimationDistribution errNormed = new ErrorEstimationDistribution(scaled);
            long[] cumulative = errNormed.getCumulative(true);
            if(cumulative.length > 1000) {
                Tuple3<Integer, Integer, long[]> cutted = cutToPercent(Tuple3.create(errNormed.min, errNormed.max, getNormedFreqs(cumulative)), 0.001);
                errNormed = new ErrorEstimationDistribution(cutted.get0(), cutted.get1(), cutted.get2(), true);

            }
            return errNormed;

        }

        public ErrorEstimationDistribution get() {

            if(rv != null)
                return rv;

            if(last == null)
                return null;

            long[] cumulative = last.get2();
            for(int i=1; i<cumulative.length; i++) {
                cumulative[i] += cumulative[i-1];
            }
            return rv = new ErrorEstimationDistribution(last.get0(), last.get1(), toCumulative(last.get2()), true);
        }

    }

    private long[] getNormed() {
        if(normed != null)
            return normed;

        return getNormedFreqs(getCumulative());
    }
    public ErrorEstimationDistribution scale(double scale) {
        return scale(0, scale);
    }

    /** scale distribution around the mean -> the SD will be get bigger / get shrinked
     *
     * if mean is 0 , the complete distribution gets normalized
     * */
    public ErrorEstimationDistribution scale(double mean, double scale) {
        long[] normed = getNormed();

        HashMap<Integer, Long> scaled = new HashMap<>();

        for(int i=0; i<normed.length; i++) {
            double fc = (min + i) * norm;

            double corrected = mean + ((fc - mean) * scale);

            MapBuilder.update(scaled, (int)(corrected * INT_FACTOR), normed[i]);

        }
        ErrorEstimationDistribution errd = new ErrorEstimationDistribution(scaled, true);
        errd.getCumulative(true);
        return errd;
    }

    public static ErrorEstimationDistribution combineAndAverage(Vector<ErrorEstimationDistribution> distributions) {
        return combine(distributions, true);
    }



    static ErrorEstimationDistribution combine(Vector<ErrorEstimationDistribution> distributions, boolean do_average) {
        if(distributions.size() == 1)
            return distributions.get(0);

        ErrorEstimationCombiner combiner = new ErrorEstimationCombiner();
        for(ErrorEstimationDistribution d : distributions) {
            combiner.add(d);
        }


        return (do_average) ? combiner.getNormed() : combiner.get();
    }

    public static ErrorEstimationDistribution combine(Vector<ErrorEstimationDistribution> distributions) {
        return combine(distributions, false);
    }

    /** for testing purposes */
    public static ErrorEstimationDistribution fromNormalDistrib(NormalDistribution nd) {
        double alpha = 0.0001;
        double min = nd.inverseCumulativeProbability(alpha);
        double max = nd.inverseCumulativeProbability(1.0 - alpha);


        HashMap<Integer, Long> intFCMap = new HashMap<>();
        for(double fc = min; fc < max; fc += norm) {

            intFCMap.put((int)(fc * INT_FACTOR), (long)(NORM_FACTOR *  (nd.cumulativeProbability(fc) - nd.cumulativeProbability(fc - norm))));
        }

        ErrorEstimationDistribution errd = new ErrorEstimationDistribution(intFCMap, true);
        errd.getCumulative(true);
        errd.getSD(nd.getMean());

        return errd;

    }

    public ErrorEstimationDistribution(HashMap<Integer, Long> vals, boolean fromLong) {
        this.vals  = vals;
        min = NumUtils.min(this.vals.keySet());
        max = NumUtils.max(this.vals.keySet());
    }

    public ErrorEstimationDistribution(HashMap<Integer, Integer> vals) {
        this.vals  = new HashMap<>();
        apply(vals.entrySet(), (_e) -> this.vals.put(_e.getKey(), _e.getValue() + 0l));
        min = NumUtils.min(this.vals.keySet());
        max = NumUtils.max(this.vals.keySet());
    }

    static long[] toCumulative(long[] simple) {
        long[] rv = new long[simple.length];
        rv[0] = simple[0];
        for(int i=1; i<simple.length; i++) {
            rv[i] = rv[i-1] + simple[i];
        }
        return rv;
    }

    static long[] getRawFreqs(long[] cumulative) {

        long[] rv = new long[cumulative.length];
        for(int i=0; i< cumulative.length; i++) {
            rv[i] = ( cumulative[i] - ((i == 0) ? 0 : cumulative[i - 1]));
        }
        return rv;

    }

    static long[] getNormedFreqs(long[] vals) {
        double NORM = NORM_FACTOR / vals[vals.length - 1];
        long[] rv = getRawFreqs(vals);
        for(int i=0; i< vals.length; i++) {
            rv[i] *= NORM;
        }
        return rv;
    }


    public static  Pair<ErrorEstimationDistribution, Integer> read(byte[] buffer, int offset) {

        int min = ByteBuffer.wrap(buffer, offset, 4).getInt();
        offset += 4;
        int max = ByteBuffer.wrap(buffer, offset, 4).getInt();
        offset += 4;
        int length = ByteBuffer.wrap(buffer, offset, 4).getInt();
        offset += 4;
        long[] cumulative = new long[length];
        for(int i=0; i<length; i++) {
            cumulative[i] = ByteBuffer.wrap(buffer, offset, 8).getLong();
            offset += 8;
        }

        ErrorEstimationDistribution errd = new ErrorEstimationDistribution(min, max);
        errd.cumulative = cumulative;

        return Pair.create(errd, offset);
    }
    public static  ErrorEstimationDistribution read(DataInputStream dis) {

        try{
           int min = dis.readInt();
           int max = dis.readInt();
            int length = dis.readInt();
            long[] cumulative = new long[length];
            for(int i=0; i<length; i++) {
                cumulative[i] = dis.readLong();
            }
            ErrorEstimationDistribution errd = new ErrorEstimationDistribution(min, max);
            errd.cumulative = cumulative;
            return errd;

        }catch (IOException ie) {
            throw new FRuntimeException("i/o error: %s while reading errdistrib", ie, ie.getMessage());
        }
    }

    public void writeToBuffer(DataOutputStream dataOutputStream) {
        try{
            dataOutputStream.writeInt(min);
            dataOutputStream.writeInt(max);

            getCumulative();
            dataOutputStream.writeInt(cumulative.length);
            for(int i=0; i<cumulative.length; i++) {
                dataOutputStream.writeLong(cumulative[i]);
            }
        } catch(IOException ie) {
            throw new FRuntimeException("i/o error: %s while writing distrib: %s", ie, ie.getMessage(), this);
        }

    }
    private ErrorEstimationDistribution(int min, int max) {
        this.min = min;
        this.max = max;
        if(Math.max(Math.abs(max), Math.abs(min)) > INT_FACTOR * 500 ) {

            log.error("got error distribution %.2f - %.2f -> this is highly unlikely, check your data - maybe you specified log data, but it is not logged??\nstack:\n%s\n",
                    getMinFC(), getMaxFC(),
                    ExceptionManager.getCallerStack());
        }
    }
    private ErrorEstimationDistribution(int min, int max, long[] normed, boolean noZ) {
        this(min, max);
        cumulative = normed;
        if(noZ) {
            getCumulative(true);
        } else {
            initZ();
        }

    }
    private ErrorEstimationDistribution(int min, int max, long[] normed) {
        this(min, max, normed, false);
    }

    private ErrorEstimationDistribution(ErrorEstimationDistribution base) {
        cumulative = Arrays.copyOf(base.cumulative, base.cumulative.length);
        zscores = (zscores == null) ? null : Arrays.copyOf(base.zscores, base.zscores.length);
        SD = base.SD;
        variance = base.variance;
        max_z = base.max_z;
        TOTAL = base.TOTAL;
        min = base.min;
        max = base.max;
    }

    public ErrorEstimationDistribution shift(double shiftFC) {
        ErrorEstimationDistribution shifted = new ErrorEstimationDistribution(this);
        int SHIFT_INT = (int)(INT_FACTOR * shiftFC);
        shifted.min = min + SHIFT_INT;
        shifted.max = max + SHIFT_INT;
        return shifted;
    }


    public ErrorEstimationDistribution(ErrorEstimationDistribution distrib, double shift) {
        int ishift = (int)(INT_FACTOR * shift);
        min = distrib.min + ishift;
        max = distrib.max + ishift;
        cumulative = distrib.cumulative;
    }

    public ErrorEstimationDistribution(Vector<ErrorEstimationDistribution> distribs, Vector<Double> shifts) {
        this(distribs, shifts, true);
    }

    public ErrorEstimationDistribution(Vector<ErrorEstimationDistribution> distribs, Vector<Double> shifts, boolean norm_each) {

        Vector<Integer> range = rangev(distribs);
        Vector<Integer> intShifts = map(shifts, (_d) -> (int)(INT_FACTOR * _d));
        min = NumUtils.min(range, (_i) -> distribs.get(_i).min + intShifts.get(_i));
        max = NumUtils.max(range, (_i) -> distribs.get(_i).max + intShifts.get(_i));

        min = Math.min(min, -1);
        max = Math.max(max, 1);

        cumulative = new long[max - min + 1];

        for(int di=0; di<distribs.size(); di++) {
            ErrorEstimationDistribution fcInfo = distribs.get(di);
            int shiftFC = intShifts.get(di);
            long[] normed_cumulative = (norm_each) ? getNormedFreqs(fcInfo.cumulative) : getRawFreqs(fcInfo.cumulative);
            for(int i=0, j = fcInfo.min + shiftFC - min; i<normed_cumulative.length; i++, j++) {
                cumulative[j] += normed_cumulative[i];
            }
        }
        initZ();
    }


    public static ErrorEstimationDistribution substract(ErrorEstimationDistribution base, ErrorEstimationDistribution tosubstract) {
        return substract(base, tosubstract, true);
    }

    public static ErrorEstimationDistribution substract(ErrorEstimationDistribution base, ErrorEstimationDistribution tosubstract, boolean zeroCentered) {
        return new ErrorEstimationDistribution(base, tosubstract, zeroCentered);
    }


    private ErrorEstimationDistribution(ErrorEstimationDistribution base, ErrorEstimationDistribution tosubstract, boolean zeroCentered) {
        min = base.min - tosubstract.max;
        max = base.max - tosubstract.min;
        long[] n1 = getNormedFreqs(base.cumulative);
        long[] n2 = getNormedFreqs(tosubstract.cumulative);

        int base_min = base.min;
        int tosubstract_min = tosubstract.min;

        if(!zeroCentered) {
            Tuple3<Integer, Integer, long[]> base_cutted =  cutToPercent(Tuple3.create(base.min, base.max, n1), 0.01);
            Tuple3<Integer, Integer, long[]> tosubtract_cutted =  cutToPercent(Tuple3.create(tosubstract.min, tosubstract.max, n2), 0.01);
            min = base_cutted.get0() - tosubtract_cutted.get1();
            max = base_cutted.get1() - tosubtract_cutted.get0();
            base_min = base_cutted.get0();
            tosubstract_min = tosubtract_cutted.get0();
            n1 = base_cutted.get2();
            n2 = tosubtract_cutted.get2();
        }

        long[] joined = new long[max - min + 1];

        for(int i=0; i< n1.length; i++) {
            int fc1 = base_min + i;
            long nfreq1 = n1[i];
            for(int j = 0; j < n2.length; j++) {
                int fc2 = tosubstract_min + j;
                long nfreq2 = n2[j];
                int fcdiff = fc1 - fc2;
                int offset = fcdiff - min;
                joined[offset] += (nfreq1 * nfreq2);
            }
        }

        if(zeroCentered) {
            this.cumulative = joined;
            initZ();
        } else {
            cumulative = toCumulative(joined);
        }


    }

    public static ErrorEstimationDistribution getDiffInfo(ErrorEstimationDistribution base, ErrorEstimationDistribution tosubstract) {
        return new ErrorEstimationDistribution(base, tosubstract, true);
    }

    public ErrorEstimationDistribution(Vector<Vector<Double>> replicate_signal_groups)
    {
        this(0, replicate_signal_groups.size(), (_i) -> replicate_signal_groups.get(_i));
    }


    public ErrorEstimationDistribution(int start, int end, Function<Integer, Vector<Double>> dataGetter) {
        this(sampleWithinReplicateLogSignalDifferences(start, end, dataGetter));
    }

    public static HashMap<Integer, Integer> sampleWithinReplicateLogSignalDifferences(int start, int end, Function<Integer, Vector<Double>> dataGetter) {
        HashMap<Integer, Integer> fcs = new HashMap<>();
        Vector<Double> logsigvals = new Vector<>();
        Logger log = LogConfig.getLogger();
        for(int i=start; i<end; i++) {

            Vector<Double> data = dataGetter.apply(i);

            if(data.size() < 2)
                continue;

            int selidx = (int)(Math.random() * data.size());
            double factor = data.get(selidx);

            logsigvals.addAll(filter_and_map(rangev(data), (_i) -> _i != selidx, (_i) -> data.get(_i) - factor) );
        }


        Collections.shuffle(logsigvals);
        fcs.clear();

        for(int i=1; i<logsigvals.size(); i++) {
            MapBuilder.update(fcs, (int) (INT_FACTOR * (0.5 * (logsigvals.get(i) - logsigvals.get(i-1)))));
        }

        int minFC = NumUtils.min(fcs.keySet());
        int maxFC = NumUtils.max(fcs.keySet());

        if(maxFC - minFC > (50 * INT_FACTOR)) {
            log.warn("STRANGE FC ESTIMATION got %.2f - %.2f logsigvals: %s\n", minFC * norm, maxFC * norm, NumUtils.getNumInfo(logsigvals).getInfoWithQ());
        }

        return fcs;
    }

    static final double SMALLEST_PVAL = Math.pow(10, -9);
    void initZ()
    {

        for(int i=1; i<cumulative.length; i++) {
            cumulative[i] += cumulative[i-1];
        }
        TOTAL = cumulative[cumulative.length - 1];

        double minP = 1.0  / ( TOTAL + 1.0);
        max_z =  Math.abs(N01.inverseCumulativeProbability(Math.max(SMALLEST_PVAL, minP)));

        zscores = new double[cumulative.length];
        int zeroPos = - min;

        if(zeroPos >= cumulative.length || zeroPos <= 0)
            throw new FRuntimeException("error getting zero pos min: %d max: %d diff: %d cumlength: %d total: %d", min, max, max - min, cumulative.length, TOTAL);

        double totalnorm = 1.0 / ((TOTAL - cumulative[zeroPos]) + 1.0);
        double zeronorm = 1.0 / (cumulative[zeroPos - 1] + 1.0);

        for(int i=0; i<cumulative.length; i++) {

            if(i == zeroPos) //logfc 0
            {
                zscores[i] = 0.0;
                continue;
            }
            long num_more_extreme = (i < zeroPos) ? ((i == 0) ? 0 : cumulative[i])
                                        : ((i == cumulative.length - 1) ? 0 : cumulative[cumulative.length - 1] - cumulative[i + 1]);

            double norm = (i < zeroPos) ? zeronorm : totalnorm;

            //long num_more_extreme = Math.max(pre, post);
            double p = 0.5 * Math.max(SMALLEST_PVAL, (num_more_extreme + 1.0) * norm);
            zscores[i] = ((i < zeroPos) ? -1.0 : 1.0) * Math.abs(N01.inverseCumulativeProbability(p));

        }

    }
    public long[] getCumulative() {
        return getCumulative(false);
    }

    public long[] getCumulative(boolean noZinit)
    {
        if(cumulative != null)
            return  cumulative;


        cumulative = new long[max - min + 1];

        //log.info("will create distrib from %d-%d num keys: %d %s\n", min, max, vals.size(), vals);


        for(Map.Entry<Integer, Long> e : vals.entrySet()) {
            cumulative[e.getKey() - min] += e.getValue();
        }
        if(!noZinit) {
            initZ();
        } else {
            for(int i=1; i<cumulative.length; i++) {
                cumulative[i] += cumulative[i-1];
            }
            TOTAL = cumulative[cumulative.length - 1];
        }
        return cumulative;
    }

    public double getMinFC() {
        return min * norm;
    }

    public double getMaxFC() {
        return max * norm;
    }

    public String toString() {
        return String.format("err: %d-%d (width: %d %.2f-%.2f) var: %.2f P: %d", min, max, max - min,  getMinFC(), getMaxFC()   , getVariance(), cumulative[cumulative.length-1]);
    }

    public double getSD() {

        return getSD(0.0);
    }

    public int getCumulativeWidth() {
        return cumulative.length;
    }

    public double getSD(double mean)
    {
        if(SD != null && Math.abs(mean) < 0.0001)
            return SD;
        getCumulative();

        double norm = 1.0 / INT_FACTOR;
        double sumsquaredError = 0.0;
        long last = 0;
        for(int i=0; i<cumulative.length; i++) {
            double fc = (i + min) * norm;

            sumsquaredError += (cumulative[i] - last) * Math.pow(fc - mean, 2.0);
            last = cumulative[i];
        }
        long total = cumulative[cumulative.length - 1];
        variance = (sumsquaredError / (total));

        return SD = Math.sqrt(variance);
    }

    public double getVariance() {
        getSD();
        return variance;
    }

    public double getCumulativeFrequencyToFoldChange(double fc) {
        int k = (int)((fc * INT_FACTOR) - min );
        //System.out.printf("cumfreq to fc in <%.2f, %.2f sd: %.2f> fc: %.2f k = %d/%d total: %d\n", min * norm, max * norm, SD, fc, k, cumulative.length, cumulative[cumulative.length-1]);
        double total_norm = 1.0 / (0.0 + cumulative[cumulative.length - 1]);
        if( k <= 0)
            return 1.0 * total_norm;

        if( k >= cumulative.length -1)
            return 1.0 - total_norm;

        return (cumulative[k]  * total_norm);
    }


    public class Peak  {
        public int peakStart;
        public int peakEnd;
        public int summit;
        long total = 0;

        Peak(long[] noncumulative, int start, long[] cumulative) {

            peakStart = start;
            summit = peakStart;


            while(summit + 1< noncumulative.length && noncumulative[summit + 1] > noncumulative[summit]) {
                summit++;
            }
            peakEnd = summit;
            while(peakEnd + 1 < noncumulative.length && noncumulative[peakEnd + 1] < noncumulative[peakEnd]) {
                peakEnd++;
            }
            total = cumulative[Math.min(cumulative.length - 1, peakEnd+1)] - cumulative[Math.max(0, peakStart - 1)];
        }

        public String toString() {
            return String.format("p:%.2f-%.2f-%.2f (%.4f%%)", (peakStart + min) * norm, (min + summit) * norm , (min + peakEnd) * norm, (100.0 * total) / TOTAL);
        }

    }

    public Vector<UPair<Double>> getCumulativeLine() {
        Vector<UPair<Double>> lineData = new Vector<>();
        for(int i=0; i<cumulative.length; i++) {
            double val = cumulative[i];
            double fc = (min + i) * norm;
            lineData.add(UPair.createU(fc, val));
        }
        return lineData;
    }
    public Vector<UPair<Double>> getNonCumulativeLine() {
        Vector<UPair<Double>> lineData = new Vector<>();
        for(int i=1; i<cumulative.length; i++) {
            double val = cumulative[i] - cumulative[i-1];
            double fc = (min + i) * norm;
            lineData.add(UPair.createU(fc, val));
        }
        return lineData;
    }

    public Peak getBestFCPeak() {
        long noncumulative[] = new long[cumulative.length];

        for(int i=1; i<noncumulative.length; i++) {
            noncumulative[i] = cumulative[i] - cumulative[i-1];
        }

        Peak bestPeak = null;
        int pos = 0;
        while(pos < cumulative.length) {
            Peak p = new Peak(noncumulative, pos, cumulative);

            pos = p.peakEnd+1;
            if(bestPeak != null && bestPeak.total >= p.total)
                continue;

            bestPeak = p;
        }

        return bestPeak;
    }
    public double getBestFCShift() {

        Peak bestPeak = getBestFCPeak();
        return (bestPeak.summit + min) * norm;
    }

    static double DEFAULT_FCWITH_FOR_PROBABLE_FC_WINDOW = 0.1;

    /** returns the fc where DEFAULT_FCWITH_FOR_PROBABLE_FC_WINDOW * full width region has the most probability mass */
    public double getMostProbableFcWindowCenter() {
        return getMostProbableFcWindowCenter(DEFAULT_FCWITH_FOR_PROBABLE_FC_WINDOW * norm * (max - min));
    }

    public double getMostProbableFcWindowCenter(double fcwidth ) {
        getCumulative();
        int step = Math.max(1, (int)(fcwidth * INT_FACTOR * 0.5));
        if(step >= cumulative.length ) {
            return norm * (min + cumulative.length >> 1);
        }
        double maxp =  0;
        double maxpPosition = -1;
        for(int i = step; i < cumulative.length - step  - 1; i++) {
            double preCum = cumulative[i - step];
            double cum = cumulative[i + step];
            if(cum - preCum <= maxp)
                continue;

            maxp = cum - preCum;
            maxpPosition = i;
        }
        return (min + maxpPosition) * norm;
    }


    public double getFoldChangeToCumulativeFrequency(double cumulativeFrequency) {

        getCumulative();
        if(cumulativeFrequency <= 0.0)
            return min * norm ;


        if(cumulativeFrequency >= 1.0)
            return max * norm;

        long val = (long)(cumulativeFrequency * cumulative[cumulative.length - 1]);


        int hitpos = Arrays.binarySearch(cumulative, val);
        if(hitpos > 0) {
            return (min + hitpos) / (INT_FACTOR + 0.0);
        }
        hitpos = Math.max(0, - hitpos - 1); //binary search returns - insertionpoint -1  =
        double halfstep_shift_factor = (hitpos == 0 || hitpos == cumulative.length -1) ? 0 : 0.5;

        return (min + hitpos + halfstep_shift_factor) * norm;

    }

    public double getConvertedZscore(double obs_fc) {
        getCumulative(); //ensure the distribution is already calculated

        if(Math.abs(obs_fc) <= MINFC)
            return 0.0;

        int k = (int)(obs_fc * INT_FACTOR);

        int rank = k - min;
        if(rank < 0)
            return - max_z;
        if(rank >= cumulative.length)
            return max_z;

        return zscores[rank];
    }

    private long getSumOccAroundZero(int radius) {
        if(radius == 0)
            return 0;

        int lower = Math.max(min, -radius);
        int upper = Math.min(max, radius);

        int startIdx = Math.min(cumulative.length - 1, Math.max(0, lower - min));
        int endIdx = Math.max(0, Math.min(cumulative.length -1 , upper - min));

        return cumulative[endIdx] - cumulative[startIdx];
    }

    public double getCenteredFCWithTargetProbabilityMass(double targetProbabilityMass) {
        assert(targetProbabilityMass > 0 && targetProbabilityMass < 1.0);

        long targetOcc = (long)(TOTAL * targetProbabilityMass);
        int maxWidth = Math.min(Math.abs(max), Math.abs(min));
        int minWidth = 0;
        int middle = 0;
        while(minWidth <= maxWidth) {
            middle = minWidth + ((maxWidth - minWidth) >> 1);

            long mass = getSumOccAroundZero(middle);
            long diff = Math.abs(mass - targetOcc);
            long preDiff = Math.abs(targetOcc - getSumOccAroundZero(middle - 1));
            long postDiff = Math.abs(targetOcc - getSumOccAroundZero(middle + 1));

            if(preDiff >= diff && postDiff >= diff)
                break;


            if(mass < targetOcc) {
                minWidth = middle + 1;
            } else {
                maxWidth = middle - 1;
            }
        }

        return middle * norm;
    }
}
