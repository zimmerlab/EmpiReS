package empires.rnaseq.simulation;

import lmu.utils.*;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.io.File;
import java.util.*;
import java.util.function.Consumer;

import static lmu.utils.IteratorUtils.rangev;
import static lmu.utils.ObjectGetter.*;

public class PositionBiasFactory {

    Logger log = LogConfig.getLogger();
    static final int MIN_READS_TO_USE_BIAS = 2000;

    NormalDistribution fragmentLengthDistrib;
    int readLength;
    boolean paired;

    TreeMap<Integer, int[]> length2biases = new TreeMap<>();
    HashMap<Integer, ScalableTranscriptInfo> transcriptLength2Model = null;
    HashMap<Integer, int[]> length2cumulative = new HashMap<>();

    int maxLength;

    public PositionBiasFactory() {
        this((HashMap<String, int[]>)null, new NormalDistribution(200, 60), 75, true);
    }

    public PositionBiasFactory(HashMap<String, int[]> transcriptBiases, NormalDistribution fragmentLengthDistrib, int readLength, boolean paired) {
        this.fragmentLengthDistrib = fragmentLengthDistrib;
        this.readLength = readLength;
        this.paired = paired;
        if(transcriptBiases == null)
            return;

        maxLength = (int)fragmentLengthDistrib.inverseCumulativeProbability(0.99);
        log.info("max frlength: %d\n", maxLength);
        HashMap<Integer, int[]>  length2mostreads = new HashMap<>();
        HashMap<Integer, Integer> length2maxreads = new HashMap<>();

        for(Map.Entry<String, int[]> e : transcriptBiases.entrySet()) {
            int nreads = NumUtils.sum(e.getValue());
            if(nreads < MIN_READS_TO_USE_BIAS ||  length2maxreads.getOrDefault(e.getValue().length, 0) >= e.getValue().length)
                continue;

            length2maxreads.put(e.getValue().length, nreads);
            length2mostreads.put(e.getValue().length, e.getValue());

        }

        transcriptLength2Model = new HashMap<>();
        length2biases.putAll(length2mostreads);
    }

    public void initLength(Set<Integer> transcriptLengths, int NUM_THREADS) {
        if(transcriptLength2Model == null)
            return;

        Vector<Integer> lengths = shuffle(toVector(transcriptLengths), false);

        int maxL = length2biases.lastKey();
        log.info("max length: %d size: %d bigger: %d/%d\n", maxL, length2biases.size(), filteredSize(lengths, (_l) -> _l > maxL), lengths.size());
        long t1 = System.currentTimeMillis();
        ThreadUtils.runLambdas(lengths,
                (_l) -> {
                    Integer key = length2biases.ceilingKey(_l);
                    if(key == null) {
                        key = length2biases.floorKey(_l);
                    }

                    log.trace("will create a model for %d with %d\n", _l, key);
                    int[] pos2freq = length2biases.get(key);
                    int[] cumulative = getCumulative(key);
                    ScalableTranscriptInfo model = new ScalableTranscriptInfo(pos2freq, cumulative, _l, readLength, maxLength);
                    synchronized (transcriptLength2Model) {
                        transcriptLength2Model.put(_l, model);
                        if(transcriptLength2Model.size() % 100 == 0) {
                            log.info("%d/%d models are ready", transcriptLength2Model.size(), transcriptLengths.size());
                        }
                    }



                }, NUM_THREADS, 100, true, -1, null);

        long t2 = System.currentTimeMillis();
        log.info("initalizing biased transcript positions ready (%d models) on %d threads took: %.2f sec", transcriptLengths.size(), NUM_THREADS,  (t2 - t1) / 1000.0);
    }

    public static double getDotScoreToUniform(int[] starts) {
        Vector<Double> unit = NumUtils.toUnitVector(NumUtils.toDoubles(starts));
        Vector<Double> uniform = NumUtils.toUnitVector(mapIndex(starts.length, (_s) -> 1.0));
        return NumUtils.getDOTScore(uniform, unit);
    }

    public static HashMap<String, int[]> readTrBiases(File f) {
        if(f == null)
            return null;

        HashMap<String, int[]> tr2startbias = new HashMap<>();

        for(String[] sp : FileUtils.getFieldSetsIterable(f, "\t")) {
            int[] starts = new int[Integer.parseInt(sp[1])];
            String[] biases = sp[2].split(",");
            for(int i = 0; i<biases.length; i++) {
                starts[i] = Integer.parseInt(biases[i]);
            }
            tr2startbias.put(sp[0], starts);
        }

        return tr2startbias;
    }
    public PositionBiasFactory(File trBiases, NormalDistribution fragmentLengthDistrib, int readLength, boolean paired) {
        this(readTrBiases(trBiases), fragmentLengthDistrib, readLength, paired);
    }

    synchronized int[] getCumulative(int key) {
        int[] cumulative = length2cumulative.get(key);
        if(cumulative != null)
            return cumulative;

        int[] pos2freq = length2biases.get(key);
        cumulative = new int[pos2freq.length];
        cumulative[0] = pos2freq[0];
        for (int i = 1; i < pos2freq.length; i++) {
            cumulative[i] = cumulative[i - 1] + pos2freq[i];
        }
        length2cumulative.put(key, cumulative);
        return cumulative;
    }


    ScalableTranscriptInfo getBiasModel(int length) {
        if(transcriptLength2Model == null)
            return null;

        ScalableTranscriptInfo model = transcriptLength2Model.get(length);
        if(model != null)
            return model;

        Integer key = length2biases.ceilingKey(length);
        if(key == null) {
            key = length2biases.floorKey(length);
        }

        log.info("will create a model for %d with %d\n", length, key);
        int[] pos2freq = length2biases.get(key);
        int[] cumulative = getCumulative(key);
        transcriptLength2Model.put(length, model = new ScalableTranscriptInfo(pos2freq, cumulative, length, readLength, maxLength));

        log.info("added model for length: %d got %d lengths so far\n", length, transcriptLength2Model.size());
        return model;
    }

    public void simulate(String gene, boolean strand, String transcriptId, HashMap<String, RegionVector> isoforms, int numReads,
                         Consumer<SimulatedRead> consumer) {

        RegionVector trVec = isoforms.get(transcriptId);
        int trLength = trVec.getCoveredLength();

        empires.rnaseq.simulation.PositionBias positionBias = new PositionBias(fragmentLengthDistrib, trLength, getBiasModel(trLength), readLength);
        log.trace("will generate %d reads for %s (l: %d)", numReads, transcriptId, trLength);
        for(int i=0; i<numReads; i++) {
            consumer.accept(new SimulatedRead(gene, strand, transcriptId, trVec, readLength, isoforms, paired, positionBias));
        }
        log.trace("generating %d reads for %s (l: %d) READY!", numReads, transcriptId, trLength);
    }


}
