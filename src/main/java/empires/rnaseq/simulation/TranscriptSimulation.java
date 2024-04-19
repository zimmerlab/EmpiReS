package nlEmpiRe.rnaseq.simulation;

import lmu.utils.*;
import nlEmpiRe.rnaseq.IsoformRegionGetter;
import nlEmpiRe.rnaseq.MultiIsoformRegion;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.Vector;
import java.util.function.Consumer;

import static lmu.utils.ObjectGetter.*;

public class TranscriptSimulation
{

    MultiIsoformRegion gene;
    HashMap<String, Integer> isoform2count;


    public TranscriptSimulation(MultiIsoformRegion gene, HashMap<String, Integer> isoform2count) {
        this.gene = gene;
        this.isoform2count = isoform2count;

    }

    public String toString() {
        return ""+ isoform2count;
    }

    public HashMap<Tuple, Double> simulateTrSetCounts(PositionBiasFactory positionBiasFactory) {
        return simulateTrSetCounts(positionBiasFactory, null);
    }


    public HashMap<Tuple, Double> simulateTrSetCounts(PositionBiasFactory positionBiasFactory, Consumer<nlEmpiRe.rnaseq.simulation.SimulatedRead> externConsumer) {

        HashMap<Tuple, Double> rv = new HashMap<>();
        simulate(positionBiasFactory,
                (_sr) ->
                {
                    if(externConsumer != null) {
                        externConsumer.accept(_sr);
                    }
                    MapBuilder.update(rv, Tuple.tupleFromCollection(toSortedVector(_sr.mapsToTranscripts, true)), 1.0);
                });
        return rv;
    }

    public Set<String> getSimulatedTranscriptIds() {
        return isoform2count.keySet();
    }

    public int getCountToSimulate(String transcriptId) {
        return isoform2count.get(transcriptId);
    }

    public void simulate(PositionBiasFactory positionBiasFactory, Consumer<nlEmpiRe.rnaseq.simulation.SimulatedRead> readConsumer) {
        for(String transcriptId : isoform2count.keySet()) {
            int count = isoform2count.get(transcriptId);
            if(count == 0)
                continue;
            positionBiasFactory.simulate(gene.id, gene.strand, transcriptId, gene.isoforms, count, readConsumer);
        }
    }


    public static HashMap<MultiIsoformRegion, TranscriptSimulation>  groupTranscriptsPerGene(Iterable<Pair<String, Integer>> transcript2count, IsoformRegionGetter annot) {
        HashMap<String, String> transcript2Gene = new HashMap<>();
        apply(annot.getRegions(), (_r) -> apply(_r.isoforms.keySet(), (_isoform) -> transcript2Gene.put(_isoform, _r.id)));
        Vector<String> missing = map_and_filter(transcript2count.iterator(), (_p) -> _p.getFirst(), (_t) -> transcript2Gene.get(_t) == null);
        if (missing.size() > 0) {
            throw new RuntimeException(String.format("unmappable transcripts: %s", missing));
        }
        HashMap<String, HashMap<String, Integer>> gene2isoforms = new HashMap<>();

        apply(transcript2count.iterator(), (_p) -> gene2isoforms.computeIfAbsent(transcript2Gene.get(_p.getFirst()), (_g) -> new HashMap<>()).put(_p.getFirst(), _p.getSecond()));

        HashMap<MultiIsoformRegion, TranscriptSimulation> rv = new HashMap<>();
        for(Map.Entry<String, HashMap<String, Integer>> e : gene2isoforms.entrySet()) {
            MultiIsoformRegion gene = annot.getRegionById(e.getKey());
            if(gene == null)
                throw new FRuntimeException("unknown gene: %s in annotation", e.getKey());

            for(String trId : e.getValue().keySet()) {
                if(!gene.isoforms.containsKey(trId))
                    throw new FRuntimeException("unknown transcript: %s to gene: %s in annotation", trId, e.getKey());

            }

            rv.put(gene, new TranscriptSimulation(gene, e.getValue()));
        }

        return rv;
    }
}
