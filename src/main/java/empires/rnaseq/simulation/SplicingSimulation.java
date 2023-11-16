package empires.rnaseq.simulation;

import lmu.utils.*;
import empires.input.RNASeqSplicingInfo;
import empires.input.ReplicateSetInfo;
import empires.rnaseq.IsoformRegionGetter;
import empires.rnaseq.MultiIsoformRegion;
import empires.test.rnaseq.SimulationConfiguration;

import java.io.File;
import java.util.*;
import java.util.function.Consumer;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static lmu.utils.ObjectGetter.*;

public class SplicingSimulation {

    Logger log = LogConfig.getLogger();
    static Pattern replicatePattern = Pattern.compile("C(\\d+)R(\\d+)");


    SimulationConfiguration configuration = new SimulationConfiguration();


    Vector<empires.rnaseq.simulation.SimulatedSplicingCondition> conditions;
    HashMap<String, MultiIsoformRegion> simulatedGenes = new HashMap<>();
    HashMap<String, SplicingSimulatedGene> simulatedGeneCounts = new HashMap<>();

    Vector<ReplicateSetInfo> replicateSetInfos = null;


    int readLength = 100;
    IsoformRegionGetter annot;
    double factor;

    public SplicingSimulation(SimulationConfiguration config, File simulation) {
        this(config, simulation, 1.0);
    }

    public SplicingSimulation(IsoformRegionGetter annot) {
        this(annot, 1.0);
    }

    public SplicingSimulation(IsoformRegionGetter annot, double factor) {
        this.annot = annot;
        this.factor = factor;
    }

    public void setConfiguration(SimulationConfiguration configuration) {
        this.configuration = configuration;
    }

    public void setConditions(Vector<String> conditionNames) {
        conditions = map(conditionNames, (_cn) -> new empires.rnaseq.simulation.SimulatedSplicingCondition(_cn));
    }

    public void addSimulation(int conditionIdx, int replicateIdx, TranscriptSimulation transcriptSimulation) {
        if(conditions == null)
            throw new FRuntimeException("invalid operation, set conditions first!");

        simulatedGenes.put(transcriptSimulation.gene.id, transcriptSimulation.gene);
        empires.rnaseq.simulation.SimulatedSplicingCondition target = conditions.get(conditionIdx);
        while(target.replicates.size() <= replicateIdx) {
            target.replicates.add(new HashMap<>());
        }
        HashMap<MultiIsoformRegion, TranscriptSimulation> gene2simul = target.replicates.get(replicateIdx);
        gene2simul.put(transcriptSimulation.gene, transcriptSimulation);


    }

    public SplicingSimulation(SimulationConfiguration config, File simulation, double factor) {
        this(config.getAnnot(), factor);
        setConfiguration(config);

        Iterator<String[]> fieldSetIterator = FileUtils.getFieldSetIterator(simulation, "\t");
        HashMap<Integer, Integer> idx2conditionIdx = new HashMap<>();
        String[] header = fieldSetIterator.next();
        for(int i=0; i<header.length; i++) {
            Matcher m = replicatePattern.matcher(header[i]);
            if(!m.matches())
                throw new FRuntimeException("unexpected replicate header in transcript simulation file: %s exptected headers match: %s", header[i], replicatePattern.pattern());

            idx2conditionIdx.put(i, RegexUtils.getInt(m, 1) - 1);

        }

        setConditions(map(toSet(idx2conditionIdx.values()), (_i) -> "cond"+(_i + 1)));


        Vector<Vector<Pair<String, Integer>>> replicate2simulation = map(idx2conditionIdx.keySet(), (_k) -> new Vector<>());

        while(fieldSetIterator.hasNext()){
            String[] sp = fieldSetIterator.next();
            String transcriptId = sp[0];
            for(int i=1; i<sp.length; i++) {
                replicate2simulation.get(i-1).add(Pair.create(transcriptId, (int)(factor * Integer.parseInt(sp[i]))));
            }
        }

        Vector<HashMap<MultiIsoformRegion, TranscriptSimulation>> replicates = map(replicate2simulation, ( _r) -> TranscriptSimulation.groupTranscriptsPerGene(_r, annot));

        for(int i=0; i<replicates.size(); i++) {

            for(MultiIsoformRegion gene : replicates.get(i).keySet()) {
                simulatedGenes.put(gene.id, gene);
            }
        }
        for(Map.Entry<Integer, Integer> e : idx2conditionIdx.entrySet()) {
            conditions.get(e.getValue()).replicates.add(replicates.get(e.getKey()));
        }



    }

    Vector<TranscriptSimulation> getSimulatedCounts(empires.rnaseq.simulation.SimulatedSplicingCondition condition, MultiIsoformRegion gene) {

        return map(condition.replicates, (_m) -> _m.get(gene));
    }

    public Vector<Vector<Integer>> getSimulatedCounts(String geneId, String trId) {
        MultiIsoformRegion gene = simulatedGenes.get(geneId);
        if(gene == null)
            return null;

        Vector<Vector<Integer>> counts = new Vector<>();
        for(empires.rnaseq.simulation.SimulatedSplicingCondition cond : conditions) {
            counts.add(map(cond.replicates, (_m) -> _m.get(gene).getCountToSimulate(trId)));
        }
        return counts;
    }


    public Set<String> getSimulatedTranscripts(String geneId) {
        MultiIsoformRegion gene = simulatedGenes.get(geneId);
        if(gene == null)
            return null;

        HashSet<String> simulatedTranscripts = new HashSet<>();
        for(empires.rnaseq.simulation.SimulatedSplicingCondition ssc : conditions) {
            for(HashMap<MultiIsoformRegion, TranscriptSimulation> m : ssc.replicates) {
                TranscriptSimulation ts = m.get(gene);
                if(ts == null)
                    continue;

                simulatedTranscripts.addAll(ts.getSimulatedTranscriptIds());
            }
        }

        return simulatedTranscripts;
    }

    public UPair<Vector<TranscriptSimulation>> getSimulatedInfo(String geneId) {
        MultiIsoformRegion gene = simulatedGenes.get(geneId);
        return UPair.createU(getSimulatedCounts(conditions.get(0), gene), getSimulatedCounts(conditions.get(1), gene));

    }


    public RNASeqSplicingInfo  simulate(int maxTrsToUsePerGene, double MIN_AVG_COUNT_ON_STRONGER_SIDE, int minCountIntMaxCountCondition) {

        HashMap<String, Vector<String>> cond2replicateNames = new HashMap<>();

        for(empires.rnaseq.simulation.SimulatedSplicingCondition simulatedSplicingCondition: conditions) {
            Vector<String> names = mapIndex(simulatedSplicingCondition.replicates.size(), (_i) -> String.format("rep%d", _i + 1));
            cond2replicateNames.put(simulatedSplicingCondition.condition, names);
        }
        RNASeqSplicingInfo rnaSeqSplicingInfo = new RNASeqSplicingInfo(maxTrsToUsePerGene, MIN_AVG_COUNT_ON_STRONGER_SIDE, cond2replicateNames);
        rnaSeqSplicingInfo.setMinCountIntMaxCountCondition(minCountIntMaxCountCondition);
        int numSimulated = 0;

        HashSet<Integer> trlengths = new HashSet<>();
        for(String g : simulatedGenes.keySet()) {
            MultiIsoformRegion gene = simulatedGenes.get(g);
            apply(getSimulatedTranscripts(g), (_t) -> trlengths.add(gene.isoforms.get(_t).getCoveredLength()));
        }
        configuration.getBiasFactory().initLength(trlengths, 20);

        for(String g : simulatedGenes.keySet()) {
            SplicingSimulatedGene ssg = simulate(g, false);
            numSimulated++;

            if(numSimulated % 200 == 0) {
                log.info("simulated %d/%d genes", numSimulated, simulatedGenes.size());
            }
            HashMap<String, Vector<HashMap<Tuple, Double>>> cond2eqClassCounts = new HashMap<>();
            for(int condIdx = 0; condIdx < conditions.size(); condIdx++) {
                String condName = conditions.get(condIdx).condition;
                cond2eqClassCounts.put(condName, ssg.condition2replicate2equivianceClasses.get(condIdx));
            }
            rnaSeqSplicingInfo.addGeneInfo(g, cond2eqClassCounts);
        }

        rnaSeqSplicingInfo.getReplicateSetInfos();
        return rnaSeqSplicingInfo;
    }

    public int getNumConditions() {
        return conditions.size();
    }

    public Set<String> getSimulatedGenes() {
        return simulatedGenes.keySet();
    }

    public empires.rnaseq.simulation.SimulatedSplicingCondition getCondition(int idx) {
        return conditions.get(idx);
    }

    public SplicingSimulatedGene simulate(String geneId, boolean cache) {
        return  simulate(geneId, cache, null);
    }

    public interface ReadInReplicateConsumer {
        public void accept(String condition, int replicateId, SimulatedRead read);
    }
    public SplicingSimulatedGene simulate(String geneId, boolean cache, ReadInReplicateConsumer readConsumer) {
        MultiIsoformRegion gene = simulatedGenes.get(geneId);
        if(gene == null)
            return null;

        SplicingSimulatedGene ssGene = new SplicingSimulatedGene();
        ssGene.geneId = geneId;

        PositionBiasFactory biasFactory = configuration.getBiasFactory();

        for(SimulatedSplicingCondition cond : conditions) {
            Vector<HashMap<Tuple, Double>> replicates = new Vector<>();
            for(int replicateId = 0; replicateId < cond.replicates.size(); replicateId++) {
                HashMap<MultiIsoformRegion, TranscriptSimulation> m = cond.replicates.get(replicateId);

                final int repId = replicateId;
                Consumer<SimulatedRead> consumer = (readConsumer == null) ? null : (_r) -> readConsumer.accept(cond.condition, repId, _r);

                TranscriptSimulation trs = m.get(gene);
                HashMap<Tuple, Double> eqClasses = (trs == null) ? new HashMap<>() : trs.simulateTrSetCounts(biasFactory, consumer);
                for(Tuple t : eqClasses.keySet()) {
                    applyIndex(t.cardinality(), (_i) -> ssGene.trsWithCounts.add(t.getAsString(_i)));
                }
                replicates.add(eqClasses);
            }
            ssGene.condition2replicate2equivianceClasses.add(replicates);

        }

        if(cache) {
            simulatedGeneCounts.put(geneId, ssGene);
        }


        return ssGene;
    }


}
