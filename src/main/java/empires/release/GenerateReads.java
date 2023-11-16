package empires.release;

import empires.rnaseq.IsoformRegionGetter;
import empires.rnaseq.simulation.SplicingSimulation;
import empires.test.rnaseq.SimulationConfiguration;
import lmu.utils.*;
import lmu.utils.tuple.Tuple3;

import java.io.File;
import java.io.PrintWriter;
import java.util.*;
import static lmu.utils.ObjectGetter.*;

public class GenerateReads {

    static void simulateReads(File trcounts, empires.test.rnaseq.SimulationConfiguration simulationConfiguration, File outdir) {
        Logger log = LogConfig.getLogger();
        SplicingSimulation simulation = new SplicingSimulation(simulationConfiguration, trcounts);

        PrintWriter pw = FileUtils.getWriter(outdir, "sample.table");
        pw.printf("label\tcondition\tbam\tfw\trw\tstrandness\n");
        for(int i=0; i<simulation.getNumConditions(); i++) {
            for(int repI = 0; repI < simulation.getCondition(i).replicates.size(); repI++) {
                String cond = simulation.getCondition(i).condition;
                pw.printf("%s\t%s\t%s\t%s\t%s\ttrue\n", empires.test.rnaseq.FastQGenerator.getReplicatePrefix(cond, repI), cond,
                        empires.test.rnaseq.FastQGenerator.getReplicatePrefix(cond, repI)+"_sorted.bam",
                        empires.test.rnaseq.FastQGenerator.getFastqFileName(cond, repI, true),
                        empires.test.rnaseq.FastQGenerator.getFastqFileName(cond, repI, false));
            }
        }
        pw.close();

        Vector<Tuple3<String, String, Integer>> gene2tr2length = new Vector<>();
        for (String geneId : simulation.getSimulatedGenes()) {
            empires.rnaseq.MultiIsoformRegion gene = simulationConfiguration.getAnnot().getRegionById(geneId);
            for(String trId : simulation.getSimulatedTranscripts(geneId)) {
                gene2tr2length.add(Tuple3.create(geneId,  trId, gene.isoforms.get(trId).getCoveredLength()));
            }
        }

        Vector<Tuple3<String, String, Integer>> too_short = filter(gene2tr2length, (_t) -> _t.get2() < simulationConfiguration.readLength);
        HashSet<String> too_short_genes = mapToSet(too_short, (_t) -> _t.get0());

        if(too_short_genes.size() > 0) {
            throw new FRuntimeException("invalid request, got %d transcripts to simulate (in %d genes) which are shorter than the requested read length: %d\n%s",
                    too_short.size(), too_short_genes.size(), simulationConfiguration.readLength,
                    StringUtils.joinObjects("\n", too_short, (_t) -> String.format("%s\t%s", _t.get0(), _t.get1())));
        }
        int numGenes = simulation.getSimulatedGenes().size();
        int numSimulated = 0;

        for (String geneId : simulation.getSimulatedGenes()) {
            Set<String> trs = simulation.getSimulatedTranscripts(geneId);

            BufferPrintWriter info = new BufferPrintWriter();
            for(String tr : trs) {
                info.printf("\t%s: counts: %s\n", tr, simulation.getSimulatedCounts(geneId, tr));
            }
            numSimulated++;
            log.info("will simulate gene %d/%d %s with counts:\n%s\n", numSimulated, numGenes, geneId, info.getBuffer());
            simulation.simulate(geneId, false, simulationConfiguration.getReadInReplicateConsumer());
            log.info("simulated gene: %d/%d %s ready", numSimulated, numGenes, geneId);

        }
        simulationConfiguration.fastQGenerator.close();
    }


    public static void main(String[] args) {
        SimpleOptionParser cmd = new SimpleOptionParser("trcounts", "gtf", "genome", "genomeidx", "od",  "mutrate",  "biaspos", "readlength", "fraglengthmean", "fraglengthsd", "nobams");
        cmd.setFile("genome", "genomeidx", "biaspos");
        cmd.setDir("od");
        cmd.setOptional("biaspos");
        cmd.setInt("readlength", "fraglengthmean");
        cmd.setDefault("readlength", "100");
        cmd.setFile("gtf");
        cmd.setFile("trcounts");
        cmd.setDouble("mutrate", "fraglengthsd");
        cmd.setDefault("mutrate", "0.01");
        cmd.setDefault("fraglengthmean", "200");
        cmd.setDefault("fraglengthsd", "60.0");
        cmd.setSwitches("nobams");



        if(!OptionParser.parseParams(args, false, false, true, true, cmd))
            return;



        File trcounts = cmd.getFile("trcounts");
        IsoformRegionGetter annot = new empires.rnaseq.GFFBasedIsoformRegionGetter(cmd.getFile("gtf"), null, null);

        File od = cmd.getFile("od");
        empires.test.rnaseq.SimulationConfiguration simulationConfiguration = new SimulationConfiguration(annot);

        simulationConfiguration.reportOutDir = od;
        simulationConfiguration.isoformRegionGetter = annot;
        simulationConfiguration.fastQGenerator = new empires.test.rnaseq.FastQGenerator(annot, cmd.getFile("genome"), cmd.getFile("genomeidx"), cmd.getDouble("mutrate"));
        simulationConfiguration.fastQGenerator.setWriteInfoFile(false);
        simulationConfiguration.fastQGenerator.setGenerateBams(!cmd.isSet("nobams"));
        simulationConfiguration.fastQGenerator.setFastQOutDir(od);
        simulationConfiguration.readLength = cmd.getInt("readlength");
        simulationConfiguration.fragmentLengthMean = cmd.getInt("fraglengthmean");
        simulationConfiguration.fragmentLengthSD = cmd.getDouble("fraglengthsd");



        simulationConfiguration.setPositionBiases(cmd.getOptionalFile("biaspos"));

        simulateReads(trcounts, simulationConfiguration, od);

    }

}
