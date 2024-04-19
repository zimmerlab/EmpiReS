package nlEmpiRe.test.rnaseq;

import lmu.utils.DataTable;
import lmu.utils.MapBuilder;
import lmu.utils.UPair;

import java.util.HashMap;
import java.util.Vector;
import java.util.function.Function;

public class SimulatedGeneInfo {
    String geneId;
    Vector<String> transcripts;
    boolean diffexp;
    boolean diffaltsplic;

    CondTrCountInfo c1CountInfo;
    CondTrCountInfo c2CountInfo;
    String cond1;
    HashMap<String, Vector<Integer>> c1_tr2count = new HashMap<>();
    HashMap<String, Vector<Integer>> c2_tr2count = new HashMap<>();

    boolean gotTrReductionError = false;

    int numTranscripts;
    int numSimulatedTranscripts;
    int numTranscriptsWithCounts;

    int numEqClasses;
    int numReducedEqClasses;

    public SimulatedGeneInfo(String geneId) {
        this.geneId = geneId;

    }

    void update(String cond, String tr, int count) {
        if(cond1 == null) {
            cond1 = cond;
        }

        HashMap<String, Vector<Integer>> target = (cond1.equals(cond)) ? c1_tr2count : c2_tr2count;
        MapBuilder.updateV(target, tr, count);
    }


    void compile() {
        c1CountInfo = new CondTrCountInfo(c1_tr2count);
        c2CountInfo = new CondTrCountInfo(c2_tr2count);
    }

    static public DataTable toTable(Vector<SimulatedGeneInfo> genes) {
        DataTable.HeaderGetterManager<SimulatedGeneInfo> hgm = DataTable.buildHeader("geneId", (SimulatedGeneInfo sg) -> sg.geneId)
                .add("num.transcripts", (sg) -> sg.numTranscripts)
                .add("num.simul.trs", (sg) -> sg.numSimulatedTranscripts)
                .add("num.tr.withcounts", (sg) -> sg.numTranscriptsWithCounts)
                .add("num.EQclasses", (sg) -> sg.numEqClasses)
                .add("num.reduced.EQclasses", (sg) -> sg.numReducedEqClasses)
                .add("lost.diffsplic", (sg) -> sg.gotTrReductionError)
                .add("simulTranscripts", (sg) -> sg.transcripts)
                .add("diffexp", (sg) -> sg.diffexp)
                .add("diffsplic", (sg) -> sg.diffaltsplic)
                ;

        for(int ci : new int[]{0,1}) {
            String cn = (ci == 0) ? "c1" : "c2";
            Function<SimulatedGeneInfo, CondTrCountInfo> g2cinfo = (s) -> (ci == 0) ?  s.c1CountInfo : s.c2CountInfo;
            for(int mi : new int[]{0,1}) {
                String mn = (mi == 0) ? "major" : "minor";
                Function<SimulatedGeneInfo, UPair<Integer>> g2mm = (s) -> (mi == 0) ? g2cinfo.apply(s).majorMinMax : g2cinfo.apply(s).minorMinMax;

                hgm.add(cn+"."+mn+".min", (s) -> g2mm.apply(s).getFirst());
                hgm.add(cn+"."+mn+".max", (s) -> g2mm.apply(s).getSecond());
            }
        }

        return DataTable.buildTable(genes, hgm);
    }
}
