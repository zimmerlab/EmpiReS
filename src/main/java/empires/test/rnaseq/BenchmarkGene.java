package empires.test.rnaseq;

import empires.SparseCumulativeDistribution;
import empires.input.RNASeqSplicingInfo;
import empires.rnaseq.simulation.SplicingSimulation;
import empires.rnaseq.simulation.TranscriptSimulation;
import lmu.utils.*;
import lmu.utils.plotting.PlotCreator;
import lmu.utils.swing.PagedDataTable;
import empires.DoubleDiffResult;
import empires.SingleFeatureDiffExp;
import lmu.utils.plotting.CachedPlotCreator;
import empires.FeatureInfo;

import java.awt.*;
import java.awt.image.BufferedImage;
import java.util.*;
import java.util.function.Function;

import static lmu.utils.ObjectGetter.*;

public class BenchmarkGene {
    public String geneId;
    public DoubleDiffResult ddr;
    boolean isTrue;
    empires.rnaseq.simulation.SplicingSimulation simulation;
    UPair<Vector<TranscriptSimulation>> simulCounts;
    empires.input.RNASeqSplicingInfo rnaSeqSplicingInfo;

    public BenchmarkGene(DoubleDiffResult ddr, Set<String> trues, empires.rnaseq.simulation.SplicingSimulation simulation, empires.input.RNASeqSplicingInfo rnaSeqSplicingInfo) {
        this.ddr = ddr;
        this.simulation = simulation;
        geneId = ddr.testName.split("\\.")[0];
        this.isTrue = (trues != null) && trues.contains(geneId);
        simulCounts = (simulation != null) ? simulation.getSimulatedInfo(geneId) :  null;
        this.rnaSeqSplicingInfo = rnaSeqSplicingInfo;
    }

    public static DataTable toExtendedTable(Vector<BenchmarkGene> v, Iterator<Pair<String, Function<BenchmarkGene, Object>>> additionalHeaders) {
        DataTable.HeaderGetterManager<BenchmarkGene> hgm = DataTable.buildHeader("gene", (BenchmarkGene bg) -> bg.geneId)
                .add("isTrue", (bg) -> bg.isTrue)
                .add("testname", (bg) -> bg.ddr.testName)
                .add("pval", bg -> bg.ddr.pval)
                .add("fdr", bg -> bg.ddr.fdr)
                .add("fc", bg -> bg.ddr.meanFC)
                .add("abs.fc", bg -> Math.abs(bg.ddr.meanFC))
                //.add("num.features1", bg -> bg.ddr.subFeatures1.size())
                //.add("num.features2", bg -> bg.ddr.subFeatures2.size())
                //.add("up-pval", bg -> bg.ddr.unpairedPval)
                //.add("up-fdr", bg -> bg.ddr.unpairedFDR)
                ;

        BenchmarkGene lead = v.firstElement();
        if(lead.simulCounts != null) {
            hgm.add("c1counts", bg -> bg.simulCounts.getFirst())
                    .add("c2counts", bg -> bg.simulCounts.getSecond());

        }

        if(additionalHeaders != null) {
            while(additionalHeaders.hasNext()) {
                Pair<String, Function<BenchmarkGene, Object>> header = additionalHeaders.next();
                hgm.add(header.getFirst(), (_b) -> header.getSecond().apply(_b));
            }
        }
        boolean gotsplicinginfo = v.size() > 0 && first(v).rnaSeqSplicingInfo != null;
        if(gotsplicinginfo) {
            hgm.add("features1", (bg) -> map(bg.ddr.subFeatures1, (_fn) -> bg.rnaSeqSplicingInfo.getEQClassWithCounts(_fn)))
                    .add("features2", (bg) -> map(bg.ddr.subFeatures2, (_fn) -> bg.rnaSeqSplicingInfo.getEQClassWithCounts(_fn)));
        } else {
            hgm.add("min.num.features", (bg) -> Math.min(bg.ddr.featureInfos1.size(), bg.ddr.featureInfos2.size()));
            hgm.add("num.features1", (bg) -> bg.ddr.featureInfos1.size());
            hgm.add("num.features2", (bg) -> bg.ddr.featureInfos2.size());


            hgm.add("features1", (bg) -> map(bg.ddr.featureInfos1, (_n) -> _n.feature));
            hgm.add("features2", (bg) -> map(bg.ddr.featureInfos2, (_n) -> _n.feature));
        }


        return DataTable.buildTable(v, hgm);
    }

    public static DataTable toTable(Vector<BenchmarkGene> v) {
        return toExtendedTable(v, new SingleIterator<>(null));
    }

    public static DataTable getTable(Vector<DoubleDiffResult> splicing, Set<String> trueGenes, empires.rnaseq.simulation.SplicingSimulation simulation, empires.input.RNASeqSplicingInfo rsi, Collection<Pair<String, Function<BenchmarkGene, Object>>> additionalHeaders) {
        NumUtils.sort(splicing, (_p) -> _p.pval, false);
        return BenchmarkGene.toExtendedTable(map(splicing, (_s) -> new BenchmarkGene(_s, trueGenes, simulation, rsi)), additionalHeaders.iterator());
    }

    public static DataTable getTable(Vector<DoubleDiffResult> splicing, Set<String> trueGenes, empires.rnaseq.simulation.SplicingSimulation simulation, empires.input.RNASeqSplicingInfo rsi, Pair<String, Function<BenchmarkGene, Object>>... additionalHeaders) {
        NumUtils.sort(splicing, (_p) -> _p.pval, false);

        return BenchmarkGene.toExtendedTable(map(splicing, (_s) -> new BenchmarkGene(_s, trueGenes, simulation, rsi)), new CollectionIterator<>(additionalHeaders));
    }

    public static PagedDataTable.MJFrame getWithDetails(DataTable dt) {
        PagedDataTable.MJFrame mjf = PagedDataTable.getDetailViewFrame(dt, null, null);


        mjf.addMenu("details", og ->
        {
            BenchmarkGene bg = (BenchmarkGene)og.getInData();
            DoubleDiffResult ddr = bg.ddr;

            double PSEUDO = 2.0;
            if(ddr.featureInfos1 != null && ddr.featureInfos2 != null) {
                PlotCreator pc = CachedPlotCreator.getPlotCreator();
                Function<Double, Double> sigfmt = (_d) -> NumUtils.logN(_d + PSEUDO, 2.0);
                Vector<BufferedImage> bims = new Vector<>();
                for(int i=0; i<2; i++) {
                    Vector<FeatureInfo> features = toVector((i == 0) ? ddr.featureInfos1 : ddr.featureInfos2);
                    double x = 0.0;


                    Vector<Pair<Double, String>> xtics = new Vector<>();
                    for(FeatureInfo fi : features) {
                        double X = x;
                        pc.scatter("", fi.cond1, (_d) -> X, (_d) -> sigfmt.apply(_d)).setColor(Color.GREEN.darker());
                        xtics.add(Pair.create(x, fi.feature));
                        double X2 = ++x;
                        pc.scatter("", fi.cond2, (_d) -> X2, (_d) -> sigfmt.apply(_d)).setColor(Color.RED);
                        x++;
                    }
                    pc.setLabels("", "log2 signal values", null);
                    pc.setTitle("EQC%d", i+1);
                    pc.setXTics(xtics);
                    bims.add(pc.getImage());
                }

                BufferedImage signals = ImageUtils.concat(bims);
                bims.clear();
                PlotCreator.BoxPlotBuilder bpb = pc.buildBoxPlot();
                for(int i=0; i<2; i++) {
                    Vector<FeatureInfo> features = toVector((i == 0) ? ddr.featureInfos1 : ddr.featureInfos2);


                    FeatureInfo firstF = first(features);

                    for(FeatureInfo ff : features) {

                        SingleFeatureDiffExp sdiffexp = new SingleFeatureDiffExp(ff.feature, ddr.diffExpManager, map(ff.cond1, (_d) -> sigfmt.apply(_d)), map(ff.cond2, (_d) -> sigfmt.apply(_d)));
                        //new SparseCumulativeDistribution(sdiffexp.scaledDistribution).drawLine(pc, ff.feature);
                        bpb.addBox(ff.feature, (ff.equals(firstF)) ? "EQC" + (i+1) : null, new SparseCumulativeDistribution(sdiffexp.scaledDistribution).quantiles,
                                (i == 0) ? Color.BLUE : Color.MAGENTA, true);

                    }


                }


                BufferedImage fcs = bpb.plot("top");
                BufferedImage combined = ImageUtils.vconcat(signals, fcs);
                ImageUtils.showImage("test", combined, false);
                return;
            }
            PlotCreator.BoxPlotBuilder bpb = ddr.drawSignalDistributionBoxPlots(CachedPlotCreator.getPlotCreator(), false);
            BufferedImage diffSignal = bpb.plot();

            bpb = ddr.drawDiffFCPlots(CachedPlotCreator.getPlotCreator(), false);
            BufferedImage diffFC = bpb.plot();

            ddr.drawDistributions(CachedPlotCreator.getPlotCreator());
            BufferedImage diffFC2 = CachedPlotCreator.getPlotCreator().getImage();
            ImageUtils.showImage(bg.geneId, ImageUtils.concat(diffSignal, diffFC, diffFC2), false);

        });
        mjf.setVisible(true);
        return mjf;
    }
    public static PagedDataTable.MJFrame showTable(Vector<DoubleDiffResult> splicing, Set<String> trueGenes, SplicingSimulation simulation, RNASeqSplicingInfo rsi) {

        DataTable dt = getTable(splicing, trueGenes, simulation, rsi);

        return getWithDetails(dt);

    }

}
