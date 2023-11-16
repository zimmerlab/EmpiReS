package empires.test.rnaseq;

import empires.rnaseq.IsoformRegionGetter;
import empires.rnaseq.simulation.PositionBiasFactory;
import empires.rnaseq.simulation.SplicingSimulation;
import lmu.utils.*;
import lmu.utils.plotting.PlotCreator;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.awt.image.BufferedImage;
import java.io.File;
import java.util.Set;
import java.util.Vector;

import static lmu.utils.ObjectGetter.*;

public class SimulationConfiguration {


    Set<String> trueSplic;
    Set<String> trueDiff;


    public File reportOutDir = null;
    public FastQGenerator fastQGenerator;
    public empires.rnaseq.IsoformRegionGetter isoformRegionGetter;
    public empires.rnaseq.simulation.PositionBiasFactory biasFactory = new empires.rnaseq.simulation.PositionBiasFactory();
    public int fragmentLengthMean = 200;
    public double fragmentLengthSD = 60;
    public int readLength = 60;
    public boolean paired = true;
    public File positionalBias = null;

    public SimulationConfiguration() {

    }


    public SimulationConfiguration(empires.rnaseq.IsoformRegionGetter isoformRegionGetter) {
        this();
        setAnnotation(isoformRegionGetter);
    }

    public empires.rnaseq.IsoformRegionGetter getAnnot() {
        return isoformRegionGetter;
    }

    public void setAnnotation(IsoformRegionGetter isoformRegionGetter) {
        this.isoformRegionGetter = isoformRegionGetter;
    }

    public SplicingSimulation.ReadInReplicateConsumer getReadInReplicateConsumer() {
        return (fastQGenerator == null) ? null : fastQGenerator.readInReplicateConsumer;
    }

    public empires.rnaseq.simulation.PositionBiasFactory getBiasFactory() {

        return (biasFactory != null) ? biasFactory : (biasFactory = new PositionBiasFactory(positionalBias, new NormalDistribution(fragmentLengthMean, fragmentLengthSD), readLength, paired));
    }

    public void setPositionBiases(File f) {
        positionalBias = f;
        biasFactory = null;
    }

    public BufferedImage writeInfos(Vector<SimulatedGeneInfo> ginfos, int MAXTRCOUNT, boolean doPlot) {

        if (reportOutDir != null) {

            SimulatedGeneInfo.toTable(ginfos).writeCSV(new File(reportOutDir, "geneInfos.tsv"));
        }


        BufferedImage bim = null;
        if(doPlot || reportOutDir != null) {

            PlotCreator pc = new PlotCreator();
            pc.setTitle("%d/%d failures target reduced size: %d", filteredSize(ginfos, (_g) -> _g.gotTrReductionError), ginfos.size(), MAXTRCOUNT);

            pc.cumhist("num.trs.withcounts", map(ginfos, (_g) -> _g.numTranscriptsWithCounts), ginfos.size(), false, false);
            pc.cumhist("num.EQs", map(ginfos, (_g) -> _g.numEqClasses), ginfos.size(), false, false);
            pc.cumhist("num.reduced.EQs", map(ginfos, (_g) -> _g.numReducedEqClasses), ginfos.size(), false, false);

            pc.setLabels("setsize", "frequency", "bottomright");
            pc.setLog(true, false);

            bim = pc.getImage();
            pc.destroy();
        }


        if(reportOutDir != null) {
            ImageUtils.saveImage(bim, new File(reportOutDir, "eqClassCumHist.png"));
        }

        return bim;

    }

}
