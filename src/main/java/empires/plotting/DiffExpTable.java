package nlEmpiRe.plotting;

import lmu.utils.*;
import lmu.utils.plotting.CachedPlotCreator;
import lmu.utils.plotting.PlotCreator;
import lmu.utils.swing.PagedDataTable;
import nlEmpiRe.*;
import nlEmpiRe.input.ReplicateSetInfo;

import java.awt.image.BufferedImage;
import java.util.HashMap;
import java.util.Vector;
import java.util.function.Function;

import static lmu.utils.IteratorUtils.rangev;
import static lmu.utils.ObjectGetter.*;

public class DiffExpTable {

    static final public double DEFAULT_FC_THRESHOLD = 1.0;
    static final public double DEFAULT_FDR_THRESHOLD = 0.05;
    static final long MAX_PC_TIMEOUT = 1000;

    Vector<DiffExpResult> expResults;
    Function<String, Boolean> labelFunc;
    DataTable table;
    PagedDataTable.MJFrame mjFrame;
    BufferedImage vulcano;
    BufferedImage pvalDistrib;
    BufferedImage errDistribs;
    double minFCThreshold;
    double maxFDRThreshold;

    NormalizedReplicateSet from;
    NormalizedReplicateSet to;
    PlotCreator plotCreator;
    long lastPlotted = System.currentTimeMillis();

    HashMap<String, DiffExpResult> lookup = null;

    public DiffExpTable(NormalizedReplicateSet from, NormalizedReplicateSet to, Vector<DiffExpResult> expResults) {
        this(from, to, expResults, DEFAULT_FC_THRESHOLD, DEFAULT_FDR_THRESHOLD, null);
    }
    public DiffExpTable(NormalizedReplicateSet from, NormalizedReplicateSet to, Vector<DiffExpResult> expResults, double minFC, double maxFDR) {
        this(from, to, expResults, minFC, maxFDR, null);
    }

    public DiffExpTable(NormalizedReplicateSet from, NormalizedReplicateSet to, Vector<DiffExpResult> expResults, double minFC, double maxFDR, Function<String, Boolean> labelFunc) {
        this.from = from;
        this.to = to;
        this.expResults = expResults;
        this.labelFunc = labelFunc;
        this.minFCThreshold = minFC;
        this.maxFDRThreshold = maxFDR;
    }


    public DiffExpResult getDiffExpById(String featureId) {
        if(lookup != null)
            return lookup.get(featureId);

        lookup = buildReverseMap(expResults, (_e) -> _e.combinedFeatureName);
        return lookup.get(featureId);
    }

    public PlotCreator getPlotCreator() {
        return CachedPlotCreator.getPlotCreator();
    }

    static int[] FCCONFINTERAL_THRESHOLDS = new int[]{90, 95, 99};

    public DataTable getTable() {
        return getTable(new Vector<>());
    }
    
    public DataTable getTable(Iterable<Pair<String, Function<DiffExpResult, Object>>> additionalHeaders) {
        if(table != null)
            return table;

        DataTable.HeaderGetterManager<DiffExpResult> hgm = DataTable.buildHeader("id", (DiffExpResult e) -> e.combinedFeatureName);
        if(labelFunc != null) {
            hgm.add("isTrue", (_e) -> labelFunc.apply(_e.combinedFeatureName));
        }
        hgm.add("numFeatures", (_e) -> _e.featureNames.size());
        if(additionalHeaders != null) {
            for(Pair<String, Function<DiffExpResult, Object>> p : additionalHeaders) {
                hgm.add(p.getFirst(), (_e) -> p.getSecond().apply(_e));
            }
        }
        hgm.add("pval", (_e) -> _e.pval);
        hgm.add("abs.log2FC", (_e) -> Math.abs(_e.estimatedFC));
        hgm.add("log2FC", (_e) -> _e.estimatedFC);
        hgm.add("fdr", (_e) -> _e.fdr);
        hgm.add("SDcorr", (_e) -> StringUtils.joinObjects(",",_e.sdCorrectionScales, (_sd) -> String.format("%.2f", _sd)));
        hgm.add("fc.pval", (_e) -> _e.fcEstimatePval);
        hgm.add("fc.fdr", (_e) -> _e.fcEstimateFDR);
        hgm.add("nonde.fcwidth", (_e) -> _e.getNonDiffExpLog2FCThresholdToDefaultPvalThreshold());

        Function<Integer, Double> interval2alpha = (fcw) -> 0.5 * ((100 - fcw) / 100.0);

        for(int fcw : FCCONFINTERAL_THRESHOLDS) {
            hgm.add("fcCI."+fcw+".start", (_e) -> (_e.combinedEmpiricalFoldChangeDistrib == null) ? Double.NEGATIVE_INFINITY :
                    _e.combinedEmpiricalFoldChangeDistrib.getFoldChangeToCumulativeFrequency(interval2alpha.apply(fcw)));

            hgm.add("fcCI."+fcw+".end", (_e) -> (_e.combinedEmpiricalFoldChangeDistrib == null) ? Double.POSITIVE_INFINITY :
                    _e.combinedEmpiricalFoldChangeDistrib.getFoldChangeToCumulativeFrequency( 1.0 - interval2alpha.apply(fcw)));

        }


        return table = DataTable.buildTable(expResults, hgm);

    }

    public BufferedImage getVulcano() {
        if(vulcano != null)
            return vulcano;

        PlotCreator pc = getPlotCreator();
        Vector<DiffExpResult> trueExp = (labelFunc == null) ? null : filter(expResults, (_e) -> labelFunc.apply(_e.combinedFeatureName));

        Function<DiffExpResult, Boolean> calledLabel = (_e) -> Math.abs(_e.estimatedFC) >= minFCThreshold && _e.fcEstimateFDR < maxFDRThreshold;
        String title = "";
        if(trueExp != null) {
            title = String.format("up: %d/%d dn: %d/%d\n(thresholds: fc: %.2f fdr: %.3f)",
                    filteredSize(trueExp, (_e) -> Math.abs(_e.estimatedFC) >= 0 && calledLabel.apply(_e)),
                    filteredSize(expResults, (_e) -> Math.abs(_e.estimatedFC) >= 0 && calledLabel.apply(_e)),

                    filteredSize(trueExp, (_e) -> Math.abs(_e.estimatedFC) <= 0 &&  calledLabel.apply(_e)),
                    filteredSize(expResults, (_e) -> Math.abs(_e.estimatedFC) <= 0 && calledLabel.apply(_e)),
                    minFCThreshold, maxFDRThreshold
                    );
        } else {
            title = String.format("up: %d dn: %d\n(thresholds: fc: %.2f fdr: %.3f)",
                    filteredSize(expResults, (_e) -> _e.estimatedFC > 0 && calledLabel.apply(_e)),

                    filteredSize(expResults, (_e) -> _e.estimatedFC < 0 && calledLabel.apply(_e)),
                    minFCThreshold, maxFDRThreshold
            );
        }

        pc.scatter(String.format("all (%d)", expResults.size()), expResults, (_e) -> _e.estimatedFC, (_e) -> -NumUtils.logN(_e.fcEstimateFDR, 10));
        if(trueExp != null) {
            pc.scatter(String.format("trues (%d)", trueExp.size()), trueExp, (_e) -> _e.estimatedFC, (_e) -> -NumUtils.logN(_e.fcEstimateFDR, 10));
        }
        pc.setTitle(title);
        pc.abline("", null, -NumUtils.logN(maxFDRThreshold, 10.0), null, null);
        pc.abline("", - minFCThreshold, null, null, null);
        pc.abline("",  minFCThreshold, null, null, null);
        pc.setLabels("log2FC", "-log10(fdr)", "bottomright");

        vulcano = pc.getImage(false);

        return  vulcano;
    }

    BufferedImage getPvalDistrib() {
        if(pvalDistrib != null)
            return pvalDistrib;

        PlotCreator pc = getPlotCreator();
        Vector<DiffExpResult> trueExp = (labelFunc == null) ? null : filter(expResults, (_e) -> labelFunc.apply(_e.combinedFeatureName));



        pc.cumhist(String.format("all(%d)", expResults.size()), expResults, (_e) -> _e.fdr, expResults.size(), false, true);
        if(trueExp != null) {
            pc.cumhist(String.format("trues(%d)", trueExp.size()), trueExp, (_e) -> _e.fdr, trueExp.size(), false, true);
        }
        pc.setLabels("fdr", "frequency", "bottomright");

        return pvalDistrib = pc.getImage(false);


    }


    BufferedImage getErrDistribs() {
        if(errDistribs != null)
            return errDistribs;

        PlotCreator pc = getPlotCreator();

        NormalizedReplicateSetPlotting pfrom = new NormalizedReplicateSetPlotting(from);
        NormalizedReplicateSetPlotting pto = new NormalizedReplicateSetPlotting(to);

        return errDistribs = ImageUtils.vconcat(pfrom.plotBackGroundDistribs(pc), pto.plotBackGroundDistribs(pc));
    }

    public PagedDataTable.MJFrame showInteractiveTable(Pair<String, Function<DiffExpResult, Object>> ... additionalHeaders) {
        return showInteractiveTable(toVector(additionalHeaders));
    }

    public PagedDataTable.MJFrame showInteractiveTable(Iterable<Pair<String, Function<DiffExpResult, Object>>> additionalHeaders) {

        if(mjFrame != null)
            return mjFrame;

        mjFrame = PagedDataTable.getDetailViewFrame(getTable(additionalHeaders), null, null);

        mjFrame.addMenu("vulcano", (og) -> ImageUtils.showImage("vulcano plot", getVulcano(), false));
        mjFrame.addMenu("fdrDistrib", (og) -> ImageUtils.showImage("fdr Distrib", getPvalDistrib(), false));

        mjFrame.addMenu("pval2fc.pval", (og) ->
        {
            PlotCreator pc = getPlotCreator();
            Vector<UPair<Double>> data = map(expResults, (_e) -> UPair.createU(-NumUtils.logN(_e.pval, 10), -NumUtils.logN(_e.fcEstimatePval, 10)));
            pc.scatter("" , data, (_p) -> _p.getFirst(), (_p) -> _p.getSecond());
            pc.setLabels("pval", "fcpval", null);
            ImageUtils.showImage("pval scatter", pc.getImage(), false);
        });
        mjFrame.addMenu("show fc estimation width plot", (og) -> {
                    DiffExpManager diffExpManager = ((DiffExpResult)og.getInData()).getDiffExpManager();
                    PlotCreator pc = getPlotCreator();
                    Vector<DiffExpResult> measured = filter(expResults, (_r) -> _r.isMeasured());
                    for(int level : toVector(50, 75, 90, 95)) {
                        Vector<Double> fcWidths = filter_and_map(expResults, (_r) -> _r.isMeasured(), (_r) -> _r.getFoldChangeWidthToConfidenceLevel(level));
                        pc.cumhist("conf="+ level, fcWidths, fcWidths.size(), false, false);
                    }
                    pc.setLabels("confidence interval width", "freq.", "bottomright");
                    ImageUtils.showImage("confidence interval widths", pc.getImage(), false);
                }
        );

        mjFrame.addMenu("show non-diff fc-thresholds", (og) -> {
                    DiffExpManager diffExpManager = ((DiffExpResult)og.getInData()).getDiffExpManager();
                    PlotCreator pc = getPlotCreator();
                    Vector<DiffExpResult> measured = filter(expResults, (_r) -> _r.isMeasured());
                    Vector<Double> pvals = map(measured, (_de) -> _de.getNonDiffExpLog2FCThresholdToDefaultPvalThreshold());
                    pc.cumhist("all measured", pvals, pvals.size(), false, false);

                    Vector<Double> nonsig_pvals = filter_and_map(measured, (_de) -> _de.fcEstimateFDR > 0.05  || Math.abs(_de.estimatedFC) < 1.0, (_de) -> _de.getNonDiffExpLog2FCThresholdToDefaultPvalThreshold());
                    pc.cumhist("non signif (fdr>0.05 || |fc|<1)", nonsig_pvals, nonsig_pvals.size(), false, false);
                    pc.setLabels("non-diff-exp FC threshold leading to non-diff-exp pval="+ DiffExpResult.DEFAULT_PVAL_THRESHOLD, "freq.", "bottomright");
                    ImageUtils.showImage("confidence interval widths", pc.getImage(), false);
                }
        );

        mjFrame.addMenu("show replicate scatters", (og) ->
        {
            DiffExpManager diffExpManager = ((DiffExpResult)og.getInData()).getDiffExpManager();
            PlotCreator pc = getPlotCreator();

            Vector<BufferedImage> condimages = new Vector<>();
            for(int i=0; i<2; i++) {
                NormalizedReplicateSet nrs = ((i == 0) ? diffExpManager.replicateSetFrom : diffExpManager.replicateSetTo);
                Vector<Double> shifts = nrs.getNormalization().getShifts();
                ReplicateSetInfo rsi = nrs.getInData();
                Vector<BufferedImage> scatters = new Vector<>();
                for(int j=0; j< nrs.getNumReplicates(); j++) {
                    for(int k=j+1; k< nrs.getNumReplicates(); k++) {
                        Vector<Double> v1 = rsi.getLog2Data().get(j);
                        Vector<Double> v2 = rsi.getLog2Data().get(k);
                        double s1 = shifts.get(j);
                        double s2 = shifts.get(k);
                        Vector<UPair<Double>> data = filter_and_map(rangev(v1.size()), (_i) -> !Double.isNaN(v1.get(_i)) && !Double.isNaN(v2.get(_i)), (_i) -> UPair.createU(v1.get(_i) + s1, v2.get(_i) + s2));

                        pc.dens_scatter2("", data, (_p) -> _p.getFirst(), (_p) -> _p.getSecond() );
                        pc.setTitle("%s %d vs %d", rsi.getReplicateSetName(), j, k);
                        pc.setLabels(j +" log2 signals", k + " log2signals", null);
                        scatters.add(pc.getImage());
                    }

                }
                condimages.add(ImageUtils.vconcat(scatters));
            }
            ImageUtils.showImage("replicate scatters", ImageUtils.concat(condimages), false);
        });

        mjFrame.addMenu("show normalization", (og) ->
        {
            DiffExpManager diffExpManager = ((DiffExpResult)og.getInData()).getDiffExpManager();
            PlotCreator pc = getPlotCreator();
            Vector<Double> diffs = diffExpManager.getMedianFCDistrib();
            pc.hist("diff", diffs, 100, false, false);


            Vector<Vector<Double>> errs = new Vector<>();

            Vector<Double> err1 = diffExpManager.replicateSetFrom.getMedianValues(diffExpManager.totalFeatures);
            Vector<Double> err2 = diffExpManager.replicateSetTo.getMedianValues(diffExpManager.totalFeatures);

            ErrorEstimationDistribution shiftErr = PairwiseMedianImpliedFCPeakErrorEstimation.getShiftError(err1, err2);
            ErrorEstimationDistribution.Peak fcpeak = shiftErr.getBestFCPeak();

            pc.setTitle("%s vs %s shift: %.2f", diffExpManager.replicateSetFrom.getInData().getReplicateSetName(),
                    diffExpManager.replicateSetTo.getInData().getReplicateSetName(), fcpeak.summit * ErrorEstimationDistribution.norm + shiftErr.getMinFC());
            pc.setLabels("log2fc", "freq", null);

            double minFC = shiftErr.getFoldChangeToCumulativeFrequency(0.01);
            double maxFC = shiftErr.getFoldChangeToCumulativeFrequency(0.99);

            pc.abline("", fcpeak.peakStart * ErrorEstimationDistribution.norm + shiftErr.getMinFC(), null, null, null);
            pc.abline("", fcpeak.summit * ErrorEstimationDistribution.norm + shiftErr.getMinFC(), null, null, null);
            pc.abline("", fcpeak.peakEnd * ErrorEstimationDistribution.norm + shiftErr.getMinFC(), null, null, null);

            pc.setLabels("between condition implied log2 fc-s", "freq", null);
            pc.setTitle("shifts: %s", diffExpManager.shifts);
            pc.abline("", 0.0, null, null, null);
            pc.setLimits(minFC, maxFC, null, null);

            BufferedImage hist = pc.getImage();

            Vector<UPair<Double>> errHist = new Vector<>();
            double step = 0.05;
            double halfstep = 0.5 * step;;
            for(double fc = minFC + halfstep; fc < maxFC; fc+=halfstep) {
                errHist.add(UPair.createU(fc, shiftErr.getCumulativeFrequencyToFoldChange(fc + halfstep) - shiftErr.getCumulativeFrequencyToFoldChange(fc - halfstep)));
            }

            pc.line("errhist", errHist, (_p) -> _p.getFirst(), (_p) -> _p.getSecond());
            pc.setLabels("normed fcdiff", "freq", null);
            pc.setTitle("fc: %.2f - %.2f shift: %.2f\n%s\npeak: %.4f", minFC, maxFC, diffExpManager.shifts.get(0), shiftErr, shiftErr.getBestFCShift());
            pc.abline("", diffExpManager.shifts.get(0), null, null, null);
            pc.setLimits(minFC, maxFC, null, null);
            BufferedImage errhist = pc.getImage();


            pc.cumhist("diff", diffs, diffs.size(), false, false);
            pc.setLabels("normed fcdiff", "freq", null);

            pc.abline("", 0.0, null, null, null);
            ImageUtils.showImage("normed fcs", ImageUtils.concat(hist, errhist, pc.getImage()), false);

        });

        mjFrame.addMenu("showErrorDistribs", (og) ->
        {
            ImageUtils.showImage("err distribs", getErrDistribs(), false);
        });


        mjFrame.addMenu("empirical combined details", (og) ->
        {

            DiffExpResult exp = (DiffExpResult) og.getInData();
            PlotCreator pc = getPlotCreator();
            BufferedImage bim = getDetailImage(pc, exp);
            String info = "detail for " + exp.combinedFeatureName;

            if(labelFunc != null) {
                info += String.format(" (%s)", labelFunc.apply(exp.combinedFeatureName));
            }
            ImageUtils.showImage(info, bim, false);

        });

        /*
        mjFrame.addMenu("new details", (og) ->
                {

                    DiffExpResult exp = (DiffExpResult) og.getInData();
                    PlotCreator pc = getPlotCreator();
                    PlotCreator.BoxPlotBuilder bpb = pc.buildBoxPlot();
                    DiffExpManager diffExpManager = exp.getDiffExpManager();
                    Vector<ErrorEstimationDistribution> perFeatureErrDistribs = exp.getFeatureErrorDistribs();
                    for (int i = 0; i < exp.featureNames.size(); i++) {
                        String fn = exp.featureNames.get(i);
                        double fc = exp.perFeatureDistribEstimates.get(i).getFirst();
                        ErrorEstimationDistribution ferr = perFeatureErrDistribs.get(i);
                        bpb.addBox(String.format("%s (%.2f)", fn, fc), new SparseCumulativeDistribution(ferr).quantiles);
                    }
                    bpb.addBox("combined", new SparseCumulativeDistribution(ErrorEstimationDistribution.fromNormalDistrib(exp.estimatedFC, exp.estimatedSD)).quantiles);
                    pc.setLabels("", "log2fc", "topright");
                    String title = String.format("%s FC: %.2f fdr: %g (%d features)", exp.combinedFeatureName, exp.estimatedFC, exp.fdr, exp.featureNames.size());
                    pc.setTitle(title);
                    BufferedImage diffFC = bpb.plot();

                    new SparseCumulativeDistribution(ErrorEstimationDistribution.fromNormalDistrib(exp.estimatedFC, exp.estimatedSD)).drawLine(pc, "est.FC distrib");

                    pc.abline("", exp.estimatedFC, null, null, null);
                    String title2 = String.format("%s FC: %.2f pval: %g (%d features)", exp.combinedFeatureName, exp.estimatedFC, exp.fdr, exp.featureNames.size());
                    pc.setTitle(title2);
                    BufferedImage diffFCFullPlot = pc.getImage();


                    String info = "detail for " + exp.combinedFeatureName;

                    if(labelFunc != null) {
                        info += String.format(" (%s)", labelFunc.apply(exp.combinedFeatureName));
                    }
                    ImageUtils.showImage(info, ImageUtils.vconcat(getSignalPlot(exp, pc), diffFC, diffFCFullPlot), false);

                });


         */
        /*mjFrame.addMenu("details", (og) ->
        {

            DiffExpResult exp = (DiffExpResult)og.getInData();
            PlotCreator pc = getPlotCreator();
            PlotCreator.BoxPlotBuilder bpb = pc.buildBoxPlot();
            DiffExpManager diffExpManager = exp.getDiffExpManager();
            for(int i=0; i<exp.featureNames.size(); i++) {
                String fn = exp.featureNames.get(i);
                double fc = exp.perFeatureDistribEstimates.get(i).getFirst();
                ErrorEstimationDistribution ferr = new ErrorEstimationDistribution(toVector(diffExpManager.getDiffError(fn)), toVector(fc));
                bpb.addBox(String.format("%s (%.2f)", fn, fc), new SparseCumulativeDistribution(ferr).quantiles);
            }
            bpb.addBox("combined", new SparseCumulativeDistribution(exp.foldChangeDistrib).quantiles);

            //new SparseCumulativeDistribution(exp.foldChangeDistrib).drawLine(pc, "diff");
            //new SparseCumulativeDistribution(exp.H0foldChangeDistrib).drawLine(pc, "H0");
            //pc.abline("", exp.estimatedFC, null, null, null);
            pc.setLabels("", "log2fc", "topright");
            String title = String.format("%s FC: %.2f fdr: %g (%d features)", exp.combinedFeatureName, exp.estimatedFC, exp.fdr, exp.featureNames.size());
            pc.setTitle(title);
            BufferedImage diffFC = bpb.plot();

            new SparseCumulativeDistribution(exp.foldChangeDistrib).drawLine(pc, "est.FC distrib");

            pc.abline("", exp.estimatedFC, null, null, null);
            String title2 = String.format("%s FC: %.2f pval: %g (%d features)", exp.combinedFeatureName, exp.estimatedFC, exp.fdr, exp.featureNames.size());
            pc.setTitle(title2);
            BufferedImage diffFCFullPlot = pc.getImage();

            //plot features
            Vector<UPair<Double>> fromData = new Vector<>();
            Vector<UPair<Double>> toData = new Vector<>();
            int xCoord = 0;
            Vector<Pair<Double, String>> xtics = new Vector<>();

            for(int i=0; i<exp.featureNames.size(); i++) {
                double x1 = xCoord++;
                double x2 = xCoord++;

                Vector<Double> n1 = from.getNormed(exp.featureNames.get(i)).replicates;
                Vector<Double> n2 = to.getNormed(exp.featureNames.get(i)).replicates;
                fromData.addAll(map(n1, (_d) -> UPair.createU(x1, _d)));
                toData.addAll(map(n2, (_d) -> UPair.createU(x2, _d)));

                xtics.add(Pair.create(0.5 * ( x1 + x2 ), exp.featureNames.get(i)));
            }

            pc.setXTics(xtics);
            pc.scatter("condition.from", fromData, (_p) -> _p.getFirst(), (_p) -> _p.getSecond());
            pc.scatter("condition.to", toData, (_p) -> _p.getFirst(), (_p) -> _p.getSecond());
            pc.setTitle(title);

            BufferedImage bim = ImageUtils.vconcat(getSignalPlot(exp, pc), diffFC, diffFCFullPlot);
            String info = "detail for " + exp.combinedFeatureName;

            if(labelFunc != null) {
                info += String.format(" (%s)", labelFunc.apply(exp.combinedFeatureName));
            }
            ImageUtils.showImage(info, bim, false);

        });
        */
        return mjFrame;
    }

    public BufferedImage getDetailImage(PlotCreator pc, DiffExpResult exp) {

        PlotCreator.BoxPlotBuilder bpb = pc.buildBoxPlot();

        Vector<ErrorEstimationDistribution> perFeatureErrDistribs = exp.perFeatureScaledDistributions;
        for (int i = 0; i < exp.featureNames.size(); i++) {
            String fn = exp.featureNames.get(i);
            double fc = exp.perFeatureDistribEstimates.get(i).getFirst();
            ErrorEstimationDistribution ferr = perFeatureErrDistribs.get(i);
            bpb.addBox(String.format("%s (%.2f)", fn, fc), new SparseCumulativeDistribution(ferr).quantiles);
        }
        bpb.addBox("combined", new SparseCumulativeDistribution(exp.combinedEmpiricalFoldChangeDistrib).quantiles);
        pc.setLabels("", "log2fc", "topright");
        String title = String.format("%s FC: %.2f fdr: %g (%d features)", exp.combinedFeatureName, exp.estimatedFC, exp.fcEstimateFDR, exp.featureNames.size());
        pc.setTitle(title);
        BufferedImage diffFC = bpb.plot();

        new SparseCumulativeDistribution(exp.combinedEmpiricalFoldChangeDistrib).drawLine(pc, "est.FC distrib");

        pc.abline("", exp.estimatedFC, null, null, null);
        String title2 = String.format("%s FC: %.2f pval: %g (%d features)", exp.combinedFeatureName, exp.estimatedFC, exp.fcEstimateFDR, exp.featureNames.size());
        pc.setTitle(title2);
        BufferedImage diffFCFullPlot = pc.getImage();


        String info = "detail for " + exp.combinedFeatureName;

        if (labelFunc != null) {
            info += String.format(" (%s)", labelFunc.apply(exp.combinedFeatureName));
        }

        return ImageUtils.vconcat(getSignalPlot(exp, pc), diffFC, diffFCFullPlot);
    }

    public BufferedImage getSignalPlot(DiffExpResult exp, PlotCreator pc) {
        //plot features
        Vector<UPair<Double>> fromData = new Vector<>();
        Vector<UPair<Double>> toData = new Vector<>();
        int xCoord = 0;
        Vector<Pair<Double, String>> xtics = new Vector<>();

        for(int i=0; i<exp.featureNames.size(); i++) {
            double x1 = xCoord++;
            double x2 = xCoord++;


            Vector<Double> n1 = from.getNormed(exp.featureNames.get(i)).replicates;
            Vector<Double> n2 = to.getNormed(exp.featureNames.get(i)).replicates;

            fromData.addAll(map(n1, (_d) -> UPair.createU(x1, _d)));
            toData.addAll(map(n2, (_d) -> UPair.createU(x2, _d)));

            xtics.add(Pair.create(0.5 * ( x1 + x2 ), exp.featureNames.get(i)));
        }

        String title = String.format("%s FC: %.2f fdr: %g (%d features)", exp.combinedFeatureName, exp.estimatedFC, exp.fdr, exp.featureNames.size());
        pc.setXTics(xtics);
        pc.scatter("condition.from", fromData, (_p) -> _p.getFirst(), (_p) -> _p.getSecond());
        pc.scatter("condition.to", toData, (_p) -> _p.getFirst(), (_p) -> _p.getSecond());
        pc.setTitle(title);

        return pc.getImage();
    }
}
