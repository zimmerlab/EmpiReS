package nlEmpiRe.test.rnaseq;


import lmu.utils.*;
import lmu.utils.fdr.PerformanceResult;
import lmu.utils.plotting.*;
import lmu.utils.swing.PagedDataTable;
import lmu.utils.tuple.Tuple3;
import nlEmpiRe.*;
import nlEmpiRe.input.DSType;
import nlEmpiRe.input.EQClassInput;
import nlEmpiRe.input.ExperimentDescriptor;
import nlEmpiRe.input.RNASeqSplicingInfo;
import lmu.utils.plotting.CachedPlotCreator;
import nlEmpiRe.rnaseq.*;
import nlEmpiRe.rnaseq.simulation.PositionBiasFactory;
import nlEmpiRe.rnaseq.simulation.SplicingSimulation;
import nlEmpiRe.rnaseq.simulation.TranscriptSimulation;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.junit.jupiter.api.Test;

import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.*;
import java.util.function.Consumer;
import java.util.function.Function;

import static lmu.utils.IteratorUtils.rangev;
import static lmu.utils.ObjectGetter.*;
import static org.junit.jupiter.api.Assertions.assertTrue;
import org.apache.logging.log4j.Logger;


class TranscriptSimulationTest {

    static final File SIMULATION_BASE = new File("/data/EmpiRe/SPLICING/RNASEQ/");
    static final File TRSIMULATION = new File(SIMULATION_BASE, "texprs.txt");
    static final File TRSIMULATION_OUTFILE = new File(SIMULATION_BASE, "gene_infos.tsv");
    static final File TRUE_SPLICING = new File(SIMULATION_BASE, "trueSplic.txt");
    static final File TRUE_DIFFEXP = new File(SIMULATION_BASE, "trueDiffexp.txt");
    static final File HUMAN_GTF = new File(SIMULATION_BASE, "Homo_sapiens.GRCh37.75.gtf");
    static IsoformRegionGetter HUMAN_IRG = null;

    public static IsoformRegionGetter getHumanAnnot() {
        if(HUMAN_IRG != null)
            return HUMAN_IRG;

        return HUMAN_IRG = new GFFBasedIsoformRegionGetter(HUMAN_GTF, null, null);
    }

    @Test
    void singleTranscriptSimulationTest() {
        IsoformRegionGetter irg = getHumanAnnot();
        HashMap<MultiIsoformRegion, nlEmpiRe.rnaseq.simulation.TranscriptSimulation> simulMap = nlEmpiRe.rnaseq.simulation.TranscriptSimulation.groupTranscriptsPerGene(toVector(Pair.create("ENST00000334267", 1000)), irg);

        PositionBiasFactory biasFactory = new PositionBiasFactory();

        Consumer<nlEmpiRe.rnaseq.simulation.SimulatedRead> consumer = (_sr) -> {
            System.out.printf("simul: %s trset: %s\n", _sr.transcriptId, _sr.mapsToTranscripts);
            assert(_sr.mapsToTranscripts.contains(_sr.transcriptId));
        };

        for(nlEmpiRe.rnaseq.simulation.TranscriptSimulation simulation : simulMap.values()) {
            simulation.simulate(biasFactory, consumer);
        }
    }

    static RegionVector convertToInner(RegionVector merged, boolean negativeStrand, RegionVector toconvert) {
        Vector<Region1D> convertedRegions = new Vector<>();
        for(Region1D r : toconvert.getRegions()) {
            int start = merged.toLocalPosition(r.getX1(), negativeStrand) + ((negativeStrand) ? 1 : 0);
            int end = merged.toLocalPosition(r.getX2() - 1 , negativeStrand) + ((negativeStrand) ? 0 : 1);
            int x1 = Math.min(start, end);
            int x2 = Math.max(start, end);
            convertedRegions.add(new Region1D(x1, x2));
        }
        return new RegionVector(toSortedVector(convertedRegions, (_x) -> _x.getX1(), true), true);
    }


    @Test
    void transcriptPairSimulationTest() {
        IsoformRegionGetter irg = getHumanAnnot();
        HashMap<MultiIsoformRegion, nlEmpiRe.rnaseq.simulation.TranscriptSimulation> simulMap = nlEmpiRe.rnaseq.simulation.TranscriptSimulation.groupTranscriptsPerGene(toVector(Pair.create("ENST00000334267", 1000)), irg);

        nlEmpiRe.rnaseq.simulation.SplicingSimulation splicingSimulation = new nlEmpiRe.rnaseq.simulation.SplicingSimulation(irg);
        splicingSimulation.setConditions(toVector("1", "2"));

        Vector<Tuple3<String, Pair<String, int[][]>, Pair<String, int[][]>>> toSimulate = toVector(
                Tuple3.create("ENSG00000103194",
                        Pair.create("ENST00000570191", new int[][]{{449, 408,588},{924, 968, 665}}),
                        Pair.create("ENST00000570053", new int[][]{{449, 408,588},{229, 234, 176}})
                )

        );

        for(Tuple3<String, Pair<String, int[][]>, Pair<String, int[][]>> gene2tpairWithCounts : toSimulate) {
            String geneId = gene2tpairWithCounts.get0();
            MultiIsoformRegion gene = irg.getRegionById(geneId);

            Pair<String, int[][]> tr1ToCounts = gene2tpairWithCounts.get1();
            Pair<String, int[][]> tr2ToCounts = gene2tpairWithCounts.get2();

            for (int cond = 0; cond < 2; cond++) {
                int[] tr1counts = tr1ToCounts.getSecond()[cond];
                int[] tr2counts = tr2ToCounts.getSecond()[cond];
                if (tr1counts.length != tr2counts.length)
                    throw new FRuntimeException("invalid input, same number of replicates expected for both transcripts (%s cond: %d)", geneId, cond + 1);

                for (int repIdx = 0; repIdx < tr1counts.length; repIdx++) {

                    HashMap<String, Integer> counts = new HashMap<>();
                    counts.put(tr1ToCounts.getFirst(), tr1counts[repIdx]);
                    counts.put(tr2ToCounts.getFirst(), tr2counts[repIdx]);
                    System.out.printf("cond: %d rep: %d counts : %s\n", cond, repIdx, counts);
                    nlEmpiRe.rnaseq.simulation.TranscriptSimulation ts = new nlEmpiRe.rnaseq.simulation.TranscriptSimulation(gene, counts);
                    splicingSimulation.addSimulation(cond, repIdx, ts);
                }
            }
        }

        for(String g : splicingSimulation.getSimulatedGenes()) {
            SimulatedGeneWithCoverage sg = new SimulatedGeneWithCoverage(irg.getRegionById(g), splicingSimulation, 5);
            sg.show();
        }

    }



    @Test
    void testResolution() {
        String[][] gene2tr = new String[][]{
                {"ENSG00000168385", "ENST00000391973", "ENST00000421717"},
                {"ENSG00000144395", "ENST00000389175", "ENST00000431807"},
                {"ENSG00000168397", "ENST00000400771", "ENST00000402096"},
                {"ENSG00000181381", "ENST00000505863", "ENST00000505890"},
                {"ENSG00000132376", "ENST00000542125", "ENST00000574561"},
                {"ENSG00000144366", "ENST00000359135", "ENST00000409805"},
                {"ENSG00000181396", "ENST00000329197", "ENST00000583897"},
                {"ENSG00000119383", "ENST00000347048", "ENST00000357197"},
                {"ENSG00000158987", "ENST00000307984", "ENST00000509018"},
                {"ENSG00000156313", "ENST00000339363", "ENST00000494707"},
                {"ENSG00000132361", "ENST00000435359", "ENST00000574426"}
        };

        IsoformRegionGetter irg = getHumanAnnot();

        NormalDistribution nd = new NormalDistribution(200, 60);

        PositionBiasFactory biasFactory = new PositionBiasFactory();

        for(int i=0; i<gene2tr.length; i++) {
            MultiIsoformRegion gene = irg.getRegionById(gene2tr[i][0]);


            HashMap<String, Integer> tr2count = new HashMap<>();
            for(int trIdx=1; trIdx<=2; trIdx++) {
                tr2count.put(gene2tr[i][trIdx], (trIdx == 1) ? 100 : 30 );
            }
            nlEmpiRe.rnaseq.simulation.TranscriptSimulation trs = new nlEmpiRe.rnaseq.simulation.TranscriptSimulation(gene, tr2count);
            HashMap<Tuple, Double> eqClasses = trs.simulateTrSetCounts(biasFactory);

            ReducedTranscriptPresentation rtp = new ReducedTranscriptPresentation(eqClasses);

            Vector<String> trKeys = toSortedVector(tr2count.keySet(), true);

            String tr1 = trKeys.get(0);
            String tr2 = trKeys.get(1);

            boolean isTested = rtp.isTestable(tr1, tr2);


            assertTrue( isTested,
                    String.format("tr: %s are not tested info: %s cluster: %s",
                            trKeys, rtp.getPartialInfos(tr1, tr2),
                            rtp.getTrCluster(tr1))
                    );


        }
    }

    static void multiBenchMark(SimulationConfiguration config, File trcounts, File trueLabels) {

        Logger log = LogConfig.getLogger();
        double RECALL_NORM = 1.0 / FileUtils.readSet(trueLabels).size();
        Vector<Integer> MAXTRS = toVector(5, 7, 10);
        Vector<Double> pseudos = toVector(1.0, 2.0, 3.0, 5.0);
        Vector<Double> factors = toVector(1.0, 2.0, 3.0);
        int NREP = 1;

        class Setup {
            int maxtr;
            double pseudo;
            double factor;
            Vector<PerformanceResult> performanceResults = new Vector<>();

            Setup(int maxtr, double pseudo, double factor) {
                this.maxtr = maxtr; this.pseudo = pseudo; this.factor = factor;
            }

            String getName() {
                return String.format("%d,%.0f,%.0f", maxtr, pseudo, factor);
            }
        }

        int numtests = MAXTRS.size() * pseudos.size() * factors.size() * NREP;
        Vector<Setup> tests = new Vector<>();
        int testidx = 0;

        for(int maxtr : MAXTRS) {
            for(double pseudo : pseudos) {
                for(double factor : factors) {
                    Setup s = new Setup(maxtr, pseudo, factor);

                    for(int i=0; i<NREP; i++) {
                        PerformanceResult pr = benchmark(false, config, trcounts, maxtr, trueLabels, pseudo, factor);
                        s.performanceResults.add(pr);
                        testidx++;
                        log.info("splictest: %d/%d : %s :: %s", testidx, numtests, s.getName(), pr);

                    }
                    tests.add(s);
                }
            }
        }

        PlotCreator pc = CachedPlotCreator.getPlotCreator();
        pc.setTitle("format: maxtr,pseudo,factor");
        PlotCreator.BoxPlotBuilder bpb = pc.buildBoxPlot();
        for(Setup s : tests) {
            bpb.addBox(s.getName(), "prec", NumUtils.getNumInfo(s.performanceResults, (_p) -> _p.PREC).quantiles, Color.RED, true);
            bpb.addBox(s.getName(), "rec", NumUtils.getNumInfo(s.performanceResults, (_p) -> _p.RECALL).quantiles, Color.BLUE, false);
            bpb.addBox(s.getName(), "recAbs", NumUtils.getNumInfo(s.performanceResults, (_p) -> _p.pred_true *RECALL_NORM ).quantiles, Color.GREEN, false);
        }


        pc.setLabels("", "prec/rec", "bottomright");
        BufferedImage bim = bpb.plot("bottomright");
        ImageUtils.showImage("simuls", bim, true);

    }

    static PerformanceResult benchmark(boolean doPlot, SimulationConfiguration config, File trcounts, int MAXTRCOUNT, File trueLabels, double RAW_PSEUDO, double MULTIPLICATIVE_FACTOR_ON_COUNTS) {
        Logger log = LogConfig.getLogger();
        log.info("init simulations");
        nlEmpiRe.rnaseq.simulation.SplicingSimulation simulation = new nlEmpiRe.rnaseq.simulation.SplicingSimulation(config, trcounts, MULTIPLICATIVE_FACTOR_ON_COUNTS);

        log.info("start simulations");
        long t1 = System.currentTimeMillis();
        RNASeqSplicingInfo rsi = simulation.simulate(MAXTRCOUNT, 1.0, 3);
        long t2 = System.currentTimeMillis();
        long simul_time = t2 - t1;
        log.info("setup splicing tests took: %.2f sec", simul_time / 1000.0);
        SplicingTest test = new SplicingTest(rsi, new AutoBackGroundContextProvider());
        Vector<String> conditions = test.getConditions();

        log.info("calc splicing conditions: %s", conditions);
        t1 = System.currentTimeMillis();
        UPair<Vector<DoubleDiffResult>> splicingPair = test.getDifferentialAlternativeSplicing(conditions.get(0), conditions.get(1));
        Vector<DoubleDiffResult> splicing = splicingPair.getFirst();
        t2 = System.currentTimeMillis();
        long DAS_test_time = t2 - t1;
        HashSet<String> trueGenes = FileUtils.readSet(trueLabels);
        Function<DoubleDiffResult, Boolean> trueLabeller = (_dr) -> trueGenes.contains(_dr.testName.split("\\.")[0]);

        PerformanceResult pr = new PerformanceResult("nlEmpire", splicing, (_dr) -> _dr.pval, trueLabeller,
                    false, (_dr) -> _dr.fdr <= 0.05, null, false);

        log.info("DAS test time: %.2f sec. simul time: %.2f sec", DAS_test_time / 1000.0, simul_time / 1000.0);
        System.out.printf("pr: %s\n", pr);
        System.out.printf("real recall: %2.2f\n", filteredSize(splicing, (_s) -> trueLabeller.apply(_s) && _s.fdr <= 0.05) / (0.0 + trueGenes.size()));

        if(doPlot) {
            NumUtils.sort(splicing, (_p) -> _p.pval, false);

            DataTable dt = BenchmarkGene.toTable(map(splicing, (_s) -> new BenchmarkGene(_s, trueGenes, simulation, rsi)));
            PagedDataTable.MJFrame mjf = PagedDataTable.getDetailViewFrame(dt, null, null);


            mjf.addMenu("details", og ->
            {
                BenchmarkGene bg = (BenchmarkGene)og.getInData();
                DoubleDiffResult ddr = bg.ddr;
                PlotCreator.BoxPlotBuilder bpb = ddr.drawSignalDistributionBoxPlots(CachedPlotCreator.getPlotCreator(), false);
                BufferedImage diffSignal = bpb.plot();

                bpb = ddr.drawDiffFCPlots(CachedPlotCreator.getPlotCreator(), false);
                BufferedImage diffFC = bpb.plot();

                ImageUtils.showImage(bg.geneId, ImageUtils.concat(diffSignal, diffFC), false);

            });
            mjf.setVisible(true);
        }


        return pr;

    }




    static SimulationConfiguration getDefaultConfig() {
        return new SimulationConfiguration(getHumanAnnot());
    }
    @Test
    void multiTranscriptSimulationTest() {
        multiTranscriptSimulationTest(false, new SimulationConfiguration(), TRSIMULATION, 5);

    }


    static void simulateReads(File trcounts, SimulationConfiguration simulationConfiguration) {
        Logger log = LogConfig.getLogger();
        nlEmpiRe.rnaseq.simulation.SplicingSimulation simulation = new nlEmpiRe.rnaseq.simulation.SplicingSimulation(simulationConfiguration, trcounts);

        Vector<Tuple3<String, String, Integer>> gene2tr2length = new Vector<>();
        for (String geneId : simulation.getSimulatedGenes()) {
            MultiIsoformRegion gene = simulationConfiguration.getAnnot().getRegionById(geneId);
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

        for (String geneId : simulation.getSimulatedGenes()) {
            MultiIsoformRegion gene = simulationConfiguration.getAnnot().getRegionById(geneId);
            simulation.simulate(geneId, false, simulationConfiguration.getReadInReplicateConsumer());
        }
        simulationConfiguration.fastQGenerator.close();
    }

    static class TestGeneInfo {

        MultiIsoformRegion gene;

        int numTranscripts = 0;
        Vector<String> transcripts;
        String errInfo = "";
        boolean isTrue;
        Vector<Vector<Integer>> tr1counts;
        Vector<Vector<Integer>> tr2counts = new Vector<>();

        int tr1MaxCount = 0;
        int tr2MaxCount = 0;
        Vector<Tuple3<String, Vector<String>, Vector<String>>> testCombis;

        int getMaxCount(Vector<Vector<Integer>> cond2rep2counts) {
            return NumUtils.max(map(cond2rep2counts, (_r) -> NumUtils.max(_r)));
        }
        public TestGeneInfo(MultiIsoformRegion gene, nlEmpiRe.rnaseq.simulation.SplicingSimulation simulation, String cond1, String cond2, RNASeqSplicingInfo rsi, Set<String> trues) {
            this.gene = gene;


            transcripts = toVector(simulation.getSimulatedTranscripts(gene.id));
            if(transcripts.size() == 0)
                return;

            tr1counts = simulation.getSimulatedCounts(gene.id, transcripts.get(0));
            isTrue = trues.contains(gene.id);
            numTranscripts = transcripts.size();
            tr1MaxCount = getMaxCount(tr1counts);
            if(numTranscripts < 2)
                return;

            tr2counts = simulation.getSimulatedCounts(gene.id, transcripts.get(1));
            tr2MaxCount = getMaxCount(tr2counts);
            StringBuffer errStr = new StringBuffer();
            testCombis = rsi.getTrPairTestFeatureCombiTests(gene.id, transcripts.get(0), transcripts.get(1), cond1, cond2, errStr);

            errInfo = errStr.toString();
        }

        public Color getColor(Set<String> trSet) {
            boolean gottr1 = trSet.contains(transcripts.get(0));
            if(transcripts.size() == 1)
                return (gottr1) ? Color.BLUE : Color.BLACK;

            boolean gottr2 = trSet.contains(transcripts.get(1));

            if(gottr1 && gottr2)
                return Color.MAGENTA;

            if(!gottr1 && !gottr2)
                return Color.BLACK;

            return (gottr1) ? Color.BLUE : Color.RED;
        }
        static public DataTable toTable(Vector<TestGeneInfo> testGeneInfos) {

            return DataTable.buildTable(
                    testGeneInfos,
                DataTable.buildHeader("gene", (TestGeneInfo tgi) -> tgi.gene.id)
                        .add("strand", t -> GenomicUtils.getStrand(t.gene.strand))

                    .add("isTrue", t -> t.isTrue)
                    .add("num.annot.trans", t -> t.gene.isoforms.size())
                    .add("nTranscripts", t -> t.numTranscripts)
                    .add("tr1maxcount", t -> t.tr1MaxCount)
                    .add("tr2maxcount", t -> t.tr2MaxCount)
                    .add("tr1counts", t -> t.tr1counts)
                    .add("tr2counts", t -> t.tr2counts)

                        .add("err", t -> t.errInfo)
                    .add("num.testcombis", t -> (t.testCombis == null) ? 0 : t.testCombis.size())
            );
        }

    }

    static void detailTestInfos(SimulationConfiguration config, File trcounts, int MAXTRCOUNT, Set<String> trues, File eqClassInputFile) {
        Logger log = LogConfig.getLogger();
        log.info("init simulations");
        nlEmpiRe.rnaseq.simulation.SplicingSimulation simulation = new nlEmpiRe.rnaseq.simulation.SplicingSimulation(config, trcounts, 1.0);



        int mincountInMaxCountCondition = 10;
        RNASeqSplicingInfo rsi = null;

        if(eqClassInputFile != null && eqClassInputFile.exists()) {
            log.info("read in EQ-classes from %s", eqClassInputFile.getAbsolutePath());
            long t1 = System.currentTimeMillis();

            HashMap<String, Vector<String>> condition2replicatenames = new HashMap<>();
            Vector<String> replicatelist = new Vector<>();
            for(int i=0; i<simulation.getNumConditions(); i++) {
                nlEmpiRe.rnaseq.simulation.SimulatedSplicingCondition ssc = simulation.getCondition(i);
                int idx = i+1;
                Vector<String> condreps = map(rangev(ssc.replicates.size()), (_i) -> ssc.condition +".rep"+ (_i+1));
                condition2replicatenames.put(ssc.condition, condreps);
                replicatelist.addAll(condreps);
            }
            System.out.printf("cond2reps: %s replist: %s\n", condition2replicatenames, replicatelist);
            rsi = new RNASeqSplicingInfo(10, 1.0,  condition2replicatenames);

            EQClassInput eqClassInput = new EQClassInput(rsi, new ExperimentDescriptor(condition2replicatenames, replicatelist, trues), DSType.READS, mincountInMaxCountCondition);
            eqClassInput.read(eqClassInputFile);
            long t2 = System.currentTimeMillis();
            long eqread_time = t2 - t1;
            log.info("reading in eq-class counts took: %.2f sec", eqread_time / 1000.0);

        }
        else {
            log.info("start simulations");
            long t1 = System.currentTimeMillis();

            rsi = simulation.simulate(MAXTRCOUNT, 1.0, mincountInMaxCountCondition);
            long t2 = System.currentTimeMillis();
            long simul_time = t2 - t1;
            log.info("setup splicing tests took: %.2f sec", simul_time / 1000.0);

        }


        SplicingTest test = new SplicingTest(rsi, new AutoBackGroundContextProvider());
        Vector<String> conditions = test.getConditions();

        String cond1 = conditions.get(0);
        String cond2 = conditions.get(1);

        NormalizedReplicateSet rs1 = test.getNormalized(cond1);
        NormalizedReplicateSet rs2 = test.getNormalized(cond2);


        DoubleDiffManager diffManager = new EmpiRe().getDoubleDiffManager(rs1, rs2);

        DiffExpManager diffExpManager = diffManager.getDiffExpManager();
        Vector<TestGeneInfo> testGeneInfos = new Vector<>();

        for(String g : rsi.getGenesForSplicingTest()) {
            TestGeneInfo tgi = new TestGeneInfo(config.getAnnot().getRegionById(g), simulation, cond1, cond2, rsi, trues);
            if(tgi.transcripts.size() == 0)
                continue;
            testGeneInfos.add(tgi);
        }

        Function<Vector<Double>, String> vfmt = (_v) -> StringUtils.joinObjects(",", _v, (_d) -> String.format("%.2f", _d));

        Function<Double, Double> logtransform = (_d) -> (Math.abs(_d) < 0.0001) ? Double.NaN : NumUtils.logN(_d, 2.0);

        PagedDataTable.getDetailViewFrame(TestGeneInfo.toTable(testGeneInfos),
                null, null)
                .addMenu("explain eqClass-reduction", (ObjectGetter og) ->
                {
                    TestGeneInfo tgi = (TestGeneInfo)og.getInData();
                    nlEmpiRe.rnaseq.simulation.SplicingSimulatedGene sgi = simulation.simulate(tgi.gene.id, false);
                    HashMap<Tuple, Double> eqClasses = new HashMap<>();

                    for (Vector<HashMap<Tuple, Double>> v : sgi.condition2replicate2equivianceClasses) {
                        for (HashMap<Tuple, Double> m : v) {
                            apply(m.entrySet(), (_e) -> MapBuilder.update(eqClasses, _e.getKey(), _e.getValue()));
                        }
                    }

                    ReducedTranscriptPresentationVisualization rtpv = new ReducedTranscriptPresentationVisualization();
                    ReducedTranscriptPresentation rtp = new ReducedTranscriptPresentation(eqClasses, 2, rtpv);
                    rtpv.showSteps();
                    return;
                })
                .addMenu("show eq-classes", (ObjectGetter og) -> {
                    TestGeneInfo tgi = (TestGeneInfo)og.getInData();
                    SimulatedGeneWithCoverage sgwc = new SimulatedGeneWithCoverage(tgi.gene, simulation, 10);
                    sgwc.show();
                })
                .addMenu("show eq-classes Image", (ObjectGetter og) -> {
                    TestGeneInfo tgi = (TestGeneInfo)og.getInData();
                    SimulatedGeneWithCoverage sgwc = new SimulatedGeneWithCoverage(tgi.gene, simulation, 10);
                    ImageUtils.showImage("rtpc image", sgwc.getImage(), false);
                })
                .addMenu("show diff signals", (ObjectGetter og) -> {
                            TestGeneInfo tgi = (TestGeneInfo) og.getInData();
                            nlEmpiRe.rnaseq.simulation.SplicingSimulatedGene ssg = simulation.simulate(tgi.gene.id, false);
                            HashMap<Tuple, Double> eqClasses = new HashMap<>();

                            for (Vector<HashMap<Tuple, Double>> v : ssg.condition2replicate2equivianceClasses) {
                                for (HashMap<Tuple, Double> m : v) {
                                    apply(m.entrySet(), (_e) -> MapBuilder.update(eqClasses, _e.getKey(), _e.getValue()));
                                }
                            }

                            System.out.printf("eqclasses: %s\n", eqClasses);

                            Vector<BufferedImage> bims = new Vector<>();

                            for (int target : toVector(10, 8, 6, 4, 2)) {
                                ReducedTranscriptPresentation rtp = new ReducedTranscriptPresentation(eqClasses, target);
                                System.out.printf("rtp target: %d got %d/%d\n", target, rtp.getNumRestrictedEQClasses(), eqClasses.size());
                                Vector<Map<Tuple,Double>> counts_cond1 = map(ssg.condition2replicate2equivianceClasses.get(0), (_m) -> rtp.reduce(_m));
                                Vector<Map<Tuple,Double>> counts_cond2 = map(ssg.condition2replicate2equivianceClasses.get(1), (_m) -> rtp.reduce(_m));

                                System.out.printf("reduced: %s\n", counts_cond1);

                                Vector<Tuple3<SingleFeatureDiffExp, Color, UPair<Vector<Double>>>> features = new Vector<>();
                                for(Tuple t : rtp.getRestrictedEQClasses()) {
                                    Vector<Double> raw1 = map(counts_cond1, (_m) -> _m.getOrDefault(t, 0.0));
                                    Vector<Double> raw2 = map(counts_cond2, (_m) -> _m.getOrDefault(t, 0.0));

                                    System.out.printf("\t%s %s %s\n", t, raw1, raw2);

                                    Vector<Double> reps1 = map(raw1, (_d) -> logtransform.apply(_d));
                                    Vector<Double> reps2 = map(raw2, (_d) -> logtransform.apply(_d));

                                    double min1 = NumUtils.min(raw1, (_d) -> _d, 0.0);
                                    double min2 = NumUtils.min(raw2, (_d) -> _d, 0.0);

                                    System.out.printf("MAX: %.2f check: %d\n", Math.max(min1, min2), mincountInMaxCountCondition);
                                    if(Math.max(min1, min2) < mincountInMaxCountCondition)
                                        continue;

                                    Set<String> trs = mapToSet(rangev(t.cardinality()), (_i) -> t.getAsString(_i));


                                    SingleFeatureDiffExp singleFeature = new SingleFeatureDiffExp(""+t, diffExpManager, reps1, reps2);
                                    if(!singleFeature.useable)
                                        continue;
                                    features.add(Tuple3.create(singleFeature, tgi.getColor(trs), UPair.createU(raw1, raw2)));

                                }

                                if(features.size() == 0)
                                    continue;

                                NumUtils.sort(features, (_p) -> _p.getFirst().getFC());
                                PlotCreator pc = CachedPlotCreator.getPlotCreator();
                                PlotCreator.BoxPlotBuilder bpb = pc.buildBoxPlot();
                                for(Tuple3<SingleFeatureDiffExp, Color, UPair<Vector<Double>>> p : features) {
                                    bpb.addBox(p.get0().feature, new SparseCumulativeDistribution(p.get0().scaledDistribution).quantiles, p.get1());
                                }

                                pc.setTitle("reduced to max " + target + " transcripts");
                                pc.setLabels("", "log2FC", null);

                                BufferedImage bim1 = bpb.plot();

                                double x = 0.0;
                                for(Tuple3<SingleFeatureDiffExp, Color, UPair<Vector<Double>>> p : features) {
                                    double X = x++;
                                    pc.scatter("", p.get2().getFirst(), (_d) -> X, (_d) -> _d).setColor(Color.DARK_GRAY);
                                    double X2 = x++;
                                    pc.scatter("", p.get2().getSecond(), (_d) -> X2, (_d) -> _d).setColor(p.get1());

                                }
                                pc.setLabels("features", "raw counts", null);
                                pc.setTitle("reduced to max " + target + " transcripts");

                                BufferedImage bim2 = pc.getImage();

                                bims.add(ImageUtils.concat(bim2, bim1));
                            }
                            if(bims.size() == 0) {
                                ImageUtils.showImage(tgi.gene.id, ImageUtils.getTextImage("no feature testable").getImage(), false);
                                return;
                            }

                            ImageUtils.showImage(tgi.gene.id, ImageUtils.vconcat(bims), false);
                        }
                    )
                .addMenu("do tests",
                        (ObjectGetter og) ->
                        {
                            TestGeneInfo tgi = (TestGeneInfo)og.getInData();
                            if(tgi.testCombis == null)
                                return;

                            for(Tuple3<String, Vector<String>, Vector<String>> t : tgi.testCombis) {
                                System.out.printf("nexts to test: %s\n", t);
                                DiffExpResult d1 = new DiffExpResult(diffExpManager, "cs1", t.get1());
                                DiffExpResult d2 = new DiffExpResult(diffExpManager, "cs2", t.get2());


                                ErrorEstimationDistribution e1 = d1.combinedEmpiricalFoldChangeDistrib;
                                ErrorEstimationDistribution e2 = d2.combinedEmpiricalFoldChangeDistrib;

                                PlotCreator pc = CachedPlotCreator.getPlotCreator();



                                Vector<Pair<String, ErrorEstimationDistribution>> toplot = toVector(
                                        Pair.create("cs1", e1), Pair.create("cs2", e2)
                                );

                                for (Pair<String, ErrorEstimationDistribution> p : toplot) {
                                    if (p.getFirst() == null)
                                        continue;

                                    new SparseCumulativeDistribution(p.getSecond()).drawLine(pc, p.getFirst());
                                }
                                BufferedImage distribs = pc.getImage();


                                PlotCreator.BoxPlotBuilder bpb = pc.buildBoxPlot();

                                for (Pair<String, ErrorEstimationDistribution> p : toplot) {
                                    if (p.getFirst() == null)
                                        continue;

                                    new SparseCumulativeDistribution(p.getSecond()).drawBox(bpb, p.getFirst());
                                }

                                pc.setLabels("", "log2fc", null);
                                BufferedImage boxes = bpb.plot();


                                ImageUtils.TextImage tim = new ImageUtils.TextImage();
                                tim.add(t.get0());
                                tim.add("cs1:" + t.get1());
                                double x = 0.0;
                                Vector<BufferedImage> signalBims = new Vector<>();
                                for (int fi = 1; fi <= 2; fi++) {
                                    Vector<String> fs = (Vector<String>) t.get(fi);

                                    for (String f : fs) {
                                        for (int c = 0; c < 2; c++) {
                                            Vector<Double> vals = (c == 0) ? rs1.getNormed(f).nonNanValues : rs2.getNormed(f).nonNanValues;
                                            double X = x++;
                                            pc.scatter(f, vals, (_x) -> X, (_x) -> _x).setColor((c == 0) ? Color.GREEN : Color.RED);
                                        }
                                        tim.add("cs"+fi+"::"+f + " c1" + vfmt.apply(rs1.getNormed(f).nonNanValues) + " c2: " + vfmt.apply(rs2.getNormed(f).nonNanValues));
                                    }

                                    pc.setLabels("", "log2vals", null);
                                    pc.setTitle("cond" + fi);
                                    signalBims.add(pc.getImage());

                                }

                                BufferedImage header = tim.getImage();
                                BufferedImage combined = ImageUtils.vconcat(header, ImageUtils.concat(signalBims),
                                        distribs, boxes);

                                ImageUtils.showImage(tgi.gene.id, combined, false);
                            }
                        }
                        );

    }

    static void eqReductionTest(SimulationConfiguration config, File trcounts) {
        nlEmpiRe.rnaseq.simulation.SplicingSimulation simulation = new nlEmpiRe.rnaseq.simulation.SplicingSimulation(config, trcounts);

        class ReduceInfo {
            int target = -1;
            Vector<Integer> sizes = new Vector<>();
            int numOk = 0;

            public ReduceInfo(int t) {
                target = t;
            }
        }
        //Vector<ReduceInfo> infos = map(toVector(-1, 20, 10, 8, 6, 5, 4, 3 , 2), (_i) -> new ReduceInfo(_i));
        Vector<ReduceInfo> infos = map(toVector(-1, 10, 5,  2), (_i) -> new ReduceInfo(_i));

        Vector<Integer> numTranscriptsWithCounts = new Vector<>();

        for(String g : simulation.getSimulatedGenes()) {
            Vector<String> trs = toVector(simulation.getSimulatedTranscripts(g));

            if(trs.size() < 2)
                continue;




            nlEmpiRe.rnaseq.simulation.SplicingSimulatedGene ssg = simulation.simulate(g, false);

            HashMap<Tuple, Double> eqClasses = new HashMap<>();

            for(Vector<HashMap<Tuple, Double>> v : ssg.condition2replicate2equivianceClasses) {
                for(HashMap<Tuple, Double> m : v) {
                    apply(m.entrySet(), (_e) -> MapBuilder.update(eqClasses,_e.getKey(), _e.getValue()));
                }
            }

            HashSet<String> trsWithCounts = new HashSet<>();
            for(Tuple t : eqClasses.keySet()) {
                applyIndex(t.cardinality(), (_i) -> trsWithCounts.add(t.getAsString(_i)));
            }
            numTranscriptsWithCounts.add(trsWithCounts.size());

            for(ReduceInfo r : infos) {
                if(r.target < 0) {
                    r.sizes.add(eqClasses.size());
                    r.numOk++;
                    continue;
                }
                ReducedTranscriptPresentation rtp = new ReducedTranscriptPresentation(eqClasses, r.target);
                r.numOk += (rtp.isTestable(trs.get(0), trs.get(1))) ? 1 : 0;
                r.sizes.add(rtp.getNumRestrictedEQClasses());
            }

        }
        int ntotal = infos.get(0).numOk;
        PlotCreator pc = new PlotCreator();

        pc.cumhist("#transcripts with counts", numTranscriptsWithCounts, numTranscriptsWithCounts.size(), false, false);

        for(ReduceInfo ri : infos) {
            pc.cumhist(String.format("target=%d lost:%.1f%%", ri.target, (100.0 * (ntotal - ri.numOk)) / ntotal), ri.sizes, ri.sizes.size(), false, false);
        }

        pc.setLabels("num features", " freq <= x",  "bottomright");
        pc.setLog(true, true);
        ImageUtils.showImage("test", pc.getImage());
    }
    static void multiTranscriptSimulationTest(boolean doPlot, SimulationConfiguration config, File trcounts, int MAXTRCOUNT) {

        nlEmpiRe.rnaseq.simulation.SplicingSimulation simulation = new SplicingSimulation(config, trcounts);

        SimulationConfiguration simulationConfiguration = getDefaultConfig();

        Vector<Integer> numTranscripts = map(simulation.getSimulatedGenes(), (_g) -> config.getAnnot().getRegionById(_g).isoforms.size());

        System.out.printf("num transcripts distrib: %s\n", NumUtils.getNumInfo(numTranscripts));


        Vector<String> toCheck = filter(simulation.getSimulatedGenes(), (_g) -> config.getAnnot().getRegionById(_g).isoforms.size() > 10);
        System.out.printf("genes with > 10 trs: %d\n", toCheck.size());


        Vector<UPair<Integer>> eq2ntr = new Vector<>();
        int numChecked = 0;


        Vector<String> errors = new Vector<>();

        Vector<SimulatedGeneInfo> simulatedGeneInfos = new Vector<>();
        Vector<Integer> trWithCounts = new Vector<>();
        Vector<Integer> eqsSizes = new Vector<>();
        Vector<Integer> restrictedEqsSizes = new Vector<>();

        for(String geneId : toCheck) {
            MultiIsoformRegion gene = config.getAnnot().getRegionById(geneId);
            nlEmpiRe.rnaseq.simulation.SplicingSimulatedGene ssg = simulation.simulate(geneId, false, simulationConfiguration.getReadInReplicateConsumer());
            HashSet<Tuple> eqs = new HashSet<>();
            HashMap<Tuple, Double> eqToCount = new HashMap<>();

            HashSet<String> simultedTranscriptIds = new HashSet<>();

            SimulatedGeneInfo simulatedGeneInfo = new SimulatedGeneInfo(gene.id);

            int minTrCount = Integer.MAX_VALUE;
            int maxTrCount = 0;
            for(int ci =0 ; ci < simulation.getNumConditions(); ci++) {
                String condname = "cond"+(ci+1);
                Vector<nlEmpiRe.rnaseq.simulation.TranscriptSimulation> repSimulation  = map(simulation.getCondition(ci).replicates, (_r) -> _r.get(gene));
                for(TranscriptSimulation trS : repSimulation) {
                    simultedTranscriptIds.addAll(trS.getSimulatedTranscriptIds());

                    for(String trid : trS.getSimulatedTranscriptIds()) {
                        simulatedGeneInfo.update("cond"+(ci+1), trid, trS.getCountToSimulate(trid));
                    }


                }
                for(HashMap<Tuple, Double> m : ssg.condition2replicate2equivianceClasses.get(ci)) {
                    eqs.addAll(m.keySet());
                    for(Map.Entry<Tuple, Double> e : m.entrySet()) {
                        MapBuilder.update(eqToCount, e.getKey(), e.getValue());
                    }

                }

            }


            eq2ntr.add(UPair.createU(eqs.size(), ssg.trsWithCounts.size()));
            if(simultedTranscriptIds.size()  < 2)
                continue;

            simulatedGeneInfos.add(simulatedGeneInfo);
            simulatedGeneInfo.compile();
            simulatedGeneInfo.numTranscripts = gene.isoforms.size();

            simulatedGeneInfo.numSimulatedTranscripts = simultedTranscriptIds.size();
            simulatedGeneInfo.transcripts = toSortedVector(simultedTranscriptIds, true);
            simulatedGeneInfo.diffaltsplic = (simulationConfiguration != null && simulationConfiguration.trueSplic.contains(gene.id));
            simulatedGeneInfo.diffexp = (simulationConfiguration != null && simulationConfiguration.trueDiff.contains(gene.id));
            simulatedGeneInfo.numEqClasses = eqToCount.size();
            simulatedGeneInfo.numTranscriptsWithCounts = ssg.trsWithCounts.size();


            ReducedTranscriptPresentation rtp = new ReducedTranscriptPresentation(eqToCount, MAXTRCOUNT);
            int numClusters = rtp.getNumClusters();

            trWithCounts.add(ssg.trsWithCounts.size());
            eqsSizes.add(eqs.size());
            restrictedEqsSizes.add(rtp.getNumRestrictedEQClasses());

            simulatedGeneInfo.numReducedEqClasses = rtp.getNumRestrictedEQClasses();

            //System.out.printf("simulated: %s , real: %d eq: %d restricted: %d clusters: %s\n", simultedTranscriptIds, ssg.trsWithCounts.size(), eqs.size(), rtp.restrictedEQCounts.size(),
            //        numClusters);



            assertTrue(numClusters <= MAXTRCOUNT, String.format("wanted to reduce: %d clusters, got %d", MAXTRCOUNT, numClusters));

            Vector<String> simTrVec = toSortedVector(simultedTranscriptIds, true);

            BufferPrintWriter bpw = new BufferPrintWriter();
            bpw.printf("test: %d error at gene: %s\n", numChecked, gene.id);

            boolean gotError = false;

            for(int i = 0; i < simTrVec.size(); i++) {
                String tr1 = simTrVec.get(i);
                for(int j = i + 1; j < simTrVec.size(); j++) {

                    String tr2= simTrVec.get(j);
                    boolean isTested = rtp.isTestable(tr1, tr2);
                    if(!isTested) {
                        gotError = true;

                        bpw.printf("\ttrs: %s,%s are not tested info: %s clusters: %d  %s score: %s\n\tremoves:\n\t%s",
                                tr1, tr2, rtp.getPartialInfos(tr1, tr2),

                                numClusters,
                                rtp.getTrCluster(tr1),
                                rtp.getPartialInfos(tr1, tr2),
                                rtp.getRemoves());
                    }



                }
            }
            simulatedGeneInfo.gotTrReductionError = gotError;
            if(gotError) {
                errors.add(bpw.getBuffer());
            }
            numChecked++;

        }

        int numErr = errors.size();


        BufferedImage bim = simulationConfiguration.writeInfos(simulatedGeneInfos, MAXTRCOUNT, doPlot);
        if(doPlot) {
            ImageUtils.showImage("multitr-test", bim);
        }


        System.out.printf("num genes with trs with counts > 10: %d/%d checks: %d\n", filteredSize(eq2ntr, (_p) -> _p.getSecond() > 10), eq2ntr.size(), numChecked);
        assertTrue( errors.size() == 0,
                String.format("got %d/%d errors\n%s\n", errors.size(), numChecked,
                StringUtils.joinObjects("\n\n", VectorUtils.slice(errors, 0, 4))));

    }

    static Vector<nlEmpiRe.rnaseq.simulation.SimulatedRead> simulateStartPos(HashMap<String, RegionVector> isoforms, int RL, double[] proportions, int start) {
        Vector<String> trids = toSortedVector(isoforms.keySet(), true);


        Vector<nlEmpiRe.rnaseq.simulation.SimulatedRead> simulatedReads = new Vector<>();

        PositionBiasFactory biasFactory = new PositionBiasFactory();
        double[] cumulative = new double[proportions.length];
        cumulative[0] = proportions[0];
        for(int i=1; i<proportions.length; i++) {
            cumulative[i] = proportions[i] + cumulative[i-1];
        }
        double TOTAL = cumulative[cumulative.length - 1];

        int[] hits = new int[cumulative.length];

        for(int i=0; i<1000; i++) {
            int hitIdx = Arrays.binarySearch(cumulative, Math.random() * TOTAL);
            if (hitIdx < 0) {
                hitIdx = -hitIdx - 1;
            }
            hits[hitIdx]++;
        }

        for(int i=0; i<hits.length; i++) {
            String transcriptId = trids.get(i);
            RegionVector src = isoforms.get(transcriptId);
            biasFactory.simulate("test",  true, transcriptId, isoforms, hits[i],
                    (_read) ->
                    {
                        simulatedReads.add(_read);
                    });
        }

        System.out.printf("after simulation: %s -> %s\n", trids, Arrays.toString(hits));
        return simulatedReads;
    }

    public static void biasPositionTest() {

        RegionVector tr1 = new RegionVector(toVector(new Region1D(1, 200), new Region1D(1000, 1500)), true);
        RegionVector tr2 = new RegionVector(toVector(new Region1D(1, 200), new Region1D(2500, 3000)), true);

        HashMap<String, RegionVector> isoforms = MapBuilder.build("tr1", tr1).add("tr2", tr2);

        int RL = 75;
        int start = 100;
        SimulatedGeneWithCoverage real = new SimulatedGeneWithCoverage("real-proportions", true, isoforms, isoforms.keySet(), 2,
                () -> simulateStartPos(isoforms, RL, new double[]{0.9, 0.1}, start));

        SimulatedGeneWithCoverage test = new SimulatedGeneWithCoverage("pos-based-simulated", true, isoforms, isoforms.keySet(), 2,
                () -> simulateStartPos(isoforms, RL, new double[]{0.3, 0.7}, start));

        ImageUtils.showImage("test", ImageUtils.vconcat(real.getImage(), test.getImage()));
    }

    public static void main(String[] args) {
        SimpleOptionParser cmd = new SimpleOptionParser("multitrtest", "trcounts", "gtf", "maxtrnum", "trues", "truediffs", "pseudo", "readcountfactor", "multibenchmark",
                "genome", "genomeidx", "fastqout", "reportout", "mutrate", "doplot", "checktrseqs", "details", "bamreference",
                "trpairvisu", "eqclasscounts", "eqredtest", "biaspos", "readlength", "fraglengthmean", "fraglengthsd");
        cmd.setFile("genome", "genomeidx","checktrseqs", "eqclasscounts", "biaspos", "bamreference");
        cmd.setDir("fastqout", "reportout");
        cmd.setSwitches("multitrtest", "multibenchmark", "doplot", "details", "trpairvisu", "eqredtest");
        cmd.setOptional("gtf", "trcounts", "truediffs", "trues", "genome", "genomeidx", "fastqout", "reportout", "checktrseqs", "eqclasscounts", "biaspos", "bamreference");
        cmd.setInt("maxtrnum", "readlength", "fraglengthmean", "fraglengthsd");
        cmd.setDefault("maxtrnum", "10");
        cmd.setDefault("readlength", "60");
        cmd.setFile("gtf", "trues");
        cmd.setFile("trcounts");
        cmd.setDouble("pseudo", "readcountfactor", "mutrate");
        cmd.setDefault("pseudo", "2.0");
        cmd.setDefault("readcountfactor", "1.0");
        cmd.setDefault("mutrate", "0.01");
        cmd.setDefault("fraglengthmean", "200");
        cmd.setDefault("fraglengthsd", "60");




        if(!OptionParser.parseParams(args, false, false, true, true, cmd))
            return;




        boolean doplot = cmd.isSet("doplot");

        File trcounts = (cmd.isOptionSet("trcounts")) ? cmd.getFile("trcounts") : TRSIMULATION;
        IsoformRegionGetter annot = (cmd.isOptionSet("gtf")) ? new GFFBasedIsoformRegionGetter(cmd.getFile("gtf"), null, null) : getHumanAnnot();

        HUMAN_IRG = annot;
        SimulationConfiguration simulationConfiguration = new SimulationConfiguration(annot);
        simulationConfiguration.reportOutDir = cmd.getOptionalFile("reportout");
        simulationConfiguration.isoformRegionGetter = annot;
        simulationConfiguration.fastQGenerator = new FastQGenerator(annot, cmd.getOptionalFile("genome"), cmd.getOptionalFile("genomeidx"), cmd.getDouble("mutrate"));
        simulationConfiguration.fastQGenerator.setFastQOutDir(cmd.getOptionalFile("fastqout"));
        simulationConfiguration.fastQGenerator.setBamfileReference(cmd.getOptionalFile("bamreference"));
        simulationConfiguration.readLength = cmd.getInt("readlength");
        simulationConfiguration.fragmentLengthMean = cmd.getInt("fraglengthmean");
        simulationConfiguration.fragmentLengthSD = cmd.getInt("fraglengthsd");



        File trueSplicFile = (cmd.isOptionSet("trues")) ? cmd.getFile("trues") : TRUE_SPLICING;
        File trueDiffFile = (cmd.isOptionSet("truediffs")) ? cmd.getFile("truediffs") : TRUE_DIFFEXP;
        simulationConfiguration.trueSplic = (trueSplicFile.exists()) ? FileUtils.readSet(trueSplicFile) : new HashSet<>();
        simulationConfiguration.trueDiff = (trueDiffFile.exists()) ? FileUtils.readSet(trueDiffFile) : new HashSet<>();
        simulationConfiguration.setPositionBiases(cmd.getOptionalFile("biaspos"));

        File trcheck = cmd.getOptionalFile("checktrseqs");

        int MAXTRNUM = cmd.getInt("maxtrnum");
        if(cmd.isSet("eqredtest")) {
            TranscriptSimulationTest.eqReductionTest(simulationConfiguration, trcounts);
            return;
        }
        if(cmd.isSet("trpairvisu")) {
            TranscriptSimulationTest test = new TranscriptSimulationTest();
            test.transcriptPairSimulationTest();
            return;
        }

        if(cmd.isSet("biasposvisu")) {
            TranscriptSimulationTest.biasPositionTest();
            return;
        }



        File eqclassCounts = cmd.getOptionalFile("eqclasscounts");

        if(cmd.isSet("details")) {
            System.out.printf("eqclasscounts: %s\n", eqclassCounts);
            detailTestInfos(simulationConfiguration, trcounts, MAXTRNUM, simulationConfiguration.trueSplic, eqclassCounts);
            return;
        }
        if(trcheck != null) {
            /*
            HashMap<String, String> trseqs = SimpleFastaReader.readAll(trcheck, (_s) -> {
                String id = _s.trim().split("\\s+")[0];
                //System.out.printf("key: %s\nid: >%s<\n", _s, id);
                return id;
            });
            SplicingSimulation simulation = new SplicingSimulation(annot, trcounts);
            int nchecked = 0;
            for(String g : simulation.getSimulatedGenes()) {
                for(String tr : simulation.getSimulatedTranscripts(g)) {
                    String trseq = simulationConfiguration.fastQGenerator.getTrSeq(g, tr);
                    String ref = trseqs.get(tr);
                    if(ref == null) {
                        System.out.printf("ref missing for %s:%s:%s\n", annot.getRegionById(g).chr, g, tr);
                        continue;
                    }
                    nchecked++;
                    if(trseq.equals(ref))
                        continue;

                    System.out.printf("trseq differs for %s:%s\nreal:%s\ngot: %s\n", g, tr, trseqs.get(tr), trseq );
                }
            }

             */
            return;
        }
        System.out.printf("got %d genes\n", toVector(annot.getRegions()).size());
        

        if(simulationConfiguration.fastQGenerator.fastqOutDir != null) {
            simulateReads(trcounts, simulationConfiguration);
            return;
        }
        if(cmd.isSet("multitrtest")) {
            multiTranscriptSimulationTest(doplot, simulationConfiguration, trcounts, MAXTRNUM);
            return;
        }

        double PSEUDO = cmd.getDouble("pseudo");
        double READCOUNTFACTOR = cmd.getDouble("readcountfactor");


        if(cmd.isSet("multibenchmark")) {
            multiBenchMark(simulationConfiguration, trcounts, trueSplicFile);
            return;
        }
        benchmark(doplot, simulationConfiguration, trcounts, MAXTRNUM, trueSplicFile, PSEUDO, READCOUNTFACTOR);



    }

}