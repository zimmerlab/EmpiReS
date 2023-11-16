package empires.test.rnaseq;

import empires.rnaseq.MultiIsoformRegion;
import empires.rnaseq.ReducedTranscriptPresentation;
import empires.rnaseq.simulation.SimulatedRead;
import empires.rnaseq.simulation.SplicingSimulation;
import lmu.utils.*;
import lmu.utils.plotting.*;
import lmu.utils.tuple.Tuple4;

import java.awt.*;
import java.awt.image.BufferedImage;
import java.util.*;
import java.util.function.Supplier;

import static lmu.utils.IteratorUtils.rangev;
import static lmu.utils.ObjectGetter.*;
import static lmu.utils.ObjectGetter.map;

public class SimulatedGeneWithCoverage {
    String id;
    boolean strand;
    RegionVector merged;
    ZoomedMultiPicPlaneFrame zoomedMultiPicPlaneFrame;

    RegionPicturePane annot;
    RegionPicturePane eqAnnot;
    int maxTranscriptToReduce;

    HashMap<String, Color> tr2color = new HashMap<>();
    static final Vector<Integer> reduceTargets = toVector(10, 8, 6, 5, 3, 2);
    static Color MIXED = Color.MAGENTA;

    static class SimulationStarter {
        empires.rnaseq.simulation.SplicingSimulation simulation;
        Supplier<Iterable<empires.rnaseq.simulation.SimulatedRead>> supplier;
        empires.rnaseq.simulation.SplicingSimulation.ReadInReplicateConsumer consumer;
        String id;

        public SimulationStarter(empires.rnaseq.simulation.SplicingSimulation simulation) {
            this.simulation = simulation;
        }


        public SimulationStarter(Supplier<Iterable<empires.rnaseq.simulation.SimulatedRead>> supplier) {
            this.supplier = supplier;
        }


        public void setConsumer(String id, empires.rnaseq.simulation.SplicingSimulation.ReadInReplicateConsumer consumer) {
            this.id = id;
            this.consumer = consumer;
        }

        public void simulate() {
            if(supplier != null) {
                for(empires.rnaseq.simulation.SimulatedRead sr : supplier.get()) {
                    consumer.accept("cond1", 0, sr);
                }
                return;
            }
            simulation.simulate(id, false, consumer);
        }

    }
    public SimulatedGeneWithCoverage(MultiIsoformRegion gene, SplicingSimulation simulation, int maxTranscriptToReduce) {
        this(gene.id, gene.strand, gene.isoforms, simulation.getSimulatedTranscripts(gene.id), maxTranscriptToReduce, new SimulationStarter(simulation));
    }

    public SimulatedGeneWithCoverage(String id, boolean strand, HashMap<String, RegionVector> isoforms, Set<String> relevant, int maxTranscriptToReduce,
                                     Supplier<Iterable<empires.rnaseq.simulation.SimulatedRead>> simulator) {
        this(id, strand, isoforms, relevant, maxTranscriptToReduce, new SimulationStarter(simulator));
    }

    private SimulatedGeneWithCoverage(String id, boolean strand, HashMap<String, RegionVector> isoforms, Set<String> relevant, int maxTranscriptToReduce,
                                     SimulationStarter simulationStarter) {

        this.id = id;
        this.strand = strand;
        merged = RegionVector.merge(isoforms.values());
        this.maxTranscriptToReduce = maxTranscriptToReduce;

        final int L = merged.getCoveredLength() + 1;

        int targetLength = 800;
        double factor = (L < targetLength) ? 1.0 : targetLength / (0.0 + L);
        ZoomedCoordinateConverter converter = new ZoomedCoordinateConverter(L, factor, null);

        zoomedMultiPicPlaneFrame = new ZoomedMultiPicPlaneFrame(L);

        zoomedMultiPicPlaneFrame.setZoomedConverter(converter);
        annot = zoomedMultiPicPlaneFrame.createRegionPane("annot");
        eqAnnot = zoomedMultiPicPlaneFrame.createRegionPane("EQannot");
        Vector<RegionPicturePane> reducedPanels = new Vector<>();
        for(int target : reduceTargets) {
            reducedPanels.add((target >= isoforms.size()) ? null : zoomedMultiPicPlaneFrame.createRegionPane("reduced EQannot maxTR=" + target));
        }



        Iterator<Color> it = filterIterator(ColorManager.getColorIterator(), (_c) -> ! _c.equals(Color.BLACK) && !_c.equals(MIXED));

        for(String tr : relevant) {
            tr2color.put(tr, it.next());
        }
        Vector<String> trs = filter(isoforms.keySet(), (_tr) -> !relevant.contains(_tr));
        trs.addAll(relevant);
        for (String tr : trs) {
            RegionVector r = isoforms.get(tr);
            RegionVector converted = empires.test.rnaseq.TranscriptSimulationTest.convertToInner(merged, !strand, r);


            if(relevant.contains(tr)) {

                Color c = tr2color.get(tr);
                annot.addRegion(tr, null, converted).setDefaultDecorator(new RegionPlot.RegionDecorator(c, true));
                eqAnnot.addRegion(tr, null, converted).setDefaultDecorator(new RegionPlot.RegionDecorator(c, true));
                for(RegionPicturePane rpp : reducedPanels) {
                    if(rpp == null)
                        continue;
                    rpp.addRegion(tr, null, converted).setDefaultDecorator(new RegionPlot.RegionDecorator(c, true));
                }


            } else {
                annot.addRegion(tr, null, converted);
            }
        }


        HashMap<Tuple, Vector<empires.rnaseq.simulation.SimulatedRead>> eq2regvecs = new HashMap<>();

        DataPicturePane dpp = zoomedMultiPicPlaneFrame.createDataPane("coverages");


        HashMap<String, HashSet<Integer>> cond2rep = new HashMap<>();

        HashMap<String, int[]> rep2coverage = new HashMap<>();
        simulationStarter.setConsumer(id,
                (String condition, int replicateId, empires.rnaseq.simulation.SimulatedRead read)
                        -> {

                    MapBuilder.update(cond2rep, condition, replicateId);
                    Tuple t = Tuple.tupleFromCollection(toSortedVector(read.mapsToTranscripts, true));
                    eq2regvecs.computeIfAbsent(t, (_t) -> new Vector<>()).add(read);
                    String repname = condition + "." + (replicateId+1);
                    int[] coverage = rep2coverage.computeIfAbsent(repname, (_n) -> new int[L]);
                    for (int fw = 0; fw < 2; fw++) {
                        RegionVector rv = (fw == 0) ? read.fw : read.rw;
                        for (Region1D r : empires.test.rnaseq.TranscriptSimulationTest.convertToInner(merged, !strand, rv).getRegions()) {
                            coverage[r.getX1()]++;
                            coverage[r.getX2()]--;
                        }
                    }


                });

        simulationStarter.simulate();


        addEQClassesToPane(eqAnnot, eq2regvecs, null);

        for(int i=0; i<reduceTargets.size(); i++) {

            int target = reduceTargets.get(i);
            if(target >= isoforms.size())
                continue;

            empires.rnaseq.ReducedTranscriptPresentation rtp = new empires.rnaseq.ReducedTranscriptPresentation(
                    buildMap(eq2regvecs.entrySet(), (_e) -> _e.getKey(), (_e) -> 0.0 + _e.getValue().size()),
                    reduceTargets.get(i)
            );

            HashMap<Tuple, Vector<empires.rnaseq.simulation.SimulatedRead>> reducedEQ2regvecs = new HashMap<>();
            for(Tuple t : eq2regvecs.keySet()) {
                Tuple reduced = rtp.reduceTuple(t);
                reducedEQ2regvecs.computeIfAbsent(reduced, (_k) -> new Vector<>()).addAll(eq2regvecs.get(t));
            }

            System.out.printf("TARGET: %d isoforms: %d\n", reduceTargets.get(i), isoforms.size());
            addEQClassesToPane(reducedPanels.get(i), reducedEQ2regvecs, rtp);

        }

        int maxVal = 0;
        Vector<String> conditions = toSortedVector(cond2rep.keySet(), true);

        for(int i = 0; i< conditions.size(); i++) {
            String cond = conditions.get(i);
            Vector<Integer> replicates = toSortedVector(cond2rep.get(cond), true);
            int nrep = replicates.size();
            for (int j = 0; j < nrep; j++) {
                String repname = cond + "." + (j + 1);
                int[] coverage = rep2coverage.get(repname);
                if (coverage == null)
                    continue;

                for (int ci = 1; ci < coverage.length; ci++) {
                    coverage[ci] += coverage[ci - 1];
                    maxVal = Math.max(coverage[ci], maxVal);
                }

                dpp.addData(repname, (i == 0) ? Color.BLUE : Color.RED).setDenseData(coverage);
            }
        }
        dpp.getHeader().setMaxValue(1.0 + maxVal);
        dpp.getHeader().setHeight(300);


    }
    public void show() {
        show(true);
    }

    void addEQClassesToPane(RegionPicturePane pane, HashMap<Tuple, Vector<empires.rnaseq.simulation.SimulatedRead>> eq2regvecs, ReducedTranscriptPresentation rtp) {

        Vector<Tuple4<Integer, String, RegionPlot.RegionDecorator,  RegionVector>> eqRegVecs = new Vector<>();

        for(Tuple t : eq2regvecs.keySet()) {
            HashSet<String> trSet = new HashSet<>();
            for(String tupleTR : map(rangev(t.cardinality()), (_i) -> t.getAsString(_i))) {
                if(rtp == null) {
                    trSet.add(tupleTR);
                    continue;
                }
                trSet.addAll(rtp.getTrCluster(tupleTR));
            }
            Vector<Color> cv  = toVector(toSet(filter_and_map(trSet, (_t) -> tr2color.containsKey(_t), (_t) -> tr2color.get(_t))));

            Vector<RegionVector> tomerge = new Vector<>();
            for(SimulatedRead sr : eq2regvecs.get(t)) {
                tomerge.add(sr.fw);
                tomerge.add(sr.rw);
            }
            RegionVector tmerged = RegionVector.merge(tomerge);

            int pcount = tomerge.size() >> 1;
            Color c = (cv.size() == 1) ? cv.get(0) : (cv.size() > 1) ? Color.MAGENTA : Color.BLACK;
            RegionPlot.RegionDecorator rd = new RegionPlot.RegionDecorator(c, true);
            eqRegVecs.add(Tuple4.create(pcount, pcount+":" + t, rd,  empires.test.rnaseq.TranscriptSimulationTest.convertToInner(merged, !strand, tmerged)));

        }

        NumUtils.sort(eqRegVecs, (_t) -> _t.get0(), true);
        for(Tuple4<Integer, String, RegionPlot.RegionDecorator, RegionVector> t : eqRegVecs) {
            pane.addRegion(t.get1(), null, t.get3()).setDefaultDecorator(t.get2());
        }


    }

    public BufferedImage getImage() {
        zoomedMultiPicPlaneFrame.update(true);
        zoomedMultiPicPlaneFrame.setAutoSize();
        return zoomedMultiPicPlaneFrame.getCombinedImage();

    }
    public void show(boolean showEqClasses) {

        zoomedMultiPicPlaneFrame.update(true);
        zoomedMultiPicPlaneFrame.setAutoSize();
        zoomedMultiPicPlaneFrame.setVisible(true);

        zoomedMultiPicPlaneFrame.setPreferredSize(new Dimension(900, 800));
        //ImageUtils.showImage("test", zoomedMultiPicPlaneFrame.getCombinedImage(), false);
    }

}
