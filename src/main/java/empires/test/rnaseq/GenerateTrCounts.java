package empires.test.rnaseq;

import empires.EmpiRe;
import empires.input.EBROWSER_INPUT;
import empires.input.RNASeqSplicingInfo;
import empires.input.ReplicateSetInfo;
import empires.rnaseq.GFFBasedIsoformRegionGetter;
import empires.rnaseq.IsoformRegionGetter;
import empires.rnaseq.MultiIsoformRegion;
import empires.rnaseq.SplicingTest;
import lmu.utils.*;
import lmu.utils.fdr.PerformanceResult;
import lmu.utils.plotting.*;
import lmu.utils.swing.PagedDataTable;
import lmu.utils.plotting.CachedPlotCreator;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.awt.*;
import java.io.File;
import java.io.PrintWriter;
import java.util.*;
import java.util.function.Function;

import static lmu.utils.IteratorUtils.rangev;
import static lmu.utils.ObjectGetter.*;

public class GenerateTrCounts {

    Logger log = LogConfig.getLogger();
    Vector<String> features = new Vector<>();
    empires.input.ReplicateSetInfo baseCondition;
    empires.NormalizedReplicateSet normalizedBaseCondition;
    HashMap<String, Vector<Integer>> tr2counts1 = new HashMap<>();
    HashMap<String, Vector<Integer>> tr2counts2 = new HashMap<>();
    Vector<Pair<Double, String>> logmean2featureName ;
    HashMap<String, Integer> baseRanks;

    Vector<Double> logmeans;
    Vector<String> logmeanFeatures;
    boolean[] masked;
    int majorIndex;
    int minorIndex;


    enum RANKTYPE {
        MAJOR,
        MINOR,
        DIFFEXP
    }

    HashMap<RANKTYPE, Integer> rank2steps = new HashMap<>();
    HashMap<RANKTYPE, Integer> rank2current = new HashMap<>();

    static class ReadData {
        Vector<String> features = new Vector<>();
        HashMap<String, Vector<Vector<Integer>>> cond2replicates = new HashMap<>();
        Vector<String> conds;

        private ReadData() {

        }

        public ReadData(File trcounts) {
            Iterator<String[]> it = FileUtils.getFieldSetIterator(trcounts, "\t");
            String[] h = it.next();
            HashMap<String, Vector<Integer>> cond2index = new HashMap<>();

            for(int i=1; i<h.length; i++) {
                String[] sp = h[i].split("\\.");
                MapBuilder.updateV(cond2index, sp[0], i);
            }

            for(Map.Entry<String, Vector<Integer>> e : cond2index.entrySet()) {
                cond2replicates.put(e.getKey(), map(e.getValue(), (_i) -> new Vector<>()));
            }
            FileUtils.Converter c = new FileUtils.Converter();

            while(it.hasNext()) {
                c.set(it.next());
                String tr = c.toString(0);
                features.add(tr);

                for(Map.Entry<String, Vector<Integer>>  e: cond2index.entrySet()) {
                    Vector<Integer> counts = map(e.getValue(), (_i) -> c.toInt(_i));
                    Vector<Vector<Integer>> target = cond2replicates.get(e.getKey());
                    for(int i=0; i<counts.size(); i++) {
                        target.get(i).add(counts.get(i));
                    }
                }

            }
            conds = toVector(cond2replicates.keySet());
            for(String cond : conds) {
                Vector<Pair<Integer, Vector<Integer>>> libsize2replicate = NumUtils.sort(map(cond2replicates.get(cond), (_v) -> Pair.create(NumUtils.sum(_v), _v)), (_p) -> _p.getFirst(), true);
                cond2replicates.put(cond, map(libsize2replicate, (_p) -> _p.getSecond()));
            }
        }

        public ReadData getDeepestCondPair(int num_replicates) {
            Vector<Pair<Integer, String>> minlibsize2cond = new Vector<>();
            for(String cond : conds) {
                Vector<Vector<Integer>> replicates = cond2replicates.get(cond);
                if(replicates.size() < num_replicates)
                    continue;

                int sum = NumUtils.sum(replicates.get(num_replicates - 1));
                minlibsize2cond.add(Pair.create(sum, cond));
            }
            if(minlibsize2cond.size() < 2)
                throw new FRuntimeException("cannot retrieve two conditions with at least %d replicates replicate nums: %s", num_replicates,
                        map(conds, (_c) -> String.format("%s: %d", _c, cond2replicates.get(_c).size())));

            NumUtils.sort(minlibsize2cond, (_p) -> _p.getFirst(), true);
            System.out.printf("condition infos sorted by the %d-th replicate library size: %s\n", num_replicates, minlibsize2cond);
            ReadData rd = new ReadData();
            rd.features = features;
            rd.conds = new Vector<>();
            for(int i=0; i<2; i++) {
                String cond = minlibsize2cond.get(i).getSecond();
                rd.conds.add(cond);
                Vector<Vector<Integer>> replicates = cond2replicates.get(cond);
                rd.cond2replicates.put(cond, VectorUtils.slice(replicates, 0, num_replicates));
            }
            return rd;
        }

        public ReadData filterLowCounts(double mincountPerReplicate) {
            ReadData rd = new ReadData();
            rd.conds = conds;

            for(Map.Entry<String, Vector<Vector<Integer>>> e : cond2replicates.entrySet()) {
                rd.cond2replicates.put(e.getKey(), map(e.getValue(), (_v) -> new Vector<>()));
            }

            int THRESHOLD = (int)(mincountPerReplicate * cond2replicates.get(conds.get(0)).size());
            for(int featureIdx = 0; featureIdx < features.size(); featureIdx++) {
                String f = features.get(featureIdx);
                int IDX = featureIdx;
                int sum = NumUtils.sum(cond2replicates.get(conds.get(0)), (_v) -> _v.get(IDX));
                if(sum < THRESHOLD)
                    continue;

                rd.features.add(f);
                for(String cond : conds) {
                    Vector<Integer>  values = map(cond2replicates.get(cond), (_v) -> _v.get(IDX));
                    Vector<Vector<Integer>> target = rd.cond2replicates.get(cond);
                    for(int i=0; i<values.size(); i++) {
                        target.get(i).add(values.get(i));
                    }
                }

            }
            return rd;
        }
    }


    public GenerateTrCounts(File trcounts, int numReps) {

        ReadData rd = new ReadData(trcounts);
        rd = rd.getDeepestCondPair(numReps).filterLowCounts(2.0);
        this.features = rd.features;
        String cond1 = rd.conds.get(0);
        String cond2 = rd.conds.get(1);
        Vector<Vector<Integer>> dataC1 = rd.cond2replicates.get(cond1);
        Vector<Vector<Integer>> dataC2 = rd.cond2replicates.get(cond2);

        Vector<Vector<Double>> data = map(dataC1, (_v) -> new Vector<>());
        for(int featureIdx = 0; featureIdx < rd.features.size(); featureIdx++) {
            int IDX = featureIdx;
            String tr = rd.features.get(featureIdx);
            Vector<Integer> c1 = map(dataC1, (_v) -> _v.get(IDX));
            for(int i=0; i<c1.size(); i++) {
                double logVal = (c1.get(i) == 0) ? Double.NaN : NumUtils.logN(c1.get(i), 2.0);
                data.get(i).add(logVal);
            }
            Vector<Integer> c2 = map(dataC2, (_v) -> _v.get(IDX));
            tr2counts1.put(tr, c1);
            tr2counts2.put(tr, c2);
        }

        baseCondition = new empires.input.ReplicateSetInfo("c1",  map(rangev(dataC1.size()), (_i) -> "C1R" + (1 + _i)), features);
        for(int i=0; i<dataC1.size(); i++) {
            baseCondition.setLog2Data(i, data.get(i));
        }

        normalizedBaseCondition = new empires.NormalizedReplicateSet(baseCondition);

        logmean2featureName = NumUtils.sort(map_and_filter(features, (_f) -> Pair.create(normalizedBaseCondition.getNormed(_f).mean, _f), (_p) -> _p.getFirst() > 0.0), (_p) -> _p.getFirst(), false);

        logmeans = map(logmean2featureName, (_p) -> _p.getFirst());
        logmeanFeatures = map(logmean2featureName, (_p) -> _p.getSecond());
        Collections.reverse(logmean2featureName);

        majorIndex = logmeanFeatures.size() - (int)(0.1 * logmeanFeatures.size());
        minorIndex = logmeanFeatures.size() - (int)(0.75 * logmeanFeatures.size());
        masked = new boolean[logmeans.size()];
        baseRanks = buildIndexMap(map(logmean2featureName, (_p) -> _p.getSecond()));

        log.info("%s log2 high expr: %.2f low: %.2f", trcounts.getName(),
                logmeans.get(majorIndex),
                logmeans.get(minorIndex)
        );


    }

    void setTargetNumSimulations(int ntotal) {
        int maxrank = (int)(logmeans.size() * 0.8);
        int step = (int)(maxrank / (0.0 + ntotal));
        rank2steps.put(RANKTYPE.MAJOR, Math.max(2, step >> 2));
        rank2current.put(RANKTYPE.MAJOR, (int) (0.95 * logmeans.size()));


        rank2steps.put(RANKTYPE.MINOR, Math.max(2, step >> 1));
        rank2current.put(RANKTYPE.MINOR, (int) (0.6 * logmeans.size()));


        rank2steps.put(RANKTYPE.DIFFEXP,  step);
        rank2current.put(RANKTYPE.DIFFEXP, (int) (0.9 * logmeans.size()));


        log.info("rankinfos %s steps: %s", rank2current, rank2steps);
    }

    public int getNextIdx(RANKTYPE ranktype) {
        int target = rank2current.get(ranktype);
        int idx = getFreeIdx(target);

        if(idx >= 0) {
            rank2current.put(ranktype, idx - rank2steps.get(ranktype));
            masked[idx] = true;
        }
        return idx;
    }
    public int getFreeIdx(int targetIdx) {
        int step = 0;
        while(targetIdx - step >= 0 ||  targetIdx + step < features.size()) {
            for(int dir = 0; dir<2; dir++) {
                int idx = targetIdx + step * ((dir == 0) ? -1 : 1);
                if(idx < 0 || idx >= features.size())
                    continue;

                if(!masked[idx]) {
                    masked[idx] = true;
                    return idx;
                }

                step++;
            }
        }
        return -1;
    }


    public int getMajorIndex() {
        while(majorIndex > 0 && masked[majorIndex])
            majorIndex--;

        if(majorIndex >= 0) {
            masked[majorIndex] = true;
        }
        return majorIndex;
    }

    public String getFeatureByIndex(int idx) {
        return logmeanFeatures.get(idx);
    }

    public double getFC(int idx1, int idx2) {
        return logmeans.get(idx2) - logmeans.get(idx1);
    }
    public int getMinorIndex() {
        while(minorIndex > 0 && masked[minorIndex])
            minorIndex--;

        if(minorIndex >= 0) {
            masked[minorIndex] = true;
        }
        return minorIndex;

    }

    public int getIndexForFoldChange(int baseIndex, double foldchange) {
        double base_meanval = logmeans.get(baseIndex);
        int step  = (foldchange > 0) ? -1 : 1;
        int hitIdx = Collections.binarySearch(logmeans, base_meanval - foldchange);
        if(hitIdx < 0) {
            hitIdx = - hitIdx - 1;
        }
        while(hitIdx >=0 && hitIdx < masked.length && masked[hitIdx]) {
            hitIdx += step;
        }
        if(hitIdx < 0 || hitIdx >= masked.length)
            return -1;

        masked[hitIdx] = true;
        return hitIdx;
    }
    public static boolean isStandardChr(String chr) {
        if(chr.equals("X") || chr.equals("Y"))
            return true;

        try {
            Integer.parseInt(chr);
            return true;
        } catch(NumberFormatException ex) {
            //silent error
        }
        return false;

    }

    static class BaseInfo {
        MultiIsoformRegion gene;
        RegionVector rv1;
        RegionVector rv2;

        RegionVector rv1_uniq;
        RegionVector rv2_uniq;
        RegionVector rv_is;
        Vector<Integer> num_bases;

        public BaseInfo(IsoformRegionGetter gtf, HashMap<String, Vector<String>> refgene2tr, String geneId) {
            Vector<String> trs = refgene2tr.get(geneId);
            gene = gtf.getRegionById(geneId);
            rv1 = gene.isoforms.get(trs.get(0));
            rv2 = gene.isoforms.get(trs.get(1));

            rv1_uniq = rv1.substract(rv2);
            rv2_uniq = rv2.substract(rv1);
            rv_is = rv1.intersect(rv2);

            num_bases = toVector(rv1_uniq.getCoveredLength(), rv2_uniq.getCoveredLength(), rv_is.getCoveredLength());
            Collections.sort(num_bases);
            Collections.reverse(num_bases);
        }

        public static PagedDataTable.MJFrame show(Vector<BaseInfo> data) {
            DataTable dt = DataTable.buildTable(data,
                    DataTable.buildHeader("gene", (BaseInfo bi) -> bi.gene.id)
                    .add("tr1length", (bi) -> bi.rv1.getCoveredLength())
                    .add("tr2length", (bi) -> bi.rv2.getCoveredLength())
                    .add("longest", (bi) -> bi.num_bases.get(0))
                    .add("second", (bi) -> bi.num_bases.get(1))
                    .add("third", (bi) -> bi.num_bases.get(2))


            );
            PagedDataTable.MJFrame mjFrame = PagedDataTable.getDetailViewFrame(dt, null, null)
                    .addMenu("visu",
                            (ObjectGetter og) -> {

                                BaseInfo bi = (BaseInfo)og.getInData();
                                int maxX = Math.max(bi.rv1.getX2(), bi.rv2.getX2());
                                int minX = Math.min(bi.rv1.getX1(), bi.rv2.getX1());

                                int maxlength = (maxX - minX + 1);
                                int targetLength = 800;
                                double factor = (maxlength < targetLength) ? 1.0 : targetLength / (0.0 + maxlength);
                                System.out.printf("maxlength: %d min: %d max: %d\n", maxlength, minX, maxX);
                                ZoomedCoordinateConverter converter = new ZoomedCoordinateConverter(maxlength, 0.1, null);
                                Object exon = new Object();
                                converter.addCategory(exon, 0.4, null);
                                for(RegionVector rv : toVector(bi.rv1, bi.rv2)) {

                                    for(Region1D r : rv.getRegions()) {
                                        converter.addRegion(new Region1D(r.getX1() - minX, r.getX2() - minX), exon);
                                    }

                                }
                                ZoomedMultiPicPlaneFrame zoomedMultiPicPlaneFrame = new ZoomedMultiPicPlaneFrame(maxlength);

                                zoomedMultiPicPlaneFrame.setZoomedConverter(converter);
                                RegionPicturePane rpp = zoomedMultiPicPlaneFrame.createRegionPane("annot");



                                RegionPlot rp = rpp.addRegion("tr1", bi.rv1, bi.rv1.translate(-minX));
                                rp.setDefaultDecorator(new RegionPlot.RegionDecorator(Color.BLUE, true));
                                rpp.addRegion("tr2", bi.rv2,bi.rv2.translate(-minX)).setDefaultDecorator(new RegionPlot.RegionDecorator(Color.RED, true));
                                rpp.addRegion("A", bi.rv1_uniq, bi.rv1_uniq.translate(-minX)).setDefaultDecorator(new RegionPlot.RegionDecorator(Color.BLUE, true));
                                rpp.addRegion("B", bi.rv2_uniq, bi.rv2_uniq.translate(-minX)).setDefaultDecorator(new RegionPlot.RegionDecorator(Color.RED, true));
                                rpp.addRegion("C", bi.rv_is, bi.rv_is.translate(-minX)).setDefaultDecorator(new RegionPlot.RegionDecorator(Color.MAGENTA, true));

                                rpp.changeRegionSettings(30, 120, 10);
                                zoomedMultiPicPlaneFrame.setTitle("test title");
                                zoomedMultiPicPlaneFrame.update(true);
                                zoomedMultiPicPlaneFrame.setAutoSize();
                                zoomedMultiPicPlaneFrame.setVisible(true);


                            });

            mjFrame.setVisible(true);
            return mjFrame;
        }
    }
    public static void main(String[] args) {


        SimpleOptionParser cmd = new SimpleOptionParser("transcriptsToSimulate", "new", "gtf", "od", "diffexp", "diffsplic", "reps", "plot", "showtable");
        cmd.setFile("transcriptsToSimulate", "new", "gtf", "diffexp", "diffsplic");
        cmd.setInt("reps");
        cmd.setDefault("reps", "3");
        cmd.setDir("od");
        cmd.setSwitches("plot", "showtable");
        if(!OptionParser.parseParams(args, true, false, true, true, cmd))
            return;


        Logger log = LogConfig.getLogger();

        File od = cmd.getFile("od");
        int reps = cmd.getInt("reps");


        GenerateTrCounts newcounts = new GenerateTrCounts(cmd.getFile("new"), reps) ;

        HashMap<String, Vector<String>> refgene2tr = new HashMap<>();
        HashMap<String, String> tr2gene = new HashMap<>();

        empires.rnaseq.GFFBasedIsoformRegionGetter gtf = new GFFBasedIsoformRegionGetter(cmd.getFile("gtf"), null, null);



        HashMap<String, Double> major2fc = new HashMap<>();
        HashMap<String, Double> minor2fc = new HashMap<>();
        HashMap<String, Double> difffcs = new HashMap<>();

        HashMap<String, String> gene2longesttr = new HashMap<>();
        Vector<Pair<String, Integer>> stdchrgene2trlength = new Vector<>();
        apply(gtf.getRegions(),
                (_gene) -> {

                    apply(_gene.isoforms.keySet(), (_tr) -> tr2gene.put(_tr, _gene.id));
                    String longesttr = Pair.convert_reverse_sorted(buildMap(_gene.isoforms.entrySet(), (_e) -> _e.getKey(), (_e) -> _e.getValue().getCoveredLength()), false).get(0).getSecond();
                    gene2longesttr.put(_gene.id, longesttr);
                    if(isStandardChr(_gene.chr)) {
                        stdchrgene2trlength.add(Pair.create(_gene.id, _gene.isoforms.get(longesttr).getCoveredLength()));
                    }
                });


        Iterator<String[]> it = FileUtils.getFieldSetIterator(cmd.getFile("transcriptsToSimulate"), "\t");
        Vector<String> reftrs  = map_and_filter(it, (_s) -> _s[0], (_tr) -> tr2gene.containsKey(_tr));

        HashSet<String> diffexp = FileUtils.readSet(cmd.getFile("diffexp"));
        HashSet<String> diffsplic = FileUtils.readSet(cmd.getFile("diffsplic"));

        Set<String> alltosimulate = SetInfo.union(diffexp, diffsplic);

        if(reftrs.size() < alltosimulate.size() * 2) {
            System.err.printf("got only %d known transcripts to simulate, less than twice as much as differential (expression + splicing) to simulate (%d)! aborting...\n",
                    reftrs.size(), alltosimulate.size());
            return;
        }

        for(String tr : reftrs) {
            MapBuilder.updateV(refgene2tr, tr2gene.get(tr), tr);
        }



        SetInfo<String> si = new SetInfo<>(diffexp, diffsplic, SetInfo.ALL_FLAGS);

        double MINIMUM_FC = 0.4;
        NormalDistribution foldChangeGetter = new NormalDistribution(0.8, 0.4);


        //every simulated entity is a swap -> all simulated values in pairs


        //for splicing we assume a major and a minor isoform -> different regimes to select source gene from
        //  then we generate a fold change the major isoform and generate a different fold change to the minor
        //  -> at the same time we generate a second alternative spliced gene with opposing changes


        // a part of the generated splicing will have one of the tanscripts unchanged, other will have both
        // changing transcript level ( which might lead to gene level differential expression - if there were such a thing

        Vector<String> diffsplicids = shuffle(toVector(diffsplic), false);

        PrintWriter info_pw = FileUtils.getWriter(od, "simul.info");

        double oneTrZeroProb = 0.3;

        Vector<String> diffvec = NumUtils.sort(filter(diffexp, (_d) -> !diffsplic.contains(_d)), (_d) -> refgene2tr.get(_d).size(), true);
        log.info("got %d/%d multi tr diffexp\n", filteredSize(diffvec, (_d) -> refgene2tr.get(_d).size() > 1), diffvec.size());



        int ntotal = diffsplic.size() + NumUtils.sum(diffvec, (_d) -> refgene2tr.get(_d).size()) >> 1;

        newcounts.setTargetNumSimulations(ntotal);


        HashSet<String> Ndiffsplic = new HashSet<>();

        HashMap<String, String> tr2feature = new HashMap<>();
        HashMap<String, String> swaps = new HashMap<>();

        HashMap<String, HashMap<String, Double>> g2tr2fc = new HashMap<>();

        int numUnchangingMajor = 0;
        for(int i=0; i<diffsplicids.size(); i+=2) {
            String g1 = diffsplicids.get(i);
            String g2 = diffsplicids.get(i+1);

            Vector<String> g1_trs = refgene2tr.get(g1);
            Vector<String> g2_trs = refgene2tr.get(g2);


            Ndiffsplic.add(g1);
            Ndiffsplic.add(g2);

            boolean first_changing = Math.random() > oneTrZeroProb;

            double baseFC = 0.0;
            for (int majorMinorIndex = 0; majorMinorIndex < 2; majorMinorIndex++) {

                String g1_tr = g1_trs.get(majorMinorIndex);
                String g2_tr = g2_trs.get(majorMinorIndex);

                int idx1 = newcounts.getNextIdx((majorMinorIndex == 0) ? RANKTYPE.MAJOR : RANKTYPE.MINOR);


                if(idx1 < 0) {
                    throw new FRuntimeException("could not get next idx: %d/%d diffsplic", i, majorMinorIndex);
                }
                String f1 = newcounts.getFeatureByIndex(idx1);
                tr2feature.put(g1_tr, f1);




                String isoformType = (majorMinorIndex == 0) ? "major" : "minor";
                boolean changing = (majorMinorIndex != 0 || first_changing);
                if (!changing) {
                    //int idx2 = (majorMinorIndex == 0) ? newcounts.getMajorIndex() : newcounts.getMajorIndex();
                    numUnchangingMajor += 2;
                    int idx2 = newcounts.getNextIdx((majorMinorIndex == 0) ? RANKTYPE.MAJOR : RANKTYPE.MINOR);
                    String f2 = newcounts.getFeatureByIndex(idx2);
                    MapBuilder.updateM(g2tr2fc, g1, g1_tr, 0.0);
                    MapBuilder.updateM(g2tr2fc, g2, g2_tr, 0.0);
                    tr2feature.put(g2_tr, f2);
                    info_pw.printf("%s\t%s\t%s\tunchanging\t%s\tDIFFSPLIC\n", g1, g1_tr, isoformType, f1);
                    info_pw.printf("%s\t%s\t%s\tunchanging\t%s\tDIFFSPLIC\n", g2, g2_tr, isoformType, f2);
                    continue;
                }
                int idx2 = -1;
                double FC = Double.NaN;
                for(int iter = 0; idx2 <0 && iter < 5; iter++) {

                    FC = baseFC + ((Math.random() < 0.5) ? 1 : -1) * (MINIMUM_FC + 0.05 + Math.abs(foldChangeGetter.sample()));

                    idx2 = newcounts.getIndexForFoldChange(idx1, FC);

                }
                if(idx2 < 0)
                    throw new FRuntimeException("could not get fc to idx: %d/%d idx: %d diffsplic FC: %.2f mean was: %.2f", i,   majorMinorIndex, idx1, FC,
                            newcounts.logmeans.get(idx1));

                String f2 = newcounts.getFeatureByIndex(idx2);
                tr2feature.put(g2_tr, f2);
                swaps.put(f1, f2);
                swaps.put(f2, f1);
                double gotFC = newcounts.getFC(idx1, idx2);

                MapBuilder.updateM(g2tr2fc, g1, g1_tr, FC);
                MapBuilder.updateM(g2tr2fc, g2, g2_tr, -FC);

                if(majorMinorIndex == 0) {
                    major2fc.put(g1_tr, FC);
                    major2fc.put(g2_tr, FC);
                } else {
                    minor2fc.put(g1_tr, FC);
                    minor2fc.put(g2_tr, FC);
                    double fcdiff = FC - baseFC;
                    difffcs.put(g1_tr, fcdiff);
                }
                System.out.printf("wanted change: %.2f got: %.2f\n", FC, gotFC);

                info_pw.printf("%s\t%s\t%s\t%.2f\t%s\tDIFFSPLIC\n", g1, g1_tr, isoformType, gotFC, f1);
                info_pw.printf("%s\t%s\t%s\t%.2f\t%s\tDIFFSPLIC\n", g2, g2_tr, isoformType, -gotFC, f2);


                baseFC = FC;
            }
        }


        log.info("after: diffsplic got %d mapping\n", tr2feature.size());

        if(cmd.isSet("plot")) {
            PlotCreator pc = CachedPlotCreator.getPlotCreator();
            Vector<Double> majorfcs = toVector(major2fc.values());
            Vector<Double> minorfcs = toVector(minor2fc.values());
            Vector<Double> diff_fcs = toVector(difffcs.values());
            pc.writeInMultipleFormats(new File(od, "fcsimul"));
            pc.cumhist("major-tr", majorfcs, majorfcs.size(), false, true);
            pc.cumhist("minor-tr", minorfcs, minorfcs.size(), false, true);
            pc.cumhist("diff", diff_fcs, diff_fcs.size(), false, true);
            pc.setLabels("log2FC", "% objects with >= x change", "bottomright");
            pc.setTitle("%d DAS genes, %d unchanging major transcripts", Ndiffsplic.size(), numUnchangingMajor);
            pc.writeImage(new File(od, "fcsimul.png"));


            Vector<Vector<Integer>> order2lengths = mapIndex(3, (_i) -> new Vector<>());
            Vector<BaseInfo> baseInfos = new Vector<>();

            for(String geneId : Ndiffsplic) {
                BaseInfo bi = new BaseInfo(gtf, refgene2tr, geneId);
                baseInfos.add(bi);
                for(int i=0; i<3; i++) {
                    order2lengths.get(i).add(bi.num_bases.get(i));
                }
            }
            pc.cumhist("longest", order2lengths.get(0), order2lengths.get(0).size(), false, true);
            pc.cumhist("second", order2lengths.get(1), order2lengths.get(1).size(), false, true);
            pc.cumhist("third", order2lengths.get(2), order2lengths.get(2).size(), false, true);
            pc.setLog(true, false);
            pc.setLabels("num bases", "% transcript pairs resulting in >= bp detectable regions", "bottomright");
            pc.setTitle("%d DAS genes", order2lengths.get(0).size());
            pc.writeInMultipleFormats(new File(od, "detectable_bases"));
            pc.writeImage(new File(od, "detectable_bases.png"));
            pc.destroy();

            if(cmd.isSet("showtable")) {
                BaseInfo.show(baseInfos);
            }
        }
        HashSet<String> Ndiffexp = new HashSet<>();
        Vector<Double> diffs = new Vector<>();

        //log.info("will simulate  %d  x 2 diff trs with step: %d", ntotal, step);
        for(int i=0; i<diffvec.size(); i+=2) {
            double FC = ((Math.random() < 0.5) ? 1 : -1) * (MINIMUM_FC + Math.abs(foldChangeGetter.sample()));
            String g1 = diffvec.get(i);
            String g2 = diffvec.get(i+1);

            Ndiffexp.add(g1);
            Ndiffexp.add(g2);

            Vector<String> g1_trs = refgene2tr.get(g1);
            Vector<String> g2_trs = refgene2tr.get(g2);

            int minNumTr = Math.min(g1_trs.size(), g2_trs.size());

            Vector<Double> real_fcs = new Vector<>();
            for(int trIdx = 0; trIdx < minNumTr; trIdx++) {
                String g1_tr = g1_trs.get(trIdx);
                String g2_tr = g2_trs.get(trIdx);
                int idx = newcounts.getNextIdx((trIdx == 0) ? RANKTYPE.MAJOR : RANKTYPE.DIFFEXP);
                int idx2 = newcounts.getIndexForFoldChange(idx, FC);
                if(idx2 < 0 ) {
                    log.warn("could not generate fc: %.2f on %d (i: %d %s:%s/%d)", FC, idx, i, g1, g1_tr, trIdx);
                }
                String fn1 = newcounts.getFeatureByIndex(idx);
                String fn2 = newcounts.getFeatureByIndex(idx2);
                tr2feature.put(g1_tr, fn1);
                tr2feature.put(g2_tr, fn2);
                swaps.put(fn1, fn2);
                swaps.put(fn2, fn1);
                double gotFC = newcounts.getFC(idx, idx2);
                if(Math.abs(gotFC) < MINIMUM_FC)
                    throw  new FRuntimeException("invlaid fc: %.2f simul: %.2f", gotFC, FC);
                real_fcs.add(gotFC);

                info_pw.printf("%s\t%s\t%s\t%.2f\t%s\tDIFFEXP\n", g1, g1_tr, "rank"+idx, gotFC, fn1);
                info_pw.printf("%s\t%s\t%s\t%.2f\t%s\tDIFFEXP\n", g2, g2_tr, "rank"+idx2, -gotFC, fn2);

            }

            if(real_fcs.size() == 1)
                continue;

            Collections.sort(real_fcs);
            double maxfcdiff = real_fcs.get(real_fcs.size() - 1) - real_fcs.get(0);

            diffs.add(maxfcdiff);

        }

        log.info("simulated diffexp maxfcdiff distrib: %s", NumUtils.getNumInfo(diffs).getInfoWithQ());
        info_pw.close();
        log.info("after: diffsplic + diffexp got %d mapping\n", tr2feature.size());

        Function<Vector<Integer>, String> fmt = (_v) -> StringUtils.joinObjects("\t", _v);


        System.out.printf("%d diffexp, %d splicing, %d both\n", diffexp.size(), diffsplic.size(), si.intersection_set.size());

        int freeIdx = 0;

        Vector<String> flist = map_and_filter(stdchrgene2trlength, (_p) -> gene2longesttr.get(_p.getFirst()), (_g) -> !refgene2tr.containsKey(_g));

        HashMap<String, Vector<String>> g2tr = new HashMap<>();

        HashSet<String> usedGenes = new HashSet<>();
        Vector<String> outFeatures = new Vector<>();
        for(int qi=0; qi<2; qi++) {
            Vector<String> srclist = (qi == 0) ? reftrs : flist;

            for(String tr : srclist) {
                if(tr2feature.containsKey(tr)) {
                    MapBuilder.updateV(g2tr, tr2gene.get(tr), tr);
                    outFeatures.add(tr);
                    continue;
                }
                String g = tr2gene.get(tr);
                if(usedGenes.contains(g))
                    continue;


                while(freeIdx < newcounts.masked.length && newcounts.masked[freeIdx]) {
                    freeIdx++;
                }
                if(freeIdx >= newcounts.masked.length)
                    break;

                tr2feature.put(tr, newcounts.logmeanFeatures.get(freeIdx++));
                outFeatures.add(tr);
                usedGenes.add(tr);
                MapBuilder.updateV(g2tr, g, tr);
            }

        }

        log.info("used: %d/%d of the features\n", tr2feature.size(), newcounts.logmeanFeatures.size());

        PrintWriter expPw = FileUtils.getWriter(od, "transcript_exprs.txt");
        expPw.printf("%s\t%s\n", StringUtils.joinObjects("\t", mapIndex(reps, (_i) -> String.format("C1R%d", 1 + _i))),
                 StringUtils.joinObjects("\t", mapIndex(reps, (_i) -> String.format("C2R%d", 1 + _i)))
        );

        Vector<String> nfeatures = new Vector<>();
        Vector<Vector<Double>> vals1 = mapIndex(reps, (_i) -> new Vector<>());
        Vector<Vector<Double>> vals2 = mapIndex(reps, (_i) -> new Vector<>());

        HashSet<String> usedRefGenes = new HashSet<>();
        //HashMap<String, Vector<String>>

        HashMap<String, Vector<Integer>> genecounts1 = new HashMap<>();
        HashMap<String, Vector<Integer>> genecounts2 = new HashMap<>();

        HashMap<String, UPair<Vector<Integer>>> simuldata = new HashMap<>();

        for(String tr : outFeatures) {
            String m1 = tr2feature.get(tr);
            if(m1 == null)
                continue;
            String m2 = swaps.getOrDefault(m1, m1);

            String gene = tr2gene.get(tr);
            usedRefGenes.add(gene);
            nfeatures.add(m1);
            Vector<Integer> c1 = newcounts.tr2counts1.get(m1);
            Vector<Integer> c2 = newcounts.tr2counts2.get(m2);

            System.out.printf("tr: %s m1: %s m2: %s c1: %s c2: %s\n", tr, m1, m2, c1, c2);
            simuldata.put(tr, UPair.createU(c1, c2));

            for(int cidx = 0; cidx < 2; cidx++) {

                HashMap<String, Vector<Integer>> geneMap = (cidx == 0) ? genecounts1 : genecounts2;
                Vector<Integer>  src = (cidx == 0) ? c1 : c2;
                Vector<Integer> old = geneMap.get(gene);
                if(old == null) {
                    geneMap.put(gene, src);
                } else {
                    Vector<Integer> sum = map(src, (_i) -> 0);
                    for(int i=0; i<src.size(); i++) {
                        sum.set(i, src.get(i) + old.get(i));
                    }
                    System.out.printf("%s,%d %s\n", gene, cidx, sum);
                    geneMap.put(gene, sum);
                }

                Vector<Vector<Double>> target = (cidx == 0) ? vals1 : vals2;


                for(int i=0; i<src.size(); i++) {
                    double v = (src.get(i) == 0) ? Double.NaN : NumUtils.logN(src.get(i), 2.0);
                    target.get(i).add(v);
                }
            }

            //System.out.printf("%s -> %s, %s %s, %s\n", m1, m2, newcounts);
            expPw.printf("%s\t%s\t%s\n", tr, fmt.apply(c1), fmt.apply(c2));

        }
        expPw.close();

        Vector<String> genefeatures = toVector(genecounts1.keySet());
        Vector<Vector<Double>> gene_vals1 = mapIndex(reps, (_i) -> new Vector<>());
        Vector<Vector<Double>> gene_vals2 = mapIndex(reps, (_i) -> new Vector<>());
        Vector<String> rep1names = mapIndex(reps, (_i) -> "C1R" + (1 + _i));
        Vector<String> rep2names = mapIndex(reps, (_i) -> "C2R" + (1 + _i));

        for(String g : genefeatures) {
            for(int cidx = 0; cidx < 2; cidx++) {
                HashMap<String, Vector<Integer>> geneMap = (cidx == 0) ? genecounts1 : genecounts2;
                Vector<Vector<Double>> target = (cidx == 0) ? gene_vals1 : gene_vals2;
                Vector<Double> logvals = map(geneMap.get(g), (_d) -> (_d == 0) ? Double.NaN : NumUtils.logN(_d, 2.0));
                applyIndex(logvals.size(), (_i) -> target.get(_i).add(logvals.get(_i)));
            }
        }

        empires.input.ReplicateSetInfo gene_rs1 = new empires.input.ReplicateSetInfo("C1", rep1names, genefeatures);
        empires.input.ReplicateSetInfo gene_rs2 = new empires.input.ReplicateSetInfo("C2", rep2names, genefeatures);

        empires.input.ReplicateSetInfo rs1 = new empires.input.ReplicateSetInfo("C1", rep1names, nfeatures);
        empires.input.ReplicateSetInfo rs2 = new empires.input.ReplicateSetInfo("C2", rep2names, nfeatures);

        HashMap<String, Vector<String>> condition2replicatenames = new HashMap<>();
        condition2replicatenames.put(rs1.getReplicateSetName(), rs1.getReplicateNames());
        condition2replicatenames.put(rs2.getReplicateSetName(), rs2.getReplicateNames());

        for(int c = 0; c<2; c++){
            empires.input.ReplicateSetInfo target = (c == 0) ? rs1 : rs2;
            Vector<Vector<Double>> data = (c == 0) ? vals1 : vals2;

            ReplicateSetInfo gene_target = (c == 0) ? gene_rs1 : gene_rs2;
            Vector<Vector<Double>> gene_data = (c == 0) ? gene_vals1 : gene_vals2;
            for(int r = 0; r < data.size(); r++) {
                target.setLog2Data(r, data.get(r));
                gene_target.setLog2Data(r, gene_data.get(r));
            }
        }

        empires.NormalizedReplicateSet g_nrs1 = new empires.NormalizedReplicateSet(gene_rs1);
        empires.NormalizedReplicateSet g_nrs2 = new empires.NormalizedReplicateSet(gene_rs2);
        empires.EmpiRe empiRe = new empires.EmpiRe();
        Vector<empires.DiffExpResult> diffExpResults = filter(empiRe.getDifferentialResults(g_nrs1, g_nrs2), (_d) -> !Ndiffsplic.contains(_d.combinedFeatureName));
        //new DiffExpTable(g_nrs1, g_nrs2, diffExpResults).showInteractiveTable(Pair.create("true", (_d) -> Ndiffexp.contains(_d.combinedFeatureName)));
        log.info("diffexp check: %s\n%s",
                new PerformanceResult("emp.old", diffExpResults, (_d) -> _d.fcEstimatePval, (_d) -> Ndiffexp.contains(_d.combinedFeatureName), false,
                (_d) -> _d.fdr <= 0.05, null),
                new PerformanceResult("emp.fcfdr", diffExpResults, (_d) -> _d.fcEstimatePval, (_d) -> Ndiffexp.contains(_d.combinedFeatureName), false,
                        (_d) -> _d.fcEstimateFDR <= 0.05, null));


        Vector<Double> fuzzy5Proportions = toVector(40., 20., 10., 5., 3.);

        empires.NormalizedReplicateSet g_nrs1_f1 = new empires.NormalizedReplicateSet(gene_rs1, new empires.ProportionBackGroundContextProvider(fuzzy5Proportions));
        empires.NormalizedReplicateSet g_nrs2_f2 = new empires.NormalizedReplicateSet(gene_rs2, new empires.ProportionBackGroundContextProvider(fuzzy5Proportions));

        empires.EmpiRe f_empiRe = new EmpiRe();
        Vector<empires.DiffExpResult> f_diffExpResults = filter(empiRe.getDifferentialResults(g_nrs1_f1, g_nrs2_f2), (_d) -> !Ndiffsplic.contains(_d.combinedFeatureName));
        //new DiffExpTable(g_nrs1, g_nrs2, diffExpResults).showInteractiveTable(Pair.create("true", (_d) -> Ndiffexp.contains(_d.combinedFeatureName)));
        log.info("fuzzy diffexp check: %s\n%s",
                new PerformanceResult("fuzzy.emp.old", f_diffExpResults, (_d) -> _d.fcEstimatePval, (_d) -> Ndiffexp.contains(_d.combinedFeatureName), false,
                        (_d) -> _d.fdr <= 0.05, null),
                new PerformanceResult("fuzzy.emp.fcfdr", f_diffExpResults, (_d) -> _d.fcEstimatePval, (_d) -> Ndiffexp.contains(_d.combinedFeatureName), false,
                        (_d) -> _d.fcEstimateFDR <= 0.05, null));



        empires.input.RNASeqSplicingInfo rsi = new RNASeqSplicingInfo(1000, 0, condition2replicatenames);
        HashMap<String, HashMap<String, Vector<HashMap<Tuple, Double>>>> simulSplicData = new HashMap<>();
        for(String g : g2tr.keySet()) {
            Vector<String> trs = g2tr.get(g);
            HashMap<String, Vector<HashMap<Tuple, Double>>> condition2replicate2EQclassCounts = new HashMap<>();
            condition2replicate2EQclassCounts.put(rs1.getReplicateSetName(), map(rs1.getReplicateNames(), (_r) -> new HashMap<>()));
            condition2replicate2EQclassCounts.put(rs2.getReplicateSetName(), map(rs2.getReplicateNames(), (_r) -> new HashMap<>()));

            for(String tr : trs) {
                String m1 = tr2feature.get(tr);
                if(m1 == null)
                    continue;
                String m2 = swaps.getOrDefault(m1, m1);
                Vector<Integer> c1 = newcounts.tr2counts1.get(m1);
                Vector<Integer> c2 = newcounts.tr2counts2.get(m2);

                Tuple key = new Tuple(tr);
                for(int i=0; i<2; i++) {
                    String cname = ((i == 0) ? rs1 : rs2).getReplicateSetName();

                    Vector<Integer> src = ( i == 0) ? c1 : c2;
                    Vector<HashMap<Tuple, Double>> target = condition2replicate2EQclassCounts.get(cname);
                    for(int r =0 ; r < src.size(); r++) {
                        target.get(r).put(key, 0.0 +src.get(r));
                    }

                }

            }
            simulSplicData.put(g, condition2replicate2EQclassCounts);
            rsi.addGeneInfo(g, condition2replicate2EQclassCounts);
        }
        empires.rnaseq.SplicingTest stest = new SplicingTest(rsi, new empires.AutoBackGroundContextProvider());

        Vector<empires.DoubleDiffResult> diffsplicResults = stest.getDifferentialAlternativeSplicing(rs1.getReplicateSetName(), rs2
                .getReplicateSetName()).getFirst();
        log.info("diff splic check: %s\n",
                new PerformanceResult("emp.splic", diffsplicResults, (_d) -> _d.pval, (_d) -> Ndiffsplic.contains(_d.testName.split("\\.")[0]), false,
                        (_d) -> _d.fdr <= 0.05, null)
        );

        PrintWriter diffexpPw = FileUtils.getWriter(od, "diffexp.trues");
        apply(Ndiffexp, (_g) -> diffexpPw.printf("%s\n", _g));
        diffexpPw.close();

        PrintWriter diffsplicPw = FileUtils.getWriter(od, "diffsplic.trues");
        apply(Ndiffsplic, (_g) -> diffsplicPw.printf("%s\n", _g));
        diffsplicPw.close();


        Vector<String> eb_genes = filter(genefeatures, (_g) -> !Ndiffsplic.contains(_g));

        PrintWriter eb_fdata = FileUtils.getWriter(empires.input.EBROWSER_INPUT.FDATA.get(od));
        apply(eb_genes, (_g) -> eb_fdata.println(_g));
        eb_fdata.close();
        PrintWriter eb_pdata = FileUtils.getWriter(empires.input.EBROWSER_INPUT.PDATA.get(od));
        apply(gene_rs1.getReplicateNames(), (_rn) -> eb_pdata.printf("%s\t0\n", _rn));
        apply(gene_rs2.getReplicateNames(), (_rn) -> eb_pdata.printf("%s\t1\n", _rn));
        eb_pdata.close();

        PrintWriter eb_trues = FileUtils.getWriter(empires.input.EBROWSER_INPUT.TRUES.get(od));
        apply(Ndiffexp, (_g) -> eb_trues.println(_g));
        eb_trues.close();

        PrintWriter eb_expr = FileUtils.getWriter(EBROWSER_INPUT.EXPRS.get(od));
        for(String g : eb_genes) {
            for(int i = 0; i< 2; i++) {
                Vector<Double> vals = ((i == 0) ? gene_rs1 : gene_rs2 ).getReplicateData(g);
                if(i > 0) {
                    eb_expr.printf("\t");
                }
                eb_expr.print(StringUtils.joinObjects("\t", vals, (_d) -> String.format("%.0f", (Double.isNaN(_d)) ? 0 : Math.pow(2.0, _d))));
            }
            eb_expr.println();
        }

        eb_expr.close();


        /*
        DataTable dt = BenchmarkGene.getTable(diffsplicResults, Ndiffsplic, null, rsi,
                toVector(Pair.create("simulInfo", (_bg) ->
                    {
                        //System.out.printf("simul: %s\n", g2tr2fc.get(_bg.geneId));
                        return ""+ g2tr2fc.get(_bg.geneId);
                    })
                ));
        BenchmarkGene.getWithDetails(dt).addMenu(
                "simdet", (_og) -> {
                    BenchmarkGene _bg = (BenchmarkGene)_og.getInData();
                    System.out.printf("simul %s: %s\n", _bg.geneId, g2tr2fc.get(_bg.geneId));
                    for(String tr : refgene2tr.get(_bg.geneId)) {
                        System.out.printf("%s:%s\t%s\n", _bg.geneId, tr, simuldata.get(tr));
                    }
                    System.out.printf("EQclass:\n%s\n\n", simulSplicData.get(_bg.geneId));
                }
        ).setTitle("test");
        */
        log.info("got %d genes", g2tr.size());
        //rsi.addGeneInfo()
        //1) diff splicing  + diff exp : both transcripts with foldchanges and different ones

        //Vector<Integer> newvals1 = newcounts.tr2counts1.get(m1);
        //Vector<Integer> newvals2 = newcounts.tr2counts2.get(m1);


        //pw.printf("%s\t%s\t%s\n", tr, fmt.apply(newvals1), fmt.apply(newvals2));

    }
}
