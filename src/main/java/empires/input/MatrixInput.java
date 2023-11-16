package empires.input;

import empires.EmpiRe;
import empires.plotting.DataOverview;
import empires.plotting.DiffExpTable;
import empires.plotting.NormalizedReplicateSetPlotting;
import lmu.utils.*;
import lmu.utils.plotting.PlotCreator;
import lmu.utils.tuple.Tuple3;
import lmu.utils.plotting.CachedPlotCreator;

import java.awt.image.BufferedImage;
import java.io.File;
import java.io.PrintWriter;
import java.util.*;
import java.util.function.Consumer;

import static lmu.utils.IteratorUtils.rangev;
import static lmu.utils.ObjectGetter.*;

public class MatrixInput {

    public static void main(String[] args) {
        empires.input.GeneralOptions generalOptions = new GeneralOptions();
        BackgroundProviderOption backgroundProviderOption = new BackgroundProviderOption();
        SimpleOptionParser cmd = new SimpleOptionParser("matrix", "labels", "nolog", "od", "condpairs", "nodiffexp", "nooverview", "minrep", "idfield", "fcthreshold",
                "fdrthreshold", "takebestn", "targetnumcpuniqfeatures", "writecenteredlog2fconcond", "noplots", "symbols2show", "Nmostsignifexamplestoplot",
                "subfeatures", "writedistribs", "gene2isoform2features", "doublediffvariant", "minfeature4splicing", "ignoreheaders");
        cmd.setFile("matrix", "labels", "condpairs", "symbols2show", "subfeatures", "gene2isoform2features", "ignoreheaders");
        cmd.setInt("minrep");
        cmd.setDefault("idfield", "id");
        cmd.setDefault("minrep", "3");
        cmd.setDouble("fcthreshold", "fdrthreshold");
        cmd.setDefault("fcthreshold", "0.9");
        cmd.setDefault("fdrthreshold", "0.05");
        cmd.setInt("takebestn", "targetnumcpuniqfeatures", "Nmostsignifexamplestoplot", "minfeature4splicing");
        cmd.setDefault("takebestn", "100");
        cmd.setDefault("minfeature4splicing", "1");
        cmd.setDefault("targetnumcpuniqfeatures", "5");
        cmd.setDefault("Nmostsignifexamplestoplot", "10");
        cmd.setDefault("doublediffvariant", empires.DoubleDiffVariant.QUICK_AND_DIRTY.getName());
        cmd.setDir("od");
        cmd.setSwitches("nolog", "nodiffexp", "nooverview", "noplots", "writedistribs");
        cmd.setOptional("condpairs", "writecenteredlog2fconcond", "symbols2show", "subfeatures", "gene2isoform2features", "ignoreheaders");


        if (!OptionParser.parseParams(args, true, false, true, true, cmd, backgroundProviderOption, generalOptions))
            return;

        generalOptions.apply();

        int minfeature4splicing = cmd.getInt("minfeature4splicing");
        empires.DoubleDiffVariant doubleDiffVariant = empires.DoubleDiffVariant.get(cmd.getValue("doublediffvariant"));
        HashMap<String, String> replicate2condition = FileUtils.readMap(cmd.getFile("labels"), "\t", 0, 1);
        HashMap<String, Vector<String>> condition2replicates = new HashMap<>();
        apply(replicate2condition.entrySet(), (_e) -> MapBuilder.updateV(condition2replicates, _e.getValue(), _e.getKey()));

        boolean noplots = cmd.isSet("noplots");
        double FCTHRESHOLD = cmd.getDouble("fcthreshold");
        double FDRTHRESHOLD = cmd.getDouble("fdrthreshold");
        int TAKEBESTN = cmd.getInt("takebestn");
        int TARGETNUM_CONDPAIR_UNIQUE_FEATURES = cmd.getInt("targetnumcpuniqfeatures");
        int NMOSTSIGNIF_TO_PLOT = cmd.getInt("Nmostsignifexamplestoplot");
        boolean nolog = cmd.isSet("nolog");

        File matrix = cmd.getFile("matrix");


        HashMap<String, Vector<String>> subfeatures = null;

        if(cmd.isOptionSet("subfeatures")) {
            HashMap<String, HashSet<String>> main2sub = new HashMap<>();
            String expformatErr = "invalid subfeature file expected format per line: mainfeaturename<tab>subfeauture";
            subfeatures = new HashMap<>();
            for(String[] sp : FileUtils.getFieldSetsIterable(cmd.getFile("subfeatures"), "\t")) {
                if(sp.length != 2) {
                    throw new FRuntimeException(expformatErr+ " record: ", toVector(sp));
                }
                MapBuilder.update(main2sub, sp[0], sp[1]);

            }
            subfeatures = buildMap(main2sub.entrySet(), (_e) -> _e.getKey(), (_e) -> toVector(_e.getValue()));
        }



        Iterator<String[]> it = FileUtils.getFieldSetIterator(matrix, "\t");
        Vector<String> header = map(it.next(), (_s) -> _s.trim());
        String idfield = cmd.getValue("idfield");

        HashSet<String> headerSet = toSet(header);
        if(filteredSize(replicate2condition.keySet(), (_n) -> !headerSet.contains(_n)) > 0) {
            throw new FRuntimeException("some replicates missing in matrix file: %s\nheaders: %s\n", SetInfo.minus(replicate2condition.keySet(), headerSet), headerSet);
        }
        Set<String> notyetmapped = SetInfo.minus(headerSet, replicate2condition.keySet());
        if(!notyetmapped.contains(idfield)) {
            System.err.printf("could not find idfield: %s in not replicated headers: %s! please provide a valid feature id field with -idfield!", idfield, notyetmapped);
            return;
        }

        /** read-in matrix data - skip idfield column */
        int idIdx = header.indexOf(idfield);
        Vector<String> features = new Vector<>();
        Vector<Vector<Double>> label2data = map(header, (_i) -> new Vector<>());

        HashSet<String> ignoreheaders = (cmd.isOptionSet("ignoreheaders")) ? FileUtils.readSet(cmd.getFile("ignoreheaders")) : new HashSet<>();
        Vector<String> alllabels = filter(header, (_h) -> !ignoreheaders.contains(_h) && !_h.equals(idfield));
        FileUtils.Converter converter = new FileUtils.Converter();
        while(it.hasNext()) {
            converter.set(it.next());
            features.add(converter.toString(idIdx));
            if(converter.length() != header.size())
                throw new FRuntimeException("invalid line: %s expected %d fields, got %d!", toVector(converter.getRaw()), converter.length(), header.size() + 1);

            for(int i=0; i<header.size(); i++) {
                if(idIdx == i || ignoreheaders.contains(header.get(i)))
                    continue;
                double v = 0.0;
                try {
                    v = converter.toDbl(i);
                } catch(Exception nfe) {
                    System.err.printf("feature:%s:%s\nvalue: %s:%s could not be converted to double! maybe a header to ignore? (-ignoreheaders)\n", header.get(idIdx), converter.toString(idIdx), header.get(i), converter.toString(i));
                    return;
                }
                v = (nolog) ? v : (v == 0.0) ? Double.NaN : NumUtils.logN(v, 2.0);
                label2data.get(i).add(v);
            }
        }

        HashMap<String, ReplicateSetInfo> conditions = new HashMap<>();

        int minrep = cmd.getInt("minrep");
        Logger log = LogConfig.getLogger();
        for(String cond : condition2replicates.keySet()) {

            Vector<String> replicates = condition2replicates.get(cond);
            if(replicates.size() < minrep) {
                log.warn("condition %s has too few replicates to analyse (%d), will be skipped!", cond, replicates.size());
                continue;
            }
            ReplicateSetInfo rsi = new ReplicateSetInfo(cond, replicates, features);
            if(subfeatures != null) {
                rsi.setCombinedFeatures(subfeatures);
            }

            for(int i=0; i<replicates.size(); i++) {
                String rep = replicates.get(i);
                int idx = first(filterIndex(header, (_s) -> rep.equals(_s)));
                //System.out.printf("add %s:%s -> %d\n", cond, rep, idx);
                rsi.setLog2Data(i, label2data.get(idx));
            }
            conditions.put(cond, rsi);
        }


        File od = cmd.getFile("od");
        PrintWriter pw = FileUtils.getWriter(od, "index.html");


        log.info("got %d conditons to check", conditions.size());

        File symbols = cmd.getOptionalFile("symbols2show");

        HashSet<String> ids2show = (symbols == null) ? new HashSet<>() : FileUtils.readSet(symbols);


        empires.BackgroundContextFuzzficationStrategyProvider backgroundContextFuzzficationStrategyProvider = backgroundProviderOption.getStrategy();
        HashMap<String, empires.NormalizedReplicateSet> allnormed = buildMap(conditions.entrySet(), (_e) -> _e.getKey(), (_e) -> new empires.NormalizedReplicateSet(_e.getValue(), backgroundContextFuzzficationStrategyProvider, true));
        HashMap<String, empires.NormalizedReplicateSet> normed = buildMap(filter(allnormed.keySet(), (_k) -> condition2replicates.get(_k).size() > 1), (_k) -> allnormed.get(_k));
        empires.plotting.DataOverview dataOverview = new DataOverview(normed);



        Vector<String> condvec = toVector(conditions.keySet());

        File condpairs = cmd.getOptionalFile("condpairs");
        Vector<UPair<String>> condpairs2test = getPairs(filter(condvec, (_c) -> normed.containsKey(_c)), true);
        if(condpairs != null) {
            condpairs2test = map(FileUtils.getFieldSetIterator(condpairs, "\t"), (_sp) -> UPair.createU(_sp[0], _sp[1]));
            HashSet<String> condset2test = new HashSet<>();
            apply(condpairs2test, (_p) -> condset2test.add(_p.getFirst()));
            apply(condpairs2test, (_p) -> condset2test.add(_p.getSecond()));
            condvec = toVector(condset2test);
        }
        NumUtils.sort(condvec, (_c) -> conditions.get(_c).getNumReplicates(), true);

        String centerCond = cmd.getOptionalValue("writecenteredlog2fconcond", condvec.get(0));

        String anchorSampleName = null;
        HashMap<String, Double> sample2shift = new HashMap<>();
        Vector<Double> anchorValues = null;

        if(centerCond != null) {
            ReplicateSetInfo rsi = conditions.get(centerCond);
            empires.Normalization norm = new empires.Normalization(rsi.getLog2Data());
            int idx = norm.getAnchorSampleIdx();
            anchorValues = rsi.getLog2Data().get(idx);

            anchorSampleName = rsi.getReplicateNames().get(idx);
            File log2fc2anchor = new File(od, String.format("log2fc2anchor_%s.tsv", anchorSampleName));
            PrintWriter logpw = FileUtils.getWriter(log2fc2anchor);
            Vector<String> samplenames = new Vector<>();
            Vector<Vector<Double>> fmatrix = new Vector<>();

            empires.ShiftedGroup ref = new empires.ShiftedGroup(0, anchorValues);

            for(ReplicateSetInfo r : conditions.values()) {

                for(int i=0; i<r.getNumReplicates(); i++) {
                    String sampleName = r.getReplicateNames().get(i);
                    if(sampleName.equals(anchorSampleName))
                        continue;

                    double shift = new empires.PairwiseMedianImpliedFCPeakErrorEstimation(ref, new empires.ShiftedGroup(1, r.getLog2Data().get(i))).shift;
                    sample2shift.put(sampleName, shift);
                    fmatrix.add(map(r.getLog2Data().get(i), (_d) -> _d + shift));
                    samplenames.add(sampleName);
                }
            }
            logpw.printf("%s\t%s\n", idfield, StringUtils.joinObjects("\t", samplenames));
            Vector<String> featureNames = rsi.getFeatureNames();
            String AS = anchorSampleName;
            for(int fIdx = 0; fIdx < featureNames.size(); fIdx++) {


                double base = anchorValues.get(fIdx);
                if(Double.isNaN(base))
                    continue;

                ReplicateSetInfo centerRsi = conditions.get(centerCond);

                String fn = featureNames.get(fIdx);
                Vector<Double> centerData = map(filter(rangev(centerRsi.getReplicateNames()), (_i) -> !centerRsi.getReplicateNames().get(_i).equals(AS)), (_i) -> centerRsi.getReplicateData(fn) .get(_i)- base);
                double centerShift = NumUtils.median(centerData);
                logpw.printf("%s", fn);

                for(int sIdx = 0; sIdx < samplenames.size(); sIdx++) {
                    logpw.printf("\t%.2f", fmatrix.get(sIdx).get(fIdx) - base - centerShift);
                }
                logpw.println();
            }

            logpw.close();

            pw.printf("<h2>anchor-normalized logfc-s</h2>\n");
            pw.printf("<a href='%s'>log2fc normalized matrix to anchor sample: %s : %s</a>\n", log2fc2anchor.getName(), centerCond, anchorSampleName);
            pw.flush();
        }


        HashMap<String, Vector<UPair<String>>> difffeature2condpairs = new HashMap<>();
        HashMap<UPair<String>, Vector<String>> cp2signifs = new HashMap<>();

        HashSet<String> alldiffs = new HashSet<>();

        double  MINFDR = Math.pow(10, -9);
        if(!cmd.isSet("nodiffexp")) {

            if(cmd.isSet("writedistribs")) {

                empires.NormalizedReplicateSet base = first(normed.values());

                Vector<Vector<Double>> normedvalues = map(features, (_d) -> new Vector<>());

                Vector<String> normedlabels = new Vector<>();

                for(empires.NormalizedReplicateSet n : normed.values()) {

                    normedlabels.addAll(n.getInData().getReplicateNames());
                    double shift = (n == base) ? 0 : new empires.DiffExpManager(base, n, features, new empires.EmpiRe()).shifts.get(1);

                    for(int i=0; i<features.size(); i++) {
                        empires.ReplicatedMeasurement rm = n.getNormed(features.get(i));
                        normedvalues.get(i).addAll(map(rm.replicates, (_v) -> _v + shift));
                    }
                    n.writeErrorDistribInfos(od);
                }
                PrintWriter npw = FileUtils.getWriter(od, "normed.tsv");
                npw.printf("feature\t%s\n", StringUtils.joinObjects("\t", normedlabels));
                for(int i=0; i<features.size(); i++) {
                    npw.printf("%s\t%s\n", features.get(i), StringUtils.joinObjects("\t", normedvalues.get(i), (_d) -> String.format("%.3f", _d)));
                }
                npw.close();
            }
            pw.printf("<h2>diffexp pairs for signif: abs(log2fc) >= %.2f  FDR <= %g</h2>\n", FCTHRESHOLD, FDRTHRESHOLD);

            HashMap<UPair<String>, Tuple3<Integer, Integer, File>> cp2sub = new HashMap<>();

            HashMap<String, Double> feature2totalsignif = new HashMap<>();

            for(UPair<String> cp : condpairs2test) {
                empires.EmpiRe emp = new EmpiRe();
                empires.NormalizedReplicateSet n1 = normed.get(cp.getFirst());
                empires.NormalizedReplicateSet n2 = normed.get(cp.getSecond());


                Vector<empires.DiffExpResult> res = emp.getDifferentialResults(n1, n2);
                if(res == null)
                    continue;


                Vector<empires.DiffExpResult> signif = filter(res, (_d) -> Math.abs(_d.estimatedFC) >= FCTHRESHOLD && _d.fcEstimateFDR <= FDRTHRESHOLD);

                NumUtils.sort(signif, (_d) -> _d.fcEstimatePval, false);


                apply(signif, (_d) ->
                    MapBuilder.update(feature2totalsignif, _d.combinedFeatureName, - NumUtils.logN(Math.max(MINFDR, _d.fcEstimateFDR), 10.0))
                );
                Vector<String> signif_totake = map(
                        ((TAKEBESTN <= 0) ? signif : VectorUtils.slice(signif, 0, Math.min(signif.size(), TAKEBESTN)))
                        , (_d) -> _d.combinedFeatureName);

                alldiffs.addAll(signif_totake);


                cp2signifs.put(cp, signif_totake);
                apply(signif_totake, (_n) -> MapBuilder.updateV(difffeature2condpairs, _n, cp));

                File cpod  = new File(od, String.format("%s_VS_%s", cp.getFirst(), cp.getSecond()));
                cp2sub.put(cp, Tuple3.create(filteredSize(signif, (_d) -> _d.estimatedFC < 0), filteredSize(signif, (_d) -> _d.estimatedFC > 0), cpod));

                cpod.mkdirs();
                empires.plotting.DiffExpTable diffExpTable = new DiffExpTable(n1, n2, res, FCTHRESHOLD, FDRTHRESHOLD);
                File table = new File(cpod, "diffexp.tsv");
                diffExpTable.getTable().writeCSV(table);

                File cp_html = new File(cpod, "index.html");

                //pw.printf("<li><a href='%s/%s'>%s vs %s</a></li>\n", cpod.getName(), cp_html.getName(), cp.getFirst(), cp.getSecond());
                pw.flush();
                PrintWriter cppw = FileUtils.getWriter(cp_html);
                cppw.printf("<h2>%s vs %s</h2>\n", cp.getFirst(), cp.getSecond());
                cppw.printf("<a href='%s'>diffexp table</a><br>\n", table.getName());

                File gene2isoform2features = cmd.getOptionalFile("gene2isoform2features");

                if(gene2isoform2features != null) {
                    File dstable = new File(cpod, "diffsplic.tsv");

                    new empires.FeatureBasedSplicingTest(n1, n2, gene2isoform2features, doubleDiffVariant, minfeature4splicing, null).getSplicingResultTable().writeCSV(dstable);
                    cppw.printf("<a href='%s'>diffsplic table</a><br>\n", dstable.getName());
                }

                if(cmd.isSet("writedistribs")) {

                    File distribout = new File(cpod, "fcdistribs.serialized");
                    new empires.NamedFoldChangeDistributionMap(res).serialize(distribout);

                    cppw.printf("<br><a href='%s'>distribs serialized</a><br>\n", distribout.getName());

                    if(subfeatures != null) {
                        empires.NamedFoldChangeDistributionMap subMap = new empires.NamedFoldChangeDistributionMap();
                        for(empires.DiffExpResult de : res) {
                            for(int i=0; i<de.featureNames.size(); i++) {

                                subMap.addDistrib(de.featureNames.get(i), de.perFeatureScaledDistributions.get(i));
                            }

                        }
                        File sub_distribout = new File(cpod, "subfeature_fcdistribs.serialized");
                        subMap.serialize(sub_distribout);

                        cppw.printf("<br><a href='%s'>subfeature distribs serialized</a><br>\n", sub_distribout.getName());

                    }
                }


                if(noplots) {
                    cppw.close();
                    continue;
                }


                File vulc = new File(cpod, "vulcano.png");

                File signal1 = new File(cpod, "signal1.png");
                File signal2 = new File(cpod, "signal2.png");

                ImageUtils.saveImage(diffExpTable.getVulcano(), vulc);

                ImageUtils.saveImage(new empires.plotting.NormalizedReplicateSetPlotting(n1).plotBackGroundDistribs(CachedPlotCreator.getPlotCreator()), signal1);
                ImageUtils.saveImage(new NormalizedReplicateSetPlotting(n2).plotBackGroundDistribs(CachedPlotCreator.getPlotCreator()), signal2);

                Vector<empires.DiffExpResult> sortedfiff = NumUtils.sort(signif, (_d) -> _d.fcEstimatePval, false);

                int totake = Math.min(NMOSTSIGNIF_TO_PLOT, sortedfiff.size());


                cppw.printf("<h2>vulcano</h2><img src='%s'><br>", vulc.getName());
                cppw.printf("<h2>signal2noise</h2><img src='%s'><br><img src='%s'><br>", signal1.getName(), signal2.getName());

                if(totake > 0) {
                    cppw.printf("<h2>%d most significant genes</h2><ol>\n", totake);
                    for(int i=0; i<totake; i++) {
                        empires.DiffExpResult diffExpResult = sortedfiff.get(i);
                        BufferedImage bim = diffExpTable.getDetailImage(CachedPlotCreator.getPlotCreator(), diffExpResult);
                        File of = new File(cpod, String.format("hit_%d.png", i));
                        ImageUtils.saveImage(bim, of);
                        cppw.printf("<li>%s<br><img src='%s'></li>\n", diffExpResult.combinedFeatureName, of.getName());
                    }

                }
                cppw.close();
            }

            log.info("processed %d condition pairs", cp2signifs.size());
            pw.flush();

            pw.printf("<h2>got %d diff features</h2>\n", difffeature2condpairs.size());

            pw.printf("<table border=1>");
            pw.printf("<tr><td></td>%s</tr>\n", StringUtils.join(" " , condvec, (_c) -> String.format("<th>%s</th>", _c)));
            for(int ci = 0; ci < condvec.size(); ci++) {
                String c1 = condvec.get(ci);
                pw.printf("<tr><th>%s</th>", c1);

                for(String c2 : condvec) {
                    UPair<String> cp = UPair.createU(c1, c2);
                    if(c1.equals(c2)) {
                        pw.printf("<td></td>");
                        continue;
                    }
                    Tuple3<Integer, Integer, File> t = cp2sub.get(cp);
                    if(t == null) {
                        t = cp2sub.get(UPair.swapU(cp));
                    }
                    if(t == null) {
                        pw.printf("<td>unchecked</td>");
                        continue;
                    }
                    pw.printf("<td><a href='%s/index.html'>up: %d dn: %d total: %d</a></td>", t.get2().getName(), t.get0(), t.get1(), t.get0() + t.get1());
                }
                pw.printf("</tr>\n");
            }

            pw.printf("</table>");
            pw.flush();

            log.info("will write most separating features");
            Vector<String> sortedFeatures = toSortedVector(difffeature2condpairs.keySet(), (_s) -> Pair.create(difffeature2condpairs.get(_s).size(), feature2totalsignif.get(_s)), false);

            HashMap<UPair<String>, Integer> max2take = buildMap(cp2signifs.entrySet(), (_e) -> _e.getKey(),
                                            (_e) -> Math.min(_e.getValue().size(), TARGETNUM_CONDPAIR_UNIQUE_FEATURES));

            HashMap<UPair<String>, Integer> cp2selected = new HashMap<>();

            String ANCHORSAMPLENAME = anchorSampleName;
            Vector<Double> ANCHORVALUES =  anchorValues;

            Vector<String> acondvec = toVector(centerCond);
            acondvec.addAll(filter(condvec, (_c) -> !_c.equals(centerCond)));

            HashMap<String, Vector<String>> subFMp = subfeatures;

            Consumer<String> plotSelected = (String f) ->
            {
                if(ANCHORSAMPLENAME == null || noplots )
                    return;



                double anchorCondShift = 0.0;

                Vector<File> imfiles = new Vector<>();

                Vector<String> fnames = (subFMp == null) ? toVector(f) : subFMp.get(f);
                for(String fname: fnames) {
                    PlotCreator.BoxPlotBuilder boxPlotBuilder = CachedPlotCreator.getPlotCreator().buildBoxPlot();
                    for(String cond : acondvec) {
                        ReplicateSetInfo rsi = conditions.get(cond);
                        int fIdx = rsi.getFeatureIdx(fname);
                        Vector<Double> cdata = new Vector<>();
                        for(int sIdx = 0; sIdx < rsi.getNumReplicates(); sIdx++) {
                            String sample = rsi.getReplicateNames().get(sIdx);
                            if(sample == ANCHORSAMPLENAME)
                                continue;

                            cdata.add(rsi.getLog2Data().get(sIdx).get(fIdx) + sample2shift.get(rsi.getReplicateNames().get(sIdx)) - ANCHORVALUES.get(fIdx));
                        }

                        if(cond.equals(centerCond)) {
                            anchorCondShift = NumUtils.median(cdata);
                        }
                        double S = anchorCondShift;
                        boxPlotBuilder.addBox(cond, NumUtils.getNumInfo(map(cdata, (_d) -> _d - S)).quantiles);
                    }

                    CachedPlotCreator.getPlotCreator().setTitle(fname);
                    CachedPlotCreator.getPlotCreator().setLabels("", String.format("log2FC to %s:%s", centerCond, ANCHORSAMPLENAME), null);
                    File bof = new File(od, "sel_"+ fname + ".png");
                    boxPlotBuilder.writePlotFile(bof);
                    imfiles.add(bof);
                }

                pw.printf("<br><table><tr>%s</tr></table><br>\n", StringUtils.joinObjects("", imfiles,
                        (_bof) -> String.format("<td><img src='%s'></td>\n", _bof.getName())));
                pw.flush();

            };
            pw.flush();



            File self = new File(od, "selected.features");
            pw.printf("<a href='%s'> most separating feature ids </a><br>\n", self.getName());
            pw.flush();

            PrintWriter selectedPw = FileUtils.getWriter(self);


            pw.printf("<ol>\n");
            for(String f : sortedFeatures) {
                Vector<UPair<String>> needed2cps = filter(difffeature2condpairs.get(f), (_cp) -> cp2selected.getOrDefault(_cp, 0) < max2take.get(_cp));
                if(needed2cps.size() == 0)
                    continue;

                pw.printf("<li>select: %s due to: %s  (signif in %s)\n", f, needed2cps, difffeature2condpairs.get(f));
                plotSelected.accept(f);

                pw.printf("</li>\n");
                selectedPw.printf("%s\n", f);
                apply(difffeature2condpairs.get(f), (_cp) -> MapBuilder.update(cp2selected, _cp));
            }
            selectedPw.close();

            log.info("most selective features written");

            pw.printf("</ol>\n");

            log.info("most separating features written\n");
            if(ids2show.size() > 0) {
                pw.printf("<h2>selected features</h2>\n");
                pw.printf("<ol>\n");
                for(String id : ids2show) {
                    pw.printf("<li>select: %s (signif in %s)\n", id, difffeature2condpairs.get(id));
                    plotSelected.accept(id);
                    pw.printf("</li>\n");

                }


                pw.printf("</ol>\n");


            }
        }


        if(!cmd.isSet("nooverview")) {
            dataOverview.write(od, pw);
        }

        pw.close();
        CachedPlotCreator.close();
    }
}
