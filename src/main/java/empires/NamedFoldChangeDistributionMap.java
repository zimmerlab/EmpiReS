package empires;

import lmu.utils.*;
import lmu.utils.fdr.BenjaminiHochberg;

import java.io.*;
import java.nio.ByteBuffer;
import java.util.*;

import static lmu.utils.ObjectGetter.*;

public class NamedFoldChangeDistributionMap {

    HashMap<String, ErrorEstimationDistribution> name2foldchangedistrib = null;

    public NamedFoldChangeDistributionMap(Vector<DiffExpResult> diffExpResuls) {

        name2foldchangedistrib = buildMap(filter(diffExpResuls, (_de) -> _de.combinedEmpiricalFoldChangeDistrib != null),
                (_de) -> _de.combinedFeatureName, (_de) -> _de.combinedEmpiricalFoldChangeDistrib);

    }
    public NamedFoldChangeDistributionMap() {
        name2foldchangedistrib = new HashMap<>();
    }

    public void addDistrib(String name, ErrorEstimationDistribution foldchangedistrib) {
        name2foldchangedistrib.put(name, foldchangedistrib);
    }

    public NamedFoldChangeDistributionMap(File f) {
        try
        {
            name2foldchangedistrib = new HashMap<>();
            DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(f)));

            while(true) {
                String name = null;
                try {
                    name = FileUtils.readString(dis);
                } catch (Exception ie) {
                    break;
                }

                name2foldchangedistrib.put(name, ErrorEstimationDistribution.read(dis));
            }

        } catch(IOException e) {
            throw new FRuntimeException("i/o error: %s while reading %s into memory", e, e.getMessage(), f.getAbsolutePath());
        }

    }
    public void serialize(File outfile) {

        try
        {
            DataOutputStream dos = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outfile)));
            for(Map.Entry<String, ErrorEstimationDistribution> dist : name2foldchangedistrib.entrySet()) {

                FileUtils.writeString(dist.getKey(), dos);
                dist.getValue().writeToBuffer(dos);
            }
            dos.close();

        } catch(IOException ie) {
            throw new FRuntimeException("i/o error: %s while writing %d named fc distribs to file: %s", ie, ie.getMessage(), name2foldchangedistrib.size(), outfile.getAbsolutePath());
        }

    }


    static class FcRangeTestDef {
        double[] range;
        String name;
    }

    static class ConfidenceIntervalTestDef {
        String name;
        int percent;
        boolean minmax = false;
        boolean median = false;
        boolean peak = false;
        double alpha;
    }

    static class TestInfos {
        String name;
        HashMap<String, UPair<Double>> fcrange2result = new HashMap();
        HashMap<String, UPair<Double>> ci2result = new HashMap<>();

        public TestInfos(String n) {
            this.name = n;
        }
    }


    public DataTable doTest(Vector<FcRangeTestDef> fctests, Vector<ConfidenceIntervalTestDef> citests) {
        apply(citests, (_t) -> _t.alpha = ((100 - _t.percent) * 0.5) / 100.0);

        HashMap<String, TestInfos> tests = buildMap(name2foldchangedistrib.keySet(), (_s) -> new TestInfos(_s));

        for(ConfidenceIntervalTestDef ctd : citests) {
            for(Map.Entry<String, ErrorEstimationDistribution> e : name2foldchangedistrib.entrySet()) {
                double start = Double.NaN;
                double end = Double.NaN;

                if(ctd.peak) {
                    start = end = e.getValue().getMostProbableFcWindowCenter();
                }
                if(ctd.median) {
                    start = end = e.getValue().getFoldChangeToCumulativeFrequency(0.5);
                }
                if(ctd.minmax) {
                    start = e.getValue().getMinFC();
                    end = e.getValue().getMaxFC();
                }
                if(Double.isNaN(start)) {
                    start = e.getValue().getFoldChangeToCumulativeFrequency(ctd.alpha);
                    end = e.getValue().getFoldChangeToCumulativeFrequency(1.0 - ctd.alpha);
                }
                tests.get(e.getKey()).ci2result.put(ctd.name, UPair.createU(start, end));
            }
        }

        Vector<String> names = toVector(name2foldchangedistrib.keySet());
        for(FcRangeTestDef fcRangeTestDef : fctests) {
            Vector<UPair<Double>> data = map(names, (_n) ->
            {
                ErrorEstimationDistribution ed = name2foldchangedistrib.get(_n);
                double pmass = ed.getCumulativeFrequencyToFoldChange(fcRangeTestDef.range[1]) - ed.getCumulativeFrequencyToFoldChange(fcRangeTestDef.range[0]);
                return UPair.createU(pmass, 1.0);
            });
            BenjaminiHochberg.adjust_pvalues(data, (_d) -> 1.0 - _d.getFirst(), (_p) -> _p.getFirst().setSecond(_p.getSecond()));

            applyIndex(names.size(), (_i) -> tests.get(names.get(_i)).fcrange2result.put(fcRangeTestDef.name, data.get(_i)));
        }

        DataTable.HeaderGetterManager<TestInfos> hgm = DataTable.buildHeader("id", (TestInfos t) -> t.name);

        for(ConfidenceIntervalTestDef ctd : citests) {
            hgm.add(ctd.name+".fcstart", (t) -> t.ci2result.get(ctd.name).getFirst());
            hgm.add(ctd.name+".fcend", (t) -> t.ci2result.get(ctd.name).getSecond());
        }
        for(FcRangeTestDef fcRangeTestDef : fctests) {
            hgm.add(fcRangeTestDef.name+".pmass", (t) -> t.fcrange2result.get(fcRangeTestDef.name).getFirst());
            hgm.add(fcRangeTestDef.name+".fdr", (t) -> t.fcrange2result.get(fcRangeTestDef.name).getSecond());
        }

        return DataTable.buildTable(toVector(tests.values()), hgm);
    }

    static final String[] FCRANGEHEADERS = new String[]{"name", "fcstart", "fcend"};

    static void checkFileHeaders(File f, String[] headers) {
        if(f == null)
            return;

        HashSet<String> gotheaders = toSet(FileUtils.getHeaders(f));
        Set<String> missing = SetInfo.minus(toSet(headers), gotheaders);
        if(missing.size() > 0)
            throw new FRuntimeException("missing headers in %s : %s headers expected: %s got: %s\n", f.getAbsolutePath(), missing, Arrays.toString(FCRANGEHEADERS), gotheaders);

    }
    static Vector<FcRangeTestDef> readFcRangeTests(File f) {
         if(f == null)
            return new Vector<>();

        checkFileHeaders(f, FCRANGEHEADERS);

        FileUtils.HeaderedReader hr = FileUtils.getHeaderedReader("\t");
        applyIndex(FCRANGEHEADERS.length, (_i) -> hr.add(FCRANGEHEADERS[_i]));

        Vector<FcRangeTestDef> rangeTestDefs = new Vector<>();

        apply(FileUtils.getHrIterator(f, hr),
                (_hr) ->
                {
                    FcRangeTestDef def = new FcRangeTestDef();
                    def.name = _hr.getString(FCRANGEHEADERS[0]);
                    System.out.printf("next: %s\n", _hr.getMap());
                    def.range = new double[]{_hr.getDbl(FCRANGEHEADERS[1]), _hr.getDbl(FCRANGEHEADERS[2])};
                    System.out.printf("got: %s\n", Arrays.toString(def.range));
                    rangeTestDefs.add(def);
                });
        return rangeTestDefs;
    }

    static final String[] CONFINTERVALHEADERS = new String[]{"name", "ciwidth"};

    static Vector<ConfidenceIntervalTestDef> readConfidenceIntervalTests(File f) {
        if(f == null)
            return new Vector<>();

        checkFileHeaders(f, CONFINTERVALHEADERS);

        FileUtils.HeaderedReader hr = FileUtils.getHeaderedReader("\t");
        applyIndex(CONFINTERVALHEADERS.length, (_i) -> hr.add(CONFINTERVALHEADERS[_i]));

        Vector<ConfidenceIntervalTestDef> confidenceIntervalTestDefs = new Vector<>();

        apply(FileUtils.getHrIterator(f, hr),
                (_hr) ->
                {
                    ConfidenceIntervalTestDef def = new ConfidenceIntervalTestDef();
                    def.name = _hr.getString(CONFINTERVALHEADERS[0]);
                    String dname = _hr.getString(CONFINTERVALHEADERS[1]);
                    def.minmax = dname.equals("minmax");
                    def.median = dname.equals("median");
                    def.peak = dname.equals("peak");

                    if(!def.minmax && !def.median && !def.peak) {
                        def.percent = _hr.getInt(CONFINTERVALHEADERS[1]);
                        if(def.percent <= 0 || def.percent >= 100)
                            throw new FRuntimeException("invalid entry: %s the %s should be 0 < < 100", _hr.getMap(), CONFINTERVALHEADERS[1]);

                    }
                    confidenceIntervalTestDefs.add(def);
                });
        return confidenceIntervalTestDefs;

    }

    public static void main(String[] args) {
        SimpleOptionParser cmd = new SimpleOptionParser("name2dists", "fcranges", "confintervals", "o");

        cmd.setFile("name2dists", "fcranges", "confintervals");
        cmd.setOutFile("o");
        cmd.setOptional("fcranges", "confintervals");

        if(!OptionParser.parseParams(args, true, false, true, true, cmd))
            return;


        Vector<ConfidenceIntervalTestDef> confTests = readConfidenceIntervalTests(cmd.getOptionalFile("confintervals"));
        Vector<FcRangeTestDef> fcRangeTests = readFcRangeTests(cmd.getOptionalFile("fcranges"));

        if(confTests.size() == 0 && fcRangeTests.size() == 0) {
            System.out.printf("no tests given (neither fcrange nor confinterval) - exiting\n");
            return;
        }

        NamedFoldChangeDistributionMap namedFoldChangeDistributionMap = new NamedFoldChangeDistributionMap(cmd.getFile("name2dists"));
        namedFoldChangeDistributionMap.doTest(fcRangeTests, confTests).writeCSV(cmd.getFile("o"));
    }
}
