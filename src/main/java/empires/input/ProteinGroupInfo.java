package empires.input;

import lmu.utils.*;
import lmu.utils.plotting.PlotCreator;

import java.awt.image.BufferedImage;
import java.io.File;

import java.util.*;

import static lmu.utils.IteratorUtils.rangev;
import static lmu.utils.ObjectGetter.*;


public class ProteinGroupInfo
{
    Logger log = LogConfig.getLogger();

    public enum IntensityType
    {
        INTENSITY("Intensity "),
        REPORTER("Reporter intensity "),
        LFQ("LFQ intensity "),
        SILAC("Ratio ")
        ;

        String prefix;

        IntensityType(String pref)
        {
            prefix = pref;
        }

        static IntensityType getType(Vector<String> headers)
        {
            for(IntensityType it : toVector(REPORTER, LFQ, SILAC, INTENSITY))
            {
                if(filteredSize(headers, (_h) -> _h.indexOf(it.prefix) >= 0) > 0)
                    return it;
            }
            throw new FRuntimeException("could not determine intensity types from headers: %s tested: %s", headers, map(IntensityType.values(), (_it) -> _it.prefix));
        }

        String getIntensityPrefix()
        {
            switch(this)
            {
                case SILAC:
                    return INTENSITY.prefix;
                default:
                    return prefix;
            }
        }

        String getPeptidePrefix()
        {
            switch(this)
            {
                case SILAC:
                case LFQ:
                    return INTENSITY.prefix;
                case REPORTER:
                    return "Reporter intensity ";
                default:
                    return prefix;
            }
        }

    }

    HashMap<String, Vector<Integer>> prot2pepids = new HashMap<>();
    Vector<UPair<String>> label2condition;

    HashSet<String> trues = null;
    HashSet<String> falses = null;

    HashMap<String, Vector<String>> cond2labels = new HashMap<>();


    public static HashMap<String, Vector<Integer>> getProtein2PeptideIndex(File pgroup) {
        HashMap<String, Vector<Integer>> prot2peps = new HashMap<>();

        apply(FileUtils.getHrIterator(pgroup, FileUtils.getHeaderedReader("\t").add("Protein IDs", "pid").add("Peptide IDs", "pepids")
                        .add("Only identified by site", "f1")
                        .add("Reverse", "f2")
                        .add("Potential contaminant", "f3")),
                (_hr) ->
                {
                    if(String.format("%s%s%s", _hr.getString("f1"), _hr.getString("f2"), _hr.getString("f3")).indexOf('+') >= 0)
                        return;

                    String id = _hr.getString("pid");
                    if(id.indexOf("REV__") >= 0|| id.indexOf("CON__") >= 0)
                        return;

                    prot2peps.put(id, map(_hr.getString("pepids").split(";"), (_s) -> Integer.parseInt(_s)));

                });

        return prot2peps;

    }

    public HashMap<String, Vector<String>> getCondition2ReplicateMap() {
        return cond2labels;
    }

    public void setTrueFalseLabels()
    {
        if(trues == null || falses == null)
            return;

        apply(id2info.values(), (_v) -> _v.setTrueFalse(trues, falses));



    }

    public void setTrues(HashSet<String> trues)
    {
        this.trues = trues;
    }

    public void setFalses(HashSet<String> falses)
    {
        this.falses = falses;
    }

    public void readPeptides(File peps)
    {
        readPeptides(peps, true);
    }

    HashMap<Integer, String> pepid2seq = new HashMap<>();

    public String getPepSeq(int id)
    {
        return pepid2seq.get(id);
    }

    HashMap<Integer, HashMap<String, Vector<Double>>> pep2cond2values = new HashMap<>();
    Vector<Pair<String, Integer>> samplelist = null;

    public Vector<Pair<String, Integer>> getSampleList()
    {
        return samplelist;
    }
    public Collection<HashMap<String, Vector<Double>>> getCond2PepValues()
    {
        return pep2cond2values.values();
    }

    public Collection<Integer> getPeptideIds()
    {
        return pep2cond2values.keySet();
    }

    public HashMap<String, Vector<Double>> getPeptideCond2Val(int pepidx) {
        return pep2cond2values.get(pepidx);
    }


    public Vector<Pair<String, Integer>> getSamplelist() {
        return samplelist;
    }


    public ReplicateSetInfo getProteinToPeps(String condition, boolean use_mod_peps_too) {

        Vector<Vector<Double>> peptide2replicates = new Vector<>();
        Vector<String> peptideNames = new Vector<>();
        HashMap<String, Integer> peptideName2Usage = new HashMap<>();
        HashMap<String, Vector<String>> protein2peptides = new HashMap<>();


        for(Id2QInfo idq : id2info.values()) {
            String proteinId = idq.id;
            Vector<Vector<Integer>> idxll = toSingleVector(idq.pepids);
            if (idq.modpepids != null && use_mod_peps_too) {
                idxll.add(idq.modpepids);
            }
            Vector<String> protein_peptide_list = new Vector<>();
            for (Vector<Integer> idxl : idxll) {
                for (int idx : idxl) {
                    if (null == pep2cond2values.get(idx))
                        continue;

                    String peptideSequence = pepid2seq.get(idx);

                    if (null != restrict_to_peptides && !restrict_to_peptides.contains(peptideSequence))
                        continue;

                    int numUsage = peptideName2Usage.getOrDefault(peptideSequence, 0);
                    peptideName2Usage.put(peptideSequence, 1 + numUsage);
                    String correctedPeptideId = (numUsage == 0) ? peptideSequence : peptideSequence + "." + (1 + numUsage);
                    peptideNames.add(correctedPeptideId);
                    protein_peptide_list.add(correctedPeptideId);
                    peptide2replicates.add(pep2cond2values.get(idx).get(condition));

                }
            }
            protein2peptides.put(proteinId, protein_peptide_list);
        }

        Vector<String> replicateNames = cond2labels.get(condition);
        System.out.printf("%s: %s\n", condition, replicateNames);

        //Vector<String> replicateNames = mapIndex(first(id2info.values()).cond2values.get(condition).size(), (_i) -> "rep"+ (_i + 1));
        ReplicateSetInfo rsi = new ReplicateSetInfo(condition, replicateNames, peptideNames);
        Double[][] logMatrix = new Double[replicateNames.size()][peptideNames.size()];

        for(int i=0; i<peptideNames.size(); i++) {

            Vector<Double> logData = map(peptide2replicates.get(i), ((_d) -> (_d == 0.0) ? Double.NaN : NumUtils.logN(_d, 2.0)));
            for(int j=0; j<logData.size(); j++) {
                logMatrix[j][i] = logData.get(j);
            }

        }

        for(int IDX : rangev(replicateNames.size())) {
            rsi.setLog2Data(IDX, map(rangev(peptideNames.size()), (_i) -> logMatrix[IDX][_i]));
        }

        rsi.setCombinedFeatures(protein2peptides);

        return rsi;

    }



    HashSet<String> restrict_to_peptides = null;

    public void restrictToPeptides(Collection<String> peps)
    {
        restrict_to_peptides = toSet(peps);
    }

    public void readPeptides(File peps, boolean do_filter)
    {
        FileUtils.HeaderedReader hr = new FileUtils.HeaderedReader("\t").add("Sequence", "seq").add("id");
        apply(label2condition, (_up) -> hr.add(type.getPeptidePrefix() + _up.getFirst(), _up.getFirst()));


        apply(FileUtils.getHrIterator(peps, hr),
                (_hr) ->
                {
                    HashMap<String, Vector<Double>> cond2values = new HashMap<>();

                    int id = _hr.getInt("id");

                    String seq = _hr.getString("seq");

                    pepid2seq.put(id, seq);

                    for(UPair<String> up : label2condition) {
                        MapBuilder.updateV(cond2values, up.getSecond(), _hr.getDbl(up.getFirst()));
                    }
                    pep2cond2values.put(id, cond2values);
                });

    }

    //can read peptide infos for it

    static class SilacConfig
    {
        double MAXALLOWED_ERROR = 0.5;

        // eg. <wt,1>,<ko,2>
        HashMap<String, UPair<Pair<String, Integer>>> silacheader2replicateinfos = new HashMap<>();

        HashMap<Pair<String, Integer>, String> rep2intensity_header = new HashMap<>();

        HashMap<UPair<Pair<String, Integer>>, Vector<Double>> reppair2errors = new HashMap<>();


        HashMap<String, Vector<Double>> readCond2ReplicateValues(FileUtils.HeaderedReader hr)
        {
            // format (cond1, rep1) / (cond2, rep2)
            //System.out.printf("intensity check : %s\n", map(rep2intensity_header.entrySet(), (_h) -> String.format("%s. %s\n", _h.getKey(), hr.getString(_h.getValue()))));
            HashMap<Pair<String, Integer>, Double> intensities = buildMap(rep2intensity_header.entrySet(),
                    (_e) -> _e.getKey(), (_e) -> hr.getDbl(_e.getValue()));

            //System.out.printf("next silac record got intensities: %s\n", intensities);
            //System.out.printf("corresponding ratios: %s\n", map(silacheader2replicateinfos.entrySet(), (_e) -> String.format("%s : %s",
            //        _e.getValue(), hr.getDbl(_e.getKey()))));

            Vector<Pair<String, Integer>> replicates = toSortedVector(rep2intensity_header.keySet(), true);
            Vector<Double> signals = map(replicates, (_r) -> hr.getDbl(rep2intensity_header.get(_r)));

            FoldChangeOptimizer fco = new FoldChangeOptimizer(replicates.size());

            HashMap<UPair<Pair<String, Integer>>, String> rpair2header = reverseMap(silacheader2replicateinfos);

            HashMap<UPair<Integer>, Double> rp2silacfc = new HashMap<>();

            for(int i=0; i<replicates.size(); i++)
            {
                Pair<String, Integer> r1 = replicates.get(i);
                if(signals.get(i) == 0.0 || Double.isNaN(signals.get(i)))
                    continue;

                for(int j=i+1; j<replicates.size(); j++)
                {
                    Pair<String, Integer> r2 = replicates.get(j);
                    if(signals.get(j) == 0.0 || Double.isNaN(signals.get(j)))
                        continue;

                    if(r1.getFirst().equals(r2.getFirst())) //same condition
                    {
                        fco.setTarget(i, j, 0.0);
                        fco.addSignals(i, j, signals.get(i), signals.get(j), 1.0);
                        continue;
                    }
                    String ratioheader1 = rpair2header.get(UPair.createU(r1, r2));
                    String ratioheader2 = rpair2header.get(UPair.createU(r2, r1));

                    if(ratioheader1 == null && ratioheader2 == null)
                        continue; //not compared

                    fco.addSignals(i, j, signals.get(i), signals.get(j), 1.0);

                    double ratio = hr.getDbl((ratioheader1 == null) ? ratioheader2 : ratioheader1);
                    if(Double.isNaN(ratio))
                        continue;

                    double fc = ((ratioheader1 == null) ? -1.0 : 1.0) * NumUtils.logN(ratio, 2.0);
                    rp2silacfc.put(UPair.createU(i,j), fc);
                    fco.setTarget(i, j, fc);

                }
            }
            if(rp2silacfc.size() == 0)
                return null;

            String id = hr.getString("ids");
            //System.out.printf("start opt for %s\n", id);
            Vector<Double> factors = fco.optimize();
            if(!fco.finishedOptimization()) {
                //System.out.printf("opt unfinished for %s\n", id);
                factors = fco.optimize(10000);
                if(!fco.finishedOptimization()) {
                    System.err.printf("cannot optimize for %s", id);
                }
            }
            //System.out.printf("end opt for %s on signals: %s fcs: %s\n", id, signals, rp2silacfc);
            if(factors.size() != replicates.size())
                return null;


            double nf = NumUtils.mean(factors);
            Vector<Double> nfactors = map(factors, (_f) -> _f / nf);

            Vector<Double> normedsignals = mapIndex(signals.size(), (_i) -> signals.get(_i) * nfactors.get(_i));

            Vector<Double> errors = new Vector<>();

            for(UPair<Integer> up : rp2silacfc.keySet()) {

                double req = rp2silacfc.get(up);

                double res = NumUtils.logN(normedsignals.get(up.getFirst()) / normedsignals.get(up.getSecond()), 2.0);

                double err = req -res;
                MapBuilder.updateV(reppair2errors, UPair.createU(replicates.get(up.getFirst()), replicates.get(up.getSecond())), err);


                errors.add(Math.abs(err));
                //System.out.printf("%s vs %s requested: %2.3f got %2.3f\n", replicates.get(up.getFirst()), replicates.get(up.getSecond()),  req, res);
            }

            if(errors.size() == 0 || NumUtils.max(errors) > MAXALLOWED_ERROR)
                return null;

            HashMap<String, Vector<Double>> rv = new HashMap<>();

            for(int i=0; i<replicates.size(); i++)
            {
                MapBuilder.updateV(rv, replicates.get(i).getFirst(), normedsignals.get(i));
            }
            return rv;
        }

    }

    SilacConfig sc = null;

    public boolean isSILAC()
    {
        return sc != null;
    }

    public BufferedImage showSILACFitErrors()
    {
        if(sc == null)
            return null;

        PlotCreator pc = new PlotCreator();

        Vector<BufferedImage> bims = new Vector<>();
        for(int abs = 0; abs < 2; abs++) {
            for (UPair<Pair<String, Integer>> up : sc.reppair2errors.keySet()) {
                String label = String.format("%s.%d vs %s.%d", up.getFirst().getFirst(), up.getFirst().getSecond(), up.getSecond().getFirst(), up.getSecond().getSecond());

                Vector<Double> v = sc.reppair2errors.get(up);
                Vector<Double> toplot = (abs == 0) ? v : map(v, (_v) -> Math.abs(_v));

                pc.cumhist(label, toplot, toplot.size(), false, true);
            }

            String x = (abs == 0) ? "FC fit error" : "absolute FC fit error";

            if(abs == 0) {
                pc.setLimits(-1.0, 1.0, null, null);
            } else
            {
                pc.setLimits(0.0, 1.0, null, null);
            }

            pc.setLabels(x, "freq", "bottomright");
            bims.add(pc.getImage());
        }
        pc.destroy();
        return ImageUtils.concat(bims);
    }


    IntensityType type = null;

    public static class Id2QInfo
    {
        public String id;
        public int pgroup_id = -1;
        public boolean is_true = false;
        public boolean only_identified_per_site = false;
        public boolean potential_contaminant = false;
        public boolean is_reverse = false;
        public boolean is_reverse_or_contaminant_id = false;
        public boolean multi_org = false;

        public static boolean checkPlus(String s)
        {
            return s != null && s.length() >0 && s.charAt(0) == '+';
        }

        HashMap<String, Vector<Double>> cond2values;
        Vector<Integer> pepids;
        Vector<Integer> modpepids;

        Vector<Integer> allpepids;

        public boolean toFilter()
        {
            return only_identified_per_site || potential_contaminant || is_reverse_or_contaminant_id || is_reverse || multi_org;
        }




        public Id2QInfo(FileUtils.HeaderedReader hr, HashMap<String, Vector<Double>> cond2values)
        {
            pgroup_id = hr.getInt("pgid");
            this.cond2values = cond2values;
            this.pepids = map(hr.getString("pepids").split(";"), (_s) -> Integer.parseInt(_s));
            this.modpepids = map(hr.getString("modpepids").split(";"), (_s) -> Integer.parseInt(_s));
            this.allpepids = toVector(toSet(VectorUtils.join(pepids, modpepids)));
            only_identified_per_site = checkPlus(hr.getString("osite"));
            potential_contaminant = checkPlus(hr.getString("pcont"));
            is_reverse = checkPlus(hr.getString("isrev"));
            id = hr.getString("ids");
            is_reverse_or_contaminant_id = id.indexOf("CON__") >= 0 || id.indexOf("REV__") >= 0;

        }

        public int getNumPepIds() {
            return pepids.size();
        }

        public boolean setTrueFalse(HashSet<String> trues, HashSet<String> falses)
        {
            Vector<String> ids = toVector(id.split(";"));
            int o1 = filteredSize(ids, (_i) -> trues.contains(_i));
            int o2 = filteredSize(ids, (_i) -> falses.contains(_i));

            multi_org = (multi_org || (o1 != ids.size() && o2 != ids.size()));
            return is_true = o1 == ids.size();

        }

        public boolean checkMultiOrg(HashSet<String> org1, HashSet<String> org2)
        {
            Vector<String> ids = toVector(id.split(";"));
            int o1 = filteredSize(ids, (_i) -> org1.contains(_i));
            int o2 = filteredSize(ids, (_i) -> org2.contains(_i));

            return multi_org = (o1 != ids.size() && o2 != ids.size());

        }
    }

    public HashMap<String, Id2QInfo> id2info = new HashMap<>();



    public ProteinGroupInfo(File pgroup, HashMap<String, String> label2condition)
    {
        this.label2condition = toSortedVector(map(label2condition.entrySet(), (_e) -> UPair.createU(_e.getKey(), _e.getValue())), true);

        HashMap<String, Vector<String>> c2labels = new HashMap<>();
        apply(this.label2condition, (_e) -> MapBuilder.updateV(c2labels, _e.getSecond(), _e.getFirst()));

        Vector<String> conds = toSortedVector(c2labels.keySet(), true);
        samplelist = new Vector<>();
        for(String c : conds) {
            int nrep = c2labels.get(c).size();
            applyIndex(nrep, (_i) -> samplelist.add(Pair.create(c, _i)));
        }


        Vector<String> headers = FileUtils.getHeaders(pgroup);

        type = IntensityType.getType(headers);

        //System.out.printf("type: %s got headers: %s\n", type, headers);
        FileUtils.HeaderedReader hr = FileUtils.getHeaderedReader("\t").add("Protein IDs", "ids").add("Peptide IDs", "pepids").add("id", "pgid");
        Vector<String> CONT = toVector("Potential contaminant", "Contaminant");
        for(String c : CONT) {
            if(null == filterOne(headers, (_h) -> _h.equals(c)))
                continue;

            hr.add(c, "pcont");

        }

        hr.add("Mod. peptide IDs", "modpepids");

        hr.add("Only identified by site", "osite");
        hr.add("Reverse", "isrev");

        HashSet<String> hset = toSet(headers);
        apply(label2condition.keySet(), (_l) -> hr.add(type.getIntensityPrefix() + _l, _l));


        if(type == IntensityType.SILAC)
        {
            System.out.printf("got SILAC\n");
            sc = new SilacConfig();
            String NORMALIZED = " normalized ";
            Vector<String> ratio_labels = filter(headers, (_h) -> _h.startsWith(type.prefix) && _h.indexOf(NORMALIZED) >= 0);

            HashMap<String,  HashMap<String, Integer>> cond2label2replidx = buildMap(toSet(label2condition.values()), (_s) -> new HashMap<>());

            for(String h : ratio_labels)
            {
                String[] info = h.substring(type.prefix.length()).split(NORMALIZED);
                String[] mass_labels = info[0].split("/");
                String sample_label = info[1];

                String ml1 = mass_labels[0];
                String ml2 = mass_labels[1];
                String label1 =  ml1 + " " + sample_label;
                String cond1 = label2condition.get(label1);
                if(cond1 == null)
                    throw new FRuntimeException("no mapping found for silac label >%s<!", label1);

                String label2 = ml2 + " " + sample_label;
                String cond2 = label2condition.get(label2);

                if(cond2 == null)
                    throw new FRuntimeException("no mapping found for silac label >%s<!", label2);


                Integer replidx1 = cond2label2replidx.get(cond1).get(label1);
                if(replidx1 == null) {
                    cond2label2replidx.get(cond1).put(label1, replidx1 = cond2label2replidx.get(cond1).size() + 1);
                    sc.rep2intensity_header.put(Pair.create(cond1, replidx1), label1);
                }

                Integer replidx2 = cond2label2replidx.get(cond2).get(label2);
                if(replidx2 == null) {
                    cond2label2replidx.get(cond2).put(label2, replidx2 = cond2label2replidx.get(cond2).size() + 1);
                    sc.rep2intensity_header.put(Pair.create(cond2, replidx2), label2);
                }
                //info[]

                hr.add(h);
                sc.silacheader2replicateinfos.put(h, UPair.createU(Pair.create(cond1, replidx1), Pair.create(cond2, replidx2)));
                //System.out.printf("%s = %s.%d / %s = %s.%d\n", label1, cond1, replidx1, label2, cond2, replidx2);
            }
        }


        apply(label2condition.entrySet(), (_e) -> MapBuilder.updateV(cond2labels, _e.getValue(), _e.getKey()));

        //condpairs = (condpairs != null) ? condpairs :   filterToSet(getPairs(toSortedVector(label2condition.values(), true), true), (_up) -> !_up.getFirst().equals(_up.getSecond()));


        Iterator<FileUtils.HeaderedReader> it = FileUtils.getHrIterator(pgroup, hr);
        while(it.hasNext())
        {
            FileUtils.HeaderedReader _hr = it.next();
            HashMap<String, Vector<Double>> cond2values = new HashMap<>();

            String id = _hr.getString("ids");

            id2info.put(id, new Id2QInfo(_hr, cond2values));


            //prot2pepids.put(id, pepids);

            if(sc != null) {
                HashMap<String, Vector<Double>> sc2v = sc.readCond2ReplicateValues(_hr);
                if(sc2v== null) {
                    System.err.printf("SILAC ratio optimizing is not possible for %s\n", _hr.getString("ids"));
                    continue;
                }
                cond2values.putAll(sc2v);
            }
            else
            {

                for(UPair<String> up : this.label2condition)
                {
                    MapBuilder.updateV(cond2values, up.getSecond(), _hr.getDbl(up.getFirst()));
                }
            }



        }
    }
}
