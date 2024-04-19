package nlEmpiRe.input;

import org.apache.logging.log4j.Logger;
import lmu.utils.*;
import lmu.utils.tuple.Tuple3;
import nlEmpiRe.rnaseq.ReducedTranscriptPresentation;
import nlEmpiRe.rnaseq.TranscriptPairInfo;

import java.util.*;

import static lmu.utils.ObjectGetter.*;

public class RNASeqSplicingInfo {

    Logger log = LogConfig.getLogger();
    int MAX_TR_NUM;
    double MIN_AVG_COUNT_ON_STRONGER_SIDE = 10;
    int minCountIntMaxCountCondition = 0;


    HashMap<String, ReducedTranscriptPresentation> gene2transcriptInfo = new HashMap<>();
    HashMap<String, Vector<TranscriptPairInfo>> gene2pairtests = new HashMap<>();
    HashMap<String, Vector<Pair<Tuple, String>>> gene2EQClassFeature = new HashMap<>();
    Vector<String> genesForSplicingTesting = new Vector<>();
    Vector<String> allGenes = new Vector<>();


    Vector<String> allGeneFeatureNames = new Vector<>();
    Vector<Tuple> eqClasses = new Vector<>();


    HashMap<String, Vector<Vector<Double>>> cond2replicate2featureLogData = new HashMap<>();
    HashMap<String, Vector<String>> condition2replicatenames;

    HashMap<String, Integer> feature2Idx;


    HashMap<String, ReplicateSetInfo> replicateSetNameLookup = null;
    Vector<ReplicateSetInfo> replicateSetInfos = null;

    private RNASeqSplicingInfo() {

    }


    public RNASeqSplicingInfo(int MAX_TR_NUM, double MIN_AVG_COUNT_ON_STRONGER_SIDE, HashMap<String, Vector<String>> condition2replicatenames) {
        this.MAX_TR_NUM = MAX_TR_NUM;
        this.MIN_AVG_COUNT_ON_STRONGER_SIDE = MIN_AVG_COUNT_ON_STRONGER_SIDE;
        this.condition2replicatenames = condition2replicatenames;

        for(String cond : condition2replicatenames.keySet()) {
            cond2replicate2featureLogData.put(cond, map(condition2replicatenames.get(cond), (_r) -> new Vector<>()));
        }
    }

    public RNASeqSplicingInfo restrictToMinCountIntMaxCountCondition(int minCountIntMaxCountCondition) {
        RNASeqSplicingInfo restricted = new RNASeqSplicingInfo();
        restricted.MAX_TR_NUM = MAX_TR_NUM;
        restricted.MIN_AVG_COUNT_ON_STRONGER_SIDE = MIN_AVG_COUNT_ON_STRONGER_SIDE;

        restricted.setMinCountIntMaxCountCondition(minCountIntMaxCountCondition);
        restricted.condition2replicatenames = condition2replicatenames;
        restricted.gene2transcriptInfo = gene2transcriptInfo;

        for(String cond : condition2replicatenames.keySet()) {
            restricted.cond2replicate2featureLogData.put(cond, map(condition2replicatenames.get(cond), (_r) -> new Vector<>()));
        }


        log.info("set maxmin: %d\n", minCountIntMaxCountCondition);
        int nskippedFeatures = 0;
        int nskippedGenes = 0;
        for(String geneId : gene2EQClassFeature.keySet()) {
            Vector<Pair<Tuple, String>> usedFeatures = gene2EQClassFeature.get(geneId);
            Vector<Pair<Tuple, String>> restrictedUsedFeatures = new Vector<>();

            for(Pair<Tuple, String> feature : usedFeatures) {

                int fIdx = feature2Idx.get(feature.getSecond());
                int maxMinCount = 0;
                for(String cond : condition2replicatenames.keySet()) {
                    Vector<Double> logdata = map(cond2replicate2featureLogData.get(cond), (_v) -> _v.get(fIdx));
                    Vector<Double> data = map(logdata, (_d) -> (Double.isNaN(_d) ? 0.0 : Math.pow(2.0, _d)));
                    maxMinCount = Math.max(maxMinCount, NumUtils.min(data).intValue());
                }

                if(maxMinCount < minCountIntMaxCountCondition) {
                    nskippedFeatures++;
                    continue;
                }


                String fn = geneId+"."+(restrictedUsedFeatures.size()+1);
                restrictedUsedFeatures.add(Pair.create(feature.getFirst(), fn));

                for(String cond : condition2replicatenames.keySet()) {
                    Vector<Vector<Double>> toWriteIn = restricted.cond2replicate2featureLogData.get(cond);
                    Vector<Double> logdata = map(cond2replicate2featureLogData.get(cond), (_v) -> _v.get(fIdx));
                    for(int repIdx = 0; repIdx < toWriteIn.size(); repIdx++) {
                        toWriteIn.get(repIdx).add(logdata.get(repIdx));
                    }
                }

            }

            if(restrictedUsedFeatures.size() == 0) {
                nskippedGenes++;
                continue;
            }


            restricted.gene2EQClassFeature.put(geneId, restrictedUsedFeatures);
            restricted.allGeneFeatureNames.addAll(map(restrictedUsedFeatures, (_p) -> _p.getSecond()));
            restricted.eqClasses.addAll(map(restrictedUsedFeatures, (_p) -> _p.getFirst()));

            Vector<TranscriptPairInfo> toTest = TranscriptPairInfo.getTrPairInfos(restrictedUsedFeatures);
            if(toTest.size() > 0 ){
                restricted.genesForSplicingTesting.add(geneId);
                restricted.gene2pairtests.put(geneId, toTest);
            }
        }

        log.info("set maxmincount: %d skipped %d/%d genes and %d/%d features\n", minCountIntMaxCountCondition,  nskippedGenes, getGenes().size(),
                nskippedFeatures, feature2Idx.size());
        return restricted;
    }

    public int getMaximalTrNumForReduction() {
        return MAX_TR_NUM;
    }


    public int getMinCountIntMaxCountCondition() {
        return minCountIntMaxCountCondition;
    }

    public void setMinCountIntMaxCountCondition(int minCountIntMaxCountCondition) {
        this.minCountIntMaxCountCondition = minCountIntMaxCountCondition;

    }

    public HashMap<String, HashMap<String, Vector<Double>>> addGeneInfo(String geneId,
                                                                        HashMap<String, Vector<HashMap<Tuple, Double>>> condition2replicate2EQclassCounts
                    )
    {
        return addGeneInfo(geneId, condition2replicate2EQclassCounts, false);
    }

        /**
         *
         * @param geneId
         * @param condition2replicate2EQclassCounts
         * @return the per condition counts for the reduced representation features
         */
    public HashMap<String, HashMap<String, Vector<Double>>> addGeneInfo(String geneId,
                HashMap<String, Vector<HashMap<Tuple, Double>>> condition2replicate2EQclassCounts,
                boolean data_already_log2
                            ) {

        boolean verbose = false; //geneId.equalsIgnoreCase("ENSG00000203279");
        if(replicateSetInfos != null)
            throw new FRuntimeException("invalid operation! addGeneInfo after invoking compile!");

        HashMap<Tuple, Double> merged = new HashMap<>();

        for(Vector<HashMap<Tuple, Double>> v : condition2replicate2EQclassCounts.values()) {
            for(int i=0; i<v.size(); i++) {
                for(Map.Entry<Tuple, Double> e : v.get(i).entrySet()) {
                    MapBuilder.update(merged, e.getKey(), e.getValue());
                }
            }
        }
        ReducedTranscriptPresentation rtp = new ReducedTranscriptPresentation(merged, MAX_TR_NUM);

        int numFeatures = rtp.getNumRestrictedEQClasses();

        HashMap<String, HashMap<String, Vector<Double>>> reducedCounts = new HashMap<>();


        Vector<Pair<Tuple, String>> usedFeatures = new Vector<>();

        for(int i=0; i<numFeatures; i++) {

            Tuple eqClass = rtp.getEqClassByIdx(i);

            HashMap<String, Vector<Double>> reducedFeatureCounts = new HashMap<>();


            int maxMinCount = 0;
            for(String cond : condition2replicatenames.keySet()) {
                Vector<Double> data = map(condition2replicate2EQclassCounts.get(cond), (_m) -> _m.getOrDefault(eqClass, 0.0));
                reducedFeatureCounts.put(cond, data);

                maxMinCount = Math.max(maxMinCount, NumUtils.min(data).intValue());
            }

            if(maxMinCount < minCountIntMaxCountCondition)
                continue;

            String fn = geneId+"."+(usedFeatures.size()+1);
            usedFeatures.add(Pair.create(eqClass, fn));



            reducedCounts.put(fn, reducedFeatureCounts);

            for(String cond : condition2replicatenames.keySet()) {
                Vector<Double> data = reducedFeatureCounts.get(cond);
                Vector<Vector<Double>> toWriteIn = cond2replicate2featureLogData.get(cond);

                Vector<Double> logdata = (data_already_log2) ? data : map(data, (_d) -> (_d == 0.0) ? Double.NaN : NumUtils.logN(_d, 2.0));


                if(verbose) {
                    log.info("%s feature: %s data: %s log: %s", geneId, fn, data, logdata);
                }
                for(int repIdx = 0; repIdx < toWriteIn.size(); repIdx++) {
                    toWriteIn.get(repIdx).add(logdata.get(repIdx));
                }

            }
        }


        if(reducedCounts.size() == 0)
            return reducedCounts;

        gene2EQClassFeature.put(geneId, usedFeatures);
        gene2transcriptInfo.put(geneId, rtp);

        allGeneFeatureNames.addAll(map(usedFeatures, (_p) -> _p.getSecond()));
        eqClasses.addAll(map(usedFeatures, (_p) -> _p.getFirst()));

        Vector<TranscriptPairInfo> toTest = TranscriptPairInfo.getTrPairInfos(usedFeatures);

        //check for min counts!!!



        if(toTest.size() > 0 ){
            genesForSplicingTesting.add(geneId);
            gene2pairtests.put(geneId, toTest);
        }


        return reducedCounts;
    }

    static Vector<Pair<Tuple, String>> EMPTYV = new Vector<>();

    public int getNumFeatures(String gene) {
        return gene2EQClassFeature.getOrDefault(gene, EMPTYV).size();
    }
    public Collection<ReplicateSetInfo> getReplicateSetInfos() {
        return compile();
    }


    public Vector<String> getFeaturesToGene(String gene) {
        return map(gene2EQClassFeature.getOrDefault(gene, new Vector<>()), (_p) -> _p.getSecond());
    }


    public Collection<String> getGenes() {
        return gene2EQClassFeature.keySet();
    }

    public Vector<String> getGenesForSplicingTest() {
        return genesForSplicingTesting;
    }

    private Collection<ReplicateSetInfo> compile() {
        if(replicateSetInfos != null)
            return replicateSetInfos;

        log.info("compile rnaseq infos. got %d features, %d genes for splicing tests", gene2EQClassFeature.size(), genesForSplicingTesting.size());

        replicateSetInfos = new Vector<>();
        replicateSetNameLookup = new HashMap<>();

        HashMap<String, Vector<String>> gene2featureIds = buildMap(gene2EQClassFeature.entrySet(), (_e) -> _e.getKey(), (_e) -> map(_e.getValue(), (_p) -> _p.getSecond()));
        for(String cond : condition2replicatenames.keySet()) {
            Vector<String> repNames = condition2replicatenames.get(cond);
            ReplicateSetInfo rsi = new ReplicateSetInfo(cond, repNames, allGeneFeatureNames);
            Vector<Vector<Double>> repLog2Data = cond2replicate2featureLogData.get(cond);

            for(int i=0; i<repLog2Data.size(); i++) {

                rsi.setLog2Data(i, repLog2Data.get(i));
            }


            rsi.setCombinedFeatures(gene2featureIds);

            replicateSetNameLookup.put(rsi.getReplicateSetName(), rsi);
            replicateSetInfos.add(rsi);
        }



        feature2Idx = buildIndexMap(allGeneFeatureNames);
        return replicateSetInfos;
    }

    public Pair<Tuple, Vector<Vector<Double>>> getEQClassWithCounts(String feature) {
        Tuple eqClass = getEQClass(feature);
        return Pair.create(eqClass, map(replicateSetInfos, (_rsi) -> _rsi.getReplicateData(feature)));
    }

    public Tuple getEQClass(String feature) {
        return eqClasses.get(feature2Idx.get(feature));
    }

    boolean isValidForTest(Vector<Double> data) {
        double avg = 0.0;
        for(double d : data) {
            if(Double.isNaN(d) || Double.isInfinite(d))
                return false;

            avg += d;
        }

        avg /= data.size();
        return (avg >= MIN_AVG_COUNT_ON_STRONGER_SIDE);
    }

    boolean testEQClassForTest(String feature, String cond1, String cond2) {

        Vector<Double> logdata1 = replicateSetNameLookup.get(cond1).getReplicateData(feature);
        Vector<Double> logdata2 = replicateSetNameLookup.get(cond2).getReplicateData(feature);

        return isValidForTest(logdata1) || isValidForTest(logdata2);

    }

    public Vector<Tuple3<String, Vector<String>, Vector<String>>> getTrPairTestFeatureCombiTests(String gene, String tr1, String tr2, String cond1, String cond2, StringBuffer err) {
        Vector<TranscriptPairInfo> pairTests = gene2pairtests.get(gene);
        if(pairTests == null || pairTests.size() == 0) {
            err.append("no pair tests found for gene");
            return null;
        }

        ReducedTranscriptPresentation rtp = gene2transcriptInfo.get(gene);

        if(rtp == null) {
            err.append("missing reduced transcript presentation");
            return null;
        }

        String leadTr1 = rtp.getLeadTr(tr1);
        String leadTr2 = rtp.getLeadTr(tr2);

        if(leadTr1 == null) {
            err.append("lead tr for " + tr1 + " is null");
            return null;
        }

        if(leadTr2 == null) {
            err.append("lead tr for " + tr2 + " is null");
            return null;
        }

        if(leadTr1.equals(leadTr2)) {
            err.append("trs not testable both have the same lead:"+ leadTr1);
            return null;
        }
        UPair<String> tocheck = UPair.createUSorted(leadTr1, leadTr2);
        TranscriptPairInfo tpi = filterOne(pairTests, (_p) -> tocheck.equals(UPair.createUSorted(_p.tr1, _p.tr2)));
        if(tpi == null) {
            err.append("no transcript pair test found");
            return null;
        }
        Vector<Pair<Tuple, String>> features = gene2EQClassFeature.get(gene);

        Vector<Tuple3<String, Vector<String>, Vector<String>>> featureCombisToTest = new Vector<>();

        HashMap<Tuple, String> eqclass2name = buildMap(features, (_p) -> _p.getFirst(), (_p) -> _p.getSecond());

        for(Tuple3<String, Vector<Tuple>, Vector<Tuple>> totest : tpi.testableCombis) {
            Vector<String> features1 = map_and_filter(totest.get1(), (_t) -> eqclass2name.get(_t), (_fn) -> testEQClassForTest(_fn, cond1, cond2));
            Vector<String> features2 = map_and_filter(totest.get2(), (_t) -> eqclass2name.get(_t), (_fn) -> testEQClassForTest(_fn, cond1, cond2));

            if(features1.size() == 0 || features2.size() == 0)
                continue;

            featureCombisToTest.add(Tuple3.create(totest.get0(), features1, features2));
        }
        return featureCombisToTest;
    }

    /** for every tuple first param: testname, second and third collection of featureIds to test against each other */
    public Vector<Tuple3<String, Vector<String>, Vector<String>>> getSplicingTestFeatureCombis(String gene, String cond1, String cond2) {

        Vector<Tuple3<String, Vector<String>, Vector<String>>> featureCombisToTest = new Vector<>();

        ReducedTranscriptPresentation rtp = gene2transcriptInfo.get(gene);
        if(rtp == null)
            return featureCombisToTest;


        Vector<TranscriptPairInfo> pairTests = gene2pairtests.get(gene);
        if(pairTests == null || pairTests.size() == 0)
            return featureCombisToTest;


        Vector<Pair<Tuple, String>> features = gene2EQClassFeature.get(gene);
        HashMap<Tuple, String> eqclass2name = buildMap(features, (_p) -> _p.getFirst(), (_p) -> _p.getSecond());


        for(TranscriptPairInfo tpi : pairTests) {
            for(Tuple3<String, Vector<Tuple>, Vector<Tuple>> totest : tpi.testableCombis) {
                    Vector<String> features1 = map_and_filter(totest.get1(), (_t) -> eqclass2name.get(_t), (_fn) -> testEQClassForTest(_fn, cond1, cond2));
                    Vector<String> features2 = map_and_filter(totest.get2(), (_t) -> eqclass2name.get(_t), (_fn) -> testEQClassForTest(_fn, cond1, cond2));

                    if(features1.size() == 0 || features2.size() == 0)
                        continue;

                    featureCombisToTest.add(Tuple3.create(totest.get0(), features1, features2));
            }
        }

        return featureCombisToTest;
    }


    public int getMinCountLevelForCoveragePercentage(String cond, int percentage) {


        ReplicateSetInfo conditionReplicates = filterOne(getReplicateSetInfos(), (_rsi) -> _rsi.getReplicateSetName().equals(cond));
        if(conditionReplicates == null)
            throw new FRuntimeException("unknown condition: %s", cond);

        HashMap<Integer, Integer> maxfeaturecounts = new HashMap<>();

        for(String g : getGenes()) {

            for(String feature: getFeaturesToGene(g)) {
                Vector<Integer> fdata = map(conditionReplicates.getReplicateData(feature), (_v) -> (int)((Double.isNaN(_v) ? 0.0 : Math.pow(2.0, _v))));
                MapBuilder.update(maxfeaturecounts, NumUtils.min(fdata));

            }
        }
        Vector<Vector<UPair<Double>>> countlevel2percentData = new Vector<>();

        Vector<Integer> keys = toVector(maxfeaturecounts.keySet());
        Collections.sort(keys);
        Collections.reverse(keys);
        double sum = 0.0 + NumUtils.sum(map(keys, (_k) -> _k * maxfeaturecounts.get(_k)));
        double covered = 0.0;
        Vector<UPair<Double>> replicateCovered = new Vector<>();
        Vector<Integer> coveragePercents = new Vector<>();
        Vector<Integer> counts = new Vector<>();

        for(int k : keys) {
            covered += k * maxfeaturecounts.get(k);
            coveragePercents.add((int)(covered * 100 / sum));
            counts.add(k);
            //System.out.printf("%s: %d count: %d covered: %.0f sum: %.0f : %.2f%%\n", cond, idx, k, covered, sum, covered * 100 / sum);
            replicateCovered.add(UPair.createU(k + 0.0, covered * 100.0 / sum));
        }
        int hitIdx = Collections.binarySearch(coveragePercents, percentage);
        if(hitIdx < 0) {
            hitIdx = - hitIdx - 1;
        }
        int minC = counts.get(hitIdx);
        double sumTest = 0.0 + NumUtils.sum(map(filter(keys, (_k) -> _k >= minC), (_k) -> _k * maxfeaturecounts.get(_k)));
        log.info("%d%% coverage at count: %d check: %.2f%%\n", percentage, minC, 100.0 * sumTest / sum);

        return minC;

    }


}
