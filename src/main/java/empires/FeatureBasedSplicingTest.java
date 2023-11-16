package empires;

import lmu.utils.*;
import lmu.utils.fdr.BenjaminiHochberg;
import empires.test.rnaseq.BenchmarkGene;

import java.io.File;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;
import java.util.Vector;
import java.util.function.Function;

import static lmu.utils.ObjectGetter.*;

public class FeatureBasedSplicingTest {

    Logger log = LogConfig.getLogger();
    static Vector<String> GPEPTIDOME_HEADERS = toVector("gene", "id", "peptide");
    static Vector<String> GENERAL_HEADERS = toVector("gene", "isoform", "feature");

    empires.NormalizedReplicateSet nrs1;
    empires.NormalizedReplicateSet nrs2;
    Set<String> sharedfeatures;
    Vector<empires.DoubleDiffResult> splicingResults;
    int num_min_features_to_test = 1;
    empires.DoubleDiffVariant doubleDiffVariant = empires.DoubleDiffVariant.QUICK_AND_DIRTY;
    SplicingTestFilter splicingTestFilter = null;

    public FeatureBasedSplicingTest(empires.NormalizedReplicateSet nrs1, empires.NormalizedReplicateSet nrs2, File mappinginfo) {
        this(nrs1, nrs2, mappinginfo, empires.DoubleDiffVariant.QUICK_AND_DIRTY, 1, null);
    }

    public FeatureBasedSplicingTest(empires.NormalizedReplicateSet nrs1, empires.NormalizedReplicateSet nrs2, File mappinginfo, DoubleDiffVariant doubleDiffVariant, int num_min_features_to_test, SplicingTestFilter splicingTestFilter) {
        this.doubleDiffVariant = doubleDiffVariant;
        this.num_min_features_to_test = num_min_features_to_test;
        this.splicingTestFilter = splicingTestFilter;
        HashSet<String> headers = toSet(FileUtils.getHeaders(mappinginfo));
        if( 0 == filteredSize(GPEPTIDOME_HEADERS, (_h) -> !headers.contains(_h))) {
            init(nrs1, nrs2, mappinginfo, GPEPTIDOME_HEADERS.get(0), GPEPTIDOME_HEADERS.get(1), GPEPTIDOME_HEADERS.get(2));
            return;
        }
        if( 0 == filteredSize(GENERAL_HEADERS, (_h) -> !headers.contains(_h))) {
            init(nrs1, nrs2, mappinginfo, GENERAL_HEADERS.get(0), GENERAL_HEADERS.get(1), GENERAL_HEADERS.get(2));
            return;
        }

        throw new FRuntimeException("no automatic recognition for genome, isoform, feature mapping from the headers: %s! current known sets are: %s, %s\nmapping file: %s\n",
                headers, GPEPTIDOME_HEADERS, GENERAL_HEADERS, mappinginfo.getAbsolutePath());
    }


    public FeatureBasedSplicingTest(empires.NormalizedReplicateSet nrs1, empires.NormalizedReplicateSet nrs2, File mappinginfo, String gene_field, String isoform_field, String feature_field) {
        init(nrs1, nrs2, mappinginfo, gene_field, isoform_field, feature_field);
    }

    void init(empires.NormalizedReplicateSet nrs1, empires.NormalizedReplicateSet nrs2, File mappinginfo,
              String gene_field, String isoform_field, String feature_field) {

        sharedfeatures = SetInfo.intersect(toSet(nrs1.inData.getFeatureNames()), toSet(nrs2.inData.getFeatureNames()));

        HashMap<String, Vector<String>> strippedFeature2Feature = new HashMap<>();

        for(String f : sharedfeatures) {
            String[] s = f.split("_");
            MapBuilder.updateV(strippedFeature2Feature, s[0], f);

        }


        HashMap<String, HashSet<String>> gene2isoforms = new HashMap<>();
        HashMap<String, HashSet<String>> gene2features = new HashMap<>();
        HashMap<String, Vector<String>> feature2isoforms = new HashMap<>();

        HashSet<String> visitedFeatures = new HashSet<>();
        HashMap<String, String> isoform2gene = new HashMap<>();

        HashSet<String> nonUniqFeatures = new HashSet<>();

        HashSet<String> nonUniqFeaturesReal = new HashSet<>();

        apply(FileUtils.getHrIterator(mappinginfo, FileUtils.getHeaderedReader("\t")
                        .add(gene_field, "gene")
                        .add(isoform_field, "isoform")
                        .add(feature_field, "feature")),
                (_hr) -> {
                    String feature = _hr.getString("feature");
                    if(!strippedFeature2Feature.containsKey(feature))
                        return;

                    if(visitedFeatures.contains(feature)) {
                        nonUniqFeatures.add(feature);
                        nonUniqFeaturesReal.addAll(strippedFeature2Feature.get(feature));
                        return;
                    }
                    visitedFeatures.add(feature);
                    Vector<String> realfeatures = strippedFeature2Feature.get(feature);

                    String gene = _hr.getString("gene");
                    Vector<String> isoforms = toSortedVector(toVector(_hr.getString("isoform").split(",")), true);
                    for(String isoform : isoforms) {
                        String oldgene = isoform2gene.getOrDefault(isoform, gene);
                        if(!oldgene.equals(gene))
                            throw new FRuntimeException("invalid mapping! isoform: %s is annotated to multiple genes to: %s and %s!", isoform2gene, oldgene, gene);

                        MapBuilder.update(gene2isoforms, gene, isoform);
                        isoform2gene.put(isoform, gene);
                    }
                    apply(realfeatures, (_rf) -> MapBuilder.update(gene2features, gene, _rf));
                    apply(realfeatures, (_rf) -> feature2isoforms.put(_rf, isoforms));

                });

        if(nonUniqFeatures.size() > 0) {
            log.warn("found %d non uniq feature definitions in %s these will be ignored. examples: %s", nonUniqFeatures.size(), mappinginfo.getAbsolutePath(),
                    VectorUtils.slice(toVector(nonUniqFeatures), 0, 10));
        }
        apply(nonUniqFeatures, (_n) -> feature2isoforms.remove(_n));

        doTests(nrs1, nrs2, gene2features, feature2isoforms);

    }

    public FeatureBasedSplicingTest(empires.NormalizedReplicateSet nrs1, empires.NormalizedReplicateSet nrs2,
                                    HashMap<String, HashSet<String>> gene2features,
                                    HashMap<String, Vector<String>> feature2isoforms


    ) {
        doTests(nrs1, nrs2, gene2features, feature2isoforms);
    }

    public interface SplicingTestFilter {
        public boolean filter(String gene, Vector<String> isoform_ids_eqclass1, Vector<String> isoform_ids_eqclass2, Vector<String> peptide_ids1, Vector<String> peptide_ids2);
    }
    void doTests(empires.NormalizedReplicateSet nrs1, NormalizedReplicateSet nrs2,
                 HashMap<String, HashSet<String>> gene2features,
                 HashMap<String, Vector<String>> feature2isoforms
                                    ) {


        this.nrs1 = nrs1;
        this.nrs2 = nrs2;

        sharedfeatures = SetInfo.intersect(toSet(nrs1.inData.getFeatureNames()), toSet(nrs2.inData.getFeatureNames()));


        DoubleDiffManager diffManager = new EmpiRe().getDoubleDiffManager(nrs1, nrs2);

        Vector<Integer> eqlassSizes = new Vector<>();

        splicingResults = new Vector<>();
        for(String gene : gene2features.keySet()) {

            Vector<String> usable_features  = filter(gene2features.get(gene), (_p) -> feature2isoforms.containsKey(_p));

            if(usable_features.size() == 0)
                continue;

            HashMap<Tuple, Vector<String>> eq2features = new HashMap<>();
            for(String feature : usable_features) {
                MapBuilder.updateV(eq2features, Tuple.tupleFromCollection(toSortedVector(feature2isoforms.get(feature), true)), feature);
            }

            Vector<Tuple> tkeys = filter(eq2features.keySet(), (_t) -> eq2features.get(_t).size() >= num_min_features_to_test);

            eqlassSizes.add(tkeys.size());
            if(tkeys.size() < 2)
                continue;



            Function<Tuple, Vector<String>> eq2isoformids = (_t) -> mapIndex(_t.cardinality(), (_i) -> _t.getAsString(_i));

            empires.DoubleDiffResult bestHit = null;
            for(UPair<Tuple> t : getPairs(tkeys, true)) {

                if(splicingTestFilter != null && !splicingTestFilter.filter(gene, eq2isoformids.apply(t.getFirst()), eq2isoformids.apply(t.getSecond()), eq2features.get(t.getFirst()), eq2features.get(t.getSecond())))
                    continue;

                empires.DoubleDiffResult ddr = diffManager.getDoubleDiffResult(gene+"."+t, eq2features.get(t.getFirst()), eq2features.get(t.getSecond()), doubleDiffVariant);

                if(bestHit == null || bestHit.pval > ddr.pval) {
                    bestHit = ddr;
                }
            }
            if(bestHit == null)
                continue;
            splicingResults.add(bestHit);
        }


        BenjaminiHochberg.adjust_pvalues(splicingResults, (_s) -> _s.pval, (_p) -> _p.getFirst().fdr = _p.getSecond());


    }

    public Vector<DoubleDiffResult> getSplicingResults() {
        return splicingResults;
    }

    public DataTable getSplicingResultTable() {
        return BenchmarkGene.getTable(splicingResults, null, null, null);
    }

    public void showSplicingResultTable() {
        BenchmarkGene.showTable(splicingResults, null, null, null);
    }
}
