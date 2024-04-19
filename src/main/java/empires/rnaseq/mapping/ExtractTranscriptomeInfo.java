package nlEmpiRe.rnaseq.mapping;

import lmu.utils.*;
import nlEmpiRe.rnaseq.*;
import java.util.*;
import org.apache.logging.log4j.Logger;

public class ExtractTranscriptomeInfo {

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


    public static void main(String[] args) {
        //gtf, genome, extract all (or only certain type of ) transcripts (from an optionally limited set of chromosomes */

        SimpleOptionParser cmd = new SimpleOptionParser("gtf", "genome", "genomeidx",  "minlength", "o");
        cmd.setFile("gtf", "genome", "genomeidx");
        cmd.setInt("minlength");
        cmd.setDefault("minlength", "100");

        if(!OptionParser.parseParams(args, true, true, true, cmd))
            return;

        int MINLENGTH = cmd.getInt("minlength");
        Logger log = LogConfig.getLogger();
        IsoformRegionGetter annot = new GFFBasedIsoformRegionGetter(cmd.getFile("gtf"), null, null);

        GenomeSequenceExtractor gex = new GenomeSequenceExtractor(cmd.getFile("genome"), cmd.getFile("genomeidx"));

        char SPACER = ' ';
        StringBuilder sb = new StringBuilder();
        List<TranscriptInfo> transcriptInfoList = new ArrayList<>();
        List<String> checkseqs = new ArrayList<>();

        HashSet<String> skipped_chrs = new HashSet<>();
        HashMap<String, Integer> biotypes = new HashMap<>();
        for(MultiIsoformRegion mir : annot.getRegionsIteratble()) {

            if(!gex.gotChr(mir.chr)) {
                if(!skipped_chrs.contains(mir.chr)) {
                    log.info("chromosome: %s is not in provided genome. skipping it.", mir.chr);
                    skipped_chrs.add(mir.chr);
                }
                continue;
            }

            //if(!isStandardChr(mir.chr))
            //    continue;

            for(String trId : mir.isoforms.keySet()) {
                if(sb.length() > 0) {
                    sb.append(SPACER);
                }

                MapBuilder.update(biotypes, mir.iso2biotype.getOrDefault(trId, mir.biotype));

                TranscriptInfo transcriptInfo = new TranscriptInfo(mir, trId, sb.length(), transcriptInfoList.size());
                if(transcriptInfo.length < MINLENGTH)
                    continue;

                String tseq = GenomicUtils.getSplicedSeq(mir.chr, mir.strand, transcriptInfo.rv, gex, true);

                transcriptInfoList.add(transcriptInfo);

                checkseqs.add(tseq);
                sb.append(tseq);
                if(transcriptInfoList.size() % 1000 == 0) {

                    log.info("after %d transcripts length: %d", transcriptInfoList.size(), sb.length());
                }

            }
        }
        log.info("ready. extracted %d transcripts length: %d", transcriptInfoList.size(), sb.length());
        //put it into  asuffix array
        TranscriptInfoBasedGenomicMapper transcriptInfoBasedGenomicMapper = new TranscriptInfoBasedGenomicMapper(transcriptInfoList, sb);
        transcriptInfoBasedGenomicMapper.serialize(cmd.getFile("o"));
    }
}
