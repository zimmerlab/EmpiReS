package empires.rnaseq.mapping;

import lmu.utils.Region1D;

public class TranscriptHit {
    empires.rnaseq.mapping.TranscriptInfo transcriptInfo;
    Region1D fragment;
    int numMismatches;

    public TranscriptHit(TranscriptInfo tinfo, int start, int end, int numMismatches) {
        this.transcriptInfo = tinfo;
        this.fragment = new Region1D(Math.min(start, end), Math.max(start,end)) ;
        this.numMismatches = numMismatches;
    }

    public static class TrHitShort {
        String gene;
        String trid;
        Region1D fragment;
        int numMismatches;

        public TrHitShort(String info) {
            String[] sp = info.split(":");
            gene = sp[0];
            trid = sp[1];
            fragment = new Region1D(Integer.parseInt(sp[2]), Integer.parseInt(sp[3]));
            numMismatches = Integer.parseInt(sp[4]);
        }


    }

    public String getShortInfo() {
        return getShortInfo(transcriptInfo.gene);
    }

    public String getRefShortInfo() {
        return getShortInfo("REF");
    }

    public String getShortInfo(String geneid) {
        return String.format("%s:%s:%d:%d:%d", geneid, transcriptInfo.transcriptId, this.fragment.getX1(), this.fragment.getX2(), numMismatches);
    }


    public String toString() {
        return String.format("hit to %s:%s %s mm: %d", transcriptInfo.gene, transcriptInfo.transcriptId, fragment, numMismatches);
    }
}
