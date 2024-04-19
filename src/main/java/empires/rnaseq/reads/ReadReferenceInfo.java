package nlEmpiRe.rnaseq.reads;

import lmu.utils.FileUtils;
import lmu.utils.Region1D;
import lmu.utils.RegionVector;

import java.io.File;
import java.util.Iterator;

public class ReadReferenceInfo {
    public int readid;
    public String chr;
    public String gene;
    public String transcript;
    public RegionVector fw;
    public RegionVector rw;

    public Region1D t_fw;
    public Region1D t_rw;

    public ReadReferenceInfo(FileUtils.HeaderedReader hr){
        this.readid = hr.getInt("id");
        this.chr = hr.getString("chr");
        this.gene = hr.getString("gene");
        this.transcript = hr.getString("transcript");
        fw = RegionVector.parseSimpleRepresentation(hr.getString("fw"));
        rw = RegionVector.parseSimpleRepresentation(hr.getString("rw"));

        t_fw = RegionVector.parseSimpleRepresentation(hr.getString("tfw")).getRegion(0);
        t_rw = RegionVector.parseSimpleRepresentation(hr.getString("trw")).getRegion(0);
    }

    public static Iterator<ReadReferenceInfo> getIterator(File f) {
        FileUtils.HeaderedReader hr = FileUtils.getHeaderedReader("\t")
                        .add("readid", "id").add("chr").add("gene").add("transcript")
                        .add("fw_regvec", "fw").add("rw_regvec", "rw")
                .add("t_fw_regvec" , "tfw").add("t_rw_regvec", "trw");


        return FileUtils.getConvertedHrIterator(FileUtils.getInputStream(f), hr,  (_hr) -> new ReadReferenceInfo(_hr), false);
    }

    public String toString() {
        return String.format("%s:%s %s-%s", gene, transcript, t_fw, t_rw);
    }
}
