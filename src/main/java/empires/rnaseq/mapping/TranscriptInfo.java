package nlEmpiRe.rnaseq.mapping;

import lmu.utils.*;
import nlEmpiRe.rnaseq.MultiIsoformRegion;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.util.*;
import java.util.function.BiConsumer;

public class TranscriptInfo {
    String chr;
    boolean strand;
    String transcriptId;
    RegionVector rv;
    String gene;
    int strStart;
    int length;
    int idx;


    public TranscriptInfo(MultiIsoformRegion gene, String transcriptId, int start, int idx) {
        this.chr = gene.chr;
        this.strand = gene.strand;
        this.transcriptId = transcriptId;
        this.rv = gene.isoforms.get(transcriptId);
        this.gene = gene.id;
        this.strStart = start;
        this.length = rv.getCoveredLength();
        this.idx = idx;
    }


    private TranscriptInfo() {

    }

    public void write(DataOutputStream dos) {
        try {
            FileUtils.writeString(chr, dos);
            dos.writeBoolean(strand);
            FileUtils.writeString(transcriptId, dos);
            String rvinfo = rv.getSimpleRepresentation();
            FileUtils.writeString(rvinfo, dos);
            FileUtils.writeString(gene, dos);
            dos.writeInt(strStart);
            dos.writeInt(length);
            dos.writeInt(idx);

        }catch (IOException ie) {
            throw new FRuntimeException("i/o error: %s while writing transcript info to data output stream", ie, ie.getMessage());
        }

    }

    public static TranscriptInfo read(DataInputStream dis) {
        try {
            TranscriptInfo ti = new TranscriptInfo();
            ti.chr = FileUtils.readString(dis);
            ti.strand = dis.readBoolean();
            ti.transcriptId = FileUtils.readString(dis);
            ti.rv = RegionVector.parseSimpleRepresentation(FileUtils.readString(dis));
            ti.gene = FileUtils.readString(dis);
            ti.strStart = dis.readInt();
            ti.length = dis.readInt();
            ti.idx = dis.readInt();
            return ti;
        }catch (IOException ie) {
            throw new FRuntimeException("error at reading from data input stream");
        }
    }
    public static Pair<TranscriptInfo, Integer> read(byte[] data, int offset) {
        TranscriptInfo ti = new TranscriptInfo();
        Pair<String, Integer> sinfo = FileUtils.readString(data, offset);
        ti.chr = sinfo.getFirst();
        offset = sinfo.getSecond();
        ti.strand = ByteBuffer.wrap(data, offset, 1).get() > 0;
        offset++;
        sinfo = FileUtils.readString(data, offset);
        ti.transcriptId = sinfo.getFirst();
        sinfo = FileUtils.readString(data, sinfo.getSecond());
        ti.rv = RegionVector.parseSimpleRepresentation(sinfo.getFirst());
        sinfo = FileUtils.readString(data, sinfo.getSecond());
        ti.gene = sinfo.getFirst();
        offset = sinfo.getSecond();
        ti.strStart = ByteBuffer.wrap(data, offset, 4).getInt();
        offset += 4;
        ti.length = ByteBuffer.wrap(data, offset, 4).getInt();
        offset += 4;
        ti.idx = ByteBuffer.wrap(data, offset, 4).getInt();
        offset += 4;

        return Pair.create(ti, offset);
    }


}
