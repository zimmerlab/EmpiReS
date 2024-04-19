package nlEmpiRe.rnaseq.simulation;

import lmu.utils.FRuntimeException;
import lmu.utils.UPair;
import org.apache.commons.math3.distribution.NormalDistribution;

public class PositionBias {

    ScalableTranscriptInfo biasModel;
    int transcriptLength;
    NormalDistribution fragmentLengthDistrib;
    int readLength;

    public PositionBias(NormalDistribution fragmentLengthDistrib, int transcriptLength, ScalableTranscriptInfo biasModel, int readLength) {
        this.biasModel = biasModel;
        this.transcriptLength = transcriptLength;
        this.fragmentLengthDistrib = fragmentLengthDistrib;
        this.readLength = readLength;
    }

    public UPair<int[]> getReads() {

        int frLength = Math.min(transcriptLength, Math.max(readLength, (int)fragmentLengthDistrib.sample()));
        int startpos = (int)(Math.random() * (transcriptLength - frLength));

        if(biasModel != null) {
            int reqLength = frLength;
            UPair<Integer> startend = biasModel.getFragment(frLength);
            startpos = startend.getFirst();
            frLength = startend.getSecond() - startpos + 1;
            if(frLength < readLength)
                throw new FRuntimeException("invalid frlength  req was: %d readlength: %d got: %d startend: %s tlength: %d", reqLength, readLength, frLength, startend, transcriptLength);
        }

        int rl = Math.min(readLength, frLength);

        return UPair.createU(
                new int[]{startpos, startpos + readLength},
                new int[]{startpos + frLength - readLength, startpos + frLength}
        );


    }


}
