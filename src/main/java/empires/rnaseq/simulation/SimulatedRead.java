package nlEmpiRe.rnaseq.simulation;

import lmu.utils.FRuntimeException;
import lmu.utils.Region1D;
import lmu.utils.RegionVector;
import lmu.utils.UPair;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.function.Consumer;

import static lmu.utils.ObjectGetter.toVector;

public class SimulatedRead {
    public String gene;
    public String transcriptId;
    public HashSet<String> mapsToTranscripts = new HashSet<>();

    public RegionVector fw;
    public RegionVector rw;

    public RegionVector tr_fw;
    public RegionVector tr_rw;

    public RegionVector merged;



    public SimulatedRead(String gene, boolean strand, String transcriptId, RegionVector trVec, int readLength, HashMap<String, RegionVector> isoforms,
                         boolean paired,
                         PositionBias bias) {

            this.gene = gene;
            this.transcriptId = transcriptId;
            UPair<int[]> fragementCoords = bias.getReads();
            int[] fwRead = fragementCoords.getFirst();
            int[] rwRead = fragementCoords.getSecond();


            tr_fw = new RegionVector(toVector(new Region1D(fwRead[0], fwRead[1])), true);
            fw = trVec.cutLocal(fwRead[0], fwRead[1], !strand);
            merged = fw;
            rw = null;

            if(paired) {
                tr_rw = new RegionVector(toVector(new Region1D(rwRead[0], rwRead[1])), true);
                rw = trVec.cutLocal(rwRead[0], rwRead[1], !strand);
            }
            if(fw == null || (paired && rw == null))
                throw new FRuntimeException("invalid region vectors for %s:%s frLength: %d/%d startpos: %d fw: %s rw: %s tr: %s %s", gene, transcriptId,
                        rwRead[1]  - fwRead[0], trVec.getCoveredLength(), fwRead[0],  fw, rw, tr_fw, tr_rw);

            merged = RegionVector.merge(fw, rw);

            for(Map.Entry<String, RegionVector> e : isoforms.entrySet())
            {
                RegionVector tr = e.getValue();
                boolean canBeProducedByTr = (merged.getNumRegions() == 1) ? tr.isSubVector3(merged, false) : tr.isSubVector3(fw, false)
                        && (rw == null || tr.isSubVector3(rw, false));

                if(!canBeProducedByTr)
                    continue;

                mapsToTranscripts.add(e.getKey());

            }

            if(!mapsToTranscripts.contains(transcriptId)) {
                throw new FRuntimeException("invalid simulation for %s:%s:%s\nfw: %s rw: %s\n", gene, transcriptId, trVec, fw, rw);
            }


    }


}
