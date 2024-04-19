package nlEmpiRe.test.rnaseq;



import lmu.utils.Region1D;
import lmu.utils.RegionVector;
import net.sf.samtools.*;
import nlEmpiRe.rnaseq.GenomeSequenceExtractor;

import java.io.File;
import java.util.Map;

public class BAMGenerator {

    final static String QUALSTRING = "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII";
    SAMFileHeader samFileHeader;
    BAMFileWriter bamWriter;


    public BAMGenerator(File referenceBamFile, GenomeSequenceExtractor gex, File outBamFile) {

        if(referenceBamFile != null) {
            SAMFileReader sam_reader =  new SAMFileReader(referenceBamFile);
            sam_reader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
            samFileHeader = sam_reader.getFileHeader().clone();
        } else {
            samFileHeader = new SAMFileHeader();
            for(Map.Entry<String, Long> e : gex.getChr2Length().entrySet()) {
                samFileHeader.addSequence(new SAMSequenceRecord(e.getKey(), e.getValue().intValue()));
            }

        }
        bamWriter  = new BAMFileWriter(outBamFile);
        bamWriter.setHeader(samFileHeader);
    }

    public void close() {
        bamWriter.close();
    }

    static Cigar toCigar(RegionVector rv) {
        Cigar cigar = new Cigar();
        Region1D last = null;
        for(Region1D r : rv.getRegions()) {

            if(last != null) {
                cigar.add(new CigarElement(r.getX1() - last.getX2(), CigarOperator.D));
            }
            last = r;
            cigar.add(new CigarElement(r.getLength(), CigarOperator.M));
        }
        return cigar;
    }


    public void writeRecords(String readId, String fw_readseq, String rw_readseq, String qualstring, String chr, boolean strand, RegionVector fw, RegionVector rw) {

        for(int i=0; i<2; i++) {
            RegionVector src = (i == 0) ? fw : rw;
            if(src == null)
                continue;

            SAMRecord sr = new SAMRecord(samFileHeader);
            sr.setReadName(readId);
            sr.setReferenceName(chr);
            sr.setReadString((i == 0) ? fw_readseq : rw_readseq);
            sr.setBaseQualityString(qualstring);
            if(rw != null) {
                sr.setProperPairFlag(true);
                sr.setNotPrimaryAlignmentFlag(false);
                sr.setMateUnmappedFlag(false);

                sr.setReadPairedFlag(true);

                sr.setFirstOfPairFlag(i == 0);
                sr.setSecondOfPairFlag(i != 0);

                sr.setMateNegativeStrandFlag((i == 0) ? strand : !strand);
            }
            sr.setReadNegativeStrandFlag((i == 0) ? !strand : strand);
            sr.setAlignmentStart(src.getX1());
            //sr.setAlignmentEnd(src.getX2());

            if(rw != null) {
                sr.setMateReferenceName(chr);
                RegionVector other = (i == 0) ? rw : fw;
                sr.setMateAlignmentStart(other.getX1());

            }

            sr.setAttribute("xP", String.format("%s:%s", chr, src.getSimpleRepresentation()));
            sr.setCigar(toCigar(src));

            bamWriter.addAlignment(sr);
        }

    }



}
