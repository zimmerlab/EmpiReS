package nlEmpiRe;

public class PairScore {
    public ShiftedGroup sg1;
    public ShiftedGroup sg2;
    public PairwiseMedianImpliedFCPeakErrorEstimation err;

    private PairScore() {

    }
    public PairScore(ShiftedGroup sg1, ShiftedGroup sg2) {
        this.sg1 = sg1; this.sg2 = sg2;
        if(sg1 == sg2)
            return;
        err = new PairwiseMedianImpliedFCPeakErrorEstimation(sg1, sg2);
    }

    public PairScore flip() {
        PairScore flipped = new PairScore();
        flipped.sg1 = sg2;
        flipped.sg2 = sg1;
        flipped.err = err;
        return flipped;

    }
}
