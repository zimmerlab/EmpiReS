package nlEmpiRe;

import lmu.utils.EnumGetter;
import nlEmpiRe.rnaseq.SplicingTest;

public enum DoubleDiffVariant {
    QUICK_AND_DIRTY("quickanddirty"),
    MISSINGVAL("missingval"),
    ALLPAIRS("allpairs")
    ;

    String name;
    DoubleDiffVariant(String n) {
        this.name = n;
        EnumGetter.add(n, this);
    }

    public String getName() {
        return name;
    }

    public static DoubleDiffVariant get(String n) {
        return EnumGetter.get(n, true);
    }
}
