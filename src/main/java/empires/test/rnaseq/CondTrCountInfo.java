package nlEmpiRe.test.rnaseq;

import lmu.utils.NumUtils;
import lmu.utils.UPair;

import java.util.HashMap;
import java.util.Vector;

import static lmu.utils.ObjectGetter.toVector;

public class CondTrCountInfo {
    String minor;
    String major;

    UPair<Integer> minorMinMax = UPair.createU(Integer.MAX_VALUE, 0);
    UPair<Integer> majorMinMax = UPair.createU(Integer.MAX_VALUE, 0);

    CondTrCountInfo(HashMap<String, Vector<Integer>> data) {
        Vector<String> trs = toVector(data.keySet());
        double mean1 = NumUtils.mean(data.get(trs.get(0)));
        double mean2 = NumUtils.mean(data.get(trs.get(1)));

        minor =  (mean1 <= mean2) ? trs.get(0) : trs.get(1);
        major =  (minor.equals(trs.get(0))) ? trs.get(1) : trs.get(0);

        minorMinMax.set(NumUtils.min(data.get(minor)), NumUtils.max(data.get(minor)));
        majorMinMax.set(NumUtils.min(data.get(major)), NumUtils.max(data.get(major)));

    }
}
