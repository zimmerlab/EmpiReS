package empires.input;

import lmu.utils.SimpleOptionParser;
import empires.DiffExpResult;
import empires.ErrorEstimationDistribution;
import empires.SingleFeatureDiffExp;

public class GeneralOptions extends  SimpleOptionParser {

    public static String COMBINEPERPAIRDISTRIBS = "combineperreppairdistribs";
    public static String NOOUTLIERCORRECTION = "nooutliercorrection";
    public static String PERCENT_ROBUST_FEATURES_USED = "percentrobustfeatures";
    public static String FCBINWIDTH = "fcbinwidth";


    public GeneralOptions() {
        super(COMBINEPERPAIRDISTRIBS, NOOUTLIERCORRECTION, PERCENT_ROBUST_FEATURES_USED, FCBINWIDTH);
        setName("EmpiRe-general options");
        setSwitches(COMBINEPERPAIRDISTRIBS, NOOUTLIERCORRECTION);
        setDouble(PERCENT_ROBUST_FEATURES_USED, FCBINWIDTH);
        setDefault(PERCENT_ROBUST_FEATURES_USED, "0.9");
        setDefault(FCBINWIDTH, "0.01");

    }


    public void apply() {
        DiffExpResult.PERCENT_ROBUST_FEATURE_USED = getDouble(PERCENT_ROBUST_FEATURES_USED);
        if (isSet(COMBINEPERPAIRDISTRIBS)) {
            SingleFeatureDiffExp.COMBINE_SHIFTED_DISTRIBS_FOR_EACH_REPLICATE_PAIR = true;
        }
        if(isSet(NOOUTLIERCORRECTION)) {
            SingleFeatureDiffExp.CORRECT_OUTLIERS = false;
        }
        ErrorEstimationDistribution.setFcWithBin(getDouble(FCBINWIDTH));
    }
}
