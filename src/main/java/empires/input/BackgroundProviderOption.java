package empires.input;

import empires.PreFuzzifiedBackGroundContextProvider;
import empires.ProportionBackGroundContextProvider;
import lmu.utils.FRuntimeException;
import lmu.utils.FileUtils;
import lmu.utils.MapBuilder;
import lmu.utils.SimpleOptionParser;
import empires.AutoBackGroundContextProvider;
import empires.BackgroundContextFuzzficationStrategyProvider;

import java.io.File;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Vector;

import static lmu.utils.ObjectGetter.map;
import static lmu.utils.ObjectGetter.toVector;

public class BackgroundProviderOption extends  SimpleOptionParser{
    public static String FUZZYPROPORTIONS = "fuzzyproportions";
    public static String FUZZYASSIGNMENTS = "fuzzyassignments";
    public static String FUZZYASSIGNMENTSID= "fuzzyassignmentsId";
    public static String MAXSDDIFF = "maxsddiff";

    public BackgroundProviderOption() {
        super(FUZZYPROPORTIONS, FUZZYASSIGNMENTS, FUZZYASSIGNMENTSID, MAXSDDIFF);
        setName("background-distribution estimate related options");
        setFile(FUZZYPROPORTIONS, FUZZYASSIGNMENTS);
        setOptional(FUZZYPROPORTIONS, FUZZYASSIGNMENTS);
        setDefault(FUZZYASSIGNMENTSID, "id");
        setDouble(MAXSDDIFF);
        setDefault(MAXSDDIFF, "0.1");
    }

    public BackgroundContextFuzzficationStrategyProvider getStrategy() {
        File proportions = getOptionalFile(FUZZYPROPORTIONS);
        if(proportions != null) {
            Vector<String> lines = FileUtils.readNLines(proportions, 1);
            String ERR = String.format("invalid file: %s expected one line containing double values separated by whitespaces", proportions.getAbsolutePath());
            if(lines.size() != 1)
                throw new FRuntimeException(ERR);
            try
            {
                return new ProportionBackGroundContextProvider(map(lines.get(0).split("\\s+"), (_s) -> Double.parseDouble(_s)));
            } catch(NumberFormatException nfe) {
                throw new FRuntimeException(ERR+" error at parsing numbers from: %s", lines.get(0));
            }
        }
        File assignments = getOptionalFile(FUZZYASSIGNMENTS);
        if(assignments != null) {
            HashMap<String, HashMap<String, String>> feature2sample2fuzzyclass = new HashMap<>();
            Vector<String> headers = FileUtils.getHeaders(assignments);
            String idfield = getValue(FUZZYASSIGNMENTSID);
            int idIdx = headers.indexOf(idfield);
            if(idIdx < 0) {
                System.err.printf("did not find id header: %s in %s! Pls. provide a valid one with -%s.\nheaders in %s: %s\n",
                        idfield, assignments.getAbsolutePath(), FUZZYASSIGNMENTSID,
                        assignments, headers);
                System.exit(-1);
            }

            Iterator<String[]> it = FileUtils.getFieldSetIterator(assignments, "\t");
            it.next();
            int ln = 0;
            while(it.hasNext()) {
                ln++;
                String[] sp = it.next();
                if(sp.length != headers.size()) {
                    System.err.printf("invalid fuzzy assignment file: %s line: %d expected %d fields, got %d!\nheaders: %s\nfields: %s\n",
                            assignments.getName(), ln, headers.size(), sp.length,
                            headers, Arrays.toString(sp));
                    System.exit(-1);
                }
                String fn = sp[idIdx];
                for(int i=0; i<sp.length; i++){
                    if(i == idIdx)
                        continue;
                    MapBuilder.updateM(feature2sample2fuzzyclass, fn, headers.get(i), sp[i]);
                }
            }
            return new PreFuzzifiedBackGroundContextProvider(feature2sample2fuzzyclass);

        }

        AutoBackGroundContextProvider.MERGE_BACKGROUNDS_WITH_MAX_SD_DIFF = getDouble(MAXSDDIFF);
        return new AutoBackGroundContextProvider();
    }


}
