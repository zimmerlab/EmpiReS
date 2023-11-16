package empires.input;

import lmu.utils.FRuntimeException;
import lmu.utils.FileUtils;
import lmu.utils.MapBuilder;

import java.io.File;
import java.net.IDN;
import java.util.*;

import static lmu.utils.ObjectGetter.*;

public class ExperimentDescriptor {

    HashMap<String, Vector<String>> condition2replicatenames;
    Vector<String> replicatelist;
    Set<String> trues;
    HashMap<String, String> rep2cond = new HashMap<>();
    HashMap<Integer, String> idx2cond = new HashMap<>();

    public ExperimentDescriptor(HashMap<String, Vector<String>> condition2replicatenames, Vector<String> samplelist, Set<String> trues) {
        this.condition2replicatenames = condition2replicatenames;
        this.replicatelist = samplelist;

        for(Map.Entry<String, Vector<String>> e : condition2replicatenames.entrySet()) {
            apply(e.getValue(), (_rep) -> rep2cond.put(_rep, e.getKey()));
        }

        applyIndex(replicatelist.size(), (_i) -> idx2cond.put(_i, rep2cond.get(replicatelist.get(_i))));
        this.trues = trues;
    }

    static Vector<String> IDNAME_VARIANTS = toVector("name", "id", "Id", "ID", "sample", "label");

    public ExperimentDescriptor(File cond2reps, File samplelist, File trues) {
        if(cond2reps != null && cond2reps.exists()) {
            condition2replicatenames = FileUtils.readMultiMap(cond2reps, "\t", 0, 1);

            for(Map.Entry<String, Vector<String>> e : condition2replicatenames.entrySet()) {
                apply(e.getValue(), (_rep) -> rep2cond.put(_rep, e.getKey()));
            }

        }


        if(samplelist != null) {
              Set<String> headers = toSet(FileUtils.getHeaders(samplelist));
              boolean gotcondition_info = headers.contains("condition");
              FileUtils.HeaderedReader hr = FileUtils.getHeaderedReader("\t");
              boolean found = false;

              for(String variant : IDNAME_VARIANTS) {
                  if(!headers.contains(variant))
                      continue;
                  found = true;
                  hr.add(variant, "name");
              }
              if(!found)
                  throw new FRuntimeException("could not find sample name column in %s headers: %s tested for: %s", samplelist.getName(), headers, IDNAME_VARIANTS);

            if(gotcondition_info) {
                hr.add("condition", "cond");
            }


            replicatelist = new Vector<>();
            if(gotcondition_info) {
                condition2replicatenames = new HashMap<>();
                rep2cond = new HashMap<>();
            }

            apply(FileUtils.getHrIterator(samplelist, hr,
                     false), (_hr) ->
                    {
                        String replicate = _hr.getString("name");
                        replicatelist.add(replicate);
                        if(gotcondition_info) {
                            String cond = _hr.getString("cond");
                            MapBuilder.updateV(condition2replicatenames, cond, replicate);
                            rep2cond.put(replicate, cond);
                        }
                    }
            );

            applyIndex(replicatelist.size(), (_i) -> idx2cond.put(_i, rep2cond.get(replicatelist.get(_i))));
         }
         if(trues != null) {
             this.trues = FileUtils.readSet(trues);
         }
    }

    public HashMap<String, Vector<String>> getCond2Reps() {
        return condition2replicatenames;
    }

    public boolean gotTrues() {
        return trues != null;
    }

    public Set<String> getTrues() {
        return trues;
    }

    public HashMap<Integer, String> getReplicateIndex2Condition() {
        return idx2cond;
    }
    public boolean isTrue(String id) {
        return trues != null && trues.contains(id);
    }
}
