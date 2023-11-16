package empires.plotting;

import lmu.utils.*;

import lmu.utils.plotting.CachedPlotCreator;
import lmu.utils.plotting.PlotCreator;
import lmu.utils.swing.PagedDataTable;
import empires.rnaseq.simulation.PositionBiasFactory;

import java.util.Vector;

import static lmu.utils.IteratorUtils.rangev;
import static lmu.utils.ObjectGetter.apply;
import static lmu.utils.ObjectGetter.filter;
import static empires.rnaseq.simulation.PositionBiasFactory.readTrBiases;

public class PositionBiasCheck {

    static class TrInfo {
        String id;
        int[] starts;
        double dotscore;
        int numReads;

        public static DataTable toTable(Vector<TrInfo> infos) {
            return DataTable.buildTable(infos,
                    DataTable.buildHeader("id", (TrInfo t) -> t.id)
                    .add("length", (t) -> t.starts.length)
                    .add("numReads", (t) -> t.numReads)
                    .add("dotscore", (t) -> t.dotscore)

            );
        }

        public void show() {
            PlotCreator pc = lmu.utils.plotting.CachedPlotCreator.getPlotCreator();
            Vector<Integer> r = rangev(starts.length);
            pc.line("", r, (_x) -> _x, (_x) -> starts[_x]);
            pc.setLabels("position", "bias", null);
            pc.setTitle(id);
            ImageUtils.showImage(String.format("%s reads: %d dot to uniform: %.2f", id, numReads, dotscore) , pc.getImage(), false);
        }
    }

    public static void main(String[] args) {
        SimpleOptionParser cmd = new SimpleOptionParser("biases");
        cmd.setFile("biases");
        if(!OptionParser.parseParams(args, true, false, true, true, cmd))
            return;


        Vector<TrInfo> trInfos = new Vector<>();

        apply(readTrBiases(cmd.getFile("biases")).entrySet(),
                (_e) ->
                {
                        TrInfo trInfo = new TrInfo();
                        trInfo.id = _e.getKey();
                        trInfo.starts = _e.getValue();
                        trInfo.numReads = NumUtils.sum(trInfo.starts);
                        trInfo.dotscore = PositionBiasFactory.getDotScoreToUniform(trInfo.starts);
                        trInfos.add(trInfo);
                });


        PlotCreator pc = CachedPlotCreator.getPlotCreator();
        Vector<TrInfo> filteredTrInfos = filter(trInfos, (_t) -> _t.numReads  > 1000);

        pc.cumhist("", filteredTrInfos, (_t) -> _t.dotscore, filteredTrInfos.size(), false, false);
        pc.setLabels("dotscore vs uniform", "frequence", null);
        ImageUtils.showImage("bias - to -uniform -dotscore", pc.getImage(), true);

        PagedDataTable.getDetailViewFrame(TrInfo.toTable(filteredTrInfos), null, null)
                .addMenu("detail", (ObjectGetter og) -> ((TrInfo)og.getInData()).show());
    }
}
