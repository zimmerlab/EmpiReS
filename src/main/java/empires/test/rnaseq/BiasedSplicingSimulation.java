package empires.test.rnaseq;

import empires.rnaseq.GFFBasedIsoformRegionGetter;
import empires.rnaseq.MultiIsoformRegion;
import lmu.utils.ImageUtils;
import lmu.utils.NumUtils;
import lmu.utils.OptionParser;
import lmu.utils.SimpleOptionParser;
import lmu.utils.plotting.PlotCreator;
import lmu.utils.plotting.CachedPlotCreator;

import java.util.Vector;
import java.util.function.Function;

import static lmu.utils.ObjectGetter.filter_and_map;
import static lmu.utils.ObjectGetter.map;

public class BiasedSplicingSimulation {
    public static void main(String[] args) {
        SimpleOptionParser cmd = new SimpleOptionParser("yeast", "human");
        cmd.setFile("yeast", "human");
        if(!OptionParser.parseParams(args, true,false, true, true, cmd))
            return;

        empires.rnaseq.GFFBasedIsoformRegionGetter yeast = new empires.rnaseq.GFFBasedIsoformRegionGetter(cmd.getFile("yeast"), null, null);
        empires.rnaseq.GFFBasedIsoformRegionGetter human = new GFFBasedIsoformRegionGetter(cmd.getFile("human"), null, null);

        Function<MultiIsoformRegion, Integer> g2longsttr = (_m) -> NumUtils.median(_m.isoforms.values(), (_t) -> _t.getCoveredLength());

        Vector<Integer> ylengths = map(yeast.getRegions(), (_m) -> g2longsttr.apply(_m));
        Vector<Integer> hlengths = filter_and_map(human.getRegions(), (_m) -> _m.isoforms.size() > 1 && _m.biotype.equals("protein_coding"),  (_m) -> g2longsttr.apply(_m));

        PlotCreator pc = CachedPlotCreator.getPlotCreator();
        pc.cumhist("yeast", ylengths, ylengths.size(), false, true);
        pc.cumhist("human", hlengths, hlengths.size(), false, true);
        pc.setLimits(0.0, 6000.0, null, null);
        pc.setLabels("glength", "rel.freq", "bottomright");
        ImageUtils.showImage("test", pc.getImage(), false);
    }
}
