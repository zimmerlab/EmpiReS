package empires.test.rnaseq;

import empires.rnaseq.ReducedTranscriptPresentation;
import lmu.utils.*;
import lmu.utils.plotting.*;
import lmu.utils.tuple.Tuple4;

import java.awt.*;
import java.awt.image.BufferedImage;
import java.util.*;

import static lmu.utils.ObjectGetter.*;

public class ReducedTranscriptPresentationVisualization implements empires.rnaseq.ReducedTranscriptPresentation.DetailProcesser {

    Vector<BufferedImage> images = new Vector<>();
    public static class EQClassDecorator extends RegionPlot.RegionDecorator
    {
        int id;

        public EQClassDecorator(int id, Color color)
        {
            super(color, true);
            this.id = id;
        }

        public java.awt.Rectangle draw(DGraphics g, int y_shift, java.awt.Rectangle box)
        {
            int x=(int)box.getX();
            int y=y_shift+(int)box.getY();
            int height=(int)box.getHeight();
            int width=(int)box.getWidth();
            java.awt.Rectangle r = super.draw(g, y_shift, box);
            g.setColor(Color.LIGHT_GRAY);

            g.drawString(String.format("%s", id), x + 2, y + height -2);

            return r;

        }
    }

    public BufferedImage plot(HashMap<Tuple, Double> restrictedEQCounts, String title, Set<String> restrictToTranscripts) {
        HashMap<String, Vector<Tuple>> tr2eqs = new HashMap<>();
        HashMap<String, Double>  tr2count = empires.rnaseq.ReducedTranscriptPresentation.getPartialCounts(restrictedEQCounts.keySet(), restrictedEQCounts, tr2eqs);
        if(restrictToTranscripts != null) {
            Set<String> toremove = SetInfo.minus(tr2count.keySet(), restrictToTranscripts);
            apply(toremove, (_k) -> tr2count.remove(_k));
        }
        int maxCount = NumUtils.max(tr2count.values()).intValue();
        Vector<Tuple> tuples = Pair.reduceSecond(Pair.convert_reverse_sorted(restrictedEQCounts, false));
        HashMap<Tuple, Integer> tuple2idx = buildIndexMap(tuples);
        Iterator<Color> cit = ColorManager.getColorIterator();
        HashMap<Tuple, Color> tuple2color = buildMap(tuples, (_t) -> cit.next());

        int spacer = 1;
        int maxlength = (NumUtils.sum(restrictedEQCounts.values())).intValue() + NumUtils.max(tr2eqs.values(), (_v) -> _v.size()) * spacer;
        int targetLength = 800;
        double factor = (maxlength < targetLength) ? 1.0 : targetLength / (0.0 + maxlength);
        ZoomedCoordinateConverter converter = new ZoomedCoordinateConverter(maxlength, factor, null);

        ZoomedMultiPicPlaneFrame zoomedMultiPicPlaneFrame = new ZoomedMultiPicPlaneFrame(maxlength);

        zoomedMultiPicPlaneFrame.setZoomedConverter(converter);
        RegionPicturePane rpp = zoomedMultiPicPlaneFrame.createRegionPane(title);


        for(String tr : Pair.reduceSecond(Pair.convert_reverse_sorted(tr2count, false))) {
            if(restrictToTranscripts != null && ! restrictToTranscripts.contains(tr))
                continue;

            System.out.printf("next tr: %s tuples: %d\n", tr, tr2eqs.get(tr).size());

            Vector<Integer> eqIndeces = toSortedVector(map(tr2eqs.get(tr), (_t) -> tuple2idx.get(_t)), true);
            Vector<Region1D> regions = new Vector<>();
            Vector<RegionPlot.RegionDecorator> decorators = new Vector<>();
            int lastWidth = 0;

            RegionPlot rp = rpp.addRegion(tr);
            rp.setDoDenseLayout(false);
            for(int eqIdx : eqIndeces) {
                Tuple t = tuples.get(eqIdx);
                int width  = restrictedEQCounts.get(t).intValue();
                Region1D r = new Region1D(lastWidth, lastWidth + width);
                lastWidth = width + lastWidth + spacer;

                r.setObject(eqIdx);
                regions.add(r);
                rp.addRegion(r, new EQClassDecorator(eqIdx, tuple2color.get(t)));

            }




        }
        zoomedMultiPicPlaneFrame.setTitle("test title");
        zoomedMultiPicPlaneFrame.update(true);
        zoomedMultiPicPlaneFrame.setAutoSize();
        return zoomedMultiPicPlaneFrame.getCombinedImage();


    }
    int iter = 0;
    empires.rnaseq.ReducedTranscriptPresentation reducer;
    HashMap<String, PriorityQueue<Pair<String, Double>>> tr2mostSim;
    HashMap<Tuple, Double> EqClassToCounts;
    HashMap<String, Double> tr2count;

    @Override
    public void init(ReducedTranscriptPresentation reducer, HashMap<String, Double> tr2Count, HashMap<String, Vector<Tuple>> tr2eqs, HashMap<Tuple, Double> EqClassToCounts) {
        this.reducer = reducer;
        this.EqClassToCounts = EqClassToCounts;
        this.tr2count = tr2Count;
        System.out.printf("init: %s\n\t%s\n", tr2Count, tr2eqs);
        images.add(plot(EqClassToCounts, "initaliziation", null));
    }

    @Override
    public void addIter(int iter, HashMap<String, Double> tr2maxInfo, String tr2remove, HashMap<String, PriorityQueue<Pair<String, Double>>> tr2mostSim) {
        System.out.printf("iter: %d\n\ttr2max: %s\ttr2mostsim: %s\n\n\n", iter, tr2maxInfo, tr2mostSim);

        this.tr2mostSim = tr2mostSim;

    }

    @Override
    public void addRemoves(int iter, Vector<Tuple4<String, String, Double, Double>> removes) {
        System.out.printf("iter: %d REMOVEs : %s\n", iter, removes);
        for(Tuple4<String, String, Double, Double> t : removes) {
            images.add(plot(EqClassToCounts, String.format("remove %s due to %s partial counts: %.0f, %.0f", t.get0(), t.get1(),
                    t.get2(), t.get3()),
                    toSet(t.get0(), t.get1())));
        }
    }

    @Override
    public void addRemove(String transcript, HashMap<Tuple, Double> restrictedEQCounts, HashMap<String, Double> restrictedTr2Count) {
        iter++;
        this.EqClassToCounts = restrictedEQCounts;
        images.add(plot(restrictedEQCounts, String.format("iter: %d after removing %s", iter, transcript), null));
        System.out.printf("remove: %s\n", transcript);
    }

    public void showSteps(){
        ImageUtils.showImage("steps", ImageUtils.vconcat(images), false);
    }
}
