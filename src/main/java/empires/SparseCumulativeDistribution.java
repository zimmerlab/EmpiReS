package nlEmpiRe;


import lmu.utils.plotting.PlotCreator;

import java.awt.*;
import java.util.Vector;
import java.util.function.BiFunction;

import static lmu.utils.IteratorUtils.rangev;
import static lmu.utils.ObjectGetter.toVector;
import static lmu.utils.ObjectGetter.map;

public class SparseCumulativeDistribution {

    public static Vector<Double> QUANTILES = toVector(0.05, 0.25, 0.5, 0.75, 0.95);

    public static class Point {
        public double statisticValue;
        public double percentDataUnder;

        Point(double val, double perc) {
            statisticValue = val;
            percentDataUnder = perc;
        }

        @Override
        public String toString() {
            return String.format("%.2f = %.2f%%", statisticValue, percentDataUnder);
        }
    }

    public Vector<Point> data = new Vector<>();
    public Vector<Double> quantiles = null;
    double estimatedValue;

    public<A> SparseCumulativeDistribution(A distribution, BiFunction<A, Double, Double> statisticValue2PercentageGetter,
                                           BiFunction<A, Double, Double> percentage2StatisticValueGetter, double valueStep) {

        double alpha = 0.01;
        double minValue = percentage2StatisticValueGetter.apply(distribution, alpha);
        double maxValue = percentage2StatisticValueGetter.apply(distribution, 1.0 - alpha);

        quantiles = map(QUANTILES, (_q) -> percentage2StatisticValueGetter.apply(distribution, _q));

        double bestW = 0.0;
        double pre = 0.0;
        for(double val = minValue; val <= maxValue; val += valueStep) {

            double perc = statisticValue2PercentageGetter.apply(distribution, val);
            if(perc - pre > bestW) {
                bestW = perc - pre;
                estimatedValue = val;
            }
            pre = perc;
            data.add(new Point(val, perc));
        }

    }

    public SparseCumulativeDistribution(ErrorEstimationDistribution fcInfos) {
        this(fcInfos, 0.01);
    }

    public SparseCumulativeDistribution(ErrorEstimationDistribution fcInfos, double fcStep) {
        this(fcInfos, (_d, _p) -> _d.getCumulativeFrequencyToFoldChange(_p ) , (_d, _p) -> _d.getFoldChangeToCumulativeFrequency(_p),  fcStep);
        estimatedValue = fcInfos.getMostProbableFcWindowCenter();
    }

    @Override
    public String toString() {
        return "" + data;
    }

    public void drawLine(PlotCreator pc, String name) {
        drawLine(pc, name, 1.0);
    }
    public void drawLine(PlotCreator pc, String name, double factor) {
        Vector<Integer> pointIndeces = rangev(1, data.size());
        pc.line(String.format("%s (%.2f)", name, estimatedValue), pointIndeces, (_i) -> (0.5 * (data.get(_i).statisticValue + data.get(_i - 1).statisticValue)),
                (_i) -> factor * (data.get(_i).percentDataUnder - data.get(_i - 1).percentDataUnder));
        pc.abline(null, estimatedValue, null, null,null);
    }

    public void drawCumulative(PlotCreator pc, String name, double factor) {
        Vector<Integer> pointIndeces = rangev(1, data.size());
        pc.line(String.format("%s (%.2f)", name, estimatedValue), pointIndeces, (_i) -> (0.5 * (data.get(_i).statisticValue + data.get(_i - 1).statisticValue)),
                (_i) -> factor * (data.get(_i).percentDataUnder));

    }

    public void drawBox(PlotCreator.BoxPlotBuilder bpb, String name) {
        bpb.addBox(name, quantiles);
    }

    public void drawBox(PlotCreator.BoxPlotBuilder bpb, String name, Color c) {
        bpb.addBox(name, quantiles, c);
    }
}
