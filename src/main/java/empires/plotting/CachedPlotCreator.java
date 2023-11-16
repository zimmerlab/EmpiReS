package empires.plotting;

import lmu.utils.plotting.PlotCreator;

public class CachedPlotCreator {
    static PlotCreator plotCreator = null;
    static long last_access = 0;

    public static PlotCreator getPlotCreator() {
        long t = System.currentTimeMillis();
        if(plotCreator == null) {
            last_access = t;
            return plotCreator = new PlotCreator();
        }
        if(t - last_access > 1000) {
            plotCreator.destroy();
            plotCreator = new PlotCreator();
        }
        last_access = t;
        return plotCreator;

    }

    public static void close() {
        if(plotCreator == null)
            return;

        plotCreator.destroy();
        plotCreator = null;
    }
}
