package empires.plotting;

import lmu.utils.ImageUtils;
import lmu.utils.UPair;
import lmu.utils.plotting.CachedPlotCreator;
import lmu.utils.plotting.PlotCreator;
import empires.NormalizedReplicateSet;

import java.awt.image.BufferedImage;
import java.io.File;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Vector;

public class DataOverview {
    HashMap<String, NormalizedReplicateSet> data;

    public DataOverview(HashMap<String, NormalizedReplicateSet> data) {
        this.data = data;
    }


    public void write(File od, PrintWriter pw) {

        PlotCreator pc = CachedPlotCreator.getPlotCreator();

        double quantile = 0.75;

        for(NormalizedReplicateSet nrs : data.values()) {
            empires.plotting.NormalizedReplicateSetPlotting nrsp = new empires.plotting.NormalizedReplicateSetPlotting(nrs);
            Vector<UPair<Double>> sig2noise = nrsp.getSignal2Error(quantile);
            pc.line(nrs.getInData().getReplicateSetName(), sig2noise, (_p) -> _p.getFirst(), (_p) -> _p.getSecond());
        }
        pc.setLabels("signal", String.format("noise at quantile: %.2f%%", 100 * quantile), "topright");
        File ov_s2n = new File(od, "overview_signal2noise.png");
        pc.writeImage(ov_s2n);

        pw.printf("<img src='%s'>\n", ov_s2n.getName());

        pw.flush();
        for(NormalizedReplicateSet nrs : data.values()) {
            String cond = nrs.getInData().getReplicateSetName();
            pw.printf("<h2>%s (%d replicates)</h2>\n", cond, nrs.getInData().getNumReplicates());
            pw.flush();
            empires.plotting.NormalizedReplicateSetPlotting nrsp = new NormalizedReplicateSetPlotting(nrs);
            pc.setWidth(800);
            pc.reset();


            File of_samples = new File(od, cond + "_samples.png");
            try {
                BufferedImage bim = nrsp.plotSamplePairHeatMap(pc);



                ImageUtils.saveImage(bim,of_samples);
                pw.printf("<img src='%s'>\n", of_samples.getName());

            } catch (Exception e) {
                pw.printf("<br><b><font color='red'>error: %s at creating sample overview</font></b><br>\n", e.getMessage());
                pc.reset();
            }



            File sig2noise = new File(od, cond+"_sig2noise.png");

            BufferedImage bim2 = nrsp.plotBackGroundDistribs(pc);
            ImageUtils.saveImage(bim2,sig2noise);

            pw.printf("\n<br>\n<img src='%s'>\n", sig2noise.getName());

            pw.flush();
        }



    }
}
