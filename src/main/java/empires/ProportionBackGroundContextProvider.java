package empires;

import lmu.utils.LogConfig;
import lmu.utils.Logger;
import lmu.utils.NumUtils;

import java.util.List;
import java.util.Vector;

import static lmu.utils.IteratorUtils.rangev;
import static lmu.utils.ObjectGetter.map;

public class ProportionBackGroundContextProvider implements BackgroundContextFuzzficationStrategyProvider {

    Logger log = LogConfig.getLogger();
    Vector<Double> quantiles = new Vector<>();

    public ProportionBackGroundContextProvider(Vector<Double> proportions) {
        assert(proportions.size() > 0);
        double sum = NumUtils.sum(proportions);
        double pre = 0.0;
        for(double d : proportions) {
            double cumulative = pre + d / sum;
            quantiles.add(cumulative);
            pre = cumulative;
        }

    }

    @Override
    public empires.BackGroundContext getContexts(String replicateSetName, Vector<String> sampleNames, Vector<ReplicatedMeasurement> normalized) {

        empires.BackGroundContext backGroundContext = new BackGroundContext();
        backGroundContext.errorBackGrounds = new Vector<>();
        backGroundContext.meanSignalValue2Error = new Vector<>();
        NumUtils.sort(normalized, (_m) -> _m.mean, false);
        int lastidx = 0;
        for(int i=0; i<quantiles.size() ; i++) {
            int idx = (i == quantiles.size()) ? normalized.size() : (int)(normalized.size() * quantiles.get(i));
            if(idx == lastidx)
                continue;
            log.info("next quantile: %d: %.2f select %d-%d/%d", i, quantiles.get(i), lastidx, idx, normalized.size());

            backGroundContext.meanSignalValue2Error.add(NumUtils.mean(map(rangev(lastidx, idx), (_i) -> normalized.get(_i).mean)));
            empires.ErrorEstimationDistribution err = new ErrorEstimationDistribution(lastidx, idx, (_i) -> normalized.get(_i).nonNanValues);

            for(int j=lastidx; j<idx; j++) {
                backGroundContext.featureIdx2ErrorBackgroundIdx.put(j, backGroundContext.errorBackGrounds.size());
            }
            lastidx = idx;

            backGroundContext.errorBackGrounds.add(err);

        }
        return backGroundContext;
    }
}
