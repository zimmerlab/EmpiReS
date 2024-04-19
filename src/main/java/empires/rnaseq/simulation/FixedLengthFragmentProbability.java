package nlEmpiRe.rnaseq.simulation;

import lmu.utils.FRuntimeException;
import lmu.utils.UPair;

import java.util.Arrays;
import java.util.Vector;

public class FixedLengthFragmentProbability {
    int length;
    double probability;

    Vector<UPair<Integer>> fragments = new Vector<>();
    Vector<Integer> frequencies = new Vector<>();
    int[] cumulativeFrequencies;

    public FixedLengthFragmentProbability(int length, double probability) {
        this.length = length;
        this.probability = probability;
    }


    public void add(UPair<Integer> fragment, int freq) {
        fragments.add(fragment);
        frequencies.add(freq);
        cumulativeFrequencies = null;
    }

    void compile() {
        if(cumulativeFrequencies != null)
            return;
        if(fragments.size() == 0)
            throw  new FRuntimeException("no fragments added yet!");
        cumulativeFrequencies = new int[frequencies.size()];
        cumulativeFrequencies[0] = frequencies.get(0);
        for(int i=1; i<cumulativeFrequencies.length; i++) {
            cumulativeFrequencies[i] = cumulativeFrequencies[i-1] + frequencies.get(i);
        }

    }
    public UPair<Integer> sample() {
        compile();
        int TOTAL = cumulativeFrequencies[cumulativeFrequencies.length - 1];
        int selected = (int)(Math.random() * TOTAL);
        int hitIdx = Arrays.binarySearch(cumulativeFrequencies, selected);
        if(hitIdx < 0) {
            hitIdx = - hitIdx - 1;
        }
        return fragments.get(Math.min(fragments.size() - 1, hitIdx));

    }
}
