package empires.rnaseq.simulation;

import lmu.utils.*;

import java.util.*;

import static lmu.utils.IteratorUtils.rangev;
import static lmu.utils.ObjectGetter.*;

public class ScalableTranscriptInfo {
    Logger log = LogConfig.getLogger();
    int length;
    TreeMap<Integer, empires.rnaseq.simulation.FixedLengthFragmentProbability> length2info = new TreeMap<>();


    ScalableTranscriptInfo(int[] input, int[] input_cumulative, int targetlength, int readLength, int maxLength) {
        this.length = targetlength;

        //remove / add empty places randomly if size does not fit
        //remove -> shrink to most popuplated area
        //add -> create a new array and select random positions to add new
        int[] inputBiases = (input.length == targetlength) ? input : null;
        if(inputBiases == null && targetlength < input.length) {
            int bestStart = 0;
            int bestVal = input_cumulative[targetlength - 1];
            for(int i = 1; i<input_cumulative.length - targetlength; i++) {
                int coverage = input_cumulative[i + targetlength] - input_cumulative[i];
                if(coverage <= bestVal)
                    continue;
                bestStart = i;
                bestVal = coverage;
            }

            inputBiases = Arrays.copyOfRange(input, bestStart, bestStart + input.length);
        }
        if(inputBiases == null) {
            //have to add empty positions
            inputBiases = new int[targetlength];
            HashSet<Integer> emptyList = toSet(shuffleN(rangev(targetlength), targetlength - input.length));
            int targetIdx = 0;
            for(int i=0; i<input.length; i++) {
                if(emptyList.contains(i)) {
                    targetIdx++;
                    continue;
                }
                inputBiases[targetIdx++] = input[i];
            }
        }

        //add pseudo value to avoid non-simulatable setup
        for(int i=0; i<inputBiases.length; i++) {
            inputBiases[i]++;
        }
        //sample start / end

        HashMap<UPair<Integer>, Integer> start2end2req = new HashMap<>();
        HashMap<Integer, Long> length2total = new HashMap<>();
        HashMap<Integer, Vector<UPair<Integer>>> length2variants = new HashMap<>();

        int[] IB = inputBiases;
        int totalNum = NumUtils.sum(input);
        int THRESHOLD = totalNum >> 1;
        Vector<Integer> nonnullPositions = filter(rangev(length), (_i) -> IB[_i] > 0);
        int initNumValids = nonnullPositions.size();
        int restrictedValids = 0;
        int totake = 200;
        int maxStart = targetlength - readLength;
        int minEnd = readLength;

        if(nonnullPositions.size() > totake * 2) {
            TreeMap<Integer, Vector<Integer>> freq2pos = new TreeMap<>();

            int prefSum = 0;
            for(int i=0; i<nonnullPositions.size(); i++) {
                int pos = nonnullPositions.get(i);
                if(pos > maxStart)
                    break;
                int freq = inputBiases[pos];

                freq2pos.computeIfAbsent(freq, (_x) -> new Vector<>()).add(pos);
                prefSum += freq;
                if(prefSum > THRESHOLD)
                    break;
            }
            HashSet<Integer> prefixPositions = new HashSet<>();
            while(freq2pos.size() > 0 && prefixPositions.size() < totake) {
                int max = freq2pos.lastKey();
                Vector<Integer> poslist = freq2pos.remove(max);
                prefixPositions.addAll(VectorUtils.slice(poslist, 0, totake - prefixPositions.size()));

            }

            int suffSum = 0;
            freq2pos.clear();
            for(int i=nonnullPositions.size() -1; i>=0; i--) {
                int pos = nonnullPositions.get(i);
                if(pos < minEnd)
                    break;
                int freq = inputBiases[pos];

                freq2pos.computeIfAbsent(freq, (_x) -> new Vector<>()).add(pos);
                suffSum += freq;
                if(suffSum > THRESHOLD)
                    break;
            }
            HashSet<Integer> suffixPositions = new HashSet<>();
            while(freq2pos.size() > 0 && suffixPositions.size() < totake) {
                int max = freq2pos.lastKey();
                Vector<Integer> poslist = freq2pos.remove(max);
                suffixPositions.addAll(VectorUtils.slice(poslist, 0, totake - suffixPositions.size()));
            }

            prefixPositions.addAll(suffixPositions);
            restrictedValids = prefixPositions.size();
            log.trace("reduce %d -> %d positions", nonnullPositions.size(), prefixPositions.size());
            nonnullPositions.clear();
            nonnullPositions.addAll(prefixPositions);
            Collections.sort(nonnullPositions);
        }

        int minPos = nonnullPositions.get(0);
        int maxPos = nonnullPositions.get(nonnullPositions.size() - 1);

        //take the 500 ost probable start pos and 500
        log.trace("will sample positions: for length: %d non-null: %d maxl: %d\n", length, nonnullPositions.size(), maxLength);
        long t1 = System.currentTimeMillis();


        long totalComb  = 0;
        for(int i = 0;  i < nonnullPositions.size(); i++) {
            int startPos = nonnullPositions.get(i);
            if(startPos > maxStart)
                break;

            int freqStart = inputBiases[startPos];

            int addedInvalid = 0;
            for(int j=i+1; addedInvalid < 3 && j<nonnullPositions.size(); j++) {
                int endPos = nonnullPositions.get(j);
                int frLength = endPos - startPos + 1;
                if (frLength < readLength)
                    continue;
                if (frLength > maxLength) {
                    addedInvalid++;
                }


                int freqEnd = inputBiases[endPos];
                int prod = freqStart * freqEnd;
                totalComb += prod;


                UPair<Integer> fr = UPair.createU(startPos, endPos);
                start2end2req.put(fr, prod);
                length2total.put(frLength, length2total.getOrDefault(frLength, 0l) + prod);
                MapBuilder.updateV(length2variants, frLength, fr);
                if (System.currentTimeMillis() - t1 > 10_000) {
                    log.info("sampling in progress... took so far: %.2f sec targetLengh: %d positions: %d in/restricted: %d/%d", (System.currentTimeMillis() - t1) / 1000.0, targetlength, nonnullPositions.size(),
                            initNumValids, restrictedValids);
                }
            }
        }
        if(length2variants.size() == 0) {
            //no valid positions found
            log.warn("NO VALID POSITION for length: %d found nonnull: %s initvalids: %d", targetlength, nonnullPositions, initNumValids);
            //add a random positions as start
            //check if there were a valid start

        }
        long t2 = System.currentTimeMillis();
        if(t2  - t1 > 10_000) {
            log.info("sampling ready took: %.2f sec targetLengh: %d positions: %d in/restricted: %d/%d", (t2 - t1) / 1000.0, targetlength, nonnullPositions.size(),
                    initNumValids, restrictedValids);
        }


        double norm = 1.0 / totalComb;

        for(int length : length2total.keySet()) {
            double prob = length2total.get(length) * norm;
            empires.rnaseq.simulation.FixedLengthFragmentProbability flp = new FixedLengthFragmentProbability(length, prob);
            for(UPair<Integer> fr : length2variants.get(length)) {
                flp.add(fr, start2end2req.get(fr));
            }
            length2info.put(length, flp);
        }
        log.trace("length2info ready");
        if(length2info.size() == 0)
            throw new FRuntimeException("no length2info! targetlength: %d readlength: %d maxlength: %d length2total: %s", targetlength, readLength,  maxLength, length2total);
    }

    public UPair<Integer> getFragment(int targetFrLength) {

        //select at most - 2 and + 2
        HashSet<Integer> validLengths = new HashSet<>();
        Integer check = targetFrLength;
        while(check != null && validLengths.size() < 3) {
            check = length2info.floorKey(check);
            if(check == null)
                break;
            validLengths.add(check);
            check--;

        }
        check = targetFrLength;
        int nadded = 0;
        while (check != null && nadded < 3) {
            check = length2info.ceilingKey(check);
            if(check == null)
                break;

            nadded++;
            validLengths.add(check);
            check++;
        }
        if(validLengths.size() == 0) {
            throw new FRuntimeException("no valid length at all for length: %d length2info: %s", targetFrLength, validLengths);
        }
        if(validLengths.size() == 1) {
            return length2info.get(first(validLengths)).sample();
        }
        Vector<Integer> keys = toVector(validLengths);
        double[] cumprob = new double[keys.size()];
        cumprob[0] = length2info.get(keys.get(0)).probability;
        for(int i=1; i<cumprob.length; i++) {
            cumprob[i] = cumprob[i-1] + length2info.get(keys.get(i)).probability;
        }
        double rnd = Math.random() * cumprob[cumprob.length - 1];
        int idx = Arrays.binarySearch(cumprob, rnd);
        if(idx < 0) {
            idx = - idx - 1;
        }
        return length2info.get(keys.get(idx)).sample();
   }
}
