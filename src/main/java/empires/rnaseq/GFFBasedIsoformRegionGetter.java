package empires.rnaseq;


import lmu.utils.*;

import java.io.File;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Vector;

import static lmu.utils.ObjectGetter.*;

/**
 * Created by csaba on 14/09/17.
 */
public class GFFBasedIsoformRegionGetter extends FileBasedIsoformRegionGetter
{

    static class GeneInfo
    {
        Logger log = LogConfig.getLogger();
        empires.rnaseq.GenomicRegion gr;
        String id;
        String name;
        String biotype;

        HashMap<String, String> prot2tr = new HashMap<>();
        HashMap<String, String> tr2prot = new HashMap<>();
        HashMap<String, String> tr2biotype = new HashMap<>();
        HashMap<String, Vector<Region1D>> trs = new HashMap<>();
        HashMap<String, Vector<Region1D>> cds = new HashMap<>();
        HashMap<String, Vector<Region1D>> tr2stop = new HashMap<>();

        public GeneInfo(empires.rnaseq.GenomicRegion gr, HashMap<String, String> infos)
        {
            this.gr = gr;
            this.id = infos.get("gene_id");
            this.biotype = infos.get("gene_biotype");
            if(this.biotype == null)
            {
                this.biotype = infos.get("gene_type");
            }
            this.name = infos.get("gene_name");

        }

        public void addExon(empires.rnaseq.GenomicRegion gr, HashMap<String, String> infos)
        {
            MapBuilder.updateV(trs, infos.get("transcript_id"), new Region1D(gr.getX1(), gr.getX2()));
        }

        public void addStopCodon(empires.rnaseq.GenomicRegion gr, HashMap<String, String> infos)
        {
            MapBuilder.updateV(tr2stop, infos.get("transcript_id"), new Region1D(gr.getX1(), gr.getX2()));

        }

        public void addCDSExon(empires.rnaseq.GenomicRegion gr, HashMap<String, String> infos, String line)
        {
            String pid = infos.get("protein_id");
            String trid = infos.get("transcript_id");
            //System.out.printf("add cds keys: >%s< pid: %s\n", infos, pid);
            pid = (pid != null) ? pid : infos.get("ccdsid");

            pid = (pid != null) ? pid : "p."+trid;
            if(pid == null)
                throw new FRuntimeException("no protein id found for cds: %s, line: >%s<", infos, line);


            MapBuilder.updateV(cds, pid, new Region1D(gr.getX1(), gr.getX2()));
            prot2tr.put(pid, trid);
            tr2prot.put(trid, pid);
        }

        static RegionVector toRV(Vector<Region1D> v)
        {
            return new RegionVector(NumUtils.sort(v, (_r) -> _r.getX1(),  false), true);
        }
        public empires.rnaseq.MultiIsoformRegion toMIR()
        {
            empires.rnaseq.MultiIsoformRegion mir = new empires.rnaseq.MultiIsoformRegion();
            mir.biotype = biotype;
            mir.id = id;
            mir.name = name;
            mir.chr  = gr.getChr();
            mir.coding2iso= prot2tr;
            mir.iso2coding = tr2prot;
            mir.strand = gr.getStrand();
            for(Map.Entry<String, Vector<Region1D>> e : trs.entrySet())
            {
                try {
                    mir.addRegionVector(e.getKey(), null, toRV(e.getValue()));
                }
                catch(Exception ex)
                {
                    throw new FRuntimeException("error: %s building transcript: %s gene: %s regions: %s", ex, ex.getMessage(), e.getKey(), id, e.getValue());
                }
            }

            for(Map.Entry<String, Vector<Region1D>> e : cds.entrySet())
            {
                String trid = prot2tr.get(e.getKey());

                Vector<Region1D> stop = tr2stop.get(trid);
                RegionVector rv = null;

                try
                {
                    //rv = toRV(toVector(toSet(e.getValue())));
                    rv = toRV(e.getValue());
                }
                catch(Exception ex)
                {
                    throw new FRuntimeException("error: %s building protein: %s gene: %s transcript: %s regions: %s", ex, ex.getMessage(),
                            e.getKey(), id, mir.coding2iso.get(e.getKey()),  e.getValue());
                }


                if(stop != null)
                {
                    //System.out.printf("%s merge: %s to %s:\n%s\n", e.getKey(), stop, rv, RegionVector.merge(rv, toRV(stop)));
                    rv = RegionVector.merge(rv, toRV(stop));
                }
                else
                {
                    RegionVector tr_rv = mir.isoforms.get(trid);
                    if(tr_rv == null) {
                        log.warn("no tr found with id: >%s< prot: %s", trid, e.getKey());
                    }
                    RegionVector _rv = rv;
                    RegionVector ending = toRV(filter(mir.isoforms.get(trid).substract(rv).getRegions(), (_r) -> (mir.strand) ? _r.getX2() > _rv.getX2() : _r.getX1() < _rv.getX1()));
                    int L = ending.getCoveredLength();
                    RegionVector stoprv = (mir.strand) ? ending.subRV(0, 3) : ending.subRV(L - 3, L);
                    //System.out.printf("%s merge from TR: %s ending %s:\nstop %s\n", e.getKey(), trid,  ending, stoprv);
                    //rv = RegionVector.merge(rv, stoprv);
                    //rv = RegionVector.merge(rv, stoprv);
                    //add from tr
                }
                //rv = (stop == null) ? rv : RegionVector.merge(rv, toRV(stop));
                RegionVector trv = mir.isoforms.get(trid);
                if(!trv.isSubVector(rv))
                {
                    mir.erroneous_isoforms.add(e.getKey());
                    //System.err.printf("gene: %s prot: %s is not a subvector of %s\n", mir.id, e.getKey(), trid);
                }
                mir.coding.put(e.getKey(), rv);
            }

            return mir;
        }
    }

    GeneInfo getOrAdd(HashMap<String, GeneInfo> genes, empires.rnaseq.GenomicRegion gr, HashMap<String, String> attribs)
    {
        String geneid = attribs.get("gene_id");
        if(geneid == null) {

            for(String alt: toVector("gene", "Parent")){

                String altv = attribs.get(alt);
                if(altv == null)
                    continue;

                log.warn("take >%s< as alternative geneid", alt);
                geneid = altv;

            }

            if(geneid == null)
                throw new FRuntimeException("no gene id given attribs: >%s< keys: %s", attribs, attribs.keySet());

        }


        GeneInfo gi = genes.get(geneid);
        if(gi != null)
            return gi;

        genes.put(geneid, gi = new GeneInfo(gr, attribs));
        return gi;
    }

    FileUtils.StringSplitter SCSP = FileUtils.getRecordSplitter(";");

    static Vector<String> splitAttribs(String record) {
        Vector<String> rv = new Vector<>();
        int start = 0;
        int lastSplitStart = 0;
        boolean inquota = false;
        while(start < record.length()) {
            int n1 = record.indexOf(';', start);
            int n2 = record.indexOf('"', start);
            if(n1 < 0 && n2 < 0) {
                String last = record.substring(start).trim();
                if(last.length() > 0) {
                    rv.add(last);
                }

                break;
            }
            int nextPos = (n1 < 0 || n2 < 0) ? (n1 < 0) ? n2 : n1 : Math.min(n1, n2);

            if(record.charAt(nextPos) == '"') {
                inquota = !inquota;
            } else {
                if(!inquota) {

                    String nrecord = record.substring(lastSplitStart, nextPos).trim();
                    rv.add(nrecord);

                    lastSplitStart = nextPos + 1;
                }
            }
            start = nextPos + 1;

        }
        return rv;
    }

    HashMap<String, String> readAttribs(String record, String line)
    {
        HashMap<String, String> attribs = new HashMap<>();
        for(String sp : splitAttribs(record))
        {

            int spidx = sp.indexOf('=');
            if(spidx >= 0) {
                log.warn("found '=' - take as attrib separator rather than ' ': %s", sp);
            } else {
                spidx = sp.indexOf(' ');
            }

            if(spidx < 0) {
                throw new FRuntimeException("invalid key attrib: >%s< in line: >%s<", sp, line);

            }


            String key = sp.substring(0, spidx);
            String value = sp.substring(spidx + 1);
            if(value.length() > 0 && value.charAt(0)=='\"')
            {
                value = value.substring(1, value.length() - 1);

            }
            attribs.put(key, value);
            //System.out.printf(">%s< add attris >%s< >%s<\n", sp, key, value);
        }
        return attribs;
    }

    interface ParseF
    {
        public void parse(GeneInfo gi, empires.rnaseq.GenomicRegion gr, HashMap<String, String> attribs, String line);
    }

    void parseLine(HashMap<String, GeneInfo> genes, FileUtils.Converter c, String line, ParseF pf)
    {
        empires.rnaseq.GenomicRegion gr = new GenomicRegion(c.toString(0), GenomicUtils.getStrand(c.toString(6)), c.toInt(3), c.toInt(4) + 1);
        HashMap<String, String>  attribs = readAttribs(c.toString(8), line);
        pf.parse(getOrAdd(genes, gr, attribs), gr,  attribs, line);

    }


    public GFFBasedIsoformRegionGetter(File f, String srcname, String version)
    {


        Iterator<String>  it = filterIterator(FileUtils.getLineIterator(f, false), (_l) -> _l.length() > 0  &&  _l.charAt(0) != '#');
        FileUtils.StringSplitter TSP = FileUtils.getRecordSplitter("\t");
        FileUtils.StringSplitter WSP = FileUtils.getRecordSplitter(" ");
        FileUtils.Converter c = new FileUtils.Converter();

        HashMap<String, GeneInfo> genes = new HashMap<>();

        int ln = 0;
        String line = null;
        while(it.hasNext())
        {
            ln++;
            c.set(TSP.split(line = it.next()));
            //c.set((line = it.next()).split("\t"));
            if(ln % 100 == 0) {
                //System.out.printf("read %d lines from %s\n", ln, f.getAbsolutePath());
            }
            if(c.length() < 7)
                throw new FRuntimeException("invalid line (too few fields (%d)) : %s", c.length(), line);

            String type = c.toString(2);

            switch(type)
            {
                case "exon":
                    parseLine(genes, c, line, (_gi, _gr, _a, _l) -> _gi.addExon(_gr, _a));
                    continue;

                case "CDS":
                    parseLine(genes, c, line, (_gi, _gr, _a, _l) -> _gi.addCDSExon(_gr, _a, _l));
                    continue;
                /*case "stop_codon":
                    genes.get(attribs.get("gene_id")).addStopCodon(gr, attribs);
                    continue;
                */
                default:
                    continue;
            }
            //System.out.printf("%s: %s %s\n", type, gr, attribs);
        }
        apply(genes.values(), (_g) ->
                {

                    MultiIsoformRegion mir = _g.toMIR();
                    regions.put(_g.id, mir);
                    process(mir);

                    //infos.add(Tuple3.create(_g.gr.getChr(), _g.gr.getStrand(), new Region1D(_g.gr.getX1(), _g.gr.getX2(), _g.id)));
                }
        );
        //System.out.printf("read %d genes\n", infos.size());

    }



}
