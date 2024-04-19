package nlEmpiRe.rnaseq;


import lmu.utils.*;
import lmu.utils.tuple.Tuple3;

import java.io.File;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Vector;

import static lmu.utils.ObjectGetter.filterOne;

public class FileBasedIsoformRegionGetter extends IsoformRegionGetter
{
    String srcname;
    String version;


    Vector<Tuple3<String,Boolean,Region1D<String>>> infos = new Vector<>();
    HashMap<String, MultiIsoformRegion> regions = new HashMap<>();

    HashMap<String, String> isoform2multi = new HashMap<>();
    HashMap<String, String> coding2multi = new HashMap<>();

    protected  FileBasedIsoformRegionGetter()
    {

    }

    void process(MultiIsoformRegion r)
    {
        ObjectGetter.apply(r.isoforms.keySet(), (_k) -> isoform2multi.put(_k, r.id));
        ObjectGetter.apply(r.coding.keySet(), (_k) -> coding2multi.put(_k, r.id));
        infos.add(Tuple3.create(r.chr, r.strand, new Region1D(r.start, r.end, r.id)));

    }

    public FileBasedIsoformRegionGetter(File f, String srcname, String version)
    {

        this.srcname = (srcname != null) ? srcname : f.getName();
        this.version = (version != null) ? version : "unknown";

        Vector<String> headers = FileUtils.getHeaders(f);
        boolean gotname = null != filterOne(headers, (_h) -> _h.equals("name"));

        FileUtils.HeaderedReader hr = FileUtils.getHeaderedReader("\t").add("regionid", "id").add("regiontype", "type")
                .add("chr", "strand", "start", "end","rvid","rvtype", "regvec");

        if(gotname) {
            hr.add("name");
        }

        Iterator<String> it = FileUtils.getLineIterator(f);
        hr.init(it.next());
        MultiIsoformRegion current = null;
        while (it.hasNext())
        {
            hr.setRecord(it.next());
            if (current == null || !current.id.equals(hr.getString("id")))
            {

                current = new MultiIsoformRegion();
                current.id = hr.getString("id");
                current.name = (!gotname) ? current.id : hr.getString("name");
                current.biotype = hr.getString("type");
                if(current.biotype.startsWith(MultiIsoformRegionWriter.PROTEIN_RV_TYPE))
                {
                    current.biotype = "protein_coding";

                }

                current.chr = hr.getString("chr");
                current.strand = GenomicUtils.getStrand(hr.getString("strand"));
                current.start = hr.getInt("start");
                current.end = hr.getInt("end");
                current.src = srcname;
                current.src_version = version;


                regions.put(current.id, current);

                infos.add(Tuple3.create(current.chr, current.strand, new Region1D(current.start, current.end, current.id)));
            }
            RegionVector rv = RegionVector.parseSimpleRepresentation(hr.getString("regvec"));
            String rvid = hr.getString("rvid");
            String type = hr.getString("rvtype");
            if(type.startsWith(MultiIsoformRegionWriter.PROTEIN_RV_TYPE))
            {
                String[] info = type.split("::");
                rv.setObject(rvid);

                current.coding.put(info[1], rv);
                current.coding2iso.put(info[1], info[2]);
                current.iso2coding.put(info[2], info[1]);
                continue;

            }


            current.addRegionVector(rvid, hr.getString("rvtype"), rv);

        }
        log.info("read %d genes from %s", regions.size(), f.getAbsolutePath());
    }

    public Pair<MultiIsoformRegion, RegionVector> getIsoform(String id)
    {
        String g = isoform2multi.get(id);
        if(g == null)
            throw new FRuntimeException("unknown isoform id: %s", id);

        MultiIsoformRegion mir = getRegionById(g);
        return Pair.create(mir, mir.isoforms.get(id));
    }

    public Pair<MultiIsoformRegion, RegionVector> getCoding(String id)
    {
        String g = coding2multi.get(id);
        if(g == null)
            throw new FRuntimeException("unknown coding sequence id: %s", id);

        MultiIsoformRegion mir = getRegionById(g);
        return Pair.create(mir, mir.coding.get(id));

    }


    public HashMap<String, MultiIsoformRegion> getRegionMap()
    {
        return regions;
    }

    public MultiIsoformRegion getRegionById(String id)
    {
        return regions.get(id);
    }

    public Iterator<MultiIsoformRegion> getRegions()
    {
        return regions.values().iterator();
    }
    /** chr, strand, region<id> */
    public Vector<Tuple3<String,Boolean,Region1D<String>>> getRegionInfos()
    {
        return infos;
    }


}

