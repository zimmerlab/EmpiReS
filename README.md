---
title: "Vignette: EmpiReS - Differential Analysis of Gene Expression and
Alternative Splicing"
---

EmpiReS is a method for differential analysis of RNA-Seq gene expression
data and for the detection of differential alternative splicing.
Furthermore, in the paper a novel benchmarking approach is presented
that does not assume a given statistical model but instead uses measured
data to simulate differential expression or differential alternative
splicing. The software package of EmpiReS provides a pipeline that
allows all these analysis - to simulate benchmarking data, to analysis
differential splicing and to analysis differential alternative splicing.
This vignette will show how to use the package for each of these tasks.

# System Requirements

All the examples given in this vignette work with the data (and jar
file) included in the EmpiReS tarball available from
<https://www.bio.ifi.lmu.de/files/download/EmpiReS/EmpiReS.tar.gz>. All
paths are relative to the directory where the tarball is extracted to.

You will need the samtools suite (version 1.10) and a java runtime
environment ($>=$ version 1.14) installed.

Before you start this vignette, you should call the `fillup.sh` script
included in the tarball, as it will generate some additional inputs
needed for the example commands in this vignette.

# Analysis of Differential Alternative Splicing {#DAS}

There are two possible inputs for EmpiReS - mapped reads as `BAM` files
or unmapped reads from `fastq` files. EmpiReS uses equivalence classes
as features and thus needs to derive counts for equivalence classes. An
equivalence class is a set of transcripts with which a set of reads is
compatible and they are, thus, constructed by deriving for each read the
set of transcripts that are compatible with it. When unmapped reads are
given the equivalence class based variant of contextmap (see EmpiReS
paper) is used to map the reads which will directly output read counts
for equivalence classes. When mapped reads from `BAM` files are used the
first step is to derive read counts for equivalence classes, so only the
first step of the pipeline differs depending on the available input and
all subsequent steps are the same.

## Derive Equivalence Class Counts

### EC-contextmap using Unmapped Reads (fastq files) {#ec-contextmap}

The EC-contextmap mapper is based on suffix arrays and LCP tables which
have to build for the transcriptome before the mapping can be performed.
To build this index the `build_index` task of the EmpiReS jar has to be
called. It needs as inputs:

-   the genomic annotation (`-gtf`) in GTF file format

-   the genome fasta file (`-genome`)

-   the index of the genome fasta file (`-genomeidx`) which can be build
    from the genome fasta file with the samtools faidx command

This task will write the index files for the given transcriptome to the
file given by the `-o` option. It needs only to be called once per
transcriptome.

To keep the tarball small the genome fasta file is not included in the
tarball. It can be downloaded from
<ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.toplevel.fa.gz>,
e.g. using wget, as is done during the `fillup.sh` script:

    wget  -O EXAMPLES/HUMAN_GENOME_GRCh37.75/Homo_sapiens.GRCh37.75.dna.toplevel.fa.gz \
    ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.toplevel.fa.gz
    gunzip EXAMPLES/HUMAN_GENOME_GRCh37.75/Homo_sapiens.GRCh37.75.dna.toplevel.fa.gz
    samtools faidx EXAMPLES/HUMAN_GENOME_GRCh37.75/Homo_sapiens.GRCh37.75.dna.toplevel.fa

This also unpacks the gzipped fasta file and builds the fasta index
using the samtools faidx command. To finally build the index for the
human transcriptome run:

    java -Xmx10G -jar empires.jar build_index \
        -gtf EXAMPLES/Homo_sapiens.GRCh37.75.gtf \
        -genome EXAMPLES/HUMAN_GENOME_GRCh37.75/Homo_sapiens.GRCh37.75.dna.toplevel.fa \
        -genomeidx EXAMPLES/HUMAN_GENOME_GRCh37.75/Homo_sapiens.GRCh37.75.dna.toplevel.fa.fai \
        -o EXAMPLES/mapping_reads_to_ECMs/human.GRCh37.65.ecm.ref

This will create the file
`EXAMPLES/mapping\_reads\_to\_ECMs/human.GRCh37.65.ecm.ref` which is a
binary index file that is used by the EC-contextmap mapper.

The mapping using the EC-contextmap mapper takes advantage of the
context of all reads, therefore all samples are mapped together. The
inputs of the `ecmapper` task are:

-   index file (`-index`) as built in the previous step

-   A table which defines which fastq files should be mapped and gives
    some additional information about the fastq files that are needed
    for the mapping (`-table`). The table contains one line for each
    sample and it has to contain at least the label of the sample that
    was measured (column \"**label**\"), its condition (column
    \"**condition**\"), paths to the fw and rw read fastq files (these
    can be gzipped or bzipped) (columns \"**fw**\" and \"**rw**\") and
    the strandness information (column \"**strandness**\". The
    **condition** column gives the condition of the sample, replicate
    measurements of the same condition will thus differ in their
    **label** but share the same value in the **condition** column. The
    **fw** and **rw** columns give the paths to the (zipped) fastq
    files. If no paired-end sequencing data but only single-read data is
    available, the **rw** column can be omitted. Finally, the
    **strandness** column contains information about the strandness of
    the measurement - a value of \"true\" indicates that the first read
    defines the strandness, while a value of \"false\" indicates that
    the second read defines the strandness and an empty string indicates
    that the measurement was strand unspecific.

        label    condition    fw                  rw                     strandness
        cond1_00    cond1       cond1_00_1.fastq.gz cond1_00_2.fastq.gz true
        cond1_01    cond1       cond1_01_1.fastq.gz cond1_01_2.fastq.gz true
        cond1_02    cond1       cond1_02_1.fastq.gz cond1_02_2.fastq.gz true
        cond2_00    cond2       cond2_00_1.fastq.gz cond2_00_2.fastq.gz true
        cond2_01    cond2       cond2_01_1.fastq.gz cond2_01_2.fastq.gz true
        cond2_02    cond2       cond2_02_1.fastq.gz cond2_02_2.fastq.gz true

-   optionally a base directory (`-basedir`) can be provided - if given
    all paths in the sample table are interpreted as relative to this
    directory

-   number of parallel threads (`-nthreads`) used for the mapping,
    default: 10.

Note that the fastq files referenced in this file are created by the
`fillup.sh` script, so you have to call it first to continue with the
vignette. The paths of the fastq files are relative to the directory\
`EXAMPLES/simulate_reads/TEST_OUTPUT/` which can be passed to EmpiReS
using the `-basedir` option. The counts for equivalence classes
resulting from the mapping of EC-contextmap are written to the path
provided by the option `-o`:

    java -Xmx10G -jar empires.jar ecmapper \
        -index EXAMPLES/mapping_reads_to_ECMs/human.GRCh37.65.ecm.ref \
        -table EXAMPLES/simulate_reads/TEST_OUTPUT/sample.table \
        -basedir EXAMPLES/simulate_reads/TEST_OUTPUT/ \
        -o EXAMPLES/stem_ecm_mapping_ecm.counts

The output file `EXAMPLES/stem_ecm_mapping_ecm.counts` directly contains
the counts for the equivalence classes:

    >ENSG00000215203    TOTAL
    reads   267.000 238.000 171.000 248.000 176.000 199.000
    >ENSG00000215203    ENST00000399770
    reads   267.000 238.000 171.000 248.000 176.000 199.000
    >ENSG00000180015    TOTAL
    reads   0.000   0.000   0.000   1.000   0.000   0.000
    >ENSG00000180015    ENST00000448646,ENST00000504731
    reads   0.000   0.000   0.000   1.000   0.000   0.000
    >ENSG00000179029    TOTAL
    reads   11.000  6.000   3.000   5.000   4.000   4.000

It contains for each gene an entry for the number of total reads and for
each equivalence class (identified by the unique set of transcripts in
the header). Each entry contains the number of reads for each sample
that was mapped in the order the samples are given in the sample table
file.

### Derive Equivalence Class Counts from Mapped Reads (BAM files) {#ecFromBAM}

If the reads were already mapped by another mapper, EmpiReS also allows
to start from BAM files. In this case the first step is to derive
equivalence classes from the annotation and BAM files and count the
reads for them using the `ec_from_bams` task of the jar file. The task
takes the following inputs:

-   a sample table (`-table`) with at least three columns: **label**
    (sample id), **bam** (path to mapping result of the sample (BAM
    file) and the **strandness** information. The table contains one
    line for each measured sample. The **strandness** column contains
    information about the strandness of the measurement - a value of
    \"true\" indicates that the first read defines the strandness, while
    a value of \"false\" indicates that the second read defines the
    strandness and an empty string indicates that the measurement was
    strand unspecific.

        label      bam                  strandness
        cond1_00    cond1_00_sorted.bam true
        cond1_01    cond1_01_sorted.bam true
        cond1_02    cond1_02_sorted.bam true
        cond2_00    cond2_00_sorted.bam true
        cond2_01    cond2_01_sorted.bam true
        cond2_02    cond2_02_sorted.bam true

-   genomic annotation (`-gtf`) in GTF format

-   optionally a base directory where all paths of the sample table are
    relative to (`-basedir`)

and writes the EC-contextmap counts to the path provided by (`-o`).

Note that the BAM files referenced in this file are created by the
`fillup.sh` script, so you have to call it first to continue with the
vignette. The paths in the sample table are relative to the directory
`EXAMPLES/simulate_reads/TEST_OUTPUT` which should be passed to the task
using the `-basedir` option.

    java -jar empires.jar ec_from_bams \
        -table EXAMPLES/simulate_reads/TEST_OUTPUT/sample.table \
        -gtf EXAMPLES/Homo_sapiens.GRCh37.75.gtf \
        -o EXAMPLES/stem_idealmapping_ecm.counts \
        -basedir EXAMPLES/simulate_reads/TEST_OUTPUT/

The output file `EXAMPLES/stem_idealmapping_ecm.counts` directly
contains the counts for the equivalence classes:

    >ENSG00000156026    TOTAL
    reads   19.00   22.00   21.00   13.00   21.00   18.00
    >ENSG00000156026    ENST00000357157,ENST00000373053,ENST00000536019,ENST00000604372,ENST00000604679,ENST00000605597
    reads   1.00    0.00    2.00    0.00    1.00    3.00
    >ENSG00000156026    ENST00000357157,ENST00000373053,ENST00000536019,ENST00000604152,ENST00000604372,ENST00000604679,ENST00000605597
    reads   1.00    0.00    1.00    0.00    0.00    1.00
    >ENSG00000156026    ENST00000373053,ENST00000536019,ENST00000604152,ENST00000604372,ENST00000604679,ENST00000605597
    reads   1.00    0.00    2.00    0.00    2.00    1.00
    >ENSG00000156026    ENST00000373053,ENST00000483185,ENST00000536019,ENST00000604152,ENST00000604372,ENST00000604679,ENST00000605597
    reads   0.00    0.00    1.00    0.00    0.00    0.00

It contains for each gene an entry for the number of total reads and for
each equivalence class (identified by the unique set of transcripts in
the header). Each entry contains the number of reads for each sample
that was mapped in the order the samples are given in the sample table
file.

## Analysis of Differential Alternative Splicing {#DE_DAS}

Finally, to identify differential alternative splicing events from the
equivalence class counts, the `diffexp_diffsplic_on_eccounts` task from
the EmpiReS jar has to be called. It needs the following inputs:

-   the equivalence class counts (`-i`) derived from BAM files (see
    section [2.1.2](#ecFromBAM){reference-type="ref"
    reference="ecFromBAM"}) or by the EC-contextmap mapper (see section
    [2.1.1](#ec-contextmap){reference-type="ref"
    reference="ec-contextmap"})

-   a sample table (`-samples`) that contains for each measurement a
    line containing the **label** of the sample and the **condition**.
    Replicate measurements should have the same value in the
    **condition** column. In our example there are two conditions with 3
    replicates each:

        label      condition
        cond1_00    cond1
        cond1_01    cond1
        cond1_02    cond1
        cond2_00    cond2
        cond2_01    cond2
        cond2_02    cond2

-   optionally the two conditions (`-cond1` and `-cond2`) that should be
    compared. If these options are not given the first two conditions
    from the samples table are used.

and writes the results of the differential expression and differential
alternative splicing analysis to the file given by the `-o` option.

The following command will analyze both differential expression and
differential alternative splicing for the given conditions:

    java -jar empires.jar diffexp_diffsplic_on_eccounts \
       -i EXAMPLES/stem_idealmapping_ecm.counts \
       -samples EXAMPLES/simulate_reads/TEST_OUTPUT/sample.table \
       -cond1 cond1 \
       -cond2 cond2 \
       -o EXAMPLES/empires_outtable_stem_simulation_on_ideal_mapping.tsv

and write the output to
`EXAMPLES/empires_outtable_stem_simulation_on_ideal_mapping.tsv`:

    gene    diffexp.fdr diffexp.log2fc  diffsplic.most.signif.test  diffsplic.fdr   diffsplic.difflog2fc
    ENSG00000013810 0,000e+00   1,650   ENSG00000013810.excl.ENST00000313288_VS_excl.ENST00000485989    1,739e-09   1,411
    ENSG00000131051 0,000e+00   -1,850  ENSG00000131051.merged.ENST00000444878.ENST00000461283_VS_excl.ENST00000461283  0,102   0,797
    ENSG00000181392 7,982e-12   -1,340  ENSG00000181392.merged.ENST00000324444.ENST00000503121_VS_excl.ENST00000324444  0,000e+00   -1,248
    ENSG00000131037 0,000e+00   1,690   ENSG00000131037.excl.ENST00000540810_VS_excl.ENST00000586329    5,368e-05   -0,798
    ENSG00000156313 2,197e-10   -1,020  ENSG00000156313.excl.ENST00000309513_VS_excl.ENST00000494707    0,000e+00   1,610
    ENSG00000132330 0,000e+00   -2,230  ENSG00000132330.excl.ENST00000254663_VS_excl.ENST00000412508    1,367e-09   1,557
    ENSG00000120334 0,000e+00   1,360   ENSG00000120334.merged.ENST00000460816.ENST00000484920_VS_excl.ENST00000484920  0,336   -0,475
    ENSG00000037897 0,000e+00   -1,580  ENSG00000037897.merged.ENST00000257848.ENST00000324871_VS_excl.ENST00000324871  0,008   0,432
    ENSG00000145604 0,000e+00   -1,630  ENSG00000145604.merged.ENST00000508514.ENST00000513263_VS_excl.ENST00000513263  0,200   -0,943

This table contains for each gene first two columns (**diffexp.fdr** and
**diffexp.log2fc**) containing the differential expression results, that
is the FDR corrected p-value and the log2 fold change of the gene.
Furthermore, there are columns for the differential alternative splicing
results: the alternative splicing event that yielded the lowest FDR
corrected p-value (**diffsplic.most.signif.test**), the FDR corrected
p-value (**diffsplic.fdr**) and the difference between the fold changes
of the two involved equivalence classes (fold change of fold changes,
**diffsplic.difflog2fc**).

# Analysis of Differential Expression

The steps described in section [2](#DAS){reference-type="ref"
reference="DAS"} will also analyze differential expression, so when
fastq or BAM files are available just follow the steps described in that
section. Additionally, EmpiReS allows to analyze differential expression
from a files that can easily be exported from an expression set in R.
For this the task `diffexp_on_eset` is used that takes the following
inputs:

-   input directory (`-inputdir`) containing three files: `exprs.txt`
    containing the feature count matrix with features as rows and
    samples as columns, `f_data.txt` containing the feature data that is
    the feature names and `p_data.txt` containing the phenotype data
    with the labels of the samples and their grouping into conditions.

    The file `exprs.txt` contains a matrix of counts without any
    headers. The rows of this matrix correspond to the features, while
    the columns correspond to the measured samples.

        239 234 152 211 202 161
        11   6    3   5     4   4
        261 277 181 290 182 213
        8     15     6    11  7   8
        328 301 258 374 279 287
        65   53  46  64  41  47
        176 192 153 186 123 123
        57   59  69  69  67  71
        108 110 90   87  65  77
        11   18  11  16  12  12

    The file `f_data.txt` contains the names of the features that are
    contained in the matrix contained in `exprs.txt`. It has the same
    number of lines as `exprs.txt`. The ordering of the rows in the
    matrix and the feature names in `f_data.txt` is the same, so that
    e.g. the first row in the matrix contains the read counts of the
    feature given in the first line of `f_data.txt`.

        ENSG00000215203
        ENSG00000179029
        ENSG00000205882
        ENSG00000131089
        ENSG00000230837
        ENSG00000228523
        ENSG00000167046
        ENSG00000168374
        ENSG00000228527
        ENSG00000025800

    The file `p_data.txt` contains a table that gives the phenotype
    information of the given data, that is how the samples that were
    measured are labeled (first column) and how they group into
    conditions (second column). Its number of rows is equal to the
    number of columns in `exprs.txt`. Again, the order of the columns in
    the matrix contained in `exprs.txt` and the rows in `p_data.txt` is
    the same, so that e.g. the first column of the matrix corresponds to
    the sample described in the first row of `p_data.txt`. The second
    column in the phenotype table contains the condition that the sample
    belongs to. Replicate measurements of the same condition share the
    same value in this column.

        C1R1    0
        C1R2    0
        C1R3    0
        C2R1    1
        C2R2    1
        C2R3    1

-   optionally the two conditions (`-cond1` and `-cond2`) that should be
    compared. If these options are not given the first two conditions
    from the `p_data.txt` table are used.

and writes the results to the output file given by the `-o` option.

When the input files are prepared in the directory
`EXAMPLES/expression_set_format/` the `diffexp_on_eset` task of the
EmpiReS jar can be called:

    java -jar empires.jar diffexp_on_eset \
       -inputdir EXAMPLES/expression_set_format/ \
       -cond1 0 \
       -cond2 1 \
       -o EXAMPLES/DE_from_eset.tsv

The output table `EXAMPLES/DE_from_eset.tsv` contains for each gene a
row that contains the log2 fold change and the FDR corrected p-value:

    gene            log2FC  fdr
    ENSG00000164818 -1,040  0,000e+00
    ENSG00000180071 -1,380  0,000e+00
    ENSG00000023228 -1,660  0,000e+00
    ENSG00000120334 1,350   0,000e+00
    ENSG00000157600 0,970   0,000e+00
    ENSG00000049883 -1,270  0,000e+00
    ENSG00000132305 -1,130  0,000e+00
    ENSG00000006625 1,480   0,000e+00
    ENSG00000037897 0,950   0,000e+00

# Benchmarking of Differential Expression and Differential Alternative Splicing

## Simulate Counts {#simCounts}

To simulate benchmarking datasets EmpiReS uses data with many replicates
which are grouped into two simulated conditions. To introduce changed
genes/transcripts the measurements of two genes/transcripts are swapped
in one of the two simulated conditions to introduce two fold changes,
without changing the overall distribution and variance between the
replicates within the conditions.

The task `simulate_counts` of the EmpiReS jar simulates counts both for
the benchmarking of differential expression as well as for differential
alternative splicing. As inputs it needs:

-   measured gene read counts with many replicates that is used as input
    data for the simulation (`-incounts`). The file should contain one
    column containing the gene id and at least three for replicates for
    both conditions, that is in total at least 7 columns. All the
    columns should contain counts of replicate measurement of the same
    real experimental condition. The column names given in the header
    define which columns should be used for which simulated condition.
    The columns should be named `<condition>.<replicate>`. The
    assignment of real replicates to to-be-simulated conditions can be
    random, or as described in the paper two technical replicates can be
    distributed to both conditions.

        gene                cond1.rep1  cond1.rep2  cond1.rep3  cond2.rep1  cond2.rep2  cond2.rep3
        ENSG00000000003  527           568         674             718         685          564
        ENSG00000000005  12         10            11                16            17              11
        ENSG00000000419  232           230             346             312          334            258
        ENSG00000000457  66         57              70              66              83              76

-   list of transcripts that should be simulated
    (`-transcriptsToSimulate`, both differential and unchanged
    transcripts). Usually, pairs of transcripts within genes inducing
    some alternative splicing. In the paper (and this example)
    transcripts inducing the longest exon skipping event for a gene are
    selected. The mapping of which transcripts belong to which gene is
    given by the GTF file.

        ENST00000589042
        ENST00000359218
        ENST00000375735
        ENST00000334267
        ENST00000424728

-   the genomic annotation (`-gtf`) in GTF file format providing the
    transcript to gene mapping

-   selected genes for differential expression (`-diffexp`) (but not
    differential splicing) is to be simulated. If the list given by
    `-transcriptsToSimulate` contains multiple transcripts for a gene in
    this file, all those transcripts will be simulated with about the
    same fold change.

        ENSG00000183354
        ENSG00000198843
        ENSG00000078319
        ENSG00000135632
        ENSG00000125834

-   selected gene ids where differential alternative splicing is to be
    simulated (`-diffsplic`). For gene ids in this file exactly two
    transcripts have to be placed in the transcripts.to.simulate file.
    The simulation will then generate different expression changes for
    these two transcripts as described in the paper.

        ENSG00000155657
        ENSG00000187240
        ENSG00000188738
        ENSG00000151150
        ENSG00000167522

The simulation will write the following files to the provided output
directory path (`-od`):

-   `transcript_exprs.txt`: the simulated number of reads for each
    simulated transcript in each sample

        C1R1    C1R2    C1R3    C2R1    C2R2    C2R3
        ENST00000589042 863 844 590 851  560     591
        ENST00000359218 45   41  59  131     111      95
        ENST00000375735 778 798 554 1793    1456    1461
        ENST00000334267 121 115 98   136     83   97

-   `simul.info`: detailed information about the simulated changes -
    there is line for each transcript that was simulated with a change
    of expression or splicing. Each line contains the gene id,
    transcript id, whether it was selected as major or minor transcript,
    the fold change implicated by the label swapping procedure, which
    gene was selected as swapping partner and whether the gene was
    simulated for differential expression or differential splicing

        ENSG00000186376 ENST00000370766 major   1,94    ENSG00000100353 DIFFSPLIC
        ENSG00000005884 ENST00000544892 major   -1,94   ENSG00000104332 DIFFSPLIC
        ENSG00000186376 ENST00000370764 minor   3,38    ENSG00000207165 DIFFSPLIC
        ENSG00000005884 ENST00000512553 minor   -3,38   ENSG00000106484 DIFFSPLIC
        ENSG00000124422 ENST00000261497 major   -1,61   ENSG00000013810 DIFFSPLIC

-   `diffexp.trues` - list of gene ids which were simulated as
    differentially expressed

        ENSG00000213025
        ENSG00000261586
        ENSG00000160888
        ENSG00000163915
        ENSG00000009307

-   `diffsplic.trues` - list of gene ids which were simulated as
    differentially alternatively spliced

        ENSG00000214357
        ENSG00000160404
        ENSG00000110148
        ENSG00000092871
        ENSG00000103855

The simulation of the counts can be called by:

    java -jar empires.jar simulate_counts \
       -transcriptsToSimulate EXAMPLES/simulate_transcript_counts/transcripts.to.simulate \
       -incounts EXAMPLES/simulate_transcript_counts/simul.stem \
       -diffexp EXAMPLES/simulate_transcript_counts/input.truediffexp \
       -diffsplic EXAMPLES/simulate_transcript_counts/input.truesplic \
       -gtf EXAMPLES/Homo_sapiens.GRCh37.75.gtf \
       -od EXAMPLES/simulate_transcript_counts/TEST_OUTPUT/

## Simulate Reads For Counts {#simReads}

Given the simulated counts for transcripts we need to simulate reads for
these counts. The task `generate_reads` generates read sequences (as
gzip compressed fastq files) along with their simulated origins as ideal
mappings (as BAM files).

It takes the following inputs:

-   the genomic annotation (`-gtf`) in GTF format

-   the genome fasta file and it's idx (`-genome` and `-genomeidx`)

-   a file containing per transcript counts for all samples to simulate
    (`-trcounts`) as provided by the output file `transcript_exprs.txt`
    of the previous step.

The task writes fastq.gz2 files, the bamfiles (unless the switch
`-nobams` is used) to the output directory provided by (`-od`).
Additionally, an overview sample table (`sample.table`) is written -
this can be then used for the subsequent evaluation steps.

    label      condition    bam                 fw                   rw                  strandness
    cond1_00    cond1       cond1_00_sorted.bam cond1_00_1.fastq.gz cond1_00_2.fastq.gz true
    cond1_01    cond1       cond1_01_sorted.bam cond1_01_1.fastq.gz cond1_01_2.fastq.gz true
    cond1_02    cond1       cond1_02_sorted.bam cond1_02_1.fastq.gz cond1_02_2.fastq.gz true
    cond2_00    cond2       cond2_00_sorted.bam cond2_00_1.fastq.gz cond2_00_2.fastq.gz true
    cond2_01    cond2       cond2_01_sorted.bam cond2_01_1.fastq.gz cond2_01_2.fastq.gz true
    cond2_02    cond2       cond2_02_sorted.bam cond2_02_1.fastq.gz cond2_02_2.fastq.gz true

There are several optional parameters to control the properties of the
reads and fragments by changing the mutation rate (`-mutrate`), read
length (`-readlength`), or the distribution of the fragment lengths
(mean: `-fraglengthmean` and standard deviation `-fraglengthsd`.
Moreover, if the optional parameter `-biaspos` is set, the simulation
will use the position specific bias derived from the . An example input
file for such (derived from yeast experiments) is provided in
`EXAMPLES/simulate_reads/yeast.tr2startfreq.out`. This table gives for
each (yeast) transcript its length and for each position the number of
mapped reads to this position:

    YFL014W 331 0,0,1,0,0,0,0,0,1,0,0,2,1,1,2,3,2,17,19,2,3,0,0,1,2,2,12,14,7,8,12,14,51,48,15,10,6,2,4,0,2,3,6,3,4,0,1,3,5,8,18,17,24,69,71,55,68,17,11,11,0,2,6,8,7,11,15,47,45,8,2,6,0,61,3,2,26,25,6,2,1,4,0,1,1,0,9,7,29,11,8,11,3,4,2,3,0,2,11,31,46,92,32,7,1,8,0,1,0,0,4,7,3,8,6,5,3,0,0,2,15,19,114,396,634,816,362,121,265,457,459,396,147,137,254,141,90,55,14,2,2,3,9,27,68,94,142,64,23,104,173,549,999,836,8,5297,3510,3034,408,150,55,199,1040,1016,1071,304,54,30,69,46,47,65,48,105,0,35,19,88,11,20,26,87,176,802,1891,2571,922,481,693,652,267,133,53,62,68,160,109,247,156,125,120,37,36,54,97,637,2128,3494,4332,4822,3631,1901,518,400,615,6,4,1186,431,3159,207,201,100,157,378,158,202,461,625,571,1761,1962,2455,2240,682,123,99,92,19,22,17,35,122,149,110,71,44,29,13,15,4,6,17,32,134,572,890,1797,1865,1516,256,340,62,78,84,36,76,150,129,510,1906,4982,7608,2400,2936,4582,1309,1922,362,45,33,28,35,99,158,101,85,182,408,1119,0,0,0,4,1,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
    YPR180W 1045    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,1,0,2,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,5,7,18,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,7,2,0,0,0,1,0,0,3,1,6,2,0,2,1,0,1,0,6,1,0,0,0,3,5,16,41,91,43,3,13,5,3,1,0,1,3,0,0,0,0,0,1,1,1,0,0,1,0,0,2,0,0,0,0,0,0,0,0,2,0,0,0,0,1,1,2,1,2,1,3,6,3,4,13,19,0,15,5,29,1,0,3,2,0,2,2,0,4,6,17,14,17,12,9,8,17,9,10,6,2,8,9,25,94,25,20,56,19,58,0,0,0,0,3,2,3,3,1,7,2,14,14,6,1,3,3,4,5,10,21,11,11,8,7,7,10,4,9,11,2,3,0,0,0,0,1,0,0,0,1,1,0,0,0,0,0,0,1,0,0,0,0,0,1,3,3,0,2,1,0,1,0,0,1,1,0,1,0,6,5,6,4,2,2,8,19,73,18,5,4,13,9,9,6,20,43,20,7,18,9,2,0,4,8,6,5,14,19,40,50,95,60,27,27,17,6,3,15,24,56,133,162,324,316,203,147,64,51,81,128,280,256,104,90,45,9,8,14,42,80,39,15,14,65,139,172,105,47,35,10,1,0,0,0,0,1,1,3,2,6,2,0,1,1,0,0,0,3,0,0,0,0,1,1,0,0,0,0,0,4,3,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,1,3,0,0,1,5,8,2,2,7,3,4,2,9,3,3,5,1,1,1,2,4,12,25,44,62,40,11,4,8,8,4,5,3,5,7,17,35,88,42,68,95,342,274,102,66,9,0,0,0,1,1,7,2,5,1,2,3,2,0,0,0,1,1,1,0,0,0,1,4,11,18,54,59,42,5,3,2,3,3,2,2,6,18,45,28,37,58,26,3,3,4,1,1,0,0,2,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,2,0,1,0,0,1,2,1,0,7,8,2,4,1,2,1,1,8,8,10,7,1,2,3,5,4,9,10,12,2,0,0,0,7,32,51,238,110,29,18,10,43,13,15,38,20,1,2,2,1,0,1,0,1,3,1,1,1,0,2,1,2,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,4,9,26,21,2,1,0,5,6,0,5,6,2,0,0,0,0,0,0,1,15,24,53,80,52,15,8,5,3,8,8,16,35,21,32,45,46,49,17,17,10,11,19,52,117,20,17,14,21,12,2,1,4,0,0,1,2,3,0,1,0,0,0,0,1,3,17,16,4,6,13,1,4,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,1,5,0,0,0,1,3,2,8,5,2,4,1,9,1,2,0,1,5,7,6,1,3,6,5,2,3,14,2,0,0,0,0,2,1,4,1,0,1,1,2,3,5,6,3,1,6,2,7,17,12,15,13,24,10,6,5,11,14,2,2,1,0,0,1,0,0,2,1,0,7,11,11,30,35,5,30,36,22,118,97,88,115,236,0,0,55,73,179,9,14,44,47,83,74,44,14,8,8,19,5,2,4,3,6,17,51,15,7,5,6,1,2,1,5,3,14,21,22,23,12,20,4,4,6,8,5,9,8,8,3,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,3,1,2,0,0,0,0,1,2,0,2,0,1,0,0,1,1,4,1,3,1,0,0,0,0,0,0,0,0,0,1,0,5,10,10,10,4,15,8,4,0,0,0,0,0,0,0,0,0,1,1,2,2,2,0,2,0,0,0,0,0,0,0,0,0,0,0,0,5,2,1,9,8,17,9,2,2,0,0,10,15,35,27,38,64,25,4,3,0,0,1,0,1,0,1,2,1,2,1,0,1,5,10,7,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0

The `generate_reads` task can be called by:

    java -jar empires.jar generate_reads \
       -gtf EXAMPLES/Homo_sapiens.GRCh37.75.gtf \
       -genome EXAMPLES/HUMAN_GENOME_GRCh37.75/Homo_sapiens.GRCh37.75.dna.toplevel.fa \
       -genomeidx EXAMPLES/HUMAN_GENOME_GRCh37.75/Homo_sapiens.GRCh37.75.dna.toplevel.fa.fai \
       -trcounts EXAMPLES/simulate_transcript_counts/TEST_OUTPUT/transcript_exprs.txt \
       -biaspos EXAMPLES/simulate_reads/yeast.tr2startfreq.out \
       -od EXAMPLES/simulate_reads/TEST_OUTPUT

The output bamfiles from the read generation are not yet coordinate
sorted, which is needed for many subsequent steps, such as to derive
EC-contextmap counts from BAM files. The sorting can easily be done with
samtools: (<http://www.htslib.org/download/>) with the command samtools
sort:

    samtools=/home/proj/software/samtools/samtools-0.1.19/samtools
    for bam in `find EXAMPLES/simulate_reads/TEST_OUTPUT/ -name "*.bam"`
    do
            prefix=`echo $bam | awk -F".bam" '{print $1}'`
            $samtools sort $bam ${prefix}_sorted
    done

The sample table written by the read generation already refers to the
paths `${prefix}_sorted.bam` generated by this snippet.

## Evaluation of Differential Expression and Differential Alternative Splicing

The steps described in sections [4.1](#simCounts){reference-type="ref"
reference="simCounts"} and [4.2](#simReads){reference-type="ref"
reference="simReads"} already produce all inputs needed to analyze
differential expression and differential alternative splicing, both from
fastq files (see section [2.1.1](#ec-contextmap){reference-type="ref"
reference="ec-contextmap"}) or from the ideal mapping BAM files
(containing the mapping exactly as simulated, see section
[2.1.2](#ecFromBAM){reference-type="ref" reference="ecFromBAM"}). To
derive the counts for the equivalence classes one can directly use the
command:

    java -Xmx10G -jar empires.jar ecmapper \
        -index EXAMPLES/mapping_reads_to_ECMs/human.GRCh37.65.ecm.ref \
        -table EXAMPLES/simulate_reads/TEST_OUTPUT/sample.table \
        -basedir EXAMPLES/simulate_reads/TEST_OUTPUT/ \
        -o EXAMPLES/stem_ecm_mapping_ecm.counts

that uses the sample table written by the `generate_reads` task. Note
that the index for the equivalence class variant of contextmap needs to
be built first, for details on how this is done see section
[2.1.1](#ec-contextmap){reference-type="ref" reference="ec-contextmap"}.
If the evaluation should be based on another mapper or on the BAM files
that contain the simulated mapping, you can instead use this command:

    java -jar empires.jar ec_from_bams \
        -table EXAMPLES/simulate_reads/TEST_OUTPUT/sample.table \
        -gtf EXAMPLES/Homo_sapiens.GRCh37.75.gtf \
        -o EXAMPLES/stem_idealmapping_ecm.counts \
        -basedir EXAMPLES/simulate_reads/TEST_OUTPUT/

This command again uses the same table created by the `generate_reads`
task which contains the paths to the ideal mapping BAM files. When the
results of other mappers should be used, the sample table together with
the `-basedir` option have to be adapted (for details see section
[2.1.2](#ecFromBAM){reference-type="ref" reference="ecFromBAM"}).

Both these commands write a count matrix
(`EXAMPLES/stem_ecm_mapping_ecm.counts` and\
`EXAMPLES/stem_idealmapping_ecm.counts`) which can be used by the task
`diffexp_diffsplic_on_eccounts` to analyze differential expression and
differential alternative splicing. If this task is additionally given
the files containing the simulated splicing events and differentially
expressed transcripts written by `simulate_counts` the evaluation
results will also be outputted to standard out and can be retrieved by a
simple grep command:

    java -jar empires.jar diffexp_diffsplic_on_eccounts \
       -i EXAMPLES/stem_idealmapping_ecm.counts \
       -samples EXAMPLES/simulate_reads/TEST_OUTPUT/sample.table \
       -cond1 cond1 \
       -cond2 cond2 \
       -o EXAMPLES/empires_outtable_stem_simulation_on_ideal_mapping.tsv \
       -truesplicing EXAMPLES/simulate_transcript_counts/TEST_OUTPUT/diffsplic.trues \
       -truediffexp EXAMPLES/simulate_transcript_counts/TEST_OUTPUT/diffexp.trues \
       | grep "PREC:" > EXAMPLES/evaluationResults.txt 

The grep-ed lines saved to `EXAMPLES/evaluationResults.txt` contain all
evaluation statistics such as precision (PREC), recall (RECALL), f1
measure (F1), area under the receiver-operator curve (AUROC), area under
the precision-recall curve (AUPR) and the precision at 80% recall (PREC
at 80,00% recall) for both differential expression and differential
alternative splicing.

    EmpiRe-diffexp 11218 pred trues: 765/789 (ntrues: 774) pred. E PREC: 96,958 RECALL: 98,837 
    F1: 97,889 AUROC: 99,989 AUPR: 99,858 PREC at 80,00% recall: 100,00
    EmpiRe-diffsplic 4843 pred trues: 891/1005 (ntrues: 1246) pred. E PREC: 88,657 RECALL: 71,509 
    F1: 79,165 AUROC: 92,662 AUPR: 88,228 PREC at 80,00% recall: 75,13