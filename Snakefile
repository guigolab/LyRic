import glob
from collections import defaultdict
import os
import itertools
import sys
import pandas as pd
from pprint import pprint
from snakemake.utils import min_version

min_version("8")


# configfile: "config.yaml"
containerized: "docker://ghcr.io/guigolab/lyric:latest"


onsuccess:
    print("LyRic workflow finished smoothly.")


onerror:
    print("LyRic workflow finished with errors, see log file.")


# prefix all shell scripts with the following commands:
shell.prefix('set -euo pipefail; export PATH="{workflow.basedir}/utils/:$PATH";')


# (it seems DRMAA does not source ~/.bashrc by itself)
include: "functions.py"


# consistent ggplot2 theme across R plots:
GGPLOT_PUB_QUALITY = "theme(axis.text= element_text(size=themeSize*1.8), axis.ticks = element_line(size=lineSize), axis.line = element_line(colour = '#595959', size=lineSize), axis.title=element_text(size = themeSize*2), panel.grid.major = element_line(colour='#d9d9d9', size=lineSize),panel.grid.minor = element_line(colour='#e6e6e6', size=minorLineSize),panel.border = element_blank(),panel.background = element_blank(), strip.background = element_rect(colour='#737373',fill='white'), legend.key.size=unit(0.5,'line'), legend.title=element_text(size=themeSize*1.2), legend.text=element_text(size=themeSize), strip.text = element_text(size = themeSize))"

# list of stat plot formats to generate (currently only one format at a time is supported):
plotFormat = [
    "png",
]
# color sizeFracs:
sizeFrac_Rpalette = "c('0+'='#b3b3b3', '0-1' ='#f765ac','1+' ='#b370f9')"

# generic long color palette for R:
long_Rpalette = "c('#8dd3c7','#ffffb3','#80ff80','#ff4dff','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd','#ccebc5','#ffed6f','#c4ff4d','#ff66a3')"

# color palette for seqPlatform's and libraryPrep's
sampleAnnot_Rpalette = "list(seqPlatform = c(ONT = '#8dd3c7', pacBioSI = '#ffffb3', pacBioSII = '#5d513a'), libraryPrep = c(CapTrap = '#ff4dff', SMARTer = '#bebada', SmartSeq2 ='#751aff', Teloprime = '#fb8072', directRNA = '#80b1d3', 'Rt' = '#8c8c8c', PcrOnt='#3333ff', R2C2='#ff9900', IsoSeq='#800080'))"

# color palette for gencode biotypes
simpleBiotypes_Rpalette = "c('lncRNA' = '#66ccff',  'nonExonic' = '#6666ff',  'protein_coding' = '#ff8c66', miRNA = '#808000', 'misc_RNA' = '#99ff99','pseudogene' = '#d98cb3', 'rRNA' = '#d9d9d9', 'ERCC' = '#8dd3c7', 'SIRV' = '#399384')"
simpleBiotypes_Rpalette_woSpikeins = "c('lncRNA' = '#66ccff',  'nonExonic' = '#6666ff',  'protein_coding' = '#ff8c66', miRNA = '#808000', 'misc_RNA' = '#99ff99','pseudogene' = '#d98cb3', 'rRNA' = '#d9d9d9')"


# which genome build corresponds to to each capture design:
CAPDESIGNTOGENOME = config["capDesignToGenome"]
# which gencode annotation GTF correspond to each capture design:
try:
    GENOMETOANNOTGTF = config["genomeToAnnotGtf"]
except KeyError:
    print("No genome annotation(s) file provided.", file=sys.stderr)

# non-overlapping targeted regions for each capture design

if config["CAPTURE"]:
    # mapping of capDesign to capDesign:
    CAPDESIGNTOCAPDESIGN = config["capDesignToCapDesign"]
# non-overlapping targeted regions for each capture design

try:
    CAPDESIGNTOTARGETSGFF = config["capDesignToTargetsGff"]
except KeyError:
    print("No primary targets file provided.", file=sys.stderr)


# which CAGE peak BED file corresponds to each genome build:
try:
    GENOMETOCAGEPEAKS = config["genomeToCAGEpeaks"]
except KeyError:
    print("No CAGE annotation file provided.", file=sys.stderr)
# which ENCODE DHS peak BED file corresponds to each genome build:
try:
    GENOMETODHSPEAKS = config["genomeToDHSpeaks"]
except KeyError:
    print("No DHS annotation file provided.", file=sys.stderr)

# URL for UCSC Track Hub data files:
if config["produceTrackHub"]:
    TRACK_HUB_DATA_URL = config["TRACK_HUB_BASE_URL"] + "dataFiles/"
# 5' and 3' end support categories for transcript models (TMs):
# "all" = all TMs (no filter applied)
# "cagePolyASupported" = CAGE (5') + PolyA (3') -supported TMs only
ENDSUPPORTcategories = ["all"]
if config.get("genomeToCAGEpeaks"):
    ENDSUPPORTcategories.append("cagePolyASupported")
# splicing status categories of TMs:
TMSPLICEDSTATUScategories = ["all", "spliced", "unspliced"]
# categories of plots wrt spike-ins:
SPIKEINcategories = ["wSpikeIns", "woSpikeIns"]
# "wSpikeIns": no filters (stats include spike-ins if present)
# "woSpikeIns": spike-ins (ERCC and SIRVs) removed from statistics

# TMPDIR to write temp files in
TMPDIR = "$TMPDIR"

MINIMUM_TMERGE_READ_SUPPORT = config["minimumTmergeReadSupport"]
READFILTERS = ["HiSS", "noFilt"]
# for tmerge:
ExonOverhangTolerance = config["exonOverhangTolerance"]

GENOMES = []
GENOMETOCAPDESIGNS = defaultdict(list)
for capD in CAPDESIGNTOGENOME:
    genome = CAPDESIGNTOGENOME[capD]
    GENOMES.append(genome)
    GENOMETOCAPDESIGNS[genome].append(capD)
GENOMES = set(GENOMES)

SIRVpresent = None
try:
    config["SIRVinfo"]
except KeyError:
    SIRVpresent = False
    print("No SIRVinfo file provided.", file=sys.stderr)
else:
    SIRVpresent = True


(TECHNAMES, CAPDESIGNS, SIZEFRACS, SAMPLEREPS) = glob_wildcards(
    "fastqs/" + "{techname, [^_/]+}_{capDesign}_{sizeFrac}_{sampleRep}.fastq.gz"
)
SAMPLEREPS = set(SAMPLEREPS)
CAPDESIGNS = set(CAPDESIGNS)
SIZEFRACS = set(SIZEFRACS)
TECHNAMES = set(TECHNAMES)


# no underscores allowed in wildcards, to avoid greedy matching since we use them as separators:
wildcard_constraints:
    capDesign="[^_/]+",
    sizeFrac="[^_/]+",
    techname="[^_/]+",
    sampleRep="[^_/]+",
    minReadSupport="|".join(MINIMUM_TMERGE_READ_SUPPORT),
    endSupport="|".join(ENDSUPPORTcategories),
    splicedStatus="|".join(TMSPLICEDSTATUScategories),
    filt="|".join(READFILTERS),
    sample_name="[^/]+",


## read sample metadata annotations into a pandas dataframe:
sampleAnnot = pd.read_table(config["SAMPLE_ANNOT"], header=0, sep="\t")

## check that samples are properly annotated:
(SAMPLES,) = glob_wildcards("fastqs/" + "{sample_name, [^/]+}.fastq.gz")
for sample_name in SAMPLES:
    sampleRow = sampleAnnot[sampleAnnot.sample_name == sample_name]
    assert len(sampleRow) < 2, (
        "Duplicate found for sample " + sample_name + " in " + config["SAMPLE_ANNOT"]
    )
    assert len(sampleRow) > 0, (
        "Sample " + sample_name + " not found in " + config["SAMPLE_ANNOT"]
    )

# drop rows of sampleAnnot whose sample_name are not in SAMPLES list:
sampleAnnot = sampleAnnot[sampleAnnot["sample_name"].isin(SAMPLES)]

##########################################################################
### Make mapping dict of sampleRep group ID (a sample, most of the times)
### to list of sampleReps belonging to this group.
### config['sampleRepGroupBy'] contains the list of attributes to group
### sampleReps on. These attributes are column names in the sample
### annotation file (config['SAMPLE_ANNOT'])
##########################################################################

# group samples contained in sampleAnnot by attributes listed in config['sampleRepGroupBy']
sampleRepGroups = sampleAnnot.groupby(config["sampleRepGroupBy"])
# initialize dict that will contain mapping of sampleRep group ID to list of sampleReps belonging to this group:
sampleRepGroupIdToSampleReps = {}
# populate sampleRepGroupIdToSampleReps:
for key, item in sampleRepGroups:
    # remove spaces from group identifiers, because they will serve as file basenames
    sampleRepGroupId = "_".join(key).replace(" ", "")
    # assign list of sample_names (dict value) to sampleRepGroupId (dict key)
    sampleRepGroupIdToSampleReps[sampleRepGroupId] = item["sample_name"].to_list()
################
### Done!
################

# convert sampleAnnot to dict to facilitate later access
sampleAnnot.set_index("sample_name", inplace=True)
sampleAnnotDict = sampleAnnot.to_dict("index")

# extract unique list of sub-projects in sampleAnnot (values in 'subProject' column) to generate separate, filtered HTML stats reports:
subProjects = set(sampleAnnot["subProject"])
# 'None' to create an HTML report without filter:
subProjects.add("ALL")


###################################################################################
### make list of authorized wildcard combinations
### inspired by
### https://stackoverflow.com/questions/41185567/how-to-use-expand-in-snakemake-when-some-particular-combinations-of-wildcards-ar
###################################################################################

AUTHORIZEDCOMBINATIONS = (
    []
)  # contains combinations of wildcards corresponding to existing FASTQs or "by*" combinations when relevant
TECHNAMESplusBY = set(TECHNAMES)
TECHNAMESplusBY.add("byTech")
SAMPLEREPSplusBY = set(SAMPLEREPS)
SAMPLEREPSplusBY.add("bySampleRep")
CAPDESIGNSplusBY = set(CAPDESIGNS)
CAPDESIGNSplusBY.add("byCapDesign")

# Populate AUTHORIZEDCOMBINATIONS:
for comb in itertools.product(
    TECHNAMESplusBY, CAPDESIGNSplusBY, SIZEFRACS, SAMPLEREPSplusBY
):
    if comb[0] == "byTech":
        filesList = glob.glob(
            "fastqs/"
            + "*"
            + "_"
            + comb[1]
            + "_"
            + comb[2]
            + "_"
            + comb[3]
            + ".fastq.gz"
        )
        if len(filesList) > 1:  # authorize "byTech" combination only if it's relevant
            authorizeComb(comb)

    elif comb[1] == "byCapDesign":
        filesList = glob.glob(
            "fastqs/"
            + comb[0]
            + "_"
            + "*"
            + "_"
            + comb[2]
            + "_"
            + comb[3]
            + ".fastq.gz"
        )
        if (
            len(filesList) > 1
        ):  # authorize "byCapDesign" combination only if it's relevant
            authorizeComb(comb)
    elif comb[3] == "bySampleRep":
        filesList = glob.glob(
            "fastqs/"
            + comb[0]
            + "_"
            + comb[1]
            + "_"
            + comb[2]
            + "_"
            + "*"
            + ".fastq.gz"
        )
        if (
            len(filesList) > 1
        ):  # authorize "bySampleRep" combination only if it's relevant
            authorizeComb(comb)
    elif os.path.isfile(
        "fastqs/" + comb[0] + "_" + comb[1] + "_" + comb[2] + ".fastq.gz"
    ) or os.path.isfile(
        "fastqs/"
        + comb[0]
        + "_"
        + comb[1]
        + "_"
        + comb[2]
        + "_"
        + comb[3]
        + ".fastq.gz"
    ):
        authorizeComb(comb)


AUTHORIZEDCOMBINATIONS = set(AUTHORIZEDCOMBINATIONS)


include: "fastqStats.smk"
include: "lrMapping.smk"
include: "srMapping.smk"
include: "polyAmapping.smk"
include: "introns.smk"
include: "processReadMappings.smk"
include: "tmClassification.smk"
include: "tmEndSupport.smk"
include: "trackHub.smk"
include: "htmlTable.smk"
include: "processInputGenome.smk"


###########################
### Transcriptome files ###
###########################
# transcriptome GTFs (per sampleRep):
gff = expand(
    "output/mappings/mergedReads/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.HiSS.tmerge.min{minReadSupport}reads.splicing_status-all.endSupport-all.gff.gz",
    filtered_product,
    techname=TECHNAMES,
    capDesign=CAPDESIGNS,
    sizeFrac=SIZEFRACS,
    sampleRep=SAMPLEREPS,
    minReadSupport=MINIMUM_TMERGE_READ_SUPPORT,
)

# transcriptome GTFs (per grouped sampleReps)
grouped_gff = expand(
    "output/mappings/mergedReads/groupedSampleReps/{groupedSampleRepBasename}.min{minReadSupport}reads.splicing_status-all.endSupport-all.gff.gz",
    groupedSampleRepBasename=sampleRepGroupIdToSampleReps.keys(),
    minReadSupport=MINIMUM_TMERGE_READ_SUPPORT,
)

# read-to-TM mapping file (per grouped sampleReps, required by LRGASP to check what each TM contains)
grouped_read_to_tm = expand(
    "output/mappings/mergedReads/groupedSampleReps/{groupedSampleRepBasename}.min{minReadSupport}reads.splicing_status-all.endSupport-all.readsToTm.tsv.gz",
    groupedSampleRepBasename=sampleRepGroupIdToSampleReps.keys(),
    minReadSupport=MINIMUM_TMERGE_READ_SUPPORT,
)

################################
### Summary statistics plots ###
################################
hiss_stats = expand(
    returnPlotFilenames(
        "output/plots/"
        + "HiSS.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.HiSS.stats"
    ),
    filtered_product,
    techname=TECHNAMESplusBY,
    capDesign=CAPDESIGNSplusBY,
    sizeFrac=SIZEFRACS,
    sampleRep=SAMPLEREPSplusBY,
)
read_length_stats = expand(
    returnPlotFilenames(
        "output/plots/"
        + "readLength.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.readLength.stats"
    ),
    filtered_product,
    techname=TECHNAMES,
    capDesign=CAPDESIGNS,
    sizeFrac=SIZEFRACS,
    sampleRep=SAMPLEREPS,
)  # facetted histograms of read length
fastq_qc = (
    expand(
        "output/fastqs/" + "qc/{techname}_{capDesign}_{sizeFrac}.{sampleRep}.dupl.txt",
        filtered_product,
        techname=TECHNAMES,
        capDesign=CAPDESIGNS,
        sizeFrac=SIZEFRACS,
        sampleRep=SAMPLEREPS,
    ),
)
lr_mapping_dupl = (
    expand(
        "output/mappings/longReadMapping/qc/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.bam.dupl.txt",
        filtered_product,
        techname=TECHNAMES,
        capDesign=CAPDESIGNS,
        sizeFrac=SIZEFRACS,
        sampleRep=SAMPLEREPS,
    ),
)
lr_mapping_basic_stats = expand(
    returnPlotFilenames(
        "output/plots/"
        + "lrMapping.basic.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.lrMapping.basic.stats"
    ),
    filtered_product,
    techname=TECHNAMESplusBY,
    capDesign=CAPDESIGNSplusBY,
    sizeFrac=SIZEFRACS,
    sampleRep=SAMPLEREPSplusBY,
)
intra_prming_stats = expand(
    returnPlotFilenames(
        "output/plots/"
        + "intraPriming.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.intraPriming.stats"
    ),
    filtered_product,
    techname=TECHNAMESplusBY,
    capDesign=CAPDESIGNSplusBY,
    sizeFrac=SIZEFRACS,
    sampleRep=SAMPLEREPSplusBY,
)
lr_mapping_spikeins_stats = expand(
    returnPlotFilenames(
        "output/plots/"
        + "lrMapping.spikeIns.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.lrMapping.spikeIns.stats"
    ),
    filtered_product,
    techname=TECHNAMESplusBY,
    capDesign=CAPDESIGNSplusBY,
    sizeFrac=SIZEFRACS,
    sampleRep=SAMPLEREPSplusBY,
)
sr_mapping_stats = expand(
    "output/plots/" + "hiSeq.mapping.stats/all.hiSeq.mapping.stats.{ext}", ext=plotFormat
)
sr_mapping_junctions_stats = expand(
    "output/plots/" + "hiSeq.SJs.stats/all.hiSeq.SJs.stats.{ext}", ext=plotFormat
)
merged_stats = expand(
    returnPlotFilenames(
        "output/plots/"
        + "merged.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.min{minReadSupport}reads.merged.stats"
    ),
    filtered_product,
    techname=TECHNAMESplusBY,
    capDesign=CAPDESIGNSplusBY,
    sizeFrac=SIZEFRACS,
    sampleRep=SAMPLEREPSplusBY,
    ext=plotFormat,
    minReadSupport=MINIMUM_TMERGE_READ_SUPPORT,
)
cage_polya_support_stats = expand(
    returnPlotFilenames(
        "output/plots/"
        + "cagePolyASupport.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.min{minReadSupport}reads.splicing_status-{splicedStatus}.cagePolyASupport.stats"
    ),
    filtered_product,
    techname=TECHNAMESplusBY,
    capDesign=CAPDESIGNSplusBY,
    sizeFrac=SIZEFRACS,
    sampleRep=SAMPLEREPSplusBY,
    minReadSupport=MINIMUM_TMERGE_READ_SUPPORT,
    splicedStatus=TMSPLICEDSTATUScategories,
)
tmerge_gencode_stats = expand(
    returnPlotFilenames(
        "output/plots/"
        + "tmerge.vs.gencode.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.splicing_status-{splicedStatus}.endSupport-{endSupport}.vs.gencode.stats"
    ),
    filtered_product,
    techname=TECHNAMESplusBY,
    capDesign=CAPDESIGNSplusBY,
    sizeFrac=SIZEFRACS,
    sampleRep=SAMPLEREPSplusBY,
    endSupport=ENDSUPPORTcategories,
    minReadSupport=MINIMUM_TMERGE_READ_SUPPORT,
    splicedStatus=TMSPLICEDSTATUScategories,
)
mature_rna_length_stats = expand(
    returnPlotFilenames(
        "output/plots/"
        + "matureRNALength.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.min{minReadSupport}reads.splicing_status-{splicedStatus}.matureRNALength.hist.stats"
    ),
    filtered_product,
    techname=TECHNAMESplusBY,
    capDesign=CAPDESIGNSplusBY,
    sizeFrac=SIZEFRACS,
    sampleRep=SAMPLEREPSplusBY,
    minReadSupport=MINIMUM_TMERGE_READ_SUPPORT,
    splicedStatus=TMSPLICEDSTATUScategories,
)
polya_reads_stats = expand(
    returnPlotFilenames(
        "output/plots/"
        + "polyAreads.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.polyAreads.stats"
    ),
    filtered_product,
    techname=TECHNAMESplusBY,
    capDesign=CAPDESIGNSplusBY,
    sizeFrac=SIZEFRACS,
    sampleRep=SAMPLEREPSplusBY,
)
tmerge_novel_loci_stats = expand(
    returnPlotFilenames(
        "output/plots/"
        + "tmerge.novelLoci.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.endSupport-{endSupport}.novelLoci.stats"
    ),
    filtered_product,
    techname=TECHNAMESplusBY,
    capDesign=CAPDESIGNSplusBY,
    sizeFrac=SIZEFRACS,
    sampleRep=SAMPLEREPSplusBY,
    endSupport=ENDSUPPORTcategories,
    minReadSupport=MINIMUM_TMERGE_READ_SUPPORT,
)
tmerge_gencode_junctions_stats = expand(
    returnPlotFilenames(
        "output/plots/"
        + "tmerge.vs.Gencode.SJs.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.vs.Gencode.SJs.stats"
    ),
    filtered_product,
    techname=TECHNAMESplusBY,
    capDesign=CAPDESIGNSplusBY,
    sizeFrac=SIZEFRACS,
    sampleRep=SAMPLEREPSplusBY,
    minReadSupport=MINIMUM_TMERGE_READ_SUPPORT,
)
tmerge_gencode_snpr_stats = expand(
    returnPlotFilenames(
        "output/plots/"
        + "tmerge.vs.gencode.SnPr.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.splicing_status-{splicedStatus}.endSupport-{endSupport}.vs.gencode.SnPr.stats"
    ),
    filtered_product,
    techname=TECHNAMESplusBY,
    capDesign=CAPDESIGNSplusBY,
    sizeFrac=SIZEFRACS,
    sampleRep=SAMPLEREPSplusBY,
    endSupport=ENDSUPPORTcategories,
    minReadSupport=MINIMUM_TMERGE_READ_SUPPORT,
    splicedStatus=TMSPLICEDSTATUScategories,
)
tmerge_gencode_length_splicing_status_stats = expand(
    returnPlotFilenames(
        "output/plots/"
        + "tmerge.vs.gencode.length.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.splicing_status-bySplicingStatus.endSupport-{endSupport}.vs.gencode.length.stats"
    ),
    filtered_product,
    techname=TECHNAMES,
    capDesign=CAPDESIGNS,
    sizeFrac=SIZEFRACS,
    sampleRep=SAMPLEREPS,
    endSupport=ENDSUPPORTcategories,
    minReadSupport=MINIMUM_TMERGE_READ_SUPPORT,
)
tmerge_gencode_length_stats = expand(
    returnPlotFilenames(
        "output/plots/"
        + "tmerge.vs.gencode.length.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.splicing_status-all.endSupport-{endSupport}.vs.gencode.length.stats"
    ),
    filtered_product,
    techname=TECHNAMES,
    capDesign=CAPDESIGNS,
    sizeFrac=SIZEFRACS,
    sampleRep=SAMPLEREPS,
    endSupport=ENDSUPPORTcategories,
    minReadSupport=MINIMUM_TMERGE_READ_SUPPORT,
)
dhs_cage_5prime_comparison_stats = expand(
    "output/plots/"
    + "dhsVsCage5primeComparison.venn.stats/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.dhsVsCage5primeComparison.venn.stats.pdf",
    filtered_product,
    techname=TECHNAMES,
    capDesign=CAPDESIGNS,
    sizeFrac=SIZEFRACS,
    sampleRep=SAMPLEREPS,
    minReadSupport=MINIMUM_TMERGE_READ_SUPPORT,
)
gene_read_coverage_stats = expand(
    returnPlotFilenames(
        "output/plots/"
        + "geneReadCoverage.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.geneReadCoverage.min{minReadSupport}reads.stats"
    ),
    filtered_product,
    techname=TECHNAMESplusBY,
    capDesign=CAPDESIGNSplusBY,
    sizeFrac=SIZEFRACS,
    sampleRep=SAMPLEREPSplusBY,
    minReadSupport=MINIMUM_TMERGE_READ_SUPPORT,
)
read_profile_density = expand(
    "output/plots/"
    + "readProfile/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.readProfile.density.png",
    filtered_product,
    techname="byTech",
    capDesign=CAPDESIGNS,
    sizeFrac=SIZEFRACS,
    sampleRep=SAMPLEREPS,
)
read_profile_heatmap = expand(
    "output/plots/"
    + "readProfile/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.readProfile.heatmap.png",
    filtered_product,
    techname="byTech",
    capDesign=CAPDESIGNS,
    sizeFrac=SIZEFRACS,
    sampleRep=SAMPLEREPS,
)
sequencing_error_stats = expand(
    returnPlotFilenames(
        "output/plots/"
        + "sequencingError.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.sequencingError.allErrors.stats"
    ),
    filtered_product,
    techname=TECHNAMESplusBY,
    capDesign=CAPDESIGNSplusBY,
    sizeFrac=SIZEFRACS,
    sampleRep=SAMPLEREPSplusBY,
)
deletions_only_stats = expand(
    returnPlotFilenames(
        "output/plots/"
        + "sequencingError.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.sequencingError.deletionsOnly.stats"
    ),
    filtered_product,
    techname=TECHNAMESplusBY,
    capDesign=CAPDESIGNSplusBY,
    sizeFrac=SIZEFRACS,
    sampleRep=SAMPLEREPSplusBY,
)
gencode_gene_detection_stats = expand(
    returnPlotFilenames(
        "output/plots/"
        + "gencode.geneDetection.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.splicing_status-{splicedStatus}.endSupport-{endSupport}.gencode.geneDetection.stats"
    ),
    filtered_product,
    techname=TECHNAMESplusBY,
    capDesign=CAPDESIGNSplusBY,
    sizeFrac=SIZEFRACS,
    sampleRep=SAMPLEREPSplusBY,
    endSupport=ENDSUPPORTcategories,
    minReadSupport=MINIMUM_TMERGE_READ_SUPPORT,
    splicedStatus=TMSPLICEDSTATUScategories,
)
read_biotypes_breakdown_stats = expand(
    returnPlotFilenames(
        "output/plots/"
        + "readToBiotypeBreakdown.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{sampleRep}_{spikeInCategories}.readToBiotypeBreakdown.stats"
    ),
    filtered_product,
    techname=TECHNAMESplusBY,
    capDesign=CAPDESIGNSplusBY,
    sizeFrac=SIZEFRACS,
    sampleRep=SAMPLEREPSplusBY,
    spikeInCategories=SPIKEINcategories,
)
tmerge_nt_coverage_genome_partition_stats = expand(
    returnPlotFilenames(
        "output/plots/"
        + "tmerge.ntCoverageByGenomePartition.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.endSupport-{endSupport}.ntCoverageByGenomePartition.stats"
    ),
    filtered_product,
    techname=TECHNAMESplusBY,
    capDesign=CAPDESIGNSplusBY,
    sizeFrac=SIZEFRACS,
    sampleRep=SAMPLEREPSplusBY,
    endSupport=ENDSUPPORTcategories,
    minReadSupport=MINIMUM_TMERGE_READ_SUPPORT,
)
tmerge_sirvs_stats = expand(
    returnPlotFilenames(
        "output/plots/"
        + "tmerge.vs.SIRVs.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.{filt}.tmerge.min{minReadSupport}reads.vs.SIRVs.stats"
    ),
    filtered_product,
    techname=TECHNAMESplusBY,
    capDesign=CAPDESIGNSplusBY,
    sizeFrac=SIZEFRACS,
    sampleRep=SAMPLEREPSplusBY,
    minReadSupport=MINIMUM_TMERGE_READ_SUPPORT,
    filt=READFILTERS,
)
tmerge_sirvs_detection_stats = expand(
    returnPlotFilenames(
        "output/plots/"
        + "tmerge.vs.SIRVs.detection.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.vs.SIRVs.detection.stats"
    ),
    filtered_product,
    techname=TECHNAMESplusBY,
    capDesign=CAPDESIGNSplusBY,
    sizeFrac=SIZEFRACS,
    sampleRep=SAMPLEREPSplusBY,
    minReadSupport=MINIMUM_TMERGE_READ_SUPPORT,
)
target_coverage_stats = expand(
    returnPlotFilenames(
        "output/plots/"
        + "targetCoverage.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.min{minReadSupport}reads.targetCoverage.stats"
    ),
    filtered_product,
    techname=TECHNAMESplusBY,
    capDesign=CAPDESIGNSplusBY,
    sizeFrac=SIZEFRACS,
    sampleRep=SAMPLEREPSplusBY,
    minReadSupport=MINIMUM_TMERGE_READ_SUPPORT,
)

################################
### HTML Summary stats table ###
################################
summary_stats = "output/html/index.html"

#################
### Track hub ###
#################
track_hub_file = "output/trackHub/" + "hub.txt"
track_hub_genomes = "output/trackHub/" + "genomes.txt"
track_hub_db = expand("output/trackHub/" + "{genome}/trackDb.txt", genome=GENOMES)


####################################
### Prepare inputs for main rule ###
####################################
inputs = []
DEFAULT_INPUTS = [gff, grouped_gff, grouped_read_to_tm]

for i in DEFAULT_INPUTS:
    inputs.append(i)

if config["produceStatPlots"]:
    inputs.append(hiss_stats)
    inputs.append(read_length_stats)
    inputs.append(fastq_qc)
    inputs.append(lr_mapping_dupl)
    inputs.append(lr_mapping_basic_stats)
    inputs.append(intra_prming_stats)
    inputs.append(lr_mapping_spikeins_stats)
    inputs.append(merged_stats)
    inputs.append(polya_reads_stats)
    inputs.append(tmerge_novel_loci_stats)
    if config.get("GENOMETOANNOTGTF"):
        inputs.append(tmerge_gencode_stats)
        inputs.append(tmerge_gencode_junctions_stats)
        inputs.append(tmerge_gencode_snpr_stats)
        inputs.append(tmerge_gencode_length_splicing_status_stats)
        inputs.append(tmerge_gencode_length_stats)
    inputs.append(gene_read_coverage_stats)
    inputs.append(read_profile_density)
    inputs.append(read_profile_heatmap)
    inputs.append(sequencing_error_stats)
    inputs.append(deletions_only_stats)
    inputs.append(gencode_gene_detection_stats)
    inputs.append(read_biotypes_breakdown_stats)
    inputs.append(tmerge_nt_coverage_genome_partition_stats)
    if config.get("genomeToCAGEpeaks"):
        inputs.append(cage_polya_support_stats)
    inputs.append(mature_rna_length_stats)
    if config.get("genomeToDHSpeaks"):
        inputs.append(dhs_cage_5prime_comparison_stats)
    if config.get("USE_MATCHED_ILLUMINA"):
        inputs.append(sr_mapping_stats)
        inputs.append(sr_mapping_junctions_stats)
    if SIRVpresent:
        inputs.append(tmerge_sirvs_stats)
        inputs.append(tmerge_sirvs_detection_stats)
    if config.get("CAPTURE"):
        inputs.append(target_coverage_stats)

if config["produceHtmlStatsTable"]:
    inputs.append(summary_stats)

if config["produceTrackHub"]:
    inputs.append(track_hub_file)
    inputs.append(track_hub_genomes)
    inputs.append(track_hub_db)


rule all:
    input:
        inputs,
