LyRic is a flexible, full-featured, automated transcriptome annotation and analysis workflow written in the [Snakemake](https://snakemake.readthedocs.io/en/stable/) language. Its core functionality is the production of:

- a set of high-quality RNA Transcript Models (TMs) mapped onto a genome sequence, based on Long-Read (LR) sequencing data. It is platform-agnostic, *i.e.* it can deal with data coming from both the ONT and PacBio platforms
- various summary statistics plots and analysis results that describe the input and output data in details
- an interactive HTML table reporting statistics for each input sample, enabling easy and intuitive sample-to-sample comparison 
- a [UCSC Track Hub](http://genome.cse.ucsc.edu/goldenPath/help/hgTrackHubHelp.html) to display output TMs, as well as various other tracks produced by LyRic.

# Dependencies


LyRic depends on the following software:

- **Conda**:  Official installation instructions [here](https://docs.anaconda.com/anaconda/install/linux/)
- **Snakemake**: Official installation instructions [here](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

Please install those as a prerequisite.


# Installation

Note that using a Snakemake-compatible HPC environment such as SGE/UGE is highly recommended.

1. `cd` to the directory where you intend to run the LyRic workflow (referred to as the **working directory** below).

2. Clone LyRic Snakefiles:

	`git clone https://github.com/julienlag/LyRic.git ./LyRic`

3. Clone relevant conda environments into `envs/`:

	`git clone https://github.com/julienlag/condaEnvExport.git envs/`

4. Clone custom software utilities repository:

	`git clone https://github.com/julienlag/utils.git utils/`


5. Customize the `*config.json` and `cluster_config.json` files according to your needs

All paths mentioned below are relative to the working directory.

# Execution


An example bash script launching LyRic in [cluster/DRMAA mode](https://snakemake.readthedocs.io/en/stable/tutorial/additional_features.html?highlight=drmaa#cluster-execution) is provided ('`run_snakemake_EXAMPLE.sh`'), together with its accompanying workflow configuration (`config_EXAMPLE.json`) and cluster configuration (`cluster_config.json`) files. Note that you will need to customize these config files manually before your first LyRic run. Once this is done, make sure you're in your working directory and issue the following command to run LyRic:

`./LyRic/run_snakemake_EXAMPLE.sh`

Please refer to Snakemake's [documentation](https://snakemake.readthedocs.io/en/stable/executing/cli.html) for more advanced usage.

# Input

## Mandatory 

- **LR FASTQ files**
	- Must be placed inside the `fastqs/` subdirectory
	- File naming scheme: **`"fastqs/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.fastq.gz"`** (no spaces or special characters allowed) where:
		- **`{techname}`** can be made up of one or more subfields, separated by a hyphen (`-`), *e.g.*: 
			- `<sequencing technology>`: *e.g.* '`ont`', '`pacBio`'
			- `<sequencing center>`: *e.g.* '`Crg`', '`Cshl`'
			- `<library preparation protocol>`: *e.g.* '`SmartSeq2`', '`CapTrap`'
		- **`{capDesign}`**: for sample names that didn't undergo targeted RNA capture, this field should match the following regex: `\S+preCap$`. For captured samples, this field should match the name of the capture design (up to the user but preferably a short string).
		- **`{sizeFrac}`**: the RNA/cDNA size fraction the file corresponds to. If no size fractionation was performed, should be `0+`.
		- **`{sampleRep}`** should match the following regex: '`(\S+)\d{2}Rep\d+`' (*e.g. `Brain01Rep1`*). The `(\S+)` prefix should match the value in the `tissue` column of the corresponding row in the sample annotation file (see below). The `\d{2}Rep\d+` suffix identifies multiple replicates of the same experiment.
		- the file basename (*i.e.* `{techname}_{capDesign}_{sizeFrac}_{sampleRep}`) should have an exact match in the `sample_name` column of the sample annotation file (see below).
- **Sample annotation file**:

	This tab-separated file contains all metadata associated to each sample/input LR FASTQ file, as well as some customizable, sample-specific LyRic run parameters. Its path is controlled by config variable `SAMPLE_ANNOT`. A mock sample annotation file, named `sample_annotations_EXAMPLE.tsv` is included in this repo. 

	The contents of each column should me mostly self-explanatory, except:

	- **`cappedSpikeIns`** (boolean): `True` if the corresponding sample contains ERCC and/or SIRV spike-ins that were artificially 5'-m7G-capped before their incorporation into the library, `False` otherwise.
	- **`cellFrac`** (string): the cell fraction, *e.g.* `Cytoplasm`, `Nucleus` or `Total` (*i.e.* whole-cell extract with no cell fractionation)
	- **`filter_SJ_Qscore`** (integer): minimum average Phred sequencing quality of read sequences +/- 3 nts around all their splice junctions for a spliced read to be considered high-confidence (see "HCGM" glossary entry below). Example values: `10` for ONT, `30` for PacBio HiFi reads.
	- **`subProject`** (string): The sub-project this dataset is part of. This is useful to produce separate interactive HTML summary stats tables (see below) for each sub-project, if desired. Can be any string except '`ALL`'.
	- **`sample_name`** (string): should match the basename of the corresponding LR FASTQ file and should be globally unique.
	- **`use_matched_HiSeq`** (boolean): controls the production of HiSS files for the corresponding LR FASTQ file (see corresponding entry in glossary below).
	
- **Genome sequences** 
	
	To map RNA sequencing reads against. One genome assembly per file, in (multi-)FASTA format. See description of config variable `GENOMESDIR` below for more details.

## Optional

- **Short-read Illumina FASTQ files**

	If present, short reads contained in these files will be used to confirm splice junctions present in the LR FASTQ files.

	Only needed if config variable `USE_MATCHED_ILLUMINA` is `True`. 


- Reference annotation GTFs in `config[genomeToAnnotGtf]`  **if required**
- TSV containing SIRV info (<transcript_id>{tab}<length>{tab}<concentration> in `config[SIRVinfo]` **if required** 
- "annotations/repeatMasker/" + CAPDESIGNTOGENOME[wildcards.capDesign] + ".repeatMasker.bed"
- `config[capDesignToTargetsGff]`: non-overlapping targeted regions (only for RNA capture samples), labelled by target type using the `gene_type` GFF attribute **if required** 

# Workflow configuration variables

The following config variables are user-customizable. These can be set either via a config file (Snakemake's `--configfile FILE` option) or directly via command line options (`--config [KEY=VALUE [KEY=VALUE ...]]]`). See [Snakemake's CLI's documentation](https://snakemake.readthedocs.io/en/stable/executing/cli.html) for more details.

## Mandatory


- **`capDesignToGenome`** (dictionary/JSON object): ***key***: capture design identifier (`{capDesign}`); ***value***: corresponding genome assembly identifier (*e.g.* `hg38`, `mm10`)

- **`CAPTURE`** (boolean): set to '`True`' if some of your data underwent RNA capture and you want to obtain related target coverage statistics etc.

- **`GENOMESDIR`** (string): The path to the directory containing the genome multifastas to map the sequencing reads against. In it, there should be one genome file per species, whose name should match values in the `capDesignToGenome{}` object (e.g. '`hg38.fa`' and '`mm10.fa`')

- **`produceHtmlStatsTable`** (boolean): set to '`True`' if you want LyRic to produce  interactive HTML summary tables containing various sample-specific statistics (see below)
- **`produceStatPlots`** (boolean): set to '`True`' if you want LyRic to produce summary statistics plots (see below)
- **`produceTrackHub`** (boolean): set to '`True`' if you want LyRic to produce a UCSC Track Hub based on your data
- **`SAMPLE_ANNOT`** (string): Path to the sample annotation file (see above)
- **`sampleRepGroupBy`** (list/JSON array): list of columns in the sample annotation file used to combine `{sampleReps}` into groups of samples. Usually, this is the combination of experimental metadata properties that form a distinct group of replicates to group by. This list is used to merge TMs from distinct `{sampleReps}`. See relevant section below for further details.

- **`PROJECT_NAME`** for track hub + job name

- **`USE_MATCHED_ILLUMINA`** (boolean)


## Optional

- **`capDesignToCapDesign`** (dictionary/JSON object): ***key***: pre- or post-capture design identifier (`{capDesign}`); ***value***: corresponding post-capture design identifier. Only needed if config variable `CAPTURE` is `True`.

- **`capDesignToTargetsGff`** (dictionary/JSON object): ***key***: capture design identifier (`{capDesign}`); ***value***: path to non-overlapping capture-targeted regions in `{capDesign}`, in GFF format. Only needed if config variable `CAPTURE` is `True`.

- **`genomeToAnnotGtf`** (dictionary/JSON object): ***key***: genome assembly identifier (*e.g.* `hg38`, `mm10`); ***value***: corresponding gene annotation file, in GTF format (*e.g.* gencode v24 GTF). Only needed if any of config variables `produceStatPlots`, `produceTrackHub` and `produceHtmlStatsTable` are `True`.

- **`genomeToCAGEpeaks`**
- **`genomeToDHSpeaks`**
- **`PROJECT_CONTACT_EMAIL`**
- **`PROJECT_LONG_NAME`**
- **`REPEATMASKER_DIR`** (string): 
- **`SIRVinfo`**
- **`TRACK_HUB_BASE_URL`**



# Output

All files output by LyRic are written under '`output/`' in the working directory. LyRic will generate various types of output files, listed below.

The production of each output type can be turned on and off using boolean Snakemake config variables, namely `config['produceStatPlots']`, `config['produceHtmlStatsTable']` and `config['produceTrackHub']`. If all these are set to false, the workflow will only produce transcriptome GTF files (this cannot be switched off).

## Transcriptome GTF files

### Sample-specific TMs

(Output directory: `output/mappings/mergedReads/`)

One GTF file per input FASTQ

### TMs merged across replicates

(Output directory: `output/mappings/mergedReads/groupedSampleReps/`)

Samples-specific TMs can be further merged across samples according to config variable `sampleRepGroupBy`. The output transcriptome files will be named based on the column values used to group by. For example, using `config_EXAMPLE.json` and `sample_annotations_EXAMPLE.tsv`:

- Samples 
	
	`pacBioSII-Uci-SmartSeq2_HpreCap_0+_Mix01Rep2` and 
	
	`pacBioSII-Uci-SmartSeq2_HpreCap_0+_Mix03Rep2` 
	
	will be merged into 
	
	`pacBioSII_SmartSeq2_Human_HpreCap_0+_Mix_Total_Maxima`

- Samples 

	`ont-Crg-CapTrap_HpreCap_0+_HEK293T01Rep2` and 
	
	`ont-Crg-CapTrap_HpreCap_0+_HEK293T01Rep3` 
	
	will be merged into 
	
	`ONT_CapTrap_Human_HpreCap_0+_HEK293T_Total_PrimeScriptII`

- All other samples will be left untouched, as they correspond to singleton replicates. Note that a file will be created in the output directory anyway, following the file naming scheme described above.

## Summary statistics plots

(Output directory:  `./output/plots/`)

Controlled by config variable '`produceStatPlots`': boolean. If `True`, multiple statistics plots in PNG format will be output inside the output directory.

## Interactive HTML summary stats table

(Output directory:  `./output/html/`)

Production of the interactive HTML summary table is controlled by config variable '`produceHtmlStatsTable`' variable. If set to `True`, detailed reports containing various per-sample statistics will be produced in the output directory. 

LyRic will produce one table per distinct `subProject` value in the input sample annotation file (as long as those correspond to actual input FASTQ files) (`./output/html/summary_table_{subProject}.html`), plus a global one containing info for all samples (`./output/html/summary_table_ALL.html`)

For each interactive HTML summary stats table, an accompanying TSV file with the same basename and the `.tsv` extension will also be produced. It contains the same data as the HTML table, in an easily parsable tab-separated format.

## UCSC Track Hub

- Controlled by config variable '`produceTrackHub`': boolean. If true, LyRic will generate a UCSC Track Hub in the `./trackHub/` subdirectory. 



# Glossary / Abbreviations

- **HCGM**: 

	(produced by snakemake rule `highConfidenceReads`).
	
	**H**igh-**C**onfidence **G**enome **M**appings. Filtered read-to-genome mappings characterized by:

	- only canonical introns (*if read is **spliced***)
	- no suspicious introns possibly arising from RT template switching ("*RT-facts*") (*if read is **spliced***)
	- minimum average sequencing quality of `filter_SJ_Qscore` on the read sequence +/- 3 nts around all their splice junctions. `filter_SJ_Qscore` is set on a per-sample basis in the corresponding column of the sample annotation file (*if read is **spliced***)
	- a detectable, clipped polyA tail on the read (*if read is **unspliced***)

- **HiSS**: 

	(produced by Snakemake rule `getHiSeqSupportedHCGMs`).

	**Hi**-**S**eq-**S**upported read mappings. Those correspond to HCGMs (see above) that have all their splice junctions supported by at least one split read in the corresponding `{capDesign}`-matched HiSeq sample, **if `use_matched_HiSeq` is set to `true`** for the corresponding sample in the sample annotation file (`config[SAMPLE_ANNOT]`). If `use_matched_HiSeq` is `false`, HiSS reads are exactly equivalent to HCGMs.

- **LR**: **L**ong sequencing **R**ead (typically produced by the PacBio and ONT platforms) 

- **TM**:

	(produced by rules `mergedReads` and `mergedReadsGroupedSampleReps` ).

	**T**ranscript **M**odel. The evidence-based model of an RNA transcript represented as the genomic coordinates of its intron-exon structure. A gene model contains a set of exon-overlapping TMs.



...


