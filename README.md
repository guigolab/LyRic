LyRic is a full-featured, automated transcriptome annotation and analysis workflow written in the [Snakemake](https://snakemake.readthedocs.io/en/stable/) language. Its core functionality is the production of:

- a set of high-quality RNA Transcript Models (TMs) mapped onto a genome sequence, based on Long-Read (LR) sequencing data (*e.g.* from the ONT or PacBio platforms)
- various summary statistics plots and analysis results that describe the input and output data in details
- an interactive HTML table reporting statistics for each input sample
- a [UCSC Track Hub](http://genome.cse.ucsc.edu/goldenPath/help/hgTrackHubHelp.html) to display output TMs as well as various other tracks produced by LyRic.

# Dependencies


LyRic depends on the following software:

- **Conda**:  Official installation instructions [here](https://docs.anaconda.com/anaconda/install/linux/)
- **Snakemake**: Official installation instructions [here](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

Please install those as a prerequisite.


# Installation

Note that using a Snakemake-compatible HPC environment such as SGE/UGE is highly recommended.

1. `cd` to the directory where you intend to run the LyRic workflow.

2. Clone LyRic Snakefiles:

	`git clone https://github.com/julienlag/LyRic.git ./`

3. Clone relevant conda environments into `envs/`:

	`git clone https://github.com/julienlag/condaEnvExport.git envs/`

4. Clone custom software utilities repository:

	`git clone https://github.com/julienlag/utils.git utils/`

	and make them executable:
	
	`chmod a+x utils/*`

5. Customize the `*config.json` and `cluster_config.json` files according to your needs

All paths mentioned below are relative to the working directory.

# Execution

# Input

- **LR FASTQ files**:
	- Must be placed inside the `fastqs/` subdirectory
	- File naming scheme: `"fastqs/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.fastq.gz"` (no spaces or special characters allowed) where:
		- `{techname}` is made up of the following subfields, separated by a hyphen (`-`): 
			- \<sequencing technology\>: *e.g.* `ont`, `pacBio`
			- \<sequencing center\>: *e.g.* `Crg`, `Cshl`
			- \<library preparation protocol\>: *e.g.* `SmartSeq2`, `CapTrap`
		- `{capDesign}`: for sample names that didn't undergo targeted RNA capture, this field should match the following regex: `\S+preCap$`. For captured samples, this field should match the name of the capture design (up to the user but preferably a short string).
		- `{sizeFrac}`: the RNA/cDNA size fraction the file corresponds to. If no size fractionation was performed, should be `0+`.
		- `{sampleRep}` should match the following regex: '`(\S+)\d{2}Rep\d+`' (*e.g. `Brain01Rep1`*). The `(\S+)` prefix should match the value in the `tissue` column of the corresponding row in the sample annotation file (see below). The `\d{2}Rep\d+` suffix identifies multiple replicates of the same experiment.

- **Sample annotation file**: This tab-separated file contains all metadata associated to each sample/input LR FASTQ file. Its path is controlled by Snakemake config value `config[SAMPLE_ANNOT]`. A mock sample annotation file, named `sample_annotations.tsv` is included in this repo. 

	The contents of each column should me mostly self-explanatory, except:

	- `sample_name` (string): should match the basename of the corresponding LR FASTQ file and should be globally unique.
	- `use_matched_HiSeq` (boolean): controls the production of HiSS files for the corresponding LR FASTQ file (see corresponding entry in glossary below).
	- `project` (string): The sub-project this dataset is part of. This is useful to produce separate interactive HTML summary stats tables (see below) for each sub-project, if desired. Can be any string except '`ALL`'.
	- `filter_SJ_Qscore` (integer): minimum average Phred sequencing quality of read sequences +/- 3 nts around all their splice junctions for a spliced read to be considered high-confidence (see "HCGM" glossary entry below). Recommended values: `10` for ONT, `30` for PacBio HiFi reads.
	
- genomes.fa
- Illumina FASTQ files in `config[HISEQ_FASTQDIR]` **if required**
- Reference annotation GTFs in `config[genomeToAnnotGtf]`  **if required**
- TSV containing SIRV info (<transcript_id>{tab}<length>{tab}<concentration> in `config[SIRVinfo]` **if required** 
- "annotations/repeatMasker/" + CAPDESIGNTOGENOME[wildcards.capDesign] + ".repeatMasker.bed"
- `config[capDesignToTargetsGff]`: non-overlapping targeted regions (only for RNA capture samples), labelled by target type using the `gene_type` GFF attribute **if required** 


# Output

All files output by LyRic are written under '`output/`' in the working directory. LyRic will generate various types of output files, listed below.

The production of each output type can be turned on and off using boolean Snakemake config variables, namely `config['produceStatPlots']`, `config['produceHtmlStatsTable']` and `config['produceTrackHub']`. If all these are set to false, the workflow will only produce transcriptome GTF files (this cannot be switched off).

## Transcriptome GTF files

## Summary statistics plots

Controlled by Snakemake's `config['produceStatPlots']` config variable: boolean. If `True`, multiple statistics plots in PNG format will be output inside the `./output/plots/` subdirectory.

## Interactive HTML summary stats table

Production of the interactive HTML summary table is controlled by Snakemake's `config['produceHtmlStatsTable']` variable. If set to `True`, detailed reports containing various per-sample statistics will be produced in the `./output/html/` subdirectory. 

LyRic will produce one table per distinct `project` value in the input sample annotation file (as long as those correspond to actual input FASTQ files) (`./output/html/summary_table_{project}.html`), plus a global one containing info for all samples (`./output/html/summary_table_ALL.html`)

For each interactive HTML summary stats table, an accompanying TSV file with the same basename and the `.tsv` extension will also be produced. It contains the same data as the HTML table, in an easily parsable tab-separated format.

## Track Hub

- `config['produceTrackHub']`: boolean. If true, LyRic will generate a UCSC Track Hub in the `./trackHub/` subdirectory. 

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


