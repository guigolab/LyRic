# Dependencies

Install Snakemake


Using a Snakemake-compatible HPC environment such as SGE/UGE is highly recommended.

* Installation

cd to the directory where you intend to run the LyRic workflow.

Clone LyRic snakefiles:
git clone https://github.com/julienlag/LyRic.git ./

Clone conda environments into envs/:

 git clone https://github.com/julienlag/condaEnvExport.git envs/

Clone custom software utilities repo:

 git clone https://github.com/julienlag/utils.git utils/
chmod a+x utils/*

Customize *config.json and cluster_config.json to your needs
 
# Input

- `config[SAMPLE_ANNOT]`
- `capDesignToTargetsGff`: non-overlapping targeted regions (only for RNA capture samples), labelled by target type using the `gene_type` GFF attribute **if required** 
- LR FASTQ files in "fastqs/" subdirectory
	- File naming scheme for LR FASTQ files: `"fastqs/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.fastq.gz"` where:
		- `{techname}` is made up of the following subfields, separated by a hyphen (`-`): 
			- \<sequencing technology\>: *e.g.* `ont`, `pacBio`
			- \<sequencing center\>: *e.g.* `Crg`, `Cshl`
			- \<library preparation protocol\>: *e.g.* `SmartSeq2`, `CapTrap`
		- `{capDesign}`: for sample names that didn't undergo targeted RNA capture, this field should match the following regex: `\S+preCap$`. For captured samples, this field should match the name of the capture design, which is up to the user.
		- `{sizeFrac}`
		- `{sampleRep}` should match the following regex: '`(\S+)\d{2}Rep\d+`' (*e.g. `Brain01Rep1`*). The `(\S+)` prefix should match the value in the 'tissue' column of the corresponding row in the `config[SAMPLE_ANNOT]` sample annotation file. `\d{2}Rep\d+` is there to identify multiple replicates of the same experiments.

	
- genomes.fa
- Illumina FASTQ files in `config[HISEQ_FASTQDIR]` **if required**
- Reference annotation GTFs in `config[genomeToAnnotGtf]`  **if required**
- TSV containing SIRV info (<transcript_id>{tab}<length>{tab}<concentration> in `config[SIRVinfo]` **if required** 
- "annotations/repeatMasker/" + CAPDESIGNTOGENOME[wildcards.capDesign] + ".repeatMasker.bed"


# Output

LyRic will produce various output files based on the following Snakemake config values: 

- config['produceStatPlots']: boolean. If true, multiple statistics plots in PNG format will be output inside the `./plots/` subdirectory. See relevant section below for more details. 
- config['produceHtmlStatsTable']: boolean. If true, a detailed HTML table will be produced in the `./html/` subdirectory. See relevant section below for more details.
- config['produceTrackHub']: boolean. If true, LyRic will generate a UCSC Track Hub in the `./trackHub/` subdirectory. See relevant section below for more details.

If all values are set to false, the workflow will only produce one transcriptome GTF file per input FASTQ file.

### Transcriptome GTF file
## Plots
## HTML Stats Table
## Track Hub

# Glossary / Abbreviations

- **HCGM**: 

	(produced by snakemake rule `highConfidenceReads`).
	
	**H**igh-**C**onfidence **G**enome **M**appings. Filtered read-to-genome mappings characterized by:

	- only canonical introns (*if read is **spliced***)
	- no suspicious introns possibly arising from RT template switching ("*RT-facts*") (*if read is **spliced***)
	- minimum average sequencing quality of `filter_SJ_Qscore` on the read sequence around all their splice junctions. `filter_SJ_Qscore` is set on a per-sample basis in the corresponding column of the sample annotation file (`config[SAMPLE_ANNOT]`) (*if read is **spliced***)
	- a detectable, clipped polyA tail on the read (*if read is **unspliced***)

- **HiSS**: 

	(produced by Snakemake rule `getHiSeqSupportedHCGMs`).

	**Hi**-**S**eq-**S**upported read mappings. Those correspond to HCGMs (see above) that have all their splice junctions supported by at least one split read in the corresponding matched HiSeq sample, **if `use_matched_HiSeq` is set to `true`** for the corresponding sample in the sample annotation file (`config[SAMPLE_ANNOT]`). If `use_matched_HiSeq` is `false`, HiSS reads are exactly equivalent to HCGMs.

- **TM**:

	(produced by rules `mergedReads` and `mergedReadsGroupedSampleReps` ).

	**T**ranscript **M**odel. The evidence-based model of an RNA transcript represented as the genomic coordinates of its intron-exon structure. A gene model contains a set of exon-overlapping TMs.



...


