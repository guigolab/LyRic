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
- LR FASTQ files in `config[LR_FASTQDIR]`
- genomes.fa
- Illumina FASTQ files in `config[HISEQ_FASTQDIR]` **if required**
- Reference annotation GTFs in `config[genomeToAnnotGtf]`  **if required**
- GTF of targeted regions in `config[TARGETSDIR] + "{capDesign}_primary_targets.exons.reduced.gene_type.segments.gtf"`  **if required** 
- TSV containing SIRV info (<transcript_id>{tab}<length>{tab}<concentration> in `config[SIRVinfo]` **if required** 
- "annotations/repeatMasker/" + CAPDESIGNTOGENOME[wildcards.capDesign] + ".repeatMasker.bed"


# Output

LyRic will produce various output files based on the following Snakemake config values: 

- config['produceStatPlots']: boolean. If true, multiple statistics plots in PNG and PDF format will be output inside the `./plots/` subdirectory. See relevant section below for more details. 
- config['produceHtmlStatsTable']: boolean. If true, a detailed HTML table will be produced in the `./html/` subdirectory. See relevant section below for more details.
- config['produceTrackHub']: boolean. If true, LyRic will generate a UCSC Track Hub in the `./trackHub/` subdirectory. See relevant section below for more details.

If all values are set to false, the workflow will only produce one transcriptome GTF file per input FASTQ file.

### Transcriptome GTF file
## Plots
## HTML Stats Table
## Track Hub
