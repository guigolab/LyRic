# LyRic

LyRic is a versatile automated transcriptome annotation and analysis workflow written in the [Snakemake](https://snakemake.readthedocs.io/en/stable/) language. Its core functionality is the production of:

1. a set of high-quality RNA **Transcript Models (TMs)** mapped onto a genome sequence, based on Long-Read (LR) RNA sequencing data.
2. various **summary statistics plots** and analysis results that describe the input and output data in details
3. an **interactive HTML table** reporting statistics for each input sample, enabling easy and intuitive sample-to-sample comparison 
4. a **[UCSC Track Hub](http://genome.cse.ucsc.edu/goldenPath/help/hgTrackHubHelp.html)** to display output TMs, as well as various other tracks produced by LyRic.

(Note that features 2, 3 and 4 can be easily switched on and off).

LyRic is platform-agnostic, *i.e.* it can deal with FASTQ data coming from both the ONT and PacBio platforms.

**Full LyRic documentation is [here](https://guigolab.github.io/LyRic/documentation.html).**


## Quickstart

### Prerequisites

* Anaconda installation (`miniconda`/`mambaforge`) 
* Snakemake v8
* Singularity

> [!NOTE]  
> It looks like Snakemake needs an installation of Anaconda even when the pipeline runs in a containerized environment

> [!TIP]
> If you use [Pixi](https://pixi.sh/) to install your conda environments you can use the provided `pixi.toml` file to setup Snakemake and other requirements. Just run the following command from the pipeline directory:
>
> ```
> pixi install
> ```

### Get the pipeline

Clone the GitHub repo to the folder you want to use as the working directory of the pipeline and move to it:

```
git clone https://github.com/guigolab/LyRic lyric_test

cd lyric_test
```

### Make a test run

The pipeline repository contains a small datasets that can be used for testing. You can run the pipeline on the test dataset with the following command:

```
snakemake
```
Execution shall take < half an hour.
