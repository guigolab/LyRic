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

* Anaconda installation (`miniconda`/`mambaforge`/`pixi`) for installing Snakemake
* Singularity

> [!NOTE]  
> It looks like Snakemake needs an installation of Anaconda even when the pipeline runs in a containerized environment

<!-- TODO: Add docs or reference on how to install Snakemake -->

[Pixi section](#setup-snakemake-with-pixi)

### Get the pipeline

Clone the GitHub repo to the folder you want to use as the working directory of the pipeline amd move to it:

```
git clone https://github.com/guigolab/lyric my_working_dir

cd my_working_dir
```

#### Setup Snakemake with Pixi

If you use [Pixi](https://pixi.sh/) to install your conda environments you can use the provided `pixi.toml` file to setup Snakemake and other requirements. Just run the following command from the pipeline directory:

```
pixi install
```

### Make a test run

From the pipeline folder extract the test dataset:

```
tar xf resources/test-data.tgz
```

Then run the pipeline on the test dataset with:

```
snakemake
```
