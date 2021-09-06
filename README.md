LyRic is a versatile automated transcriptome annotation and analysis workflow written in the [Snakemake](https://snakemake.readthedocs.io/en/stable/) language. Its core functionality is the production of:

- a set of high-quality RNA Transcript Models (TMs) mapped onto a genome sequence, based on Long-Read (LR) RNA sequencing data.
- various summary statistics plots and analysis results that describe the input and output data in details
- an interactive HTML table reporting statistics for each input sample, enabling easy and intuitive sample-to-sample comparison 
- a [UCSC Track Hub](http://genome.cse.ucsc.edu/goldenPath/help/hgTrackHubHelp.html) to display output TMs, as well as various other tracks produced by the pipeline.

 It is platform-agnostic, *i.e.* it can deal with data coming from both the ONT and PacBio platforms.

Full LyRic documentation is [here](documentation.md).


