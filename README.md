* Dependencies

Install Snakemake

gtfToGenePred, genePredToBed


# SQANTI3 is required only in 'full' mode
You should follow these instructions to install SQANTI3: https://github.com/ConesaLab/SQANTI3/wiki/SQANTI3-dependencies-and-installation 
SQANTI and cDNA_Cupcake should be installed under ~/bin/SQANTI3/ and ~/bin/SQANTI3/cDNA_Cupcake/, respectively.

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
 
