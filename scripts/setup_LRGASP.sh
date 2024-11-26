# cd to your working directory, then:

mkdir fastqs
mkdir genomes
cd genomes
#download genomes
cd genomes
synapse get syn25683364  
synapse get syn25683365
zcat lrgasp_grch38_sirvs.fasta.gz >hg38.fa
zcat lrgasp_grcm39_sirvs.fasta.gz > mm10.fa
cd ..

cd fastqs 
# name fastqs according to sample_name values in sample_annotations_LRGASP.tsv

cd ..

# follow further instructions in https://github.com/julienlag/LyRic to install and execute LyRic

