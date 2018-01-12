import glob
from collections import defaultdict
import os.path
from itertools import product


# # path to dropbox folder where to sync output R plots:
DROPBOX_PLOTS=config["DROPBOX_PLOTSDIR"]
GGPLOT_PUB_QUALITY=config["GGPLOT_PUB_QUALITY"]
#TECHNAME=config["TECHNAME"]
CAPDESIGNTOGENOME=config["capDesignToGenome"]
#mappingDir=config["PB_MAPPINGS"]
sizeFrac_Rpalette=config["SIZEFRACTION_RPALETTE"]
#print(CAPDESIGNTOGENOME)

# no underscores allowed in wildcards, to avoid greedy matching since we use them as separators
wildcard_constraints:
 	capDesign = "[^_/]+",
 	sizeFrac = "[^_/]+",
 	techname = "[^_/]+"

# get CAPDESIGNS (capture designs, i.e. Hv1, Hv2, Mv1) and SIZEFRACS (size fractions) variables from FASTQ file names (warning: this will generate duplicate entries):
(TECHNAMES, CAPDESIGNS, SIZEFRACS) = glob_wildcards(config["FQPATH"] + "{techname}_{capDesign}_{sizeFrac}.fastq")
# remove duplicate entries:
CAPDESIGNS=set(CAPDESIGNS)
SIZEFRACS=set(SIZEFRACS)
TECHNAMES=set(TECHNAMES)



adaptersTSV = "demultiplexing/all_adapters.tsv"
f = open(adaptersTSV, 'r')
BARCODES = []
BARCODESUNDETER = []
CAPDESIGNTOBARCODES = defaultdict(list)
for line in f:
	columns = line.split("\t")
	barcodeId = columns[0].split("_")
	if columns[0] != 'UP':
		BARCODES.append(columns[0])
		BARCODESUNDETER.append(columns[0])
		CAPDESIGNTOBARCODES[barcodeId[0]].append(columns[0])

BARCODESUNDETER.append("Undeter")
BARCODES=set(BARCODES)
BARCODESUNDETER=set(BARCODESUNDETER)

include: "demultiplex.snakefile.py"

#to avoid AmbiguousRuleException: (and in that order, otherwise {capDesign}_{sizeFrac}.Undeter.fastq will be generated using rule demultiplexFastqs and not getUndeterminedReads, which will always give empty Undeter files :
ruleorder: getUndeterminedReads > demultiplexFastqs

#pseudo-rule specifying the target files we ultimately want.
rule all:
	input:
		expand(config["PLOTSDIR"] + "{techname}_{capDesign}_all.readlength.{ext}", techname=TECHNAMES, capDesign=CAPDESIGNS, ext=config["PLOTFORMATS"]), # facetted histograms of read length
		expand(config["STATSDATADIR"] + "{techname}.fastq.readCounts.tsv", techname=TECHNAMES), #read counts per fastq
		expand(config["PLOTSDIR"] + "{techname}.fastq.UP.stats.{ext}", techname=TECHNAMES, ext=config["PLOTFORMATS"]), # UP reads plots
		expand(config["PLOTSDIR"] + "{techname}.fastq.BC.stats.{ext}", techname=TECHNAMES, ext=config["PLOTFORMATS"]), # barcode reads plots
		expand(config["PLOTSDIR"] + "{techname}.fastq.foreignBC.stats.{ext}", techname=TECHNAMES, ext=config["PLOTFORMATS"]), #foreign barcode reads plots

		expand ("mappings/" + "{techname}_{capDesign}_{sizeFrac}.{barcodesU}.bam", techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodesU=BARCODESUNDETER),  # STAR-mapped reads
		expand(config["PLOTSDIR"] + "{techname}.ambiguousBarcodes.reads.stats.{ext}", techname=TECHNAMES, ext=config["PLOTFORMATS"]), # ambiguous barcodes plots
		expand(config["PLOTSDIR"] + "{techname}_{capDesign}.adapters.location.stats.{ext}",techname=TECHNAMES, capDesign=CAPDESIGNS, ext=config["PLOTFORMATS"]), #location of adapters over reads
		expand(config["PLOTSDIR"] + "{techname}.chimeric.reads.stats.{ext}", techname=TECHNAMES, ext=config["PLOTFORMATS"]), # stats on chimeric reads
		expand(config["PLOTSDIR"] + "{techname}.finalDemul.reads.stats.{ext}", techname=TECHNAMES, ext=config["PLOTFORMATS"]), #final demultiplexing stats
		expand(config["DEMULTIPLEX_DIR"] + "qc/{techname}_{capDesign}_{sizeFrac}.demul.QC1.txt", techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS), # QC on demultiplexing (checks that there is only one barcode assigned per read
		expand(config["DEMULTIPLEX_DIR"] + "{techname}_{capDesign}_{sizeFrac}.{barcodes}.fastq",techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES), #get demultiplexed FASTQs
		expand(config["DEMULTIPLEX_DIR"] + "{techname}_{capDesign}_{sizeFrac}.Undeter.fastq", techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS), # get Undetermined (non-demultiplexed) reads
		expand(config["DEMULTIPLEX_DIR"] + "qc/{techname}_{capDesign}_{sizeFrac}.demul.QC2.txt", techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS), # QC on demultiplexing (checks that there is only one barcode assigned per read
		expand(config["DEMULTIPLEX_DIR"] + "qc/{techname}_{capDesign}_{sizeFrac}.{barcodes}.demul.QC3.txt", techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES),
		expand(config["PLOTSDIR"] + "{techname}.demultiplexing.perSample.stats.{ext}", techname=TECHNAMES, ext=config["PLOTFORMATS"]),
		expand(config["PLOTSDIR"] + "{techname}.mapping.perSample.perFraction.stats.{ext}", techname=TECHNAMES, ext=config["PLOTFORMATS"])
#		expand("{capDesign}_{sizeFrac}.{barcodesU}.test", capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodesU=BARCODESUNDETER)


rule getFastqReadCounts:
	input: config["FQPATH"] + "{techname}_{capDesign}_{sizeFrac}.fastq"
	output: temp(config["STATSDATADIR"] + "{techname}_{capDesign}_{sizeFrac}.fastq.readCounts.tsv")
	shell:
		'''
let total=$(cat {input} | wc -l)/4
echo -e "$(basename {input})\t$total" > {output}
		'''

rule aggFastqReadCounts:
	input:  expand (config["STATSDATADIR"] + "{techname}_{capDesign}_{sizeFrac}.fastq.readCounts.tsv", techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS)
	output: config["STATSDATADIR"] + "{techname}.fastq.readCounts.tsv"
	shell: "cat {input} | sort > {output}"

#get read lengths for all FASTQ files:
rule getReadLength:
	input: config["FQPATH"] + "{techname}_{capDesign}_{sizeFrac}.fastq"
	output: temp("{techname}_{capDesign}_{sizeFrac}.readlength.tsv")
	shell: "fastq2fasta.pl {input} | FastaToTbl | awk -v s={wildcards.sizeFrac} '{{print s\"\\t\"length($2)}}' > {output}"

#aggregate read length data over all fractions of a given capDesign:
rule aggReadLength:
	input: expand("{{techname}}_{{capDesign}}_{sizeFrac}.readlength.tsv", sizeFrac=SIZEFRACS)
	#input: glob.glob(os.path.join("{capDesign}_*.readlength.tsv"))
	output: config["STATSDATADIR"] + "{techname}_{capDesign}_all.readlength.tsv"
	shell: "cat {input} > {output}"

# plot histograms with R:
rule plotReadLength:
	input: config["STATSDATADIR"] + "{techname}_{capDesign}_all.readlength.tsv"
	output: config["PLOTSDIR"] + "{techname}_{capDesign}_all.readlength.{ext}"
	shell:
		'''
echo "library(ggplot2)
library(plyr)
dat <- read.table('{input}', header=F, as.is=T, sep='\\t')
colnames(dat)<-c('fraction','length')
medians <- ddply(dat, 'fraction', summarise, grp.median=median(length))
sizes <- ddply(dat, 'fraction', summarise, grp.length=length(length))
ggplot(dat, aes(x=length)) +
geom_histogram(binwidth=200, aes(fill=fraction)) +
scale_fill_manual(values={sizeFrac_Rpalette}) +
facet_grid( . ~ fraction) +
annotate(geom='text', label=c(paste('n= ', sizes[,'grp.length'])), hjust=0, vjust=c(3.8), x=10, y=c(Inf), color=c('black'), size=6) +
annotate(geom='text', label=c(paste('Median (nt)= ', medians[,'grp.median'])), hjust=0, vjust=c(5.8), x=10, y=c(Inf), color=c('black'), size=6) +
coord_cartesian(xlim=c(0, 5000)) +
theme_bw(base_size=17) +
{GGPLOT_PUB_QUALITY} + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('{output}', width=10, height=3)
" > {output}.r
cat {output}.r | R --slave
dropbox_uploader.sh upload {output} {DROPBOX_PLOTS};

	 	'''



# # STAR mapping of reads:
# rule readMapping:
# #	wildcard_constraints:
# #		 barcodesU = lambda wildcards: {wildcards.capDesign} + "_.+"
# 	input:
# #		reads = returnCapDesignBarcodesFastqs,
# 		reads = config["DEMULTIPLEX_DIR"] + "{capDesign}_{sizeFrac}.{barcodesU}.fastq",
# 		reference=  lambda wildcards: "/users/rg/buszczynska/Projects/RNA.capture/pacbio.postCapture.20150317/mapping/reference/" + CAPDESIGNTOGENOME[wildcards.capDesign] + "/star.index/"
# 	output:
# 		config["PB_MAPPINGS"] + "{capDesign}_{sizeFrac}.{barcodesU}.bam"
# 	shell:
# 		'''
# # (Recommended STAR parameters from https://github.com/PacificBiosciences/cDNA_primer/wiki/Bioinfx-study:-Optimizing-STAR-aligner-for-Iso-Seq-data )

# sameBarcodeCapDesign=$(echo "{wildcards.capDesign} {wildcards.barcodesU}" | perl -slane 'if($F[1]=~/$F[0]/ || $F[1] eq "Undeter"){{print "1"}} else {{print "0"}}')

# if [ {TECHNAME} == "ont" && {wildcards.barcodesU} == "Undeter" ]; then
# run=0
# else
# run=1
# fi


# if [ $sameBarcodeCapDesign == 1 && $run == 1 ]; then
# /users/rg/jlagarde/bin/STAR-2.5.3a/bin/Linux_x86_64/STARlong \
# --runMode alignReads \
# --outSAMattributes NH HI NM MD \
# --readNameSeparator space \
# --outFilterMultimapScoreRange 1 \
# --outFilterMismatchNmax 2000 \
# --scoreGapNoncan -20 \
# --scoreGapGCAG -4 \
# --scoreGapATAC -8 \
# --scoreDelOpen -1 \
# --scoreDelBase -1 \
# --scoreInsOpen -1 \
# --scoreInsBase -1 \
# --alignEndsType Local \
# --seedSearchStartLmax 50 \
# --seedPerReadNmax 100000 \
# --seedPerWindowNmax 1000 \
# --alignTranscriptsPerReadNmax 100000 \
# --alignTranscriptsPerWindowNmax 10000 \
# --readFilesIn {input.reads} \
# --genomeDir {input.reference} \
# --outFileNamePrefix {pacBioMappingDir}/{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.barcodesU} \
# --outStd SAM \
# | samtools view -b -u -S - \
# | samtools sort -m 20000000000 - >{output}
# samtools index {output}
# else
# touch {output}
# fi
# 		'''

# GMAP mapping of reads:
rule readMapping:
#	wildcard_constraints:
#		 barcodesU = lambda wildcards: {wildcards.capDesign} + "_.+"
	input:
#		reads = returnCapDesignBarcodesFastqs,
		reads = config["DEMULTIPLEX_DIR"] + "{techname}_{capDesign}_{sizeFrac}.{barcodesU}.fastq",
		genome = lambda wildcards: "/users/rg/jlagarde/genomes/" + CAPDESIGNTOGENOME[wildcards.capDesign] + ".fa"
#	params:
#		reference=  lambda wildcards: CAPDESIGNTOGENOME[wildcards.capDesign]
	threads: 12
	output:
		"mappings/" + "{techname}_{capDesign}_{sizeFrac}.{barcodesU}.bam"
	shell:
		'''
sameBarcodeCapDesign=$(echo "{wildcards.capDesign} {wildcards.barcodesU}" | perl -slane 'if($F[1]=~/$F[0]/ || $F[1] eq "Undeter"){{print "1"}} else {{print "0"}}')
run=1
if [[ $sameBarcodeCapDesign == 1 && $run == 1 ]]; then
echoerr "Mapping"
minimap2 -t {threads} -L -ax splice {input.genome} {input.reads} > {output}.tmp
echoerr "Mapping done"
echoerr "Creating/sorting BAM"
cat {output}.tmp | samtools view -b -u -S - | samtools sort --threads {threads} -T $TMPDIR -m 5G - >{output}
echoerr "Done creating/sorting BAM"
echoerr "Indexing BAM"
samtools index {output}
rm {output}.tmp
else
touch {output}
fi
		'''

# gmap -D /users/rg/jlagarde/genomes/gmapdb/ -d {params.reference} -f samse -n 0 -t {threads} {input.reads} > {output}.tmp


rule getMappingStats:
	input:
		bams = "mappings/" + "{capDesign}_{sizeFrac}.{barcodesU}.bam",
		fastqs = config["DEMULTIPLEX_DIR"] + "{techname}_{capDesign}_{sizeFrac}.{barcodesU}.fastq"
	output: temp(config["STATSDATADIR"] + "{techname}_{capDesign}_{sizeFrac}.{barcodesU}.mapping.perSample.perFraction.stats.tsv")
	shell:
		'''
sameBarcodeCapDesign=$(echo "{wildcards.capDesign} {wildcards.barcodesU}" | perl -slane 'if($F[1]=~/$F[0]/ || $F[1] eq "Undeter"){{print "1"}} else {{print "0"}}')
if [ $sameBarcodeCapDesign == 1 ]; then
totalReads=$(cat {input.fastqs} | fastq2tsv.pl | wc -l)
mappedReads=$(samtools view  -F 4 {input.bams}|cut -f1|sort|uniq|wc -l)
echo -e "{wildcards.capDesign}\t{wildcards.sizeFrac}\t{wildcards.barcodesU}\t$totalReads\t$mappedReads" | awk '{{print $0"\t"$5/$4}}' > {output}
else
touch {output}
fi
		'''
rule aggMappingStats:
	input:expand(config["STATSDATADIR"] + "{techname}_{capDesign}_{sizeFrac}.{barcodesU}.mapping.perSample.perFraction.stats.tsv", techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodesU=BARCODESUNDETER)
	output: config["STATSDATADIR"] + "{techname}.mapping.perSample.perFraction.stats.tsv"
	shell:
		'''
cat {input} | sort > {output}
		'''

rule plotMappingStats:
	input: config["STATSDATADIR"] + "{techname}.mapping.perSample.perFraction.stats.tsv"
	output: config["PLOTSDIR"] + "{techname}.mapping.perSample.perFraction.stats.{ext}"
	shell:
		'''
echo "library(ggplot2)
library(plyr)
library(scales)
dat <- read.table('{input}', header=F, as.is=T, sep='\\t')
colnames(dat)<-c('capDesign', 'sizeFraction','barcode','totalReads', 'mappedReads', 'percentMappedReads')
ggplot(dat, aes(x=barcode, y=percentMappedReads, fill=sizeFraction)) +
geom_bar(width=0.75,stat='identity', position=position_dodge(width=0.9)) +
scale_fill_manual(values={sizeFrac_Rpalette}) +
geom_hline(aes(yintercept=1), linetype='dashed', alpha=0.7)+
facet_grid(sizeFraction ~ capDesign, scales='free') +
geom_text(aes(group=sizeFraction, y=0.01, label = paste(sep='',percent(percentMappedReads),' / ','(',comma(mappedReads),')')), angle=90, size=5, hjust=0, vjust=0.5, position = position_dodge(width=0.9)) +
scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
xlab ('Sample (barcode)') +
theme_bw(base_size=17) +
{GGPLOT_PUB_QUALITY} + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('{output}', width=13, height=9)
" > {output}.r
cat {output}.r | R --slave
dropbox_uploader.sh upload {output} {DROPBOX_PLOTS};

		'''

# rule test:
# 	input: returnCapDesignBarcodesFastqs
# #config["DEMULTIPLEX_DIR"] + "{capDesign}_{sizeFrac}.{barcodesU}.fastq"
# 	output: "{capDesign}_{sizeFrac}.{barcodesU}.test" if returnCapDesignBarcodesFastqs else "./tmp/{capDesign}_{sizeFrac}.{barcodesU}.test"
# 	shell:
# 		'''
# head {input} > {output}
# 		'''






###########################
##### run BAMQC
###############################
