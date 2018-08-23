rule checkNoDuplicateReadIds:
	input: FQ_CORR_PATH + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.fastq.gz" if config["DEMULTIPLEX"] else  FQ_CORR_PATH + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.{barcodes}.fastq.gz"
	output: FQ_CORR_PATH + "qc/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.dupl.txt" if config["DEMULTIPLEX"] else  FQ_CORR_PATH + "qc/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.{barcodes}.dupl.txt"
	shell:
		'''
zcat {input} | fastq2tsv.pl | cut -f1 | sort| uniq -dc > {output}
count=$(cat {output} | wc -l)
if [ $count -gt 0 ]; then echo "$count duplicate read IDs found"; mv {output} {output}.tmp; exit 1; fi
		'''


rule getFastqReadCounts:
	input: FQ_CORR_PATH + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.fastq.gz" if config["DEMULTIPLEX"] else FQ_CORR_PATH + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.{barcodes}.fastq.gz"
	output: config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.fastq.readCounts.tsv" if config["DEMULTIPLEX"] else config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.{barcodes}.fastq.readCounts.tsv"
	shell:
		'''
let total=$(zcat {input} | wc -l)/4 || true # "true" serves to avoid "let" exiting with status > 0 when its output is = 0
echo -e "$(basename {input})\t$total" > {output}
		'''

rule aggFastqReadCounts:
	input: lambda wildcards: expand (config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.fastq.readCounts.tsv", filtered_product, techname=wildcards.techname, corrLevel={wildcards.corrLevel}, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS) if config["DEMULTIPLEX"] else expand (config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.{barcodes}.fastq.readCounts.tsv", filtered_product, techname=wildcards.techname, corrLevel={wildcards.corrLevel}, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=wildcards.barcodes)
	output: config["STATSDATADIR"] + "{techname}Corr{corrLevel}.fastq.readCounts.tsv" if config["DEMULTIPLEX"] else config["STATSDATADIR"] + "{techname}Corr{corrLevel}.{barcodes}.fastq.readCounts.tsv"
	shell: "cat {input} | sort > {output}"

#get read lengths for all FASTQ files:
rule getReadLength:
	input: FQ_CORR_PATH + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.fastq.gz" if config["DEMULTIPLEX"] else FQ_CORR_PATH + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.{barcodes}.fastq.gz"
	output: config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.readlength.tsv" if config["DEMULTIPLEX"] else config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.{barcodes}.readlength.tsv"
	shell: "zcat {input} | fastq2tsv.pl | awk -v s={wildcards.sizeFrac} '{{print s\"\\t\"length($2)}}' > {output}"

#aggregate read length data over all fractions of a given capDesign:
rule aggReadLength:
	input: lambda wildcards: expand(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.readlength.tsv", filtered_product, techname=wildcards.techname, corrLevel={wildcards.corrLevel}, capDesign=wildcards.capDesign, sizeFrac=SIZEFRACS) if config["DEMULTIPLEX"] else expand(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.{barcodes}.readlength.tsv", filtered_product, techname=wildcards.techname, corrLevel={wildcards.corrLevel}, capDesign=wildcards.capDesign, sizeFrac=SIZEFRACS, barcodes=wildcards.barcodes)
	#input: glob.glob(os.path.join("{capDesign}_*.readlength.tsv"))
	output: config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_all.readlength.tsv" if config["DEMULTIPLEX"] else config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}.{barcodes}_all.readlength.tsv"
	shell: "cat {input} > {output}"

# plot histograms with R:
rule plotReadLength:
	input: config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_all.readlength.tsv" if config["DEMULTIPLEX"] else config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}.{barcodes}_all.readlength.tsv"
	output: config["PLOTSDIR"] + "{techname}Corr{corrLevel}_{capDesign}_all.readlength.{ext}" if config["DEMULTIPLEX"] else config["PLOTSDIR"] + "{techname}Corr{corrLevel}_{capDesign}.{barcodes}_all.readlength.{ext}"
	shell:
		'''
echo "library(ggplot2)
library(plyr)
dat <- read.table('{input}', header=F, as.is=T, sep='\\t')
colnames(dat)<-c('fraction','length')
dat\$fraction_f=factor(dat\$fraction, levels=names({sizeFrac_Rpalette}), ordered=TRUE)
medians <- ddply(dat, 'fraction_f', summarise, grp.median=median(length))
sizes <- ddply(dat, 'fraction_f', summarise, grp.length=length(length))
ggplot(dat, aes(x=length)) +
geom_histogram(binwidth=200, aes(fill=fraction_f)) +
scale_fill_manual(values={sizeFrac_Rpalette}) +
#facet_grid( . ~ names({sizeFrac_Rpalette})) +
facet_grid( . ~ fraction_f) +
annotate(geom='text', label=c(paste('n= ', sizes[,'grp.length'])), hjust=0, vjust=c(3.8), x=10, y=c(Inf), color=c('black'), size=6) +
annotate(geom='text', label=c(paste('Median (nt)= ', medians[,'grp.median'])), hjust=0, vjust=c(5.8), x=10, y=c(Inf), color=c('black'), size=6) +
coord_cartesian(xlim=c(0, 5000)) +
theme_bw(base_size=17) +
{GGPLOT_PUB_QUALITY} + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('{output}', width=10, height=3)
" > {output}.r
cat {output}.r | R --slave

	 	'''

