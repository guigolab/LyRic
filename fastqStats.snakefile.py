rule checkNoDuplicateReadIds:
	input: config["FQPATH"] + "{techname}_{capDesign}_{sizeFrac}.fastq"
	output: config["FQPATH"] + "qc/{techname}_{capDesign}_{sizeFrac}.dupl.txt"
	shell:
		'''
cat {input} | fastq2tsv.pl | cut -f1 | sort| uniq -dc > {output}
count=$(cat {output} | wc -l)
if [ $count -gt 0 ]; then echo "$count duplicate read IDs found"; mv {output} {output}.tmp; exit 1; fi
		'''


rule getFastqReadCounts:
	input: config["FQPATH"] + "{techname}_{capDesign}_{sizeFrac}.fastq"
	output: temp(config["STATSDATADIR"] + "{techname}_{capDesign}_{sizeFrac}.fastq.readCounts.tsv")
	shell:
		'''
let total=$(cat {input} | wc -l)/4 || true # "true" serves to avoid "let" exiting with status > 0 when its output is = 0
echo -e "$(basename {input})\t$total" > {output}
		'''

rule aggFastqReadCounts:
	input: lambda wildcards: expand (config["STATSDATADIR"] + "{techname}_{capDesign}_{sizeFrac}.fastq.readCounts.tsv", filtered_product, techname=wildcards.techname, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS)
#	input: lambda wildcards: expand("{techname}_{capDesign}_{sizeFrac}.readlength.tsv", filtered_product, techname=wildcards.techname, capDesign=wildcards.capDesign, sizeFrac=SIZEFRACS)
	output: config["STATSDATADIR"] + "{techname}.fastq.readCounts.tsv"
	shell: "cat {input} | sort > {output}"

#get read lengths for all FASTQ files:
rule getReadLength:
	input: config["FQPATH"] + "{techname}_{capDesign}_{sizeFrac}.fastq"
	output: temp("{techname}_{capDesign}_{sizeFrac}.readlength.tsv")
	shell: "fastq2tsv.pl {input} | awk -v s={wildcards.sizeFrac} '{{print s\"\\t\"length($2)}}' > {output}"

#aggregate read length data over all fractions of a given capDesign:
rule aggReadLength:
	input: lambda wildcards: expand("{techname}_{capDesign}_{sizeFrac}.readlength.tsv", filtered_product, techname=wildcards.techname, capDesign=wildcards.capDesign, sizeFrac=SIZEFRACS)
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

