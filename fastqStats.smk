rule checkNoDuplicateReadIds:
	input: FQ_CORR_PATH + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.fastq.gz" if config["DEMULTIPLEX"] else  FQ_CORR_PATH + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.{barcodes}.fastq.gz"
	output: FQ_CORR_PATH + "qc/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.dupl.txt" if config["DEMULTIPLEX"] else  FQ_CORR_PATH + "qc/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.{barcodes}.dupl.txt"
	shell:
		'''
zcat {input} | fastq2tsv.pl | cut -f1 | sort -T {config[TMPDIR]} | uniq -dc > {output}
count=$(cat {output} | wc -l)
if [ $count -gt 0 ]; then echo "$count duplicate read IDs found"; mv {output} {output}.tmp; exit 1; fi
		'''



#get read lengths for all FASTQ files:
rule getReadLength:
	input: FQ_CORR_PATH + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.fastq.gz" if config["DEMULTIPLEX"] else FQ_CORR_PATH + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.{barcodes}.fastq.gz"
	output: temp(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.readlength.tsv") if config["DEMULTIPLEX"] else temp(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.{barcodes}.readlength.tsv")
	params:
		bc=lambda wildcards: 'allTissues' if config["DEMULTIPLEX"] else wildcards.barcodes
	shell:
		'''
zcat {input} | fastq2tsv.pl | awk -v t={wildcards.techname}Corr{wildcards.corrLevel} -v c={wildcards.capDesign} -v si={wildcards.sizeFrac} -v b={params.bc} '{{print t\"\\t\"c\"\\t\"si\"\\t\"b\"\\t\"length($2)}}'| sed 's/Corr0/\tNo/' | sed 's/Corr{lastK}/\tYes/' > {output}
		'''

#aggregate read length data over all fractions of a given capDesign:
rule aggReadLength:
	input: lambda wildcards: expand(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.readlength.tsv", filtered_product, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS) if config["DEMULTIPLEX"] else expand(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.{barcodes}.readlength.tsv", filtered_product, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES)
	#input: glob.glob(os.path.join("{capDesign}_*.readlength.tsv"))
	output: config["STATSDATADIR"] + "all.readlength.tsv"
	shell:
		'''
echo -e "seqTech\tcorrectionLevel\tcapDesign\tsizeFrac\ttissue\tlength" > {output}
cat {input} |sort -T {config[TMPDIR]}  >> {output}

		'''
#rule addAllReadLength:
#	input: lambda wildcards:


# plot histograms with R:
rule plotReadLength:
	input: config["STATSDATADIR"] + "all.readlength.tsv"
	output: config["PLOTSDIR"] + "readLength.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.readLength.stats.{ext}"  if config["DEMULTIPLEX"] else config["PLOTSDIR"] + "readLength.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.readLength.stats.{ext}"
	params:
		filterDat=lambda wildcards: merge_figures_params(wildcards.capDesign, wildcards.sizeFrac, 'allTissues', wildcards.corrLevel, wildcards.techname) if config["DEMULTIPLEX"] else merge_figures_params(wildcards.capDesign, wildcards.sizeFrac, wildcards.barcodes, wildcards.corrLevel, wildcards.techname)
	shell:
		'''
echo "library(ggplot2)
library(dplyr)
library(scales)
library(data.table)
dat<-fread('{input}', header=T, sep='\\t')
{params.filterDat[10]}
{params.filterDat[0]}
{params.filterDat[1]}
{params.filterDat[2]}
{params.filterDat[3]}
{params.filterDat[4]}
{params.filterDat[5]}
{params.filterDat[8]}

dat\$sizeFrac_f=factor(dat\$sizeFrac, levels=names({sizeFrac_Rpalette}), ordered=TRUE)
summaryStats = transform(summarise(group_by(dat, seqTech, sizeFrac_f, capDesign, tissue), Label = paste0('N= ', comma(length(length)), '\\n', 'Median= ', comma(median(length)))))

ggplot(dat, aes(x=length)) +
geom_histogram(aes(y=..density..,fill=sizeFrac_f), binwidth=200) +
scale_fill_manual(values={sizeFrac_Rpalette}) +
facet_grid( seqTech + sizeFrac_f ~ capDesign + tissue, scales='free_y') +
geom_text(data = summaryStats, aes(label = Label, x = 10, y = Inf), hjust=0, vjust=3.8,  size=geom_textSize) +

coord_cartesian(xlim=c(0, 3500)) +
scale_y_continuous(labels=scientific)+
scale_x_continuous(labels=comma)+
#theme_bw(base_size=17) +
{GGPLOT_PUB_QUALITY} + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('{output}', width=plotWidth, height=plotHeight)
" > {output}.r
cat {output}.r | R --slave

	 	'''
