rule basicFASTQqc:
	input: FQ_CORR_PATH + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.fastq.gz" if config["DEMULTIPLEX"] else  FQ_CORR_PATH + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.fastq.gz"
	output: FQ_CORR_PATH + "qc/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.dupl.txt" if config["DEMULTIPLEX"] else  FQ_CORR_PATH + "qc/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.{barcodes}.dupl.txt"
	shell:
		'''

# check that read IDs don't contain spaces (frequent with ONT)
zcat {input} | perl -lane 'if ($F[3]) {{die}}' 

# check that there are no read ID duplicates
zcat {input} | fastq2tsv.pl | cut -f1 | sort -T {config[TMPDIR]} | uniq -dc > {output}
count=$(cat {output} | wc -l)
if [ $count -gt 0 ]; then echo "$count duplicate read IDs found"; mv {output} {output}.tmp; exit 1; fi
		'''


#get read lengths for all FASTQ files:
rule getReadLength:
	input: FQ_CORR_PATH + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.fastq.gz" if config["DEMULTIPLEX"] else FQ_CORR_PATH + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.fastq.gz"
	output: config["STATSDATADIR"] + "tmp/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.readlength.tsv" if config["DEMULTIPLEX"] else config["STATSDATADIR"] + "tmp/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.{barcodes}.readlength.tsv"
	params:
		bc=lambda wildcards: 'allTissues' if config["DEMULTIPLEX"] else wildcards.barcodes
	shell:
		'''
uuidTmpOut=$(uuidgen)
zcat {input} | fastq2tsv.pl | awk -v t={wildcards.techname}Corr{wildcards.corrLevel} -v c={wildcards.capDesign} -v si={wildcards.sizeFrac} -v b={params.bc} '{{print t\"\\t\"c\"\\t\"si\"\\t\"b\"\\t\"length($2)}}'| sed 's/Corr0/\tNo/' | sed 's/Corr{lastK}/\tYes/' > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
		'''

#aggregate read length data over all fractions of a given capDesign:
rule aggReadLength:
	input: lambda wildcards: expand(config["STATSDATADIR"] + "tmp/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.readlength.tsv", filtered_product, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS) if config["DEMULTIPLEX"] else expand(config["STATSDATADIR"] + "tmp/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.{barcodes}.readlength.tsv", filtered_product, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=wildcards.capDesign, sizeFrac=SIZEFRACS, barcodes=BARCODES)
	#input: glob.glob(os.path.join("{capDesign}_*.readlength.tsv"))
	output: config["STATSDATADIR"] + "all.{capDesign}.readlength.tsv"
	shell:
		'''
uuidTmpOut=$(uuidgen)
echo -e "seqTech\tcorrectionLevel\tcapDesign\tsizeFrac\ttissue\tlength" > {config[TMPDIR]}/$uuidTmpOut
cat {input} |sort -T {config[TMPDIR]}  >> {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}

		'''
#rule addAllReadLength:
#	input: lambda wildcards:


# plot histograms with R:
rule plotReadLength:
	input: config["STATSDATADIR"] + "all.{capDesign}.readlength.tsv"
	output: returnPlotFilenames(config["PLOTSDIR"] + "readLength.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.readLength.stats")  if config["DEMULTIPLEX"] else returnPlotFilenames(config["PLOTSDIR"] + "readLength.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.readLength.stats")
	params:
		filterDat=lambda wildcards: merge_figures_params(wildcards.capDesign, wildcards.sizeFrac, 'allTissues', wildcards.corrLevel, wildcards.techname) if config["DEMULTIPLEX"] else merge_figures_params(wildcards.capDesign, wildcards.sizeFrac, wildcards.barcodes, wildcards.corrLevel, wildcards.techname)
	shell:
		'''
echo "
library(ggplot2)
library(cowplot)
library(scales)
library(gridExtra)
library(grid)
library(ggplotify)
library(dplyr)
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
dat %>%
  group_by(seqTech, sizeFrac_f, capDesign, tissue) %>%
  summarise(n=n(), med=median(length)) -> datSumm


summaryStats = transform(datSumm, LabelN = paste0('N= ', comma(n)), LabelM = paste0( 'Median= ', comma(med))) 

plotBase <- \\"p <- ggplot(dat, aes(x=length)) +
geom_histogram(aes(y=..density..,fill=sizeFrac_f), binwidth=100) +
geom_vline(data = summaryStats, aes(xintercept=med), color='#ff0055', linetype='solid', size=1.2) +
scale_fill_manual(values={sizeFrac_Rpalette}) +
geom_text(data = summaryStats, aes(label = LabelN, x = Inf, y = Inf), hjust=1, vjust=1,  size=geom_textSize, fontface = 'bold') +
geom_text(data = summaryStats, aes(label = LabelM, x = med, y = Inf), hjust=-0.1, vjust=2.5,  size=geom_textSize, fontface = 'bold', color='#ff0055') +
coord_cartesian(xlim=c(0, 3500)) +
#scale_y_continuous(labels=scientific)+
scale_x_continuous(labels=comma, name='Read length (nts)')+
{GGPLOT_PUB_QUALITY} + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + \\"

{params.filterDat[12]}

save_plot('{output[0]}', legendOnly, base_width=wLegendOnly, base_height=hLegendOnly)
save_plot('{output[1]}', legendOnly, base_width=wLegendOnly, base_height=hLegendOnly)

save_plot('{output[2]}', pXy, base_width=wXyPlot, base_height=hXyPlot)
save_plot('{output[3]}', pXy, base_width=wXyPlot, base_height=hXyPlot)

save_plot('{output[4]}', pXyNoLegend, base_width=wXyNoLegendPlot, base_height=hXyNoLegendPlot)
save_plot('{output[5]}', pXyNoLegend, base_width=wXyNoLegendPlot, base_height=hXyNoLegendPlot)

save_plot('{output[6]}', pYx, base_width=wYxPlot, base_height=hYxPlot)
save_plot('{output[7]}', pYx, base_width=wYxPlot, base_height=hYxPlot)

save_plot('{output[8]}', pYxNoLegend, base_width=wYxNoLegendPlot, base_height=hYxNoLegendPlot)
save_plot('{output[9]}', pYxNoLegend, base_width=wYxNoLegendPlot, base_height=hYxNoLegendPlot)


" > $(dirname {output[0]})/$(basename {output[0]} .legendOnly.png).r
 set +eu
conda activate R_env
set -eu

cat $(dirname {output[0]})/$(basename {output[0]} .legendOnly.png).r | R --slave


	 	'''
