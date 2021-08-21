rule basicFASTQqc:
	input: "fastqs/" + "{techname}_{capDesign}_{sizeFrac}_{barcodes}.fastq.gz"
	output: "output/fastqs/" + "qc/{techname}_{capDesign}_{sizeFrac}.{barcodes}.dupl.txt"
	shell:
		'''

# check that read IDs don't contain spaces (frequent with ONT)
#zcat {input} | perl -lane 'if ($F[3]) {{die}}' 

# check that there are no read ID duplicates
zcat {input} | fastq2tsv.pl | awk '{{print $1}}' | sort -T {config[TMPDIR]} | uniq -dc > {output}
count=$(cat {output} | wc -l)
if [ $count -gt 0 ]; then echo "$count duplicate read IDs found"; mv {output} {output}.tmp; exit 1; fi
		'''

rule fastqTimestamps:
	input: expand("fastqs/" + "{techname}_{capDesign}_{sizeFrac}_{barcodes}.fastq.gz", filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES)

	output: "output/statsFiles/" + "all.fastq.timestamps.tsv"
	shell:
		'''
uuid=$(uuidgen)
echo -e "sample_name\\tFASTQ_modified" > {config[TMPDIR]}/$uuid
for file in $(echo {input}); do
echo -e "$(basename $file .fastq.gz)\t$(stat -L --printf='%y' $file  | awk '{{print $1}}')"
done | sort >> {config[TMPDIR]}/$uuid
mv {config[TMPDIR]}/$uuid {output}

		'''

#get read lengths for all FASTQ files:
rule getReadLengthSummary:
	input: "fastqs/" + "{techname}_{capDesign}_{sizeFrac}_{barcodes}.fastq.gz"
	output: 
		reads="output/statsFiles/" + "tmp/{techname}_{capDesign}_{sizeFrac}.{barcodes}.readlength.tsv.gz",
		summ="output/statsFiles/" + "tmp/{techname}_{capDesign}_{sizeFrac}.{barcodes}.readlengthSummary.tsv"
	conda: "envs/R_env.yml"
	params:
		bc=lambda wildcards: wildcards.barcodes
	shell:
		'''
uuidTmpOut=$(uuidgen)
echo -e "seqTech\\tcapDesign\\tsizeFrac\\ttissue\\tlength" |gzip > {config[TMPDIR]}/$uuidTmpOut.gz

zcat {input} | fastq2tsv.pl | perl -F"\\t" -slane '$F[0]=~s/^(\S+).*/$1/; print join("\\t", @F)' | awk -v t={wildcards.techname} -v c={wildcards.capDesign} -v si={wildcards.sizeFrac} -v b={params.bc} '{{print t\"\\t\"c\"\\t\"si\"\\t\"b\"\\t\"length($2)}}'| gzip >> {config[TMPDIR]}/$uuidTmpOut.gz

echo "
library(data.table)
library(tidyverse)
dat<-fread('{config[TMPDIR]}/$uuidTmpOut.gz', header=T, sep='\\t')
dat %>%
  group_by(seqTech, sizeFrac, capDesign, tissue) %>%
  summarise(n=n(), median=median(length), mean=mean(length), max=max(length)) -> datSumm

write_tsv(datSumm, '{config[TMPDIR]}/$uuidTmpOut.1')

" | R --slave

mv {config[TMPDIR]}/$uuidTmpOut.gz {output.reads}
mv {config[TMPDIR]}/$uuidTmpOut.1 {output.summ}
		'''

rule aggReadLengthSummary:
	input: lambda wildcards: expand("output/statsFiles/" + "tmp/{techname}_{capDesign}_{sizeFrac}.{barcodes}.readlengthSummary.tsv", filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES)
	output: 
		summary="output/statsFiles/" + "all.readlength.summary.tsv"
	shell:
		'''
uuid=$(uuidgen)
head -n1 {input[0]} > {config[TMPDIR]}/$uuid
tail -q -n+2 {input} |sort >> {config[TMPDIR]}/$uuid
mv {config[TMPDIR]}/$uuid {output}

		'''

# plot histograms with R:
rule plotReadLength:
	input: "output/statsFiles/" + "tmp/{techname}_{capDesign}_{sizeFrac}.{barcodes}.readlength.tsv.gz"
	output: returnPlotFilenames("output/plots/" + "readLength.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{barcodes}.readLength.stats")
	conda: "envs/R_env.yml"
	params:
		filterDat=lambda wildcards: multi_figures(wildcards.capDesign, wildcards.sizeFrac, wildcards.barcodes, wildcards.techname)
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
{params.filterDat[technameFilterString]}
{params.filterDat[capDesignFilterString]}

{params.filterDat[sizeFracFilterString]}
{params.filterDat[tissueFilterString]}
{params.filterDat[substSeqTechString]}
{params.filterDat[substTissueString]}
{params.filterDat[graphDimensions]}

wXyPlot = wXyPlot * 1.2

dat\$sizeFrac_f=factor(dat\$sizeFrac, levels=names({sizeFrac_Rpalette}), ordered=TRUE)
dat %>%
  group_by(seqTech, sizeFrac_f, capDesign, tissue) %>%
  summarise(n=n(), med=median(length)) -> datSumm


summaryStats = transform(datSumm, LabelN = paste0('N= ', comma(n)), LabelM = paste0( 'Median= ', comma(med))) 

plotBase <- \\"p <- ggplot(dat, aes(x=length)) +
geom_histogram(aes(y=..density..,fill=sizeFrac_f), binwidth=100) +
geom_vline(data = summaryStats, aes(xintercept=med), color='#ff0055', linetype='solid', size=lineSize) +
scale_fill_manual(values={sizeFrac_Rpalette}) +
geom_text(data = summaryStats, aes(label = LabelN, x = Inf, y = Inf), hjust=1, vjust=1,  size=geom_textSize, fontface = 'bold') +
geom_text(data = summaryStats, aes(label = LabelM, x = med, y = Inf), hjust=-0.1, vjust=2.5,  size=geom_textSize, fontface = 'bold', color='#ff0055') +
coord_cartesian(xlim=c(0, 3500)) +
#scale_y_continuous(labels=scientific)+
scale_x_continuous(labels=comma, name='Read length (nts)')+
{GGPLOT_PUB_QUALITY} + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + \\"

{params.filterDat[facetPlotSetup]}

wYxPlot = wYxPlot * 1.2
wYxNoLegendPlot<- wYxPlot - wLegendOnly

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

cat $(dirname {output[0]})/$(basename {output[0]} .legendOnly.png).r | R --slave


	 	'''
