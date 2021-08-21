rule makeIntrons:
	input: "output/mappings/readBedToGff/{techname}_{capDesign}_{sizeFrac}_{barcodes}.gff.gz"
	output: "output/mappings/makeIntrons/{techname}_{capDesign}_{sizeFrac}_{barcodes}.introns.gff.gz"
	shell:
		'''
uuidTmpOut=$(uuidgen)
zcat {input} | makeIntrons.pl - | sort -T {config[TMPDIR]}  -k1,1 -k4,4n -k5,5n  |gzip> {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
		'''

rule getIntronMotif:
	input:
		introns = "output/mappings/makeIntrons/{techname}_{capDesign}_{sizeFrac}_{barcodes}.introns.gff.gz",
		genome = lambda wildcards: config["GENOMESDIR"] + CAPDESIGNTOGENOME[wildcards.capDesign] + ".sorted.fa"
	output:
		gff = "output/mappings/getIntronMotif/{techname}_{capDesign}_{sizeFrac}_{barcodes}.introns.gff.gz",
		tsv = "output/mappings/getIntronMotif/{techname}_{capDesign}_{sizeFrac}_{barcodes}.transcripts.tsv.gz"
	conda: "envs/perl_env.yml"
	shell:
		'''
uuid=$(uuidgen)
mkdir -p {config[TMPDIR]}/$uuid

zcat {input.introns} | grep -vP "^ERCC"| extract_intron_strand_motif.pl - {input.genome} {config[TMPDIR]}/$uuid/$(basename {output.gff} .introns.gff.gz)

gzip {config[TMPDIR]}/$uuid/*
mv {config[TMPDIR]}/$uuid/* $(dirname {output.gff})
		'''
		
rule getGencodeSpliceJunctions:
	input: lambda wildcards: GENOMETOANNOTGTF[CAPDESIGNTOGENOME[wildcards.capDesign]]
	output: "output/annotations/spliceJunctions/{capDesign}.gencode.spliceJunctions.list"
	shell:
		'''
uuidTmpOut=$(uuidgen)
cat {input} | awk '$3=="exon"' |sort -T {config[TMPDIR]}  -k1,1 -k4,4n -k5,5n | makeIntrons.pl - | awk '{{print $1"_"$4"_"$5"_"$7}}' |sort -T {config[TMPDIR]} |uniq > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}

		'''

rule getClsSpliceJunctions:
	input:"output/mappings/mergedReads/{techname}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.tmerge.min{minReadSupport}reads.splicing_status:spliced.endSupport:all.gff.gz"
	output: "output/mappings/mergedReads/spliceJunctions/{techname}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.spliceJunctions.list"
	shell:
		'''
uuidTmpOut=$(uuidgen)
zcat {input} | awk '$3=="exon"' |sort -T {config[TMPDIR]}  -k1,1 -k4,4n -k5,5n | makeIntrons.pl - | awk '{{print $1"_"$4"_"$5"_"$7}}' |sort -T {config[TMPDIR]} |uniq > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}

		'''

rule getCompareClsGencodeSJsStats:
	input:
		gencodeSJs="output/annotations/spliceJunctions/{capDesign}.gencode.spliceJunctions.list",
		clsSJs="output/mappings/mergedReads/spliceJunctions/{techname}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.spliceJunctions.list"
	output: "output/statsFiles/" + "tmp/{techname}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.vs.Gencode.SJs.stats.tsv"
	shell:
		'''
uuidTmpOut=$(uuidgen)
#annSJs=$(cat {input.gencodeSJs} | wc -l)
clsSJs=$(cat {input.clsSJs} | wc -l)
commonSJs=$(comm -1 -2 {input.gencodeSJs} {input.clsSJs} | wc -l)
novelSJs=$(comm -1 -3 {input.gencodeSJs} {input.clsSJs} | wc -l)
echo -e "{wildcards.techname}\t{wildcards.capDesign}\t{wildcards.sizeFrac}\t{wildcards.barcodes}\t$clsSJs\t$commonSJs\t$novelSJs"  > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}

		'''

rule aggCompareClsGencodeSJsStats:
	input: lambda wildcards: expand("output/statsFiles/" +"tmp/{techname}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.vs.Gencode.SJs.stats.tsv",filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES, minReadSupport=wildcards.minReadSupport)
	output: "output/statsFiles/" + "all.tmerge.min{minReadSupport}reads.vs.Gencode.SJs.stats.tsv"
	shell:
		'''
uuidTmpOut=$(uuidgen)
echo -e "seqTech\tcapDesign\tsizeFrac\ttissue\tcategory\tcount\tpercent" > {config[TMPDIR]}/$uuidTmpOut
cat {input} | awk '{{ print $1"\\t"$2"\\t"$3"\\t"$4"\\tcommon\\t"$6"\t"$6/$5"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\tnovel\\t"$7"\t"$7/$5}}' | sort -T {config[TMPDIR]}  >> {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}


		'''

rule plotCompareClsGencodeSJsStats:
	input: "output/statsFiles/" + "all.tmerge.min{minReadSupport}reads.vs.Gencode.SJs.stats.tsv"
	output: returnPlotFilenames("output/plots/" + "tmerge.vs.Gencode.SJs.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.vs.Gencode.SJs.stats")
	conda: "envs/R_env.yml"
	params:
		filterDat=lambda wildcards: multi_figures(wildcards.capDesign, wildcards.sizeFrac, wildcards.barcodes, wildcards.techname)
	shell:
		'''
echo "
library(ggplot2)

library(cowplot)
library(plyr)
library(scales)
library(gridExtra)
library(grid)
library(ggplotify)
library(data.table)

dat <- read.table('{input}', header=T, as.is=T, sep='\\t')
{params.filterDat[technameFilterString]}
{params.filterDat[capDesignFilterString]}

{params.filterDat[sizeFracFilterString]}
{params.filterDat[tissueFilterString]}
{params.filterDat[substSeqTechString]}
{params.filterDat[substTissueString]}
{params.filterDat[graphDimensions]}

dat\$category<-factor(dat\$category, ordered=TRUE, levels=rev(c('common', 'novel')))

plotBase <- \\"p <- ggplot(dat[order(dat\$category), ], aes(x=1, y=count, fill=category)) +
geom_bar(stat='identity') +
scale_fill_manual (values=c(annOnly='#7D96A2',common='#83A458', novel='#B8CF7E'), labels=c(annOnly='Only in GENCODE', common='In sample+GENCODE', novel='Only in sample')) +
ylab('# Splice Junctions')+
geom_text(position = 'stack', size=geom_textSize, aes(y = count, label = paste(sep='',percent(round(percent, digits=2)),' / ','(',comma(count),')'), hjust = 0.5, vjust = 1))+
{params.filterDat[hideXaxisLabels]}
{GGPLOT_PUB_QUALITY} + 
theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
\\"

{params.filterDat[facetPlotSetup]}

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
