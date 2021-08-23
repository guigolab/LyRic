rule extractPooledTmsFivepEnds:
	input: "output/mappings/mergedReads/{techname}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.splicing_status:all.endSupport:all.bed"
	output: "output/mappings/mergedReads/5pEnds/{techname}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.5pEnds.bed"
	shell:
		'''
uuidTmpOut=$(uuidgen)
cat {input} |extractTranscriptEndsFromBed12.pl 5 |sort -T {config[TMPDIR]}  -k1,1 -k2,2n -k3,3n > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
		'''

rule cageSupportedfivepEnds:
	input:
		fivePends="output/mappings/mergedReads/5pEnds/{techname}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.5pEnds.bed",
		tms="output/mappings/mergedReads/{techname}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.splicing_status:all.endSupport:all.bed",
		cagePeaks=lambda wildcards: GENOMETOCAGEPEAKS[CAPDESIGNTOGENOME[wildcards.capDesign]],
		genome = lambda wildcards: config["GENOMESDIR"] + CAPDESIGNTOGENOME[wildcards.capDesign] + ".sorted.genome"
	output: "output/mappings/mergedReads/cageSupported/{techname}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.cageSupported5pEnds.bed"
	conda: "envs/xtools_env.yml"
	shell:
		'''
uuidTmpOut=$(uuidgen)


cat {input.fivePends} | sort -T {config[TMPDIR]}  -k1,1 -k2,2n -k3,3n  | bedtools slop -s -l 50 -r 50 -i stdin -g {input.genome} | bedtools intersect -u -s -a stdin -b {input.cagePeaks} | cut -f4 | fgrep -w -f - {input.tms} > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
		'''

rule dhsSupportedfivepEnds:
	input:
		fivePends="output/mappings/mergedReads/5pEnds/{techname}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.5pEnds.bed",
		tms="output/mappings/mergedReads/{techname}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.splicing_status:all.endSupport:all.bed",
		dhsPeaks=lambda wildcards: GENOMETODHSPEAKS[CAPDESIGNTOGENOME[wildcards.capDesign]],
	output: "output/mappings/mergedReads/dhsSupported/{techname}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.dhsSupported5pEnds.bed"
	conda: "envs/xtools_env.yml"
	shell:
		'''
uuidTmpOut=$(uuidgen)


cat {input.fivePends} | sort -T {config[TMPDIR]}  -k1,1 -k2,2n -k3,3n  | bedtools intersect -u -a stdin -b {input.dhsPeaks} | cut -f4 | fgrep -w -f - {input.tms} > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
		'''

rule plotDhsVsCage5primeComparisonStats:
	input: 
		dhs="output/mappings/mergedReads/dhsSupported/{techname}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.dhsSupported5pEnds.bed",
		cage="output/mappings/mergedReads/cageSupported/{techname}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.cageSupported5pEnds.bed"
	output: "output/plots/" + "dhsVsCage5primeComparison.venn.stats/{techname}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.dhsVsCage5primeComparison.venn.stats.pdf"
	conda: "envs/R_env.yml"
	shell:
		'''
uuidTmpCage=$(uuidgen)
uuidTmpDhs=$(uuidgen)
cat {input.cage} | cut -f4 |sort|uniq > {config[TMPDIR]}/$uuidTmpCage
cat {input.dhs} | cut -f4 |sort|uniq > {config[TMPDIR]}/$uuidTmpDhs

area1=$(cat {config[TMPDIR]}/$uuidTmpCage | wc -l)
area2=$(cat {config[TMPDIR]}/$uuidTmpDhs | wc -l)
oneItwo=$(comm -1 -2 {config[TMPDIR]}/$uuidTmpCage {config[TMPDIR]}/$uuidTmpDhs |wc -l)

echo "
pdf(file='{output}', bg = 'white')
library(VennDiagram)
venn.plot <- draw.pairwise.venn(
area1=$area1,
area2=$area2,
cross.area=$oneItwo,
fill = c('#ff8000', '#0086b3'),
cat.pos=c(-20,20),
col='black',
alpha=0.6,
fontfamily = rep('Helvetica', 3),
cat.fontfamily = rep('Helvetica', 2),
cat.fontface = rep('bold.italic', 2),
euler.d=TRUE,
scaled=TRUE,
cex=3,
cat.cex=2,
cat.dist = c(0.05, 0.05),
cat.col= c('#ff8000', '#0086b3'),
lty = rep('blank', 2),
category=c('CAGE-supported', 'DHS-supported')
)
dev.off()
" > {output}.r

cat {output}.r | R --slave

		'''



rule extractPooledTmsThreepEnds:
	input: "output/mappings/mergedReads/{techname}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.splicing_status:all.endSupport:all.bed"
	output: "output/mappings/mergedReads/3pEnds/{techname}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.3pEnds.bed"
	shell:
		'''
uuidTmpOut=$(uuidgen)
cat {input} |extractTranscriptEndsFromBed12.pl 3 |sort -T {config[TMPDIR]}  -k1,1 -k2,2n -k3,3n > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
		'''

rule polyASupportedthreepEnds:
	input:
		threePends="output/mappings/mergedReads/3pEnds/{techname}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.3pEnds.bed",
		tms="output/mappings/mergedReads/{techname}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.splicing_status:all.endSupport:all.bed",
		polyAsites="output/mappings/removePolyAERCCs/{techname}_{capDesign}_{sizeFrac}_{barcodes}.polyAsitesNoErcc.bed",
		genome = lambda wildcards: config["GENOMESDIR"] + CAPDESIGNTOGENOME[wildcards.capDesign] + ".sorted.genome"
	conda: "envs/xtools_env.yml"
	output: "output/mappings/mergedReads/polyASupported/{techname}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.polyASupported3pEnds.bed"
	shell:
		'''
uuid=$(uuidgen)
uuidTmpOut=$(uuidgen)
cat {input.polyAsites} |sort -T {config[TMPDIR]}  -k1,1 -k2,2n -k3,3n  > {config[TMPDIR]}/$uuid.polyAsites.bed


cat {input.threePends} | sort -T {config[TMPDIR]}  -k1,1 -k2,2n -k3,3n  | bedtools slop -s -l 5 -r 5 -i stdin -g {input.genome} | bedtools intersect -u -s -a stdin -b {config[TMPDIR]}/$uuid.polyAsites.bed | cut -f4 | tgrep -F -w -f - {input.tms} > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
		'''

rule getCagePolyASupport:
	input:
		polyA="output/mappings/mergedReads/polyASupported/{techname}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.polyASupported3pEnds.bed",
		cage="output/mappings/mergedReads/cageSupported/{techname}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.cageSupported5pEnds.bed",
		tms="output/mappings/mergedReads/{techname}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.tmerge.min{minReadSupport}reads.splicing_status:{splicedStatus}.endSupport:all.gff.gz"
	output:
		stats="output/statsFiles/" + "tmp/{techname}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.splicing_status:{splicedStatus}.cagePolyASupport.stats.tsv",
		FLbed="output/mappings/mergedReads/{techname}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.splicing_status:{splicedStatus}.endSupport:cagePolyASupported.bed",
		cageOnlyBed="output/mappings/mergedReads/{techname}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.splicing_status:{splicedStatus}.endSupport:cageOnlySupported.bed",
		polyAOnlyBed="output/mappings/mergedReads/{techname}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.splicing_status:{splicedStatus}.endSupport:polyAOnlySupported.bed",
		noCageNoPolyABed="output/mappings/mergedReads/{techname}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.splicing_status:{splicedStatus}.endSupport:noCageNoPolyASupported.bed",
		FLgff="output/mappings/mergedReads/{techname}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.tmerge.min{minReadSupport}reads.splicing_status:{splicedStatus}.endSupport:cagePolyASupported.gff.gz"

	shell:
		'''

uuid=$(uuidgen)
uuidTmpOutS=$(uuidgen)
uuidTmpOutB=$(uuidgen)
uuidTmpOutG=$(uuidgen)
zcat {input.tms} > {config[TMPDIR]}/$uuid.tms
cat {config[TMPDIR]}/$uuid.tms |extractGffAttributeValue.pl transcript_id | sort -T {config[TMPDIR]} |uniq > {config[TMPDIR]}/$uuid.all.list

cat {input.polyA} | tgrep -F -w -f {config[TMPDIR]}/$uuid.all.list | cut -f4 | sort -T {config[TMPDIR]} |uniq > {config[TMPDIR]}/$uuid.polyA.list
cat {input.cage} | tgrep -F -w -f {config[TMPDIR]}/$uuid.all.list | cut -f4 | sort -T {config[TMPDIR]} |uniq > {config[TMPDIR]}/$uuid.cage.list
cat {config[TMPDIR]}/$uuid.polyA.list {config[TMPDIR]}/$uuid.cage.list |sort -T {config[TMPDIR]} |uniq > {config[TMPDIR]}/$uuid.cageOrPolyA.list
comm -1 -2 {config[TMPDIR]}/$uuid.polyA.list {config[TMPDIR]}/$uuid.cage.list |sort -T {config[TMPDIR]} |uniq > {config[TMPDIR]}/$uuid.cage+PolyA.list

comm -2 -3 {config[TMPDIR]}/$uuid.all.list {config[TMPDIR]}/$uuid.cageOrPolyA.list > {config[TMPDIR]}/$uuid.noCageNoPolyA.list
noCageNoPolyA=$(cat {config[TMPDIR]}/$uuid.noCageNoPolyA.list |wc -l)
tgrep -F -w -f {config[TMPDIR]}/$uuid.noCageNoPolyA.list {config[TMPDIR]}/$uuid.tms | gff2bed_full.pl - | sort -T {config[TMPDIR]}  -k1,1 -k2,2n -k3,3n  > {config[TMPDIR]}/$uuid.noCageNoPolyA.bed
mv {config[TMPDIR]}/$uuid.noCageNoPolyA.bed {output.noCageNoPolyABed}

comm -2 -3 {config[TMPDIR]}/$uuid.cage.list {config[TMPDIR]}/$uuid.polyA.list > {config[TMPDIR]}/$uuid.cageOnly.list
cageOnly=$(cat {config[TMPDIR]}/$uuid.cageOnly.list |wc -l)
tgrep -F -w -f {config[TMPDIR]}/$uuid.cageOnly.list {config[TMPDIR]}/$uuid.tms | gff2bed_full.pl - | sort -T {config[TMPDIR]}  -k1,1 -k2,2n -k3,3n  > {config[TMPDIR]}/$uuid.cageOnly.bed
mv {config[TMPDIR]}/$uuid.cageOnly.bed {output.cageOnlyBed}

comm -2 -3 {config[TMPDIR]}/$uuid.polyA.list {config[TMPDIR]}/$uuid.cage.list > {config[TMPDIR]}/$uuid.polyAOnly.list
polyAOnly=$(cat {config[TMPDIR]}/$uuid.polyAOnly.list |wc -l)
tgrep -F -w -f {config[TMPDIR]}/$uuid.polyAOnly.list {config[TMPDIR]}/$uuid.tms | gff2bed_full.pl - | sort -T {config[TMPDIR]}  -k1,1 -k2,2n -k3,3n  > {config[TMPDIR]}/$uuid.polyAOnly.bed
mv {config[TMPDIR]}/$uuid.polyAOnly.bed {output.polyAOnlyBed}

cageAndPolyA=$(cat {config[TMPDIR]}/$uuid.cage+PolyA.list | wc -l)
let total=$noCageNoPolyA+$cageOnly+$polyAOnly+$cageAndPolyA
tgrep -F -w -f {config[TMPDIR]}/$uuid.cage+PolyA.list {config[TMPDIR]}/$uuid.tms |sort -T {config[TMPDIR]}  -k1,1 -k4,4n -k5,5n  > {config[TMPDIR]}/$uuidTmpOutG
cat {config[TMPDIR]}/$uuidTmpOutG | gzip > {output.FLgff}
cat {config[TMPDIR]}/$uuidTmpOutG |gff2bed_full.pl - | sort -T {config[TMPDIR]}  -k1,1 -k2,2n -k3,3n  > {config[TMPDIR]}/$uuidTmpOutB
mv {config[TMPDIR]}/$uuidTmpOutB {output.FLbed}
echo -e "{wildcards.techname}\t{wildcards.capDesign}\t{wildcards.sizeFrac}\t{wildcards.barcodes}\t$total\t$cageOnly\t$cageAndPolyA\t$polyAOnly\t$noCageNoPolyA"  > {config[TMPDIR]}/$uuidTmpOutS
mv {config[TMPDIR]}/$uuidTmpOutS {output.stats}
		'''


rule aggCagePolyAStats:
	input: lambda wildcards: expand("output/statsFiles/" + "tmp/{techname}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.splicing_status:{splicedStatus}.cagePolyASupport.stats.tsv",filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES, minReadSupport=wildcards.minReadSupport, splicedStatus=wildcards.splicedStatus)
	output: "output/statsFiles/" + "all.min{minReadSupport}reads.splicing_status:{splicedStatus}.cagePolyASupport.stats.tsv"
	shell:
		'''
uuidTmpOut=$(uuidgen)
echo -e "seqTech\tcapDesign\tsizeFrac\ttissue\tcategory\tcount\tpercent" > {config[TMPDIR]}/$uuidTmpOut

cat {input} | awk '{{print $1"\\t"$2"\\t"$3"\\t"$4"\\tcageOnly\\t"$6"\\t"$6/$5"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\tcageAndPolyA\\t"$7"\\t"$7/$5"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\tpolyAOnly\\t"$8"\\t"$8/$5"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\tnoCageNoPolyA\\t"$9"\\t"$9/$5}}'  | sort -T {config[TMPDIR]}  >> {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
		'''

rule plotCagePolyAStats:
	input: "output/statsFiles/" + "all.min{minReadSupport}reads.splicing_status:{splicedStatus}.cagePolyASupport.stats.tsv"
	output: returnPlotFilenames("output/plots/" + "cagePolyASupport.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.splicing_status:{splicedStatus}.cagePolyASupport.stats")
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

dat <- read.table('{input}', header=T, as.is=T, sep='\\t')
{params.filterDat[technameFilterString]}
{params.filterDat[capDesignFilterString]}

{params.filterDat[sizeFracFilterString]}
{params.filterDat[tissueFilterString]}
{params.filterDat[substSeqTechString]}
{params.filterDat[substTissueString]}
{params.filterDat[graphDimensions]}

dat\$category<-factor(dat\$category, ordered=TRUE, levels=rev(c('cageOnly', 'cageAndPolyA', 'polyAOnly', 'noCageNoPolyA')))
plotBase <- \\"p <- ggplot(dat[order(dat\$category), ], aes(x=1, y=count, fill=category)) +
geom_bar(stat='identity') +
ylab('# TMs') +
scale_y_continuous(labels=comma)+
scale_fill_manual (values=c(cageOnly='#e5b3e5', cageAndPolyA='#C453C4', polyAOnly = '#b3e0ff', noCageNoPolyA='#a6a6a6'), labels = c(cageOnly = '5´-complete only', cageAndPolyA = '5´+3´-complete', polyAOnly = '3´-complete only', noCageNoPolyA = '5´+3´-incomplete'))+
xlab('') +
guides(fill = guide_legend(title='Category'))+
#geom_text(position = 'stack', size=geom_textSize, aes( y = count, label = paste(sep='',percent(round(percent, digits=2)),' / ','(',comma(count),')'), hjust = 0.5, vjust = 1))+

{GGPLOT_PUB_QUALITY} + 
theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
\\"

{params.filterDat[facetPlotSetup]}

save_plot('{output[0]}', legendOnly, base_width=wLegendOnly, base_height=hLegendOnly)
save_plot('{output[1]}', legendOnly, base_width=wLegendOnly, base_height=hLegendOnly)

save_plot('{output[2]}', pXy, base_width=wXyPlot, base_height=hXyPlot)
save_plot('{output[3]}', pXy, base_width=wXyPlot, base_height=hXyPlot)

save_plot('{output[4]}', pXyNoLegend, base_width=wXyNoLegendPlot, base_height=hXyNoLegendPlot)


" > $(dirname {output[0]})/$(basename {output[0]} .legendOnly.png).r

cat $(dirname {output[0]})/$(basename {output[0]} .legendOnly.png).r | R --slave

		'''


