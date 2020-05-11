rule extractPooledTmsFivepEnds:
	input: "mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.splicing_status:all.endSupport:all.bed"
	output: "mappings/nonAnchoredMergeReads/5pEnds/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.5pEnds.bed"
	shell:
		'''
uuidTmpOut=$(uuidgen)
cat {input} |extractTranscriptEndsFromBed12.pl 5 |sort -T {config[TMPDIR]}  -k1,1 -k2,2n -k3,3n > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
		'''

rule cageSupportedfivepEnds:
	input:
		fivePends="mappings/nonAnchoredMergeReads/5pEnds/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.5pEnds.bed",
		tms="mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.splicing_status:all.endSupport:all.bed",
		cagePeaks=lambda wildcards: GENOMETOCAGEPEAKS[CAPDESIGNTOGENOME[wildcards.capDesign]],
		genome = lambda wildcards: config["GENOMESDIR"] + CAPDESIGNTOGENOME[wildcards.capDesign] + ".genome"
	output: "mappings/nonAnchoredMergeReads/cageSupported/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.cageSupported5pEnds.bed"
	shell:
		'''
uuidTmpOut=$(uuidgen)
cat {input.fivePends} | sort -T {config[TMPDIR]}  -k1,1 -k2,2n -k3,3n  | bedtools slop -s -l 50 -r 50 -i stdin -g {input.genome} | bedtools intersect -u -s -a stdin -b {input.cagePeaks} | cut -f4 | fgrep -w -f - {input.tms} > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
		'''

rule dhsSupportedfivepEnds:
	input:
		fivePends="mappings/nonAnchoredMergeReads/5pEnds/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.5pEnds.bed",
		tms="mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.splicing_status:all.endSupport:all.bed",
		dhsPeaks=lambda wildcards: GENOMETODHSPEAKS[CAPDESIGNTOGENOME[wildcards.capDesign]],
	output: "mappings/nonAnchoredMergeReads/dhsSupported/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.dhsSupported5pEnds.bed"
	shell:
		'''
uuidTmpOut=$(uuidgen)
cat {input.fivePends} | sort -T {config[TMPDIR]}  -k1,1 -k2,2n -k3,3n  | bedtools intersect -u -a stdin -b {input.dhsPeaks} | cut -f4 | fgrep -w -f - {input.tms} > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
		'''

rule plotDhsVsCage5primeComparisonStats:
	input: 
		dhs="mappings/nonAnchoredMergeReads/dhsSupported/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.dhsSupported5pEnds.bed",
		cage="mappings/nonAnchoredMergeReads/cageSupported/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.cageSupported5pEnds.bed"
	output: config["PLOTSDIR"] + "dhsVsCage5primeComparison.venn.stats/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.dhsVsCage5primeComparison.venn.stats.pdf"
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
	input: "mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.splicing_status:all.endSupport:all.bed"
	output: "mappings/nonAnchoredMergeReads/3pEnds/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.3pEnds.bed"
	shell:
		'''
uuidTmpOut=$(uuidgen)
cat {input} |extractTranscriptEndsFromBed12.pl 3 |sort -T {config[TMPDIR]}  -k1,1 -k2,2n -k3,3n > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
		'''

rule polyASupportedthreepEnds:
	input:
		threePends="mappings/nonAnchoredMergeReads/3pEnds/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.3pEnds.bed",
		tms="mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.splicing_status:all.endSupport:all.bed",
		polyAsites="mappings/removePolyAERCCs/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.polyAsitesNoErcc.bed",
		genome = lambda wildcards: config["GENOMESDIR"] + CAPDESIGNTOGENOME[wildcards.capDesign] + ".genome"
	output: "mappings/nonAnchoredMergeReads/polyASupported/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.polyASupported3pEnds.bed"
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
		polyA="mappings/nonAnchoredMergeReads/polyASupported/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.polyASupported3pEnds.bed",
		cage="mappings/nonAnchoredMergeReads/cageSupported/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.cageSupported5pEnds.bed",
		tms="mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.tmerge.min{minReadSupport}reads.splicing_status:{splicedStatus}.endSupport:all.gff"
	output:
		stats=config["STATSDATADIR"] + "tmp/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.splicing_status:{splicedStatus}.cagePolyASupport.stats.tsv",
		FLbed="mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.splicing_status:{splicedStatus}.endSupport:cagePolyASupported.bed",
		FLgff="mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.tmerge.min{minReadSupport}reads.splicing_status:{splicedStatus}.endSupport:cagePolyASupported.gff"



	shell:
		'''

uuid=$(uuidgen)
uuidTmpOutS=$(uuidgen)
uuidTmpOutB=$(uuidgen)
uuidTmpOutG=$(uuidgen)
cat {input.tms} |extractGffAttributeValue.pl transcript_id | sort -T {config[TMPDIR]} |uniq > {config[TMPDIR]}/$uuid.all.list

cat {input.polyA} | fgrep -w -f {config[TMPDIR]}/$uuid.all.list | cut -f4 | sort -T {config[TMPDIR]} |uniq > {config[TMPDIR]}/$uuid.polyA.list
cat {input.cage} | fgrep -w -f {config[TMPDIR]}/$uuid.all.list | cut -f4 | sort -T {config[TMPDIR]} |uniq > {config[TMPDIR]}/$uuid.cage.list
cat {config[TMPDIR]}/$uuid.polyA.list {config[TMPDIR]}/$uuid.cage.list |sort -T {config[TMPDIR]} |uniq > {config[TMPDIR]}/$uuid.cageOrPolyA.list
comm -1 -2 {config[TMPDIR]}/$uuid.polyA.list {config[TMPDIR]}/$uuid.cage.list |sort -T {config[TMPDIR]} |uniq > {config[TMPDIR]}/$uuid.cage+PolyA.list
noCageNoPolyA=$(comm -2 -3 {config[TMPDIR]}/$uuid.all.list {config[TMPDIR]}/$uuid.cageOrPolyA.list |wc -l)
cageOnly=$(comm -2 -3 {config[TMPDIR]}/$uuid.cage.list {config[TMPDIR]}/$uuid.polyA.list |wc -l)
polyAOnly=$(comm -2 -3 {config[TMPDIR]}/$uuid.polyA.list {config[TMPDIR]}/$uuid.cage.list |wc -l)
cageAndPolyA=$(cat {config[TMPDIR]}/$uuid.cage+PolyA.list | wc -l)
let total=$noCageNoPolyA+$cageOnly+$polyAOnly+$cageAndPolyA
fgrep -w -f {config[TMPDIR]}/$uuid.cage+PolyA.list {input.tms} |sort -T {config[TMPDIR]}  -k1,1 -k4,4n -k5,5n  > {config[TMPDIR]}/$uuidTmpOutG
mv {config[TMPDIR]}/$uuidTmpOutG {output.FLgff}
cat {output.FLgff} |gff2bed_full.pl - | sort -T {config[TMPDIR]}  -k1,1 -k2,2n -k3,3n  > {config[TMPDIR]}/$uuidTmpOutB
mv {config[TMPDIR]}/$uuidTmpOutB {output.FLbed}
echo -e "{wildcards.techname}Corr{wildcards.corrLevel}\t{wildcards.capDesign}\t{wildcards.sizeFrac}\t{wildcards.barcodes}\t$total\t$cageOnly\t$cageAndPolyA\t$polyAOnly\t$noCageNoPolyA"  > {config[TMPDIR]}/$uuidTmpOutS
mv {config[TMPDIR]}/$uuidTmpOutS {output.stats}
		'''


rule aggCagePolyAStats:
	input: lambda wildcards: expand(config["STATSDATADIR"] + "tmp/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.splicing_status:{splicedStatus}.cagePolyASupport.stats.tsv",filtered_product_merge, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNSplusMERGED, sizeFrac=SIZEFRACS, barcodes=BARCODESpluSMERGED, minReadSupport=wildcards.minReadSupport, splicedStatus=wildcards.splicedStatus)
	output: config["STATSDATADIR"] + "all.min{minReadSupport}reads.splicing_status:{splicedStatus}.cagePolyASupport.stats.tsv"
	shell:
		'''
uuidTmpOut=$(uuidgen)
echo -e "seqTech\tcorrectionLevel\tcapDesign\tsizeFrac\ttissue\tcategory\tcount\tpercent" > {config[TMPDIR]}/$uuidTmpOut

cat {input} | awk '{{print $1"\\t"$2"\\t"$3"\\t"$4"\\tcageOnly\\t"$6"\\t"$6/$5"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\tcageAndPolyA\\t"$7"\\t"$7/$5"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\tpolyAOnly\\t"$8"\\t"$8/$5"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\tnoCageNoPolyA\\t"$9"\\t"$9/$5}}' | sed 's/Corr0/\tNo/' | sed 's/Corr{lastK}/\tYes/' | sort -T {config[TMPDIR]}  >> {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
		'''

rule plotCagePolyAStats:
	input: config["STATSDATADIR"] + "all.min{minReadSupport}reads.splicing_status:{splicedStatus}.cagePolyASupport.stats.tsv"
	output: returnPlotFilenames(config["PLOTSDIR"] + "cagePolyASupport.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.splicing_status:{splicedStatus}.cagePolyASupport.stats")
	params:
		filterDat=lambda wildcards: merge_figures_params(wildcards.capDesign, wildcards.sizeFrac, wildcards.barcodes, wildcards.corrLevel, wildcards.techname)
	shell:
		'''
echo "
library(cowplot)
library(plyr)
library(scales)
library(gridExtra)
library(grid)
library(ggplotify)

dat <- read.table('{input}', header=T, as.is=T, sep='\\t')
{params.filterDat[10]}
{params.filterDat[0]}
{params.filterDat[1]}
{params.filterDat[2]}
{params.filterDat[3]}
{params.filterDat[4]}
{params.filterDat[5]}
{params.filterDat[8]}

dat\$category<-factor(dat\$category, ordered=TRUE, levels=rev(c('cageOnly', 'cageAndPolyA', 'polyAOnly', 'noCageNoPolyA')))
plotBase <- \\"ggplot(dat[order(dat\$category), ], aes(x=factor(correctionLevel), y=count, fill=category)) +
geom_bar(stat='identity') +
ylab('# TMs') +
scale_y_continuous(labels=comma)+
scale_fill_manual (values=c(cageOnly='#e5b3e5', cageAndPolyA='#C453C4', polyAOnly = '#b3e0ff', noCageNoPolyA='#a6a6a6'))+
xlab('{params.filterDat[6]}') +
guides(fill = guide_legend(title='Category'))+
#geom_text(position = 'stack', size=geom_textSize, aes(x = factor(correctionLevel), y = count, label = paste(sep='',percent(round(percent, digits=2)),' / ','(',comma(count),')'), hjust = 0.5, vjust = 1))+
{params.filterDat[7]}
{GGPLOT_PUB_QUALITY} + \\"

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
cat $(dirname {output[0]})/$(basename {output[0]} .legendOnly.png).r | R --slave

		'''



rule plotMetaTmEndsStats:
	input: config["STATSDATADIR"] + "all.{capDesign}.min{minReadSupport}.endSupport:{endSupport}.TmStats.stats.tsv"
	output: returnPlotFilenames(config["PLOTSDIR"] + "TmEndsStats.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.splicing_status:{splicedStatus}.endSupport:{endSupport}.TmEndsStats.meta.stats"),

	params:
		filterDat=lambda wildcards: merge_figures_params(wildcards.capDesign, wildcards.sizeFrac, wildcards.barcodes, wildcards.corrLevel, wildcards.techname, wildcards.splicedStatus)
	shell:
		'''
echo "
library(cowplot)
library(plyr)
library(scales)
library(gridExtra)
library(grid)
library(ggplotify)
library(data.table)
dat <- fread('{input}', header=T, sep='\\t')
{params.filterDat[10]}
{params.filterDat[0]}
{params.filterDat[1]}
{params.filterDat[2]}
{params.filterDat[3]}
{params.filterDat[4]}
{params.filterDat[5]}
{params.filterDat[8]}
{params.filterDat[11]}

dat\$end <- factor(dat\$end)

plotBase <- \\"ggplot(dat, aes(x=normDistance, color=end)) +
stat_density(geom='line', adjust=3) +
scale_color_manual(values=c('5' = '#009900', '3' = '#800000')) +
xlab('Normalized distance of read ends to TM\\'s TSS') +
{GGPLOT_PUB_QUALITY}+ \\"

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
cat $(dirname {output[0]})/$(basename {output[0]} .legendOnly.png).r | R --slave



		'''


rule plotAbsFiveTmEndsStats:
	input: config["STATSDATADIR"] + "all.{capDesign}.min{minReadSupport}.endSupport:{endSupport}.TmStats.stats.tsv"
	output:
		five=returnPlotFilenames(config["PLOTSDIR"] + "TmEndsStats.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.splicing_status:{splicedStatus}.endSupport:{endSupport}.TmEndsStats.5p.abs.stats"),
	params:
		filterDat=lambda wildcards: merge_figures_params(wildcards.capDesign, wildcards.sizeFrac, wildcards.barcodes, wildcards.corrLevel, wildcards.techname, wildcards.splicedStatus)
	shell:
		'''

echo "
library(cowplot)
library(scales)
library(gridExtra)
library(grid)
library(ggplotify)
library(data.table)
library(dplyr)

dat <- fread('{input}', header=T, sep='\\t')
{params.filterDat[10]}
{params.filterDat[0]}
{params.filterDat[1]}
{params.filterDat[2]}
{params.filterDat[3]}
{params.filterDat[4]}
{params.filterDat[5]}
{params.filterDat[8]}
{params.filterDat[11]}

dat <- subset(dat, end==5)

dat\$end <- factor(dat\$end)
summaryStatsFive = transform(summarise(group_by(dat,seqTech, sizeFrac, tissue), Label = paste0('N= ', comma(length(distance)), '\\n', 'Median= ', comma(median(distance)))))


plotBase <- \\"ggplot(dat, aes(x=distance)) +
geom_histogram(aes(y=..density..,fill=end), binwidth=50) +
scale_fill_manual(values=c('5' = '#009900'), name='End', labels =c('5' = '5´')) +
xlab('Distance of read ends\\nto TM´s TSS (mature RNA nts)') +
geom_text(data = summaryStatsFive, aes(label = Label, x = 0, y = Inf), hjust=0, vjust=1,  size=geom_textSize) +
coord_cartesian(xlim=c(0, 2000)) +
{GGPLOT_PUB_QUALITY}+ \\"

{params.filterDat[12]}

save_plot('{output.five[0]}', legendOnly, base_width=wLegendOnly, base_height=hLegendOnly)
save_plot('{output.five[1]}', legendOnly, base_width=wLegendOnly, base_height=hLegendOnly)

save_plot('{output.five[2]}', pXy, base_width=wXyPlot, base_height=hXyPlot)
save_plot('{output.five[3]}', pXy, base_width=wXyPlot, base_height=hXyPlot)

save_plot('{output.five[4]}', pXyNoLegend, base_width=wXyNoLegendPlot, base_height=hXyNoLegendPlot)
save_plot('{output.five[5]}', pXyNoLegend, base_width=wXyNoLegendPlot, base_height=hXyNoLegendPlot)

save_plot('{output.five[6]}', pYx, base_width=wYxPlot, base_height=hYxPlot)
save_plot('{output.five[7]}', pYx, base_width=wYxPlot, base_height=hYxPlot)

save_plot('{output.five[8]}', pYxNoLegend, base_width=wYxNoLegendPlot, base_height=hYxNoLegendPlot)
save_plot('{output.five[9]}', pYxNoLegend, base_width=wYxNoLegendPlot, base_height=hYxNoLegendPlot)

" > $(dirname {output.five[0]})/$(basename {output.five[0]} .5p.abs.stats.legendOnly.png).r
cat $(dirname {output.five[0]})/$(basename {output.five[0]} .5p.abs.stats.legendOnly.png).r | R --slave


		'''


rule plotAbsThreeTmEndsStats:
	input: config["STATSDATADIR"] + "all.{capDesign}.min{minReadSupport}.endSupport:{endSupport}.TmStats.stats.tsv"
	output:
		three=returnPlotFilenames(config["PLOTSDIR"] + "TmEndsStats.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.splicing_status:{splicedStatus}.endSupport:{endSupport}.TmEndsStats.3p.abs.stats"),
	params:
		filterDat=lambda wildcards: merge_figures_params(wildcards.capDesign, wildcards.sizeFrac, wildcards.barcodes, wildcards.corrLevel, wildcards.techname, wildcards.splicedStatus)
	shell:
		'''

echo "
library(cowplot)
library(scales)
library(gridExtra)
library(grid)
library(ggplotify)
library(data.table)
library(dplyr)

dat <- fread('{input}', header=T, sep='\\t')
{params.filterDat[10]}
{params.filterDat[0]}
{params.filterDat[1]}
{params.filterDat[2]}
{params.filterDat[3]}
{params.filterDat[4]}
{params.filterDat[5]}
{params.filterDat[8]}
{params.filterDat[11]}

dat <- subset(dat, end==3)
dat\$end <- factor(dat\$end)

summaryStatsThree = transform(summarise(group_by(dat,seqTech, sizeFrac, tissue), Label = paste0('N= ', comma(length(distance)), '\\n', 'Median= ', comma(median(distance)))))


plotBase <- \\"ggplot(dat, aes(x=distance)) +
geom_histogram(aes(y=..density..,fill=end), binwidth=50) +
scale_fill_manual(values=c('3' = '#800000'), name='End', labels =c('3' = '3´')) +
xlab('Distance of read ends\\n to TM´s 3´ end (mature RNA nts)') +
geom_text(data = summaryStatsThree, aes(label = Label, x = -Inf, y = Inf), hjust=0, vjust=1,  size=geom_textSize) +
coord_cartesian(xlim=c(-2000, 0)) +
{GGPLOT_PUB_QUALITY}+ \\"

{params.filterDat[12]}

save_plot('{output.three[0]}', legendOnly, base_width=wLegendOnly, base_height=hLegendOnly)
save_plot('{output.three[1]}', legendOnly, base_width=wLegendOnly, base_height=hLegendOnly)

save_plot('{output.three[2]}', pXy, base_width=wXyPlot, base_height=hXyPlot)
save_plot('{output.three[3]}', pXy, base_width=wXyPlot, base_height=hXyPlot)

save_plot('{output.three[4]}', pXyNoLegend, base_width=wXyNoLegendPlot, base_height=hXyNoLegendPlot)
save_plot('{output.three[5]}', pXyNoLegend, base_width=wXyNoLegendPlot, base_height=hXyNoLegendPlot)

save_plot('{output.three[6]}', pYx, base_width=wYxPlot, base_height=hYxPlot)
save_plot('{output.three[7]}', pYx, base_width=wYxPlot, base_height=hYxPlot)

save_plot('{output.three[8]}', pYxNoLegend, base_width=wYxNoLegendPlot, base_height=hYxNoLegendPlot)
save_plot('{output.three[9]}', pYxNoLegend, base_width=wYxNoLegendPlot, base_height=hYxNoLegendPlot)

" > $(dirname {output.three[0]})/$(basename {output.three[0]} .5p.abs.stats.legendOnly.png).r
cat $(dirname {output.three[0]})/$(basename {output.three[0]} .5p.abs.stats.legendOnly.png).r | R --slave


		'''


