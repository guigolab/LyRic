rule extractPooledTmsFivepEnds:
	input: "mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.all.bed"
	output: temp("mappings/nonAnchoredMergeReads/5pEnds/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.5pEnds.bed")
	shell:
		'''
cat {input} |extractTranscriptEndsFromBed12.pl 5 |sort -T {config[TMPDIR]}  -k1,1 -k2,2n -k3,3n > {output}
		'''

rule cageSupportedfivepEnds:
	input:
		fivePends="mappings/nonAnchoredMergeReads/5pEnds/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.5pEnds.bed",
		tms="mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.all.bed",
		cagePeaks=lambda wildcards: CAPDESIGNTOCAGEPEAKS[wildcards.capDesign],
		genome = lambda wildcards: config["GENOMESDIR"] + CAPDESIGNTOGENOME[wildcards.capDesign] + ".genome"
	output: temp("mappings/nonAnchoredMergeReads/cageSupported/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.cageSupported5pEnds.bed")
	shell:
		'''
cat {input.fivePends} | sort -T {config[TMPDIR]}  -k1,1 -k2,2n -k3,3n  | bedtools slop -s -l 50 -r 50 -i stdin -g {input.genome} | bedtools intersect -u -s -a stdin -b {input.cagePeaks} | cut -f4 | fgrep -w -f - {input.tms} > {output}
		'''


rule extractPooledTmsThreepEnds:
	input: "mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.all.bed"
	output: temp("mappings/nonAnchoredMergeReads/3pEnds/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.3pEnds.bed")
	shell:
		'''
cat {input} |extractTranscriptEndsFromBed12.pl 3 |sort -T {config[TMPDIR]}  -k1,1 -k2,2n -k3,3n > {output}
		'''

rule polyASupportedthreepEnds:
	input:
		threePends="mappings/nonAnchoredMergeReads/3pEnds/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.3pEnds.bed",
		tms="mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.all.bed",
		polyAsites="mappings/removePolyAERCCs/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.polyAsitesNoErcc.bed",
		genome = lambda wildcards: config["GENOMESDIR"] + CAPDESIGNTOGENOME[wildcards.capDesign] + ".genome"
	output: temp("mappings/nonAnchoredMergeReads/polyASupported/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.polyASupported3pEnds.bed")
	shell:
		'''
uuid=$(uuidgen)
cat {input.polyAsites} |sort -T {config[TMPDIR]}  -k1,1 -k2,2n -k3,3n  > {config[TMPDIR]}/$uuid.polyAsites.bed
cat {input.threePends} | sort -T {config[TMPDIR]}  -k1,1 -k2,2n -k3,3n  | bedtools slop -s -l 5 -r 5 -i stdin -g {input.genome} | bedtools intersect -u -s -a stdin -b {config[TMPDIR]}/$uuid.polyAsites.bed | cut -f4 | tgrep -F -w -f - {input.tms} > {output}
		'''

rule getCagePolyASupport:
	input:
		polyA="mappings/nonAnchoredMergeReads/polyASupported/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.polyASupported3pEnds.bed",
		cage="mappings/nonAnchoredMergeReads/cageSupported/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.cageSupported5pEnds.bed",
		tms="mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.all.gff"
	output:
		stats=temp(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.splicing_status:{splicedStatus}.cagePolyASupport.stats.tsv"),
		FLbed="mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.splicing_status:{splicedStatus}.cage+polyASupported.bed",
		FLgff="mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.splicing_status:{splicedStatus}.cage+polyASupported.gff"
	params:
		grepSpliced = lambda wildcards: '| fgrep \'spliced \"1\"\'' if wildcards.splicedStatus == "spliced" else '| fgrep \'spliced \"0\"\'' if wildcards.splicedStatus == "unspliced" else ''


	shell:
		'''

uuid=$(uuidgen)
cat {input.tms} {params.grepSpliced} |extractGffAttributeValue.pl transcript_id | sort -T {config[TMPDIR]} |uniq > {config[TMPDIR]}/$uuid.all.list

cat {input.polyA} | fgrep -w -f {config[TMPDIR]}/$uuid.all.list | cut -f4 | sort -T {config[TMPDIR]} |uniq > {config[TMPDIR]}/$uuid.polyA.list
cat {input.cage} | fgrep -w -f {config[TMPDIR]}/$uuid.all.list | cut -f4 | sort -T {config[TMPDIR]} |uniq > {config[TMPDIR]}/$uuid.cage.list
cat {config[TMPDIR]}/$uuid.polyA.list {config[TMPDIR]}/$uuid.cage.list |sort -T {config[TMPDIR]} |uniq > {config[TMPDIR]}/$uuid.cageOrPolyA.list
comm -1 -2 {config[TMPDIR]}/$uuid.polyA.list {config[TMPDIR]}/$uuid.cage.list |sort -T {config[TMPDIR]} |uniq > {config[TMPDIR]}/$uuid.cage+PolyA.list
noCageNoPolyA=$(comm -2 -3 {config[TMPDIR]}/$uuid.all.list {config[TMPDIR]}/$uuid.cageOrPolyA.list |wc -l)
cageOnly=$(comm -2 -3 {config[TMPDIR]}/$uuid.cage.list {config[TMPDIR]}/$uuid.polyA.list |wc -l)
polyAOnly=$(comm -2 -3 {config[TMPDIR]}/$uuid.polyA.list {config[TMPDIR]}/$uuid.cage.list |wc -l)
cageAndPolyA=$(cat {config[TMPDIR]}/$uuid.cage+PolyA.list | wc -l)
let total=$noCageNoPolyA+$cageOnly+$polyAOnly+$cageAndPolyA
fgrep -w -f {config[TMPDIR]}/$uuid.cage+PolyA.list {input.tms} {params.grepSpliced} |sort -T {config[TMPDIR]}  -k1,1 -k4,4n -k5,5n  > {output.FLgff}
cat {output.FLgff} |gff2bed_full.pl - | sort -T {config[TMPDIR]}  -k1,1 -k2,2n -k3,3n  > {output.FLbed}
echo -e "{wildcards.techname}Corr{wildcards.corrLevel}\t{wildcards.capDesign}\t{wildcards.sizeFrac}\t{wildcards.barcodes}\t$total\t$cageOnly\t$cageAndPolyA\t$polyAOnly\t$noCageNoPolyA"  > {output.stats}
		'''


rule aggCagePolyAStats:
	input: lambda wildcards: expand(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.splicing_status:{splicedStatus}.cagePolyASupport.stats.tsv",filtered_product_merge, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNSplusMERGED, sizeFrac=SIZEFRACSpluSMERGED, barcodes=BARCODESpluSMERGED, minReadSupport=wildcards.minReadSupport, splicedStatus=wildcards.splicedStatus)
	output: config["STATSDATADIR"] + "all.min{minReadSupport}reads.splicing_status:{splicedStatus}.cagePolyASupport.stats.tsv"
	shell:
		'''
echo -e "seqTech\tcorrectionLevel\tcapDesign\tsizeFrac\ttissue\tcategory\tcount\tpercent" > {output}

cat {input} | awk '{{print $1"\\t"$2"\\t"$3"\\t"$4"\\tcageOnly\\t"$6"\\t"$6/$5"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\tcageAndPolyA\\t"$7"\\t"$7/$5"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\tpolyAOnly\\t"$8"\\t"$8/$5"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\tnoCageNoPolyA\\t"$9"\\t"$9/$5}}' | sed 's/Corr0/\tNo/' | sed 's/Corr{lastK}/\tYes/' | sort -T {config[TMPDIR]}  >> {output}
		'''

rule plotCagePolyAStats:
	input: config["STATSDATADIR"] + "all.min{minReadSupport}reads.splicing_status:{splicedStatus}.cagePolyASupport.stats.tsv"
	output: config["PLOTSDIR"] + "cagePolyASupport.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.splicing_status:{splicedStatus}.cagePolyASupport.stats.{ext}"
	params:
		filterDat=lambda wildcards: merge_figures_params(wildcards.capDesign, wildcards.sizeFrac, wildcards.barcodes, wildcards.corrLevel, wildcards.techname)
	shell:
		'''
echo "library(ggplot2)
library(plyr)
library(scales)
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
ggplot(dat[order(dat\$category), ], aes(x=factor(correctionLevel), y=count, fill=category)) +
geom_bar(stat='identity') +
ylab('# CLS TMs') +
scale_y_continuous(labels=comma)+
#scale_fill_manual (values=c(cageOnly='#66B366', cageAndPolyA='#82865f', polyAOnly = '#D49090', noCageNoPolyA='#a6a6a6'))+
scale_fill_manual (values=c(cageOnly='#98cd98', cageAndPolyA='#C453C4', polyAOnly = '#b3e0ff', noCageNoPolyA='#a6a6a6'))+

facet_grid( seqTech + sizeFrac ~ capDesign + tissue, scales='free_y')+
xlab('{params.filterDat[6]}') +
guides(fill = guide_legend(title='Category'))+
geom_text(position = 'stack', size=geom_textSize, aes(x = factor(correctionLevel), y = count, label = paste(sep='',percent(round(percent, digits=2)),' / ','(',comma(count),')'), hjust = 0.5, vjust = 1))+
{params.filterDat[7]}
{GGPLOT_PUB_QUALITY}
ggsave('{output}', width=plotWidth, height=plotHeight)
" > {output}.r
cat {output}.r | R --slave

		'''



rule plotMetaTmEndsStats:
	input: config["STATSDATADIR"] + "all.min{minReadSupport}.{endSupport}.TmStats.stats.tsv"
	output: config["PLOTSDIR"] + "TmEndsStats.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.splicing_status:{splicedStatus}.{endSupport}.TmEndsStats.meta.stats.{ext}",

	params:
		filterDat=lambda wildcards: merge_figures_params(wildcards.capDesign, wildcards.sizeFrac, wildcards.barcodes, wildcards.corrLevel, wildcards.techname, wildcards.splicedStatus)
	shell:
		'''
echo "library(ggplot2)
library(data.table)
library(plyr)
library(scales)
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

ggplot(dat, aes(x=normDistance, color=end)) +
stat_density(geom='line', adjust=3) +
scale_color_manual(values=c('5p' = '#009900', '3p' = '#800000')) +
facet_grid( seqTech + sizeFrac ~ capDesign + tissue, scales='free_y')+
xlab('Normalized distance of read ends to TM\\'s TSS') +
{GGPLOT_PUB_QUALITY}
ggsave('{output}', width=plotWidth, height=plotHeight)
" > {output}.r
cat {output}.r | R --slave


		'''


rule plotAbsTmEndsStats:
	input: config["STATSDATADIR"] + "all.min{minReadSupport}.{endSupport}.TmStats.stats.tsv"
	output:
		five=config["PLOTSDIR"] + "TmEndsStats.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.splicing_status:{splicedStatus}.{endSupport}.TmEndsStats.5p.abs.stats.{ext}",
		three=config["PLOTSDIR"] + "TmEndsStats.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.splicing_status:{splicedStatus}.{endSupport}.TmEndsStats.3p.abs.stats.{ext}",
	params:
		filterDat=lambda wildcards: merge_figures_params(wildcards.capDesign, wildcards.sizeFrac, wildcards.barcodes, wildcards.corrLevel, wildcards.techname, wildcards.splicedStatus)
	shell:
		'''
#R --version
#module list
#ldd /nfs/users2/rg/jlagarde/R_libs/scales/libs/scales.so

echo "library(ggplot2)
library(data.table)
library(dplyr)
library(scales)
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

datFive <- subset(dat, end=='5p')
summaryStatsFive = transform(summarise(group_by(datFive,seqTech, sizeFrac, capDesign, tissue), Label = paste0('N= ', comma(length(distance)), '\\n', 'Median= ', comma(median(distance)))))

datThree <- subset(dat, end=='3p')
summaryStatsThree = transform(summarise(group_by(datThree,seqTech, sizeFrac, capDesign, tissue), Label = paste0('N= ', comma(length(distance)), '\\n', 'Median= ', comma(median(distance)))))



ggplot(datFive, aes(x=distance)) +
geom_histogram(aes(y=..density..,fill=end), binwidth=50) +
scale_fill_manual(values=c('5p' = '#009900'), name='End', labels =c('5p' = '5\\'')) +
facet_grid( seqTech + sizeFrac ~ capDesign + tissue, scales='free_y')+
xlab('Distance of read ends\\nto TM\\'s TSS (mature RNA nts)') +
geom_text(data = summaryStatsFive, aes(label = Label, x = 1000, y = Inf), hjust=0, vjust=3.8,  size=geom_textSize) +
coord_cartesian(xlim=c(0, 2000)) +
{GGPLOT_PUB_QUALITY}
ggsave('{output.five}', width=plotWidth, height=plotHeight)

ggplot(datThree, aes(x=distance)) +
geom_histogram(aes(y=..density..,fill=end), binwidth=50) +
scale_fill_manual(values=c('3p' = '#800000'), name='End', labels =c('3p' = '3\\'')) +
facet_grid( seqTech + sizeFrac ~ capDesign + tissue, scales='free_y')+
xlab('Distance of read ends\\n to TM\\'s 3\\' end (mature RNA nts)') +
geom_text(data = summaryStatsThree, aes(label = Label, x = -1500, y = Inf), hjust=0, vjust=3.8,  size=geom_textSize) +
coord_cartesian(xlim=c(-2000, 0)) +
{GGPLOT_PUB_QUALITY}
ggsave('{output.three}', width=plotWidth, height=plotHeight)


" > {output.five}.r
cat {output.five}.r | R --slave

		'''
