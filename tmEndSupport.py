rule extractPooledTmsFivepEnds:
	input: "mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.bed"
	output: "mappings/nonAnchoredMergeReads/5pEnds/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.5pEnds.bed"
	shell:
		'''
cat {input} |extractTranscriptEndsFromBed12.pl 5 |sortbed> {output}
		'''

rule cageSupportedfivepEnds:
	input:
		fivePends="mappings/nonAnchoredMergeReads/5pEnds/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.5pEnds.bed",
		tms="mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.bed",
		cagePeaks=lambda wildcards: CAPDESIGNTOCAGEPEAKS[wildcards.capDesign],
		genome = lambda wildcards: config["GENOMESDIR"] + CAPDESIGNTOGENOME[wildcards.capDesign] + ".genome"
	output: "mappings/nonAnchoredMergeReads/cageSupported/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.cageSupported.bed"
	shell:
		'''
cat {input.fivePends} | sortbed | bedtools slop -s -l 50 -r 50 -i stdin -g {input.genome} | bedtools intersect -u -s -a stdin -b {input.cagePeaks} | cut -f4 | fgrep -w -f - {input.tms} > {output}
		'''


rule extractPooledTmsThreepEnds:
	input: "mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.bed"
	output: "mappings/nonAnchoredMergeReads/5pEnds/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.3pEnds.bed"
	shell:
		'''
cat {input} |extractTranscriptEndsFromBed12.pl 3 |sortbed> {output}
		'''

rule polyASupportedthreepEnds:
	input:
		threePends="mappings/nonAnchoredMergeReads/5pEnds/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.3pEnds.bed",
		tms="mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.bed",
		polyAsites="mappings/removePolyAERCCs/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.polyAsitesNoErcc.bed",
		genome = lambda wildcards: config["GENOMESDIR"] + CAPDESIGNTOGENOME[wildcards.capDesign] + ".genome"
	output: "mappings/nonAnchoredMergeReads/polyASupported/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.polyASupported.bed"
	shell:
		'''
cat {input.polyAsites} |sortbed > $TMPDIR/polyAsites.bed
cat {input.threePends} | sortbed | bedtools slop -s -l 5 -r 5 -i stdin -g {input.genome} | bedtools intersect -u -s -a stdin -b $TMPDIR/polyAsites.bed | cut -f4 | fgrep -w -f - {input.tms} > {output}
		'''

rule getCagePolyASupport:
	input:
		polyA="mappings/nonAnchoredMergeReads/polyASupported/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.polyASupported.bed",
		cage="mappings/nonAnchoredMergeReads/cageSupported/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.cageSupported.bed",
		tms="mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.bed"
	output:
		stats=temp(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.cagePolyASupport.stats.tsv"),
		FLbed="mappings/nonAnchoredMergeReads/cage+polyASupported/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.cage+polyASupported.bed"
	shell:
		'''
cat {input.polyA} | cut -f4 | sort|uniq > $TMPDIR/polyA.list
cat {input.cage} | cut -f4 | sort|uniq > $TMPDIR/cage.list
cat {input.tms} | cut -f4 | sort|uniq > $TMPDIR/all.list
cat $TMPDIR/polyA.list $TMPDIR/cage.list |sort|uniq > $TMPDIR/cageOrPolyA.list
comm -1 -2 $TMPDIR/polyA.list $TMPDIR/cage.list |sort|uniq > $TMPDIR/cage+PolyA.list
noCageNoPolyA=$(comm -2 -3 $TMPDIR/all.list $TMPDIR/cageOrPolyA.list |wc -l)
cageOnly=$(comm -2 -3 $TMPDIR/cage.list $TMPDIR/polyA.list |wc -l)
polyAOnly=$(comm -2 -3 $TMPDIR/polyA.list $TMPDIR/cage.list |wc -l)
cageAndPolyA=$(cat $TMPDIR/cage+PolyA.list | wc -l)
let total=$noCageNoPolyA+$cageOnly+$polyAOnly+$cageAndPolyA
fgrep -w -f $TMPDIR/cage+PolyA.list {input.tms} > {output.FLbed}

echo -e "{wildcards.techname}Corr{wildcards.corrLevel}\t{wildcards.capDesign}\t{wildcards.sizeFrac}\t{wildcards.barcodes}\t$total\t$cageOnly\t$cageAndPolyA\t$polyAOnly\t$noCageNoPolyA"  > {output.stats}
		'''


rule aggCagePolyAStats:
	input: expand(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.cagePolyASupport.stats.tsv",filtered_product_merge, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACSpluSMERGED, barcodes=BARCODESpluSMERGED)
	output: config["STATSDATADIR"] + "all.cagePolyASupport.stats.tsv"
	shell:
		'''
echo -e "seqTech\tcorrectionLevel\tcapDesign\tsizeFrac\ttissue\tcategory\tcount\tpercent" > {output}

cat {input} | awk '{{print $1"\\t"$2"\\t"$3"\\t"$4"\\tcageOnly\\t"$6"\\t"$6/$5"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\tcageAndPolyA\\t"$7"\\t"$7/$5"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\tpolyAOnly\\t"$8"\\t"$8/$5"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\tnoCageNoPolyA\\t"$9"\\t"$9/$5}}' | sed 's/Corr0/\tNo/' | sed 's/Corr{lastK}/\tYes/' | sort >> {output}
		'''

rule plotCagePolyAStats:
	input: config["STATSDATADIR"] + "all.cagePolyASupport.stats.tsv"
	output: config["PLOTSDIR"] + "cagePolyASupport.stats/{capDesign}_Corr{corrLevel}_{sizeFrac}_{barcodes}.cagePolyASupport.stats.{ext}"
	params:
		filterDat=lambda wildcards: merge_figures_params(wildcards.capDesign, wildcards.sizeFrac, wildcards.barcodes, wildcards.corrLevel)
	shell:
		'''
echo "library(ggplot2)
library(plyr)
library(scales)
dat <- read.table('{input}', header=T, as.is=T, sep='\\t')
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
scale_fill_manual (values=c(cageOnly='#66B366', cageAndPolyA='#82865f', polyAOnly = '#D49090', noCageNoPolyA='#a6a6a6'))+
facet_grid( seqTech + sizeFrac ~ capDesign + tissue)+
xlab('{params.filterDat[6]}') +
guides(fill = guide_legend(title='Category'))+
geom_text(position = 'stack', size=geom_textSize, aes(x = factor(correctionLevel), y = count, ymax=count, label = paste(sep='',percent(round(percent, digits=2)),' / ','(',comma(count),')'), hjust = 0.5, vjust = 1))+
{params.filterDat[7]}
{GGPLOT_PUB_QUALITY}
ggsave('{output}', width=plotWidth, height=plotHeight)
" > {output}.r
cat {output}.r | R --slave

		'''
