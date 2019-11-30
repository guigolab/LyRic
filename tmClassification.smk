rule compareTargetsToTms:
	input:
		tms= "mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.splicing_status:all.endSupport:all.gff",
		targetedSegments=config["TARGETSDIR"] + "{capDesign}_primary_targets.exons.reduced.gene_type.segments.gtf"
	output: "mappings/nonAnchoredMergeReads/vsTargets/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.gfftsv.gz"
	shell:
		'''
uuidTmpOut=$(uuidgen)
bedtools intersect -wao -a {input.targetedSegments} -b {input.tms} |gzip > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}

		'''

rule getTargetCoverageStats:
	input: "mappings/nonAnchoredMergeReads/vsTargets/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.gfftsv.gz"
	output: temp(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.targetCoverage.stats.tsv")
	shell:
		'''
uuidTmpOut=$(uuidgen)
for type in `zcat {input} | extractGffAttributeValue.pl gene_type | sort -T {config[TMPDIR]} |uniq`; do
all=$(zcat {input} | fgrep "gene_type \\"$type\\";" | extractGffAttributeValue.pl transcript_id | sort -T {config[TMPDIR]} |uniq|wc -l)
detected=$(zcat {input} | fgrep "gene_type \\"$type\\";" | awk '$NF>0' | extractGffAttributeValue.pl transcript_id | sort -T {config[TMPDIR]} |uniq|wc -l)
let undetected=$all-$detected || true
echo -e "{wildcards.techname}Corr{wildcards.corrLevel}\t{wildcards.capDesign}\t{wildcards.sizeFrac}\t{wildcards.barcodes}\t$type\t$all\t$detected" | awk '{{print $0"\t"$7/$6}}'
done > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
		'''

rule aggTargetCoverageStats:
	input: lambda wildcards: expand(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.targetCoverage.stats.tsv",filtered_product_merge,  techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODESpluSMERGED, minReadSupport=wildcards.minReadSupport)
#lambda wildcards: expand(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_pooled.targetCoverage.stats.tsv", techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS)
	output: config["STATSDATADIR"] + "all.min{minReadSupport}reads.targetCoverage.stats.tsv"
	shell:
		'''
uuidTmpOut=$(uuidgen)
echo -e "seqTech\tcorrectionLevel\tcapDesign\tsizeFrac\ttissue\ttargetType\ttotalTargets\tdetectedTargets\tpercentDetectedTargets" > {config[TMPDIR]}/$uuidTmpOut
cat {input} | grep -v erccSpikein | sed 's/Corr0/\tNo/' | sed 's/Corr{lastK}/\tYes/'| sort -T {config[TMPDIR]}  >> {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
		'''

rule plotTargetCoverageStats:
	input: config["STATSDATADIR"] + "all.min{minReadSupport}reads.targetCoverage.stats.tsv"
	output: returnPlotFilenames(config["PLOTSDIR"] + "targetCoverage.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.targetCoverage.stats")
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


plotWidth = plotWidth + 2
plotHeight = plotHeight + 2

plotBase <- \\"ggplot(dat, aes(x=factor(correctionLevel), y=percentDetectedTargets, fill=targetType)) +
geom_bar(width=0.75,stat='identity', position=position_dodge(width=0.9)) +
scale_fill_manual(values={long_Rpalette}) +
geom_hline(aes(yintercept=1), linetype='dashed', alpha=0.7) +
geom_text(size=geom_textSize, aes(group=targetType, y=0.01, label = paste(sep='',percent(percentDetectedTargets),' / ','(',comma(detectedTargets),')')), angle=90, size=2.5, hjust=0, vjust=0.5, position = position_dodge(width=0.9)) +
ylab('% targeted regions detected') +
xlab('{params.filterDat[6]}') +
scale_y_continuous(limits = c(0, 1), labels = scales::percent)+
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


rule gffcompareToAnnotation:
	input:
		annot=lambda wildcards: CAPDESIGNTOANNOTGTF[wildcards.capDesign],
		tm="mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.splicing_status:all.endSupport:{endSupport}.gff"
	output: "mappings/nonAnchoredMergeReads/gffcompare/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.endSupport:{endSupport}.vs.gencode.simple.tsv"
	shell:
		'''
pref=$(basename {output} .simple.tsv)
annotFullPath=$(fullpath {input.annot})
tmFullPath=$(fullpath {input.tm})
cd $(dirname {output})
gffcompare -o $pref -r $annotFullPath $tmFullPath
cat $pref.tracking | simplifyGffCompareClasses.pl - > $(basename {output})

		'''

if SIRVpresent:
	rule gffcompareToSirvAnnotation:
		input:
			annot=config["SIRVgff"],
			tm="mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.splicing_status:all.endSupport:all.gff"
		output: "mappings/nonAnchoredMergeReads/gffcompare/SIRVs/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.vs.SIRVs.simple.tsv"
		shell:
			'''
pref=$(basename {output} .simple.tsv)
annotFullPath=$(fullpath {input.annot})
cat {input.tm} | awk '$1=="SIRVome_isoforms" ' > $(dirname {output})/$(basename {input.tm})
cd $(dirname {output})
gffcompare -o $pref -r $annotFullPath $(basename {input.tm})
cat $pref.tracking | simplifyGffCompareClasses.pl - > $(basename {output})

			'''

rule getGffCompareSirvStats:
	input:"mappings/nonAnchoredMergeReads/gffcompare/SIRVs/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.vs.SIRVs.simple.tsv"
	output: temp(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.vs.SIRVs.stats.tsv")
	shell:
		'''
uuidTmpOut=$(uuidgen)
file=$(dirname {input})/$(basename {input} .simple.tsv)
for level in `echo Baselevel Exonlevel Intronchainlevel Intronlevel Locuslevel Transcriptlevel`; do
Sn=`cat $file |grep "level:" |sed 's/ //g'| sed 's/:/\\t/'|sed 's/|$//'|sed 's/|/\\t/g' | awk -v l=$level '$1==l' |cut -f2`
Sp=`cat $file |grep "level:" |sed 's/ //g'| sed 's/:/\\t/'|sed 's/|$//'|sed 's/|/\\t/g' | awk -v l=$level '$1==l' |cut -f3`
echo -e "{wildcards.techname}Corr{wildcards.corrLevel}\t{wildcards.capDesign}\t{wildcards.sizeFrac}\t{wildcards.barcodes}\t$level\t$Sn\t$Sp";
done |sed 's/level//g' > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}

		'''

rule aggGffCompareSirvStats:
	input: lambda wildcards:expand(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.vs.SIRVs.stats.tsv", filtered_product_merge, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODESpluSMERGED, minReadSupport=wildcards.minReadSupport)
	output: config["STATSDATADIR"] + "all.tmerge.min{minReadSupport}reads.vs.SIRVs.stats.tsv"
	shell:
		'''
uuidTmpOut=$(uuidgen)
echo -e "seqTech\tcorrectionLevel\tcapDesign\tsizeFrac\ttissue\tlevel\tmetric\tvalue" > {config[TMPDIR]}/$uuidTmpOut
cat {input} | awk '{{print $1"\\t"$2"\\t"$3"\\t"$4"\\t"$5"\\tSn\\t"$6"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\t"$5"\\tPr\\t"$7}}'| sed 's/Corr0/\tNo/' | sed 's/Corr{lastK}/\tYes/' | sort -T {config[TMPDIR]}  >> {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}

		'''

rule plotGffCompareSirvStats:
	input:config["STATSDATADIR"] + "all.tmerge.min{minReadSupport}reads.vs.SIRVs.stats.tsv"
	output: returnPlotFilenames(config["PLOTSDIR"] + "tmerge.vs.SIRVs.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.vs.SIRVs.stats")
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
cbPalette <- c('Sn'='#cc6600', 'Pr'='#2d8659')
dat <- read.table('{input}', header=T, as.is=T, sep='\\t')
{params.filterDat[10]}
{params.filterDat[0]}
{params.filterDat[1]}
{params.filterDat[2]}
{params.filterDat[3]}
{params.filterDat[4]}
{params.filterDat[5]}
{params.filterDat[8]}

plotHeight = plotHeight +1
plotWidth = plotWidth +1
plotBase <- \\"ggplot(dat, aes(x=level, y=value)) +
geom_point(aes(color=metric, shape=correctionLevel), size=3, alpha=0.8) +
scale_colour_manual (values=cbPalette, name='Metric', breaks=c('Sn', 'Pr'))+
scale_shape_manual(values=c(16,21), name='Error correction') +
ylab('Sn | Pr (%)') +
xlab('Evaluation level') +
scale_y_continuous() +
expand_limits(y=c(0,100))+
theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
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

if SIRVpresent:
	rule getSirvDetectionStats:
		input:
			gffC="mappings/nonAnchoredMergeReads/gffcompare/SIRVs/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.vs.SIRVs.simple.tsv",
			tm="mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.splicing_status:all.endSupport:all.gff",
			sirvInfo=config["SIRVinfo"]
		output:temp(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.vs.SIRVs.detection.stats.tsv")
		shell:
			'''
uuid=$(uuidgen)
uuidTmpOut=$(uuidgen)
sirvDetectionStats.pl {input.sirvInfo} $(dirname {input.gffC})/$(basename {input.gffC} .simple.tsv).$(basename {input.tm}).refmap > {config[TMPDIR]}/$uuid
cat {config[TMPDIR]}/$uuid | while read id l c ca; do echo -e "{wildcards.techname}Corr{wildcards.corrLevel}\t{wildcards.capDesign}\t{wildcards.sizeFrac}\t{wildcards.barcodes}\t$id\t$l\t$c\t$ca"; done > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}

			'''

rule aggSirvDetectionStats:
	input: lambda wildcards:expand(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.vs.SIRVs.detection.stats.tsv", filtered_product_merge, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODESpluSMERGED, minReadSupport=wildcards.minReadSupport)
	output: config["STATSDATADIR"] + "all.tmerge.min{minReadSupport}reads.vs.SIRVs.detection.stats.tsv"
	shell:
		'''
uuidTmpOut=$(uuidgen)
echo -e "seqTech\tcorrectionLevel\tcapDesign\tsizeFrac\ttissue\tSIRVid\tlength\tconcentration\tdetectionStatus" > {config[TMPDIR]}/$uuidTmpOut
cat {input} | sed 's/Corr0/\tNo/' | sed 's/Corr{lastK}/\tYes/' | sort -T {config[TMPDIR]}  >> {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}

		'''

rule plotSirvDetectionStats:
	input:config["STATSDATADIR"] + "all.tmerge.min{minReadSupport}reads.vs.SIRVs.detection.stats.tsv"
	output: returnPlotFilenames(config["PLOTSDIR"] + "tmerge.vs.SIRVs.detection.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.vs.SIRVs.detection.stats")
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

palette <- c('end-to-end' = '#00e600', 'absent' = '#666666', 'partial' = '#ff0066')
dat <- read.table('{input}', header=T, as.is=T, sep='\\t')

{params.filterDat[10]}
{params.filterDat[0]}
{params.filterDat[1]}
{params.filterDat[2]}
{params.filterDat[3]}
{params.filterDat[4]}
{params.filterDat[5]}
{params.filterDat[8]}


plotBase <- \\"ggplot(dat, aes(x=concentration, y=length, color=detectionStatus, shape=correctionLevel)) + geom_point(alpha=0.93) +
coord_trans(x='log2') +
scale_color_manual(values=palette) +
scale_shape_manual(values=c('Yes'=16, 'No'=17), name='Error correction') +
xlab('SIRV molarity (fmol/uL)') +
ylab('SIRV length (nt)') +
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


rule colorBedAccordingToGffCompare:
	input:
		classes="mappings/nonAnchoredMergeReads/gffcompare/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.endSupport:{endSupport}.vs.gencode.simple.tsv",
		tm="mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.splicing_status:all.endSupport:{endSupport}.bed"
	output: "mappings/nonAnchoredMergeReads/colored/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.endSupport:{endSupport}.bed"
	shell:
			'''
uuidTmpOut=$(uuidgen)
colorNovelTxBed.pl {input.classes} {input.tm} > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}

			'''


rule getGffCompareStats:
	input: "mappings/nonAnchoredMergeReads/gffcompare/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.endSupport:{endSupport}.vs.gencode.simple.tsv"
	output: temp(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.endSupport:{endSupport}.vs.gencode.stats.tsv")
	shell:
		'''
uuidTmpOut=$(uuidgen)
cat {input} |cut -f2 | sort -T {config[TMPDIR]} |uniq -c | awk -v s={wildcards.techname}Corr{wildcards.corrLevel} -v c={wildcards.capDesign} -v si={wildcards.sizeFrac} -v b={wildcards.barcodes} '{{print s"\t"c"\t"si"\t"b"\t"$2"\t"$1}}' | sed 's/Corr0/\tNo/' | sed 's/Corr{lastK}/\tYes/' > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
		'''

rule getGffCompareGencodeSnPrStats:
	input:"mappings/nonAnchoredMergeReads/gffcompare/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.endSupport:{endSupport}.vs.gencode.simple.tsv"
	output: temp(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.endSupport:{endSupport}.vs.gencode.SnPr.stats.tsv")
	shell:
		'''
uuidTmpOut=$(uuidgen)
file=$(dirname {input})/$(basename {input} .simple.tsv)
for level in `echo Baselevel Exonlevel Intronchainlevel Intronlevel Locuslevel Transcriptlevel`; do
Sn=`cat $file |grep "level:" |sed 's/ //g'| sed 's/:/\\t/'|sed 's/|$//'|sed 's/|/\\t/g' | awk -v l=$level '$1==l' |cut -f2` 
Sp=`cat $file |grep "level:" |sed 's/ //g'| sed 's/:/\\t/'|sed 's/|$//'|sed 's/|/\\t/g' | awk -v l=$level '$1==l' |cut -f3` 
echo -e "{wildcards.techname}Corr{wildcards.corrLevel}\t{wildcards.capDesign}\t{wildcards.sizeFrac}\t{wildcards.barcodes}\t$level\t$Sn\t$Sp";
done |sed 's/level//g' > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}

		'''

rule aggGffCompareGencodeSnPrStats:
	input: lambda wildcards:expand(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.endSupport:{endSupport}.vs.gencode.SnPr.stats.tsv", filtered_product_merge, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODESpluSMERGED, endSupport=wildcards.endSupport, minReadSupport=wildcards.minReadSupport)
	output: config["STATSDATADIR"] + "all.tmerge.min{minReadSupport}reads.endSupport:{endSupport}.vs.gencode.SnPr.stats.tsv"
	shell:
		'''
uuidTmpOut=$(uuidgen)
echo -e "seqTech\tcorrectionLevel\tcapDesign\tsizeFrac\ttissue\tlevel\tmetric\tvalue" > {config[TMPDIR]}/$uuidTmpOut
cat {input} | awk '{{print $1"\\t"$2"\\t"$3"\\t"$4"\\t"$5"\\tSn\\t"$6"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\t"$5"\\tPr\\t"$7}}'| sed 's/Corr0/\tNo/' | sed 's/Corr{lastK}/\tYes/' | sort -T {config[TMPDIR]}  >> {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}

		'''

rule plotGffCompareGencodeSnPrStats:
	input:config["STATSDATADIR"] + "all.tmerge.min{minReadSupport}reads.endSupport:{endSupport}.vs.gencode.SnPr.stats.tsv"
	output: returnPlotFilenames(config["PLOTSDIR"] + "tmerge.vs.gencode.SnPr.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.endSupport:{endSupport}.vs.gencode.SnPr.stats")
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

cbPalette <- c('Sn'='#cc6600', 'Pr'='#2d8659')
dat <- read.table('{input}', header=T, as.is=T, sep='\\t')
{params.filterDat[10]}
{params.filterDat[0]}
{params.filterDat[1]}
{params.filterDat[2]}
{params.filterDat[3]}
{params.filterDat[4]}
{params.filterDat[5]}
{params.filterDat[8]}

plotHeight = plotHeight +1
plotWidth = plotWidth +1

plotBase <- \\"ggplot(dat, aes(x=level, y=value)) +
geom_point(aes(color=metric, shape=correctionLevel), size=3, alpha=0.8) +
scale_colour_manual (values=cbPalette, name='Metric', breaks=c('Sn', 'Pr'))+
scale_shape_manual(values=c(16,21), name='Error correction') +
ylab('Sn | Pr (%)') +
xlab('Evaluation level') +
scale_y_continuous() +
expand_limits(y=c(0,100))+
theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
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



rule aggGffCompareStats:
	input: lambda wildcards: expand(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.endSupport:{endSupport}.vs.gencode.stats.tsv",filtered_product_merge, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNSplusMERGED, sizeFrac=SIZEFRACS, barcodes=BARCODESpluSMERGED, endSupport=wildcards.endSupport, minReadSupport=wildcards.minReadSupport)
	output: config["STATSDATADIR"] + "all.tmerge.min{minReadSupport}reads.endSupport:{endSupport}.vs.gencode.stats.tsv"
	shell:
		'''
uuidTmpOut=$(uuidgen)
echo -e "seqTech\tcorrectionLevel\tcapDesign\tsizeFrac\ttissue\tcategory\tcount" > {config[TMPDIR]}/$uuidTmpOut
cat {input} >> {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
		'''

rule plotGffCompareStats:
	input: config["STATSDATADIR"] + "all.tmerge.min{minReadSupport}reads.endSupport:{endSupport}.vs.gencode.stats.tsv"
	output: returnPlotFilenames(config["PLOTSDIR"] + "tmerge.vs.gencode.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.endSupport:{endSupport}.vs.gencode.stats")
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


dat\$category<-factor(dat\$category, ordered=TRUE, levels=rev(c('Intergenic', 'Extends', 'Intronic', 'Overlaps', 'Antisense', 'Equal', 'Included')))
palette <- c('Intergenic' = '#0099cc', 'Extends' ='#00bfff', 'Intronic' = '#4dd2ff', 'Overlaps' = '#80dfff', 'Antisense' = '#ccf2ff', 'Equal' = '#c65353', 'Included' ='#d98c8c')

plotBase <- \\"ggplot(dat[order(dat\$category), ], aes(x=factor(correctionLevel), y=count, fill=category)) +
geom_bar(stat='identity') +
scale_fill_manual(values=palette) +
ylab('# TMs') +
xlab('{params.filterDat[6]}') +
guides(fill = guide_legend(title='Category'))+
scale_y_continuous(labels=scientific)+
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

rule simplifyGencode:
	input: lambda wildcards: CAPDESIGNTOANNOTGTF[wildcards.capDesign]
	output: "annotations/simplified/{capDesign}.gencode.simplified_biotypes.gtf"
	shell:
		'''
uuidTmpOut=$(uuidgen)
cat {input}  | simplifyGencodeGeneTypes.pl - | sort -T {config[TMPDIR]}  -k1,1 -k4,4n -k5,5n  > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
		'''

rule mergeTmsWithGencode:
	input:
		annot="annotations/simplified/{capDesign}.gencode.simplified_biotypes.gtf",
		tm="mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.splicing_status:all.endSupport:{endSupport}.gff"
	output: "mappings/nonAnchoredMergeReads/gencodeMerge/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.endSupport:{endSupport}.+gencode.gff.gz"
	threads:1
	shell:
		'''
uuidTmpOut=$(uuidgen)
cat {input.annot} {input.tm}  | skipcomments | sort -T {config[TMPDIR]}  -k1,1 -k4,4n -k5,5n  | tmerge --exonOverhangTolerance {config[exonOverhangTolerance]} - |sort -T {config[TMPDIR]}  -k1,1 -k4,4n -k5,5n  |gzip > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
		'''

rule makeClsGencodeLoci:
	input: "mappings/nonAnchoredMergeReads/gencodeMerge/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.endSupport:{endSupport}.+gencode.gff.gz"
	params: locusPrefix=config["PROJECT_NAME"]
	output: temp("mappings/nonAnchoredMergeReads/gencodeLociMerge/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.endSupport:{endSupport}.+gencode.loci.gff.gz")
	shell:
		'''
uuid=$(uuidgen)
uuidTmpOut=$(uuidgen)
zcat {input} > {config[TMPDIR]}/$uuid
bedtools intersect -s -wao -a {config[TMPDIR]}/$uuid -b {config[TMPDIR]}/$uuid |fgrep -v ERCC| buildLoci.pl --locPrefix {params.locusPrefix}: - |sort -T {config[TMPDIR]}  -k1,1 -k4,4n -k5,5n  | gzip> {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}

		'''


rule mergeWithRef:
	input:
		clsGencode="mappings/nonAnchoredMergeReads/gencodeLociMerge/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.endSupport:{endSupport}.+gencode.loci.gff.gz",
		gencode="annotations/simplified/{capDesign}.gencode.simplified_biotypes.gtf"
	output: "mappings/nonAnchoredMergeReads/mergeToRef/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.endSupport:{endSupport}.+gencode.loci.refmerged.gff.gz"
	shell:
		'''
uuid=$(uuidgen)
uuidTmpOut=$(uuidgen)
zcat  {input.clsGencode} > {config[TMPDIR]}/$uuid
mergeToRef.pl {input.gencode} {config[TMPDIR]}/$uuid | sort -T {config[TMPDIR]}  -k1,1 -k4,4n -k5,5n  |gzip > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
		'''

rule getNovelIntergenicLoci:
	input:
		gencode="annotations/simplified/{capDesign}.gencode.simplified_biotypes.gtf",
		tmergeGencode="mappings/nonAnchoredMergeReads/mergeToRef/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.endSupport:{endSupport}.+gencode.loci.refmerged.gff.gz"
	output:"mappings/nonAnchoredMergeReads/mergeToRef/novelLoci/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.endSupport:{endSupport}.novelLoci.gff.gz"
	shell:
		'''
uuid1=$(uuidgen)
uuid2=$(uuidgen)
uuid3=$(uuidgen)
uuid4=$(uuidgen)
uuid5=$(uuidgen)
uuidTmpOut=$(uuidgen)
cat {input.gencode} |awk '$3=="exon"' | extract_locus_coords.pl -| sort -T {config[TMPDIR]}  -k1,1 -k2,2n -k3,3n  > {config[TMPDIR]}/$uuid1
zcat {input.tmergeGencode} | tgrep -F 'gene_ref_status "novel";' > {config[TMPDIR]}/$uuid4
cat {config[TMPDIR]}/$uuid4 | extract_locus_coords.pl - | sort -T {config[TMPDIR]}  -k1,1 -k2,2n -k3,3n  > {config[TMPDIR]}/$uuid2
bedtools intersect -v -a {config[TMPDIR]}/$uuid2 -b {config[TMPDIR]}/$uuid1 > {config[TMPDIR]}/$uuid5
cat {config[TMPDIR]}/$uuid5 |tgrep -F -v ERCC |cut -f4 | sort -T {config[TMPDIR]} |uniq > {config[TMPDIR]}/$uuid3
zcat {input.tmergeGencode}| tgrep -F -w -f {config[TMPDIR]}/$uuid3 - |gzip > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
		'''

rule getNovelIntergenicLociStats:
	input:
		tmergeGencode="mappings/nonAnchoredMergeReads/mergeToRef/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.endSupport:{endSupport}.+gencode.loci.refmerged.gff.gz",
		intergenic="mappings/nonAnchoredMergeReads/mergeToRef/novelLoci/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.endSupport:{endSupport}.novelLoci.gff.gz"
	output: temp(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.endSupport:{endSupport}.novelLoci.stats.tsv")
	shell:
		'''
uuidTmpOut=$(uuidgen)
totalNovel=$(zcat {input.tmergeGencode} | tgrep -F 'gene_ref_status "novel";' |extractGffAttributeValue.pl gene_id | sort -T {config[TMPDIR]} | uniq | wc -l)
interg=$(zcat {input.intergenic} | extractGffAttributeValue.pl gene_id | sort -T {config[TMPDIR]} | uniq | wc -l)
echo -e "{wildcards.techname}Corr{wildcards.corrLevel}\t{wildcards.capDesign}\t{wildcards.sizeFrac}\t{wildcards.barcodes}\t$totalNovel\t$interg"  > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}

		'''

rule aggNovelIntergenicLociStats:
	input: lambda wildcards: expand(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.endSupport:{endSupport}.novelLoci.stats.tsv",filtered_product_merge, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNSplusMERGED, sizeFrac=SIZEFRACS, barcodes=BARCODESpluSMERGED, endSupport=wildcards.endSupport, minReadSupport=wildcards.minReadSupport)
	output: config["STATSDATADIR"] + "all.tmerge.min{minReadSupport}reads.endSupport:{endSupport}.novelLoci.stats.tsv"
	shell:
		'''
uuidTmpOut=$(uuidgen)
echo -e "seqTech\tcorrectionLevel\tcapDesign\tsizeFrac\ttissue\tcategory\tcount\tpercent" > {config[TMPDIR]}/$uuidTmpOut
cat {input} | awk '{{if ($5!=0) print $1"\\t"$2"\\t"$3"\\t"$4"\\tintergenic\\t"$6"\\t"$6/$5"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\tintronic\\t"$5-$6"\\t"($5-$6)/$5; else print $1"\\t"$2"\\t"$3"\\t"$4"\\tintergenic\\t"$6"\\t0\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\tintronic\\t"$5-$6"\\t0"}}'| sed 's/Corr0/\tNo/' | sed 's/Corr{lastK}/\tYes/' | sort -T {config[TMPDIR]}  >> {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
		'''


rule plotNovelIntergenicLociStats:
	input: config["STATSDATADIR"] + "all.tmerge.min{minReadSupport}reads.endSupport:{endSupport}.novelLoci.stats.tsv"
	output: returnPlotFilenames(config["PLOTSDIR"] + "tmerge.novelLoci.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.endSupport:{endSupport}.novelLoci.stats")
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

dat\$category<-factor(dat\$category, ordered=TRUE, levels=rev(c('intronic', 'intergenic')))
plotBase <- \\"ggplot(dat[order(dat\$category), ], aes(x=factor(correctionLevel), y=count, fill=category)) +
geom_bar(stat='identity') +
ylab('# Novel CLS loci') +
scale_y_continuous(labels=comma)+
scale_fill_manual (values=c(intronic='#d98c8c', intergenic='#33ccff'))+
xlab('{params.filterDat[6]}') +
guides(fill = guide_legend(title='Category\\n(w.r.t. GENCODE)'))+
geom_text(position = 'stack', size=geom_textSize, aes(x = factor(correctionLevel), y = count, label = paste(sep='',percent(round(percent, digits=2)),' / ','(',comma(count),')'), hjust = 0.5, vjust = 1))+
{params.filterDat[7]}
{GGPLOT_PUB_QUALITY}  + \\"

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



rule tmergeAll:
	input:
		tm=lambda wildcards: expand("mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.splicing_status:{splicedStatus}.endSupport:{endSupport}.gff", filtered_product_merge, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=wildcards.capDesign, sizeFrac=SIZEFRACSnoSIZESELECTONLY, barcodes=BARCODES, endSupport=wildcards.endSupport,  minReadSupport=wildcards.minReadSupport, splicedStatus=wildcards.splicedStatus)
		#gencode="annotations/simplified/{capDesign}.gencode.simplified_biotypes.gtf",
	output: "mappings/nonAnchoredMergeReads/tmergeAll/{capDesign}_min{minReadSupport}reads.splicing_status:{splicedStatus}.endSupport:{endSupport}.tmerge.gff"
	shell:
		'''
uuid=$(uuidgen)
uuidTmpOut=$(uuidgen)
for file in `echo {input.tm}`; do
bn=$(basename $file .tmerge.min{wildcards.minReadSupport}reads.{wildcards.endSupport}.gff)
cat $file | perl -sne '$_=~s/transcript_id \"(\S+)\"/transcript_id \"=$var=$1\"/g; print' -- -var=$bn
done > {config[TMPDIR]}/$uuid

cat {config[TMPDIR]}/$uuid  | skipcomments | sort -T {config[TMPDIR]} -k1,1 -k4,4n -k5,5n | tmerge --exonOverhangTolerance {config[exonOverhangTolerance]} - |sort -T {config[TMPDIR]}  -k1,1 -k4,4n -k5,5n  > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}

		'''

rule getSampleComparisonStats:
	input: "mappings/nonAnchoredMergeReads/tmergeAll/{capDesign}_min{minReadSupport}reads.splicing_status:{splicedStatus}.endSupport:{endSupport}.tmerge.gff"
	output: 
		fullMatrix=config["STATSDATADIR"] + "all.sampleComparison.{capDesign}_min{minReadSupport}reads.splicing_status:{splicedStatus}.endSupport:{endSupport}.stats.tsv",
		simpsonMatrix=config["STATSDATADIR"] + "all.sampleComparison.{capDesign}_min{minReadSupport}reads.splicing_status:{splicedStatus}.endSupport:{endSupport}.overlap_coeff.tsv",
		jaccardMatrix=config["STATSDATADIR"] + "all.sampleComparison.{capDesign}_min{minReadSupport}reads.splicing_status:{splicedStatus}.endSupport:{endSupport}.jaccard_ind.tsv",
		OneMinusSimpsonMatrix=config["STATSDATADIR"] + "all.sampleComparison.{capDesign}_min{minReadSupport}reads.splicing_status:{splicedStatus}.endSupport:{endSupport}.oneMinusSimpson_coeff.tsv",
		OneMinusJaccardMatrix=config["STATSDATADIR"] + "all.sampleComparison.{capDesign}_min{minReadSupport}reads.splicing_status:{splicedStatus}.endSupport:{endSupport}.oneMinusJaccard_coeff.tsv"
	shell:
		'''
cat {input} | tmergeToBinaryMatrix.pl - $(dirname {output.simpsonMatrix})/$(basename {output.simpsonMatrix} .overlap_coeff.tsv) |perl -ne '$_=~s/(\S+):Corr\d+_(\S+)/$1 $2/g; print' > {output.fullMatrix}

#convert simpson matrix to dissimilarity matrix (1-simpson)
 cat {output.simpsonMatrix} |perl -ne 'chomp; @line=split "\\t"; for($i=0; $i<=$#line;$i++){{$t=$line[$i];if ($t=~/^-?(?:\d+\.?|\.\d)\d*\z/ ){{$line[$i]=1-$line[$i]}}}}; print join("\\t", @line)."\\n"' > {output.OneMinusSimpsonMatrix}
 cat {output.jaccardMatrix} |perl -ne 'chomp; @line=split "\\t"; for($i=0; $i<=$#line;$i++){{$t=$line[$i];if ($t=~/^-?(?:\d+\.?|\.\d)\d*\z/ ){{$line[$i]=1-$line[$i]}}}}; print join("\\t", @line)."\\n"' > {output.OneMinusJaccardMatrix}

		'''

rule plotSampleComparisonStats:
	input: 
		simpson=config["STATSDATADIR"] + "all.sampleComparison.{capDesign}_min{minReadSupport}reads.splicing_status:{splicedStatus}.endSupport:{endSupport}.oneMinusSimpson_coeff.tsv",
		jaccard=config["STATSDATADIR"] + "all.sampleComparison.{capDesign}_min{minReadSupport}reads.splicing_status:{splicedStatus}.endSupport:{endSupport}.oneMinusJaccard_coeff.tsv",

	output: 
		simpson=config["PLOTSDIR"] + "sampleComparison.stats/{capDesign}/{capDesign}_min{minReadSupport}reads.splicing_status:{splicedStatus}.endSupport:{endSupport}.dendrogram.sampleComparison.simpson.png",
		jaccard=config["PLOTSDIR"] + "sampleComparison.stats/{capDesign}/{capDesign}_min{minReadSupport}reads.splicing_status:{splicedStatus}.endSupport:{endSupport}.dendrogram.sampleComparison.jaccard.png"
	shell:
		'''
echo "

library(cluster)
library(dendextend)
dat <- read.table('{input.simpson}', header=T, as.is=T, sep='\t', row.names=1)
dend <- as.dendrogram(hclust(as.dist(dat), method ='ward.D'))
lib_type <- rep('Other', length(rownames(dat)))
is_x <- grepl('CapTrap', rownames(dat))
lib_type[is_x] <- 'CapTrap'
is_x <- grepl('Telop', rownames(dat))
lib_type[is_x] <- 'Telop'
is_x <- grepl('Smarter', rownames(dat))
lib_type[is_x] <- 'Smarter'
is_x <- grepl('dRNA', rownames(dat))
lib_type[is_x] <- 'dRNA'
lib_type <- factor(lib_type)
n_lib_types <- length(unique(lib_type))
cols_lib <- colorspace::rainbow_hcl(n_lib_types, c = 100, l  = 70, start = -360, end =-180)
col_lib_type <- cols_lib[lib_type]

seqPlatform <- rep('Other', length(rownames(dat)))
is_x <- grepl('ont', rownames(dat))
seqPlatform[is_x] <- 'ont'
is_x <- grepl('pacBio', rownames(dat))
seqPlatform[is_x] <- 'pacBio'
seqPlatform <- factor(seqPlatform)
n_seqPlatforms <- length(unique(seqPlatform))
cols_seqPlatform <- colorspace::rainbow_hcl(n_seqPlatforms, c = 50, l  = 60, start = -350, end =200)
col_seqPlatform <- cols_seqPlatform[seqPlatform]

tissue <- rep('Other', length(rownames(dat)))
is_x <- grepl('Brain', rownames(dat))
tissue[is_x] <- 'Brain'
is_x <- grepl('Heart', rownames(dat))
tissue[is_x] <- 'Heart'
tissue <- factor(tissue)
n_tissues <- length(unique(tissue))
cols_tissue <- colorspace::rainbow_hcl(n_tissues, c = 120, l  = 100, start = -300, end =180)
col_tissue <- cols_tissue[tissue]

png(file = '{output.simpson}', width=1500, height=1500)
par(mar = c(10,1,1,25), cex=3)
plot(dend, axes=F, horiz=TRUE)
colored_bars(cbind(col_lib_type, col_seqPlatform, col_tissue), dend, rowLabels = c('Library type', 'Seq platform', 'Biosample'), horiz=TRUE, text_shift=16)
dev.off()

" > {output.simpson}.r
cat {output.simpson}.r | R --slave

echo "

library(cluster)
library(dendextend)
dat <- read.table('{input.jaccard}', header=T, as.is=T, sep='\t', row.names=1)
dend <- as.dendrogram(hclust(as.dist(dat)))
lib_type <- rep('Other', length(rownames(dat)))
is_x <- grepl('CapTrap', rownames(dat))
lib_type[is_x] <- 'CapTrap'
is_x <- grepl('Telop', rownames(dat))
lib_type[is_x] <- 'Telop'
is_x <- grepl('Smarter', rownames(dat))
lib_type[is_x] <- 'Smarter'
is_x <- grepl('dRNA', rownames(dat))
lib_type[is_x] <- 'dRNA'
lib_type <- factor(lib_type)
n_lib_types <- length(unique(lib_type))
cols_lib <- colorspace::rainbow_hcl(n_lib_types, c = 100, l  = 70, start = -360, end =-180)
col_lib_type <- cols_lib[lib_type]

seqPlatform <- rep('Other', length(rownames(dat)))
is_x <- grepl('ont', rownames(dat))
seqPlatform[is_x] <- 'ont'
is_x <- grepl('pacBio', rownames(dat))
seqPlatform[is_x] <- 'pacBio'
seqPlatform <- factor(seqPlatform)
n_seqPlatforms <- length(unique(seqPlatform))
cols_seqPlatform <- colorspace::rainbow_hcl(n_seqPlatforms, c = 50, l  = 60, start = -350, end =200)
col_seqPlatform <- cols_seqPlatform[seqPlatform]

tissue <- rep('Other', length(rownames(dat)))
is_x <- grepl('Brain', rownames(dat))
tissue[is_x] <- 'Brain'
is_x <- grepl('Heart', rownames(dat))
tissue[is_x] <- 'Heart'
tissue <- factor(tissue)
n_tissues <- length(unique(tissue))
cols_tissue <- colorspace::rainbow_hcl(n_tissues, c = 120, l  = 100, start = -300, end =180)
col_tissue <- cols_tissue[tissue]

png(file = '{output.jaccard}', width=1500, height=1500)
par(mar = c(10,1,1,25), cex=3)
plot(dend, axes=F, horiz=TRUE)
colored_bars(cbind(col_lib_type, col_seqPlatform, col_tissue), dend, rowLabels = c('Library type', 'Seq platform', 'Biosample'), horiz=TRUE, text_shift=16)
dev.off()

" > {output.jaccard}.r
cat {output.jaccard}.r | R --slave

		'''
