rule compareTargetsToTms:
	input:
		tms= "mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.gff",
		targetedSegments=config["TARGETSDIR"] + "{capDesign}_primary_targets.exons.reduced.gene_type.segments.gtf"
	output: "mappings/nonAnchoredMergeReads/vsTargets/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.gfftsv.gz"
	shell:
		'''
bedtools intersect -wao -a {input.targetedSegments} -b {input.tms} |gzip > {output}

		'''

rule getTargetCoverageStats:
	input: "mappings/nonAnchoredMergeReads/vsTargets/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.gfftsv.gz"
	output: temp(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.targetCoverage.stats.tsv")
	shell:
		'''
for type in `zcat {input} | extractGffAttributeValue.pl gene_type | sort|uniq`; do
all=$(zcat {input} | fgrep "gene_type \\"$type\\";" | extractGffAttributeValue.pl transcript_id | sort|uniq|wc -l)
detected=$(zcat {input} | fgrep "gene_type \\"$type\\";" | awk '$NF>0' | extractGffAttributeValue.pl transcript_id | sort|uniq|wc -l)
let undetected=$all-$detected || true
echo -e "{wildcards.techname}Corr{wildcards.corrLevel}\t{wildcards.capDesign}\t{wildcards.sizeFrac}\t{wildcards.barcodes}\t$type\t$all\t$detected" | awk '{{print $0"\t"$7/$6}}'
done > {output}
		'''

rule aggTargetCoverageStats:
	input: expand(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.targetCoverage.stats.tsv",filtered_product_merge, techname=TECHNAMESplusMERGED, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACSpluSMERGED, barcodes=BARCODESpluSMERGED)
#lambda wildcards: expand(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_pooled.targetCoverage.stats.tsv", techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS)
	output: config["STATSDATADIR"] + "all.targetCoverage.stats.tsv"
	shell:
		'''
echo -e "seqTech\tcorrectionLevel\tcapDesign\tsizeFrac\ttissue\ttargetType\ttotalTargets\tdetectedTargets\tpercentDetectedTargets" > {output}
cat {input} | grep -v erccSpikein | sed 's/Corr0/\tNo/' | sed 's/Corr{lastK}/\tYes/'| sort >> {output}
		'''

rule plotTargetCoverageStats:
	input: config["STATSDATADIR"] + "all.targetCoverage.stats.tsv"
	output: config["PLOTSDIR"] + "targetCoverage.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.targetCoverage.stats.{ext}"
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


plotWidth = plotWidth + 1
plotHeight = plotHeight + 1

ggplot(dat, aes(x=factor(correctionLevel), y=percentDetectedTargets, fill=targetType)) +
geom_bar(width=0.75,stat='identity', position=position_dodge(width=0.9)) +
scale_fill_manual(values={long_Rpalette}) +
facet_grid( seqTech +sizeFrac ~ capDesign + tissue) +
geom_hline(aes(yintercept=1), linetype='dashed', alpha=0.7) +
geom_text(size=geom_textSize, aes(group=targetType, y=0.01, label = paste(sep='',percent(percentDetectedTargets),' / ','(',comma(detectedTargets),')')), angle=90, size=2.5, hjust=0, vjust=0.5, position = position_dodge(width=0.9)) +
ylab('% targeted regions detected') +
xlab('{params.filterDat[6]}') +
scale_y_continuous(limits = c(0, 1), labels = scales::percent)+
{params.filterDat[7]}
{GGPLOT_PUB_QUALITY}
ggsave('{output}', width=plotWidth, height=plotHeight)
" > {output}.r
cat {output}.r | R --slave


		'''


rule gffcompareToAnnotation:
	input:
		annot=lambda wildcards: CAPDESIGNTOANNOTGTF[wildcards.capDesign],
		tm="mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.gff"
	output: "mappings/nonAnchoredMergeReads/gffcompare/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.vs.gencode.simple.tsv"
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
			tm="mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.gff"
		output: "mappings/nonAnchoredMergeReads/gffcompare/SIRVs/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.vs.SIRVs.simple.tsv"
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
	input:"mappings/nonAnchoredMergeReads/gffcompare/SIRVs/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.vs.SIRVs.simple.tsv"
	output: temp(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.vs.SIRVs.stats.tsv")
	shell:
		'''
file=$(dirname {input})/$(basename {input} .simple.tsv)
for level in `echo Baselevel Exonlevel Intronchainlevel Intronlevel Locuslevel Transcriptlevel`; do
Sn=`cat $file |grep "level:" |sed 's/ //g'| sed 's/:/\\t/'|sed 's/|$//'|sed 's/|/\\t/g' | awk -v l=$level '$1==l' |cut -f2` || Sn='NaN'
Sp=`cat $file |grep "level:" |sed 's/ //g'| sed 's/:/\\t/'|sed 's/|$//'|sed 's/|/\\t/g' | awk -v l=$level '$1==l' |cut -f3` || Sp='NaN'
echo -e "{wildcards.techname}Corr{wildcards.corrLevel}\t{wildcards.capDesign}\t{wildcards.sizeFrac}\t{wildcards.barcodes}\t$level\t$Sn\t$Sp";
done |sed 's/level//g' > {output}

		'''

rule aggGffCompareSirvStats:
	input:expand(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.vs.SIRVs.stats.tsv", filtered_product_merge, techname=TECHNAMESplusMERGED, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACSpluSMERGED, barcodes=BARCODESpluSMERGED)
	output: config["STATSDATADIR"] + "all.tmerge.vs.SIRVs.stats.tsv"
	shell:
		'''
echo -e "seqTech\tcorrectionLevel\tcapDesign\tsizeFrac\ttissue\tlevel\tmetric\tvalue" > {output}
cat {input} | awk '{{print $1"\\t"$2"\\t"$3"\\t"$4"\\t"$5"\\tSn\\t"$6"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\t"$5"\\tPr\\t"$7}}'| sed 's/Corr0/\tNo/' | sed 's/Corr{lastK}/\tYes/' | sort >> {output}

		'''

rule plotGffCompareSirvStats:
	input:config["STATSDATADIR"] + "all.tmerge.vs.SIRVs.stats.tsv"
	output: config["PLOTSDIR"] + "tmerge.vs.SIRVs.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.vs.SIRVs.stats.{ext}"
	params:
		filterDat=lambda wildcards: merge_figures_params(wildcards.capDesign, wildcards.sizeFrac, wildcards.barcodes, wildcards.corrLevel, wildcards.techname)
	shell:
		'''
echo "library(ggplot2)
library(plyr)
library(scales)
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

dat\$levelCorrlevel <- paste(sep='', dat\$level, ' (Corr: ', dat\$correctionLevel, ')')
plotHeight = plotHeight +1
ggplot(dat, aes(x=levelCorrlevel, y=value)) +
geom_point(aes(color=metric, shape=correctionLevel), size=4, alpha=0.8) +
scale_colour_manual (values=cbPalette, name='Metric', breaks=c('Sn', 'Pr'))+
scale_shape_manual(values=c(16,21), name='Error correction') +
ylab('Sn | Pr (%)') +
xlab('Evaluation level / (Error correction)') +
scale_y_continuous() +
expand_limits(y=c(0,100))+
theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
facet_grid( seqTech + sizeFrac ~ capDesign + tissue)+
{GGPLOT_PUB_QUALITY}
ggsave('{output}', width=plotWidth, height=plotHeight)
" > {output}.r
cat {output}.r | R --slave
		'''

if SIRVpresent:
	rule getSirvDetectionStats:
		input:
			gffC="mappings/nonAnchoredMergeReads/gffcompare/SIRVs/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.vs.SIRVs.simple.tsv",
			tm="mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.gff",
			sirvInfo=config["SIRVinfo"]
		output:temp(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.vs.SIRVs.detection.stats.tsv")
		shell:
			'''
uuid=$(uuidgen)
sirvDetectionStats.pl {input.sirvInfo} $(dirname {input.gffC})/$(basename {input.gffC} .simple.tsv).$(basename {input.tm}).refmap > $TMPDIR/$uuid
cat $TMPDIR/$uuid | while read id l c ca; do echo -e "{wildcards.techname}Corr{wildcards.corrLevel}\t{wildcards.capDesign}\t{wildcards.sizeFrac}\t{wildcards.barcodes}\t$id\t$l\t$c\t$ca"; done > {output}

			'''

rule aggSirvDetectionStats:
	input:expand(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.vs.SIRVs.detection.stats.tsv", filtered_product_merge, techname=TECHNAMESplusMERGED, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACSpluSMERGED, barcodes=BARCODESpluSMERGED)
	output: config["STATSDATADIR"] + "all.tmerge.vs.SIRVs.detection.stats.tsv"
	shell:
		'''
echo -e "seqTech\tcorrectionLevel\tcapDesign\tsizeFrac\ttissue\tSIRVid\tlength\tconcentration\tdetectionStatus" > {output}
cat {input} | sed 's/Corr0/\tNo/' | sed 's/Corr{lastK}/\tYes/' | sort >> {output}

		'''

rule plotSirvDetectionStats:
	input:config["STATSDATADIR"] + "all.tmerge.vs.SIRVs.detection.stats.tsv"
	output: config["PLOTSDIR"] + "tmerge.vs.SIRVs.detection.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.vs.SIRVs.detection.stats.{ext}"
	params:
		filterDat=lambda wildcards: merge_figures_params(wildcards.capDesign, wildcards.sizeFrac, wildcards.barcodes, wildcards.corrLevel, wildcards.techname)
	shell:
		'''
echo "library(ggplot2)
library(plyr)
library(scales)
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


ggplot(dat, aes(x=concentration, y=length, color=detectionStatus, shape=correctionLevel)) + geom_point(alpha=0.23) +
coord_trans(x='log2') +
scale_color_manual(values=palette) +
scale_shape_manual(values=c('Yes'=16, 'No'=17), name='Error correction') +
xlab('SIRV molarity (fmol/uL)') +
ylab('SIRV length (nt)') +
facet_grid( seqTech + sizeFrac ~ capDesign + tissue)+
{GGPLOT_PUB_QUALITY}
ggsave('{output}', width=plotWidth, height=plotHeight)
" > {output}.r
cat {output}.r | R --slave


		'''


rule colorBedAccordingToGffCompare:
	input:
		classes="mappings/nonAnchoredMergeReads/gffcompare/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.vs.gencode.simple.tsv",
		tm="mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.bed"
	output: "mappings/nonAnchoredMergeReads/colored/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.bed"
	shell:
			'''
colorNovelTxBed.pl {input.classes} {input.tm} > {output}

			'''

rule colorFLBedAccordingToGffCompare:
	input:
		classes="mappings/nonAnchoredMergeReads/gffcompare/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.vs.gencode.simple.tsv",
		tm="mappings/nonAnchoredMergeReads/cage+polyASupported/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.cage+polyASupported.bed"
	output: "mappings/nonAnchoredMergeReads/cage+polyASupported/colored/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}_FL.tmerge.bed"
	shell:
			'''
colorNovelTxBed.pl {input.classes} {input.tm} > {output}

			'''



rule getGffCompareStats:
	input: "mappings/nonAnchoredMergeReads/gffcompare/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.vs.gencode.simple.tsv"
	output: temp(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.vs.gencode.stats.tsv")
	shell:
		'''
cat {input} |cut -f2 | sort|uniq -c | awk -v s={wildcards.techname}Corr{wildcards.corrLevel} -v c={wildcards.capDesign} -v si={wildcards.sizeFrac} -v b={wildcards.barcodes} '{{print s"\t"c"\t"si"\t"b"\t"$2"\t"$1}}' | sed 's/Corr0/\tNo/' | sed 's/Corr{lastK}/\tYes/'>> {output}
		'''


#echo -e "seqTech\tcorrectionLevel\tcapDesign\tcategory\tcount" > {output}

rule aggGffCompareStats:
	input: expand(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.vs.gencode.stats.tsv",filtered_product_merge, techname=TECHNAMESplusMERGED, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACSpluSMERGED, barcodes=BARCODESpluSMERGED)
	output: config["STATSDATADIR"] + "all.tmerge.vs.gencode.stats.tsv"
	shell:
		'''
echo -e "seqTech\tcorrectionLevel\tcapDesign\tsizeFrac\ttissue\tcategory\tcount" > {output}
cat {input} >> {output}
		'''

rule plotGffCompareStats:
	input: config["STATSDATADIR"] + "all.tmerge.vs.gencode.stats.tsv"
	output: config["PLOTSDIR"] + "tmerge.vs.gencode.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.vs.gencode.stats.{ext}"
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


dat\$category<-factor(dat\$category, ordered=TRUE, levels=rev(c('Intergenic', 'Extends', 'Intronic', 'Overlaps', 'Antisense', 'Equal', 'Included')))
palette <- c('Intergenic' = '#0099cc', 'Extends' ='#00bfff', 'Intronic' = '#4dd2ff', 'Overlaps' = '#80dfff', 'Antisense' = '#ccf2ff', 'Equal' = '#c65353', 'Included' ='#d98c8c')

ggplot(dat[order(dat\$category), ], aes(x=factor(correctionLevel), y=count, fill=category)) +
geom_bar(stat='identity') +
scale_fill_manual(values=palette) +
facet_grid( seqTech + sizeFrac ~ capDesign + tissue)+
 ylab('# TMs') +
xlab('{params.filterDat[6]}') +
guides(fill = guide_legend(title='Category'))+
scale_y_continuous(labels=scientific)+
{params.filterDat[7]}
{GGPLOT_PUB_QUALITY}
ggsave('{output}', width=plotWidth, height=plotHeight)
" > {output}.r
cat {output}.r | R --slave

		'''

rule simplifyGencode:
	input: lambda wildcards: CAPDESIGNTOANNOTGTF[wildcards.capDesign]
	output: "annotations/simplified/{capDesign}.gencode.simplified_biotypes.gtf"
	shell:
		'''
cat {input}  | simplifyGencodeGeneTypes.pl - | sortgff > {output}
		'''

rule mergeTmsWithGencode:
	input:
		annot="annotations/simplified/{capDesign}.gencode.simplified_biotypes.gtf",
		tm="mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.gff"
	output: temp("mappings/nonAnchoredMergeReads/gencodeMerge/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge+gencode.gff.gz")
	threads:8
	shell:
		'''
cat {input.annot} {input.tm}  | skipcomments | sortgff | tmerge --cpu {threads} --exonOverhangTolerance 25 - |sortgff |gzip > {output}
		'''

rule makeClsGencodeLoci:
	input: "mappings/nonAnchoredMergeReads/gencodeMerge/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge+gencode.gff.gz"
	params: locusPrefix=config["PROJECT_NAME"]
	output: temp("mappings/nonAnchoredMergeReads/gencodeLociMerge/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge+gencode.loci.gff.gz")
	shell:
		'''
uuid=$(uuidgen)
zcat {input} > $TMPDIR/$uuid
bedtools intersect -s -wao -a $TMPDIR/$uuid -b $TMPDIR/$uuid |fgrep -v ERCC| buildLoci.pl --locPrefix {params.locusPrefix}: - |sortgff | gzip> {output}

		'''


rule mergeWithRef:
	input:
		clsGencode="mappings/nonAnchoredMergeReads/gencodeLociMerge/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge+gencode.loci.gff.gz",
		gencode="annotations/simplified/{capDesign}.gencode.simplified_biotypes.gtf"
	output: "mappings/nonAnchoredMergeReads/mergeToRef/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge+gencode.loci.refmerged.gff.gz"
	shell:
		'''
uuid=$(uuidgen)
zcat  {input.clsGencode} > $TMPDIR/$uuid
mergeToRef.pl {input.gencode} $TMPDIR/$uuid | sortgff |gzip > {output}
		'''

rule getNovelIntergenicLoci:
	input:
		gencode="annotations/simplified/{capDesign}.gencode.simplified_biotypes.gtf",
		tmergeGencode="mappings/nonAnchoredMergeReads/mergeToRef/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge+gencode.loci.refmerged.gff.gz"
	output:"mappings/nonAnchoredMergeReads/mergeToRef/novelLoci/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.novelLoci.gff.gz"
	shell:
		'''
uuid1=$(uuidgen)
uuid2=$(uuidgen)
uuid3=$(uuidgen)
cat {input.gencode} |awk '$3=="exon"' | extract_locus_coords.pl -| sortbed > $TMPDIR/$uuid1
zcat {input.tmergeGencode} | fgrep 'gene_ref_status "novel";' | extract_locus_coords.pl - | sortbed > $TMPDIR/$uuid2
bedtools intersect -v -a $TMPDIR/$uuid2 -b $TMPDIR/$uuid1 |fgrep -v ERCC |cut -f4 | sort|uniq > $TMPDIR/$uuid3
zcat {input.tmergeGencode}| fgrep -w -f $TMPDIR/$uuid3 - |gzip > {output}
		'''

rule getNovelIntergenicLociStats:
	input:
		tmergeGencode="mappings/nonAnchoredMergeReads/mergeToRef/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge+gencode.loci.refmerged.gff.gz",
		intergenic="mappings/nonAnchoredMergeReads/mergeToRef/novelLoci/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.novelLoci.gff.gz"
	output: temp(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.novelLoci.stats.tsv")
	shell:
		'''
totalNovel=$(zcat {input.tmergeGencode} | fgrep 'gene_ref_status "novel";' |extractGffAttributeValue.pl gene_id | sort| uniq | wc -l)
interg=$(zcat {input.intergenic} | extractGffAttributeValue.pl gene_id | sort| uniq | wc -l)
echo -e "{wildcards.techname}Corr{wildcards.corrLevel}\t{wildcards.capDesign}\t{wildcards.sizeFrac}\t{wildcards.barcodes}\t$totalNovel\t$interg"  > {output}

		'''

rule aggNovelIntergenicLociStats:
	input: expand(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.novelLoci.stats.tsv",filtered_product_merge, techname=TECHNAMESplusMERGED, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACSpluSMERGED, barcodes=BARCODESpluSMERGED)
	output: config["STATSDATADIR"] + "all.tmerge.novelLoci.stats.tsv"
	shell:
		'''
echo -e "seqTech\tcorrectionLevel\tcapDesign\tsizeFrac\ttissue\tcategory\tcount\tpercent" > {output}
cat {input} | awk '{{print $1"\\t"$2"\\t"$3"\\t"$4"\\tintergenic\\t"$6"\\t"$6/$5"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\tintronic\\t"$5-$6"\\t"($5-$6)/$5}}'| sed 's/Corr0/\tNo/' | sed 's/Corr{lastK}/\tYes/' | sort >> {output}
		'''


rule plotNovelIntergenicLociStats:
	input: config["STATSDATADIR"] + "all.tmerge.novelLoci.stats.tsv"
	output: config["PLOTSDIR"] + "tmerge.novelLoci.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.novelLoci.stats.{ext}"
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

dat\$category<-factor(dat\$category, ordered=TRUE, levels=rev(c('intronic', 'intergenic')))
ggplot(dat[order(dat\$category), ], aes(x=factor(correctionLevel), y=count, fill=category)) +
geom_bar(stat='identity') +
ylab('# Novel CLS loci') +
scale_y_continuous(labels=comma)+
scale_fill_manual (values=c(intronic='#d98c8c', intergenic='#33ccff'))+
facet_grid( seqTech + sizeFrac ~ capDesign + tissue)+
xlab('{params.filterDat[6]}') +
guides(fill = guide_legend(title='Category\\n(w.r.t. GENCODE)'))+
geom_text(position = 'stack', size=geom_textSize, aes(x = factor(correctionLevel), y = count, ymax=count, label = paste(sep='',percent(round(percent, digits=2)),' / ','(',comma(count),')'), hjust = 0.5, vjust = 1))+
{params.filterDat[7]}
{GGPLOT_PUB_QUALITY}
ggsave('{output}', width=plotWidth, height=plotHeight)
" > {output}.r
cat {output}.r | R --slave

		'''

