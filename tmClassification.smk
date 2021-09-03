rule compareTargetsToTms:
	input:
		tms= "output/mappings/mergedReads/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.HiSS.tmerge.min{minReadSupport}reads.splicing_status-all.endSupport-all.gff.gz",
		targetedSegments=lambda wildcards: CAPDESIGNTOTARGETSGFF[CAPDESIGNTOCAPDESIGN[wildcards.capDesign]]
	conda: "envs/xtools_env.yml"
	output: "output/mappings/mergedReads/vsTargets/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.HiSS.tmerge.min{minReadSupport}reads.gfftsv.gz"
	shell:
		'''
uuidTmpOut=$(uuidgen)

zcat {input.tms} > {TMPDIR}/$uuidTmpOut.1
bedtools intersect -wao -a {input.targetedSegments} -b {TMPDIR}/$uuidTmpOut.1 |gzip > {TMPDIR}/$uuidTmpOut
mv {TMPDIR}/$uuidTmpOut {output}

		'''

rule getTargetCoverageStats:
	input: "output/mappings/mergedReads/vsTargets/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.HiSS.tmerge.min{minReadSupport}reads.gfftsv.gz"
	output: "output/statsFiles/" + "tmp/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.min{minReadSupport}reads.targetCoverage.stats.tsv"
	shell:
		'''
uuidTmpOut=$(uuidgen)
for type in `zcat {input} | extractGffAttributeValue.pl gene_type | sort -T {TMPDIR} |uniq`; do
all=$(zcat {input} | fgrep "gene_type \\"$type\\";" | extractGffAttributeValue.pl transcript_id | sort -T {TMPDIR} |uniq|wc -l)
detected=$(zcat {input} | fgrep "gene_type \\"$type\\";" | awk '$NF>0' | extractGffAttributeValue.pl transcript_id | sort -T {TMPDIR} |uniq|wc -l)
let undetected=$all-$detected || true
echo -e "{wildcards.techname}\t{wildcards.capDesign}\t{wildcards.sizeFrac}\t{wildcards.sampleRep}\t$type\t$all\t$detected" | awk '{{print $0"\t"$7/$6}}'
done > {TMPDIR}/$uuidTmpOut
mv {TMPDIR}/$uuidTmpOut {output}
		'''

rule aggTargetCoverageStats:
	input: lambda wildcards: expand("output/statsFiles/" + "tmp/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.min{minReadSupport}reads.targetCoverage.stats.tsv",filtered_product,  techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, sampleRep=SAMPLEREPS, minReadSupport=wildcards.minReadSupport)
	output: "output/statsFiles/" + "all.min{minReadSupport}reads.targetCoverage.stats.tsv"
	shell:
		'''
uuidTmpOut=$(uuidgen)
echo -e "seqTech\tcapDesign\tsizeFrac\tsampleRep\ttargetType\ttotalTargets\tdetectedTargets\tpercentDetectedTargets" > {TMPDIR}/$uuidTmpOut
cat {input} | grep -v erccSpikein | sort -T {TMPDIR}  >> {TMPDIR}/$uuidTmpOut
mv {TMPDIR}/$uuidTmpOut {output}
		'''

rule plotTargetCoverageStats:
	input: "output/statsFiles/" + "all.min{minReadSupport}reads.targetCoverage.stats.tsv"
	output: returnPlotFilenames("output/plots/" + "targetCoverage.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.min{minReadSupport}reads.targetCoverage.stats")
	conda: "envs/R_env.yml"
	params:
		filterDat=lambda wildcards: multi_figures(wildcards.capDesign, wildcards.sizeFrac, wildcards.sampleRep, wildcards.techname)
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
{params.filterDat[sampleRepFilterString]}
{params.filterDat[substSeqTechString]}
{params.filterDat[substSampleRepString]}
{params.filterDat[graphDimensions]}


wXyPlot = wXyPlot * 1.2
hXyPlot = hXyPlot * 1.2

plotBase <- \\"p <- ggplot(dat, aes(x=1, y=percentDetectedTargets, fill=targetType)) +
geom_bar(width=0.75,stat='identity', position=position_dodge(width=0.9)) +
scale_fill_manual(values={long_Rpalette}) +
geom_hline(aes(yintercept=1), linetype='dashed', alpha=0.7, size=lineSize) +
geom_text(size=geom_textSize, aes(group=targetType, y=0.01, label = paste(sep='',percent(percentDetectedTargets),' / ','(',comma(detectedTargets),')')), angle=90, size=2.5, hjust=0, vjust=0.5, position = position_dodge(width=0.9)) +
ylab('% targeted regions detected') +
xlab('') +
scale_y_continuous(limits = c(0, 1), labels = scales::percent)+
{params.filterDat[hideXaxisLabels]}
{GGPLOT_PUB_QUALITY} + \\"

{params.filterDat[facetPlotSetup]}

save_plot('{output[0]}', legendOnly, base_width=wLegendOnly, base_height=hLegendOnly)
save_plot('{output[1]}', pXy, base_width=wXyPlot, base_height=hXyPlot)

save_plot('{output[2]}', pXyNoLegend, base_width=wXyNoLegendPlot, base_height=hXyNoLegendPlot)
save_plot('{output[3]}', pYx, base_width=wYxPlot, base_height=hYxPlot)

save_plot('{output[4]}', pYxNoLegend, base_width=wYxNoLegendPlot, base_height=hYxNoLegendPlot)


" > $(dirname {output[0]})/$(basename {output[0]} .legendOnly.png).r

cat $(dirname {output[0]})/$(basename {output[0]} .legendOnly.png).r | R --slave


		'''


rule gffcompareToAnnotation:
	input:
		annot=lambda wildcards: GENOMETOANNOTGTF[CAPDESIGNTOGENOME[wildcards.capDesign]],
		tm="output/mappings/mergedReads/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.HiSS.tmerge.min{minReadSupport}reads.splicing_status-{splicedStatus}.endSupport-{endSupport}.gff.gz"
	output: 
		standard="output/mappings/mergedReads/gffcompare/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.splicing_status-{splicedStatus}.endSupport-{endSupport}.vs.gencode.simple.tsv",
		adjustedSn="output/mappings/mergedReads/gffcompare/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.splicing_status-{splicedStatus}.endSupport-{endSupport}.vs.gencode.adj.simple.tsv"

	shell:
		'''
pref=$(basename {output.standard} .simple.tsv)
annotFullPath=$(fullpath {input.annot})
uuid=$(uuidgen)
outdir="$PWD/$(dirname {output.standard})"
zcat {input.tm} > {TMPDIR}/$uuid
cd {TMPDIR}
gffcompare -o ${{uuid}}PREF -r $annotFullPath $uuid
cat ${{uuid}}PREF.tracking | simplifyGffCompareClasses.pl - > ${{uuid}}PREF.simple.tsv

mv ${{uuid}}PREF.simple.tsv $outdir/$(basename {output.standard})
mv ${{uuid}}PREF.loci $outdir/$pref.loci
mv ${{uuid}}PREF.stats $outdir/$pref
mv ${{uuid}}PREF.tracking $outdir/$pref.tracking
mv ${{uuid}}PREF.$uuid.refmap $outdir/$pref.refmap
mv ${{uuid}}PREF.$uuid.tmap $outdir/$pref.tmap

##### ADJUSTED METRICS (gffcompare -R)

pref=$(basename {output.adjustedSn} .simple.tsv)
gffcompare -o ${{uuid}}PREF -r $annotFullPath -R $uuid
cat ${{uuid}}PREF.tracking | simplifyGffCompareClasses.pl - > ${{uuid}}PREF.simple.tsv
mv ${{uuid}}PREF.simple.tsv $outdir/$(basename {output.adjustedSn})
mv ${{uuid}}PREF.loci $outdir/$pref.loci
mv ${{uuid}}PREF.stats $outdir/$pref
mv ${{uuid}}PREF.tracking $outdir/$pref.tracking
mv ${{uuid}}PREF.$uuid.refmap $outdir/$pref.refmap
mv ${{uuid}}PREF.$uuid.tmap $outdir/$pref.tmap

		'''

if SIRVpresent:
	rule gffcompareToSirvAnnotation:
		input:
			annot="output/annotations/{capDesign}.SIRVs.gff",
			tm="output/mappings/mergedReads/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.{filt}.tmerge.min{minReadSupport}reads.splicing_status-all.endSupport-all.gff.gz"
		output: "output/mappings/mergedReads/gffcompare/SIRVs/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.{filt}.tmerge.min{minReadSupport}reads.vs.SIRVs.simple.tsv"
		shell:
			'''
pref=$(basename {output} .simple.tsv)
annotFullPath=$(fullpath {input.annot})
uuid=$(uuidgen)
outdir="$PWD/$(dirname {output})"
zcat {input.tm} | awk '$1 ~ /SIRV/ ' > {TMPDIR}/$uuid
cd {TMPDIR}
gffcompare -o ${{uuid}}PREF -r $annotFullPath $uuid
cat ${{uuid}}PREF.tracking | simplifyGffCompareClasses.pl - > ${{uuid}}PREF.simple.tsv
mv ${{uuid}}PREF.simple.tsv $outdir/$(basename {output})
mv ${{uuid}}PREF.loci $outdir/$pref.loci
mv ${{uuid}}PREF.stats $outdir/$pref
mv ${{uuid}}PREF.tracking $outdir/$pref.tracking
mv ${{uuid}}PREF.$uuid.refmap $outdir/$pref.refmap
mv ${{uuid}}PREF.$uuid.tmap $outdir/$pref.tmap


			'''

rule getGffCompareSirvStats:
	input:"output/mappings/mergedReads/gffcompare/SIRVs/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.{filt}.tmerge.min{minReadSupport}reads.vs.SIRVs.simple.tsv"
	output: "output/statsFiles/" + "tmp/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.{filt}.tmerge.min{minReadSupport}reads.vs.SIRVs.stats.tsv"
	shell:
		'''
uuidTmpOut=$(uuidgen)
uuid=$(uuidgen)
file=$(dirname {input})/$(basename {input} .simple.tsv)
for level in `echo Baselevel Exonlevel Intronchainlevel Intronlevel Locuslevel Transcriptlevel`; do

cat $file |grep "level:" > {TMPDIR}/$uuid || true
SnDEFAULT=0
SpDEFAULT='NA'

Sn=`cat {TMPDIR}/$uuid |sed 's/ //g'| sed 's/:/\\t/'|sed 's/|$//'|sed 's/|/\\t/g' | awk -v l=$level '$1==l' |cut -f2`
Sn=${{Sn:-$SnDEFAULT}} #see https://vaneyckt.io/posts/safer_bash_scripts_with_set_euxo_pipefail/ default variable assignments


Sp=`cat {TMPDIR}/$uuid |sed 's/ //g'| sed 's/:/\\t/'|sed 's/|$//'|sed 's/|/\\t/g' | awk -v l=$level '$1==l' |cut -f3`
Sp=${{Sp:-$SpDEFAULT}} #see https://vaneyckt.io/posts/safer_bash_scripts_with_set_euxo_pipefail/ default variable assignments


echo -e "{wildcards.techname}\t{wildcards.capDesign}\t{wildcards.sizeFrac}\t{wildcards.sampleRep}\t$level\t$Sn\t$Sp";
done |sed 's/level//g' > {TMPDIR}/$uuidTmpOut
mv {TMPDIR}/$uuidTmpOut {output}

		'''

rule aggGffCompareSirvStats:
	input: lambda wildcards:expand("output/statsFiles/" + "tmp/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.{filt}.tmerge.min{minReadSupport}reads.vs.SIRVs.stats.tsv", filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, sampleRep=SAMPLEREPS, minReadSupport=wildcards.minReadSupport, filt=wildcards.filt)
	output: "output/statsFiles/" + "all.{filt}.tmerge.min{minReadSupport}reads.vs.SIRVs.stats.tsv"
	shell:
		'''
uuidTmpOut=$(uuidgen)
echo -e "seqTech\tcapDesign\tsizeFrac\tsampleRep\tlevel\tmetric\tvalue" > {TMPDIR}/$uuidTmpOut
cat {input} | awk '{{print $1"\\t"$2"\\t"$3"\\t"$4"\\t"$5"\\tSn\\t"$6"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\t"$5"\\tPr\\t"$7}}' | sort -T {TMPDIR}  >> {TMPDIR}/$uuidTmpOut
mv {TMPDIR}/$uuidTmpOut {output}

		'''

rule plotGffCompareSirvStats:
	input:"output/statsFiles/" + "all.{filt}.tmerge.min{minReadSupport}reads.vs.SIRVs.stats.tsv"
	output: returnPlotFilenames("output/plots/" + "tmerge.vs.SIRVs.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.{filt}.tmerge.min{minReadSupport}reads.vs.SIRVs.stats")
	conda: "envs/R_env.yml"
	params:
		filterDat=lambda wildcards: multi_figures(wildcards.capDesign, wildcards.sizeFrac, wildcards.sampleRep, wildcards.techname)
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
library(ggforce)

cbPalette <- c('Sn'='#ffb366', 'Pr'='#2d8659')
dat <- read.table('{input}', header=T, as.is=T, sep='\\t')
{params.filterDat[technameFilterString]}
{params.filterDat[capDesignFilterString]}

{params.filterDat[sizeFracFilterString]}
{params.filterDat[sampleRepFilterString]}
{params.filterDat[substSeqTechString]}
{params.filterDat[substSampleRepString]}
{params.filterDat[graphDimensions]}

plotBase <- \\"p <- ggplot(dat, aes(x=level, y=value)) +
#geom_mark_rect(aes(filter = level == 'Transcript' & metric == 'Pr'), size=2, expand = unit(7, 'mm'), radius = unit(4, 'mm'), color='#4d4d4d', fill='#ffff00') +
geom_point(aes(color=metric), shape=18, size=lineSize*15, alpha=0.7) +
scale_colour_manual (values=cbPalette, name='Metric', breaks=c('Sn', 'Pr'))+
ylab('Sn | Pr (%)') +
xlab('Evaluation level') +
scale_x_discrete(limits=c('Base', 'Locus', 'Exon', 'Intron', 'Intronchain', 'Transcript')) +
scale_y_continuous() +
expand_limits(y=c(0,100))+
theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
{GGPLOT_PUB_QUALITY} + \\"

{params.filterDat[facetPlotSetup]}

save_plot('{output[0]}', legendOnly, base_width=wLegendOnly, base_height=hLegendOnly)
save_plot('{output[1]}', pXy, base_width=wXyPlot, base_height=hXyPlot)

save_plot('{output[2]}', pXyNoLegend, base_width=wXyNoLegendPlot, base_height=hXyNoLegendPlot)
save_plot('{output[3]}', pYx, base_width=wYxPlot, base_height=hYxPlot)

save_plot('{output[4]}', pYxNoLegend, base_width=wYxNoLegendPlot, base_height=hYxNoLegendPlot)


" > $(dirname {output[0]})/$(basename {output[0]} .legendOnly.png).r

cat $(dirname {output[0]})/$(basename {output[0]} .legendOnly.png).r | R --slave
		'''

if SIRVpresent:
	rule getSirvDetectionStats:
		input:
			gffC="output/mappings/mergedReads/gffcompare/SIRVs/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.HiSS.tmerge.min{minReadSupport}reads.vs.SIRVs.simple.tsv",
			sirvInfo=config["SIRVinfo"]
		output:"output/statsFiles/" + "tmp/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.vs.SIRVs.detection.stats.tsv"
		shell:
			'''
uuid=$(uuidgen)
uuidTmpOut=$(uuidgen)
sirvDetectionStats.pl {input.sirvInfo} $(dirname {input.gffC})/$(basename {input.gffC} .simple.tsv).refmap > {TMPDIR}/$uuid
cat {TMPDIR}/$uuid | while read id l c ca; do echo -e "{wildcards.techname}\t{wildcards.capDesign}\t{wildcards.sizeFrac}\t{wildcards.sampleRep}\t$id\t$l\t$c\t$ca"; done > {TMPDIR}/$uuidTmpOut
mv {TMPDIR}/$uuidTmpOut {output}

			'''

rule aggSirvDetectionStats:
	input: lambda wildcards:expand("output/statsFiles/" + "tmp/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.vs.SIRVs.detection.stats.tsv", filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, sampleRep=SAMPLEREPS, minReadSupport=wildcards.minReadSupport)
	output: "output/statsFiles/" + "all.tmerge.min{minReadSupport}reads.vs.SIRVs.detection.stats.tsv"
	shell:
		'''
uuidTmpOut=$(uuidgen)
echo -e "seqTech\tcapDesign\tsizeFrac\tsampleRep\tSIRVid\tlength\tconcentration\tdetectionStatus" > {TMPDIR}/$uuidTmpOut
cat {input} | sort -T {TMPDIR}  >> {TMPDIR}/$uuidTmpOut
mv {TMPDIR}/$uuidTmpOut {output}

		'''

rule plotSirvDetectionStats:
	input:"output/statsFiles/" + "all.tmerge.min{minReadSupport}reads.vs.SIRVs.detection.stats.tsv"
	output: returnPlotFilenames("output/plots/" + "tmerge.vs.SIRVs.detection.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.vs.SIRVs.detection.stats")
	conda: "envs/R_env.yml"
	params:
		filterDat=lambda wildcards: multi_figures(wildcards.capDesign, wildcards.sizeFrac, wildcards.sampleRep, wildcards.techname)
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

palette <- c('end-to-end' = '#00e600', 'absent' = '#666666', 'partial' = '#ff0066')
dat <- read.table('{input}', header=T, as.is=T, sep='\\t')

{params.filterDat[technameFilterString]}
{params.filterDat[capDesignFilterString]}

{params.filterDat[sizeFracFilterString]}
{params.filterDat[sampleRepFilterString]}
{params.filterDat[substSeqTechString]}
{params.filterDat[substSampleRepString]}
{params.filterDat[graphDimensions]}


plotBase <- \\"p <- ggplot(dat, aes(x=concentration, y=length, color=detectionStatus)) + geom_point(alpha=0.8, shape=18) +
coord_trans(x='log2') +
scale_color_manual(values=palette) +
xlab('SIRV molarity (fmol/uL)') +
ylab('SIRV length (nt)') +
{GGPLOT_PUB_QUALITY} + \\"

{params.filterDat[facetPlotSetup]}

save_plot('{output[0]}', legendOnly, base_width=wLegendOnly, base_height=hLegendOnly)
save_plot('{output[1]}', pXy, base_width=wXyPlot, base_height=hXyPlot)

save_plot('{output[2]}', pXyNoLegend, base_width=wXyNoLegendPlot, base_height=hXyNoLegendPlot)
save_plot('{output[3]}', pYx, base_width=wYxPlot, base_height=hYxPlot)

save_plot('{output[4]}', pYxNoLegend, base_width=wYxNoLegendPlot, base_height=hYxNoLegendPlot)


" > $(dirname {output[0]})/$(basename {output[0]} .legendOnly.png).r

cat $(dirname {output[0]})/$(basename {output[0]} .legendOnly.png).r | R --slave



		'''


rule colorBedAccordingToGffCompare:
	input:
		classes="output/mappings/mergedReads/gffcompare/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.splicing_status-all.endSupport-{endSupport}.vs.gencode.simple.tsv",
		tm="output/mappings/mergedReads/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.splicing_status-all.endSupport-{endSupport}.bed"
	output: "output/mappings/mergedReads/colored/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.splicing_status-all.endSupport-{endSupport}.bed"
	shell:
			'''
uuidTmpOut=$(uuidgen)
colorNovelTxBed.pl {input.classes} {input.tm} > {TMPDIR}/$uuidTmpOut
mv {TMPDIR}/$uuidTmpOut {output}

			'''


rule getGffCompareStats:
	input: "output/mappings/mergedReads/gffcompare/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.splicing_status-{splicedStatus}.endSupport-{endSupport}.vs.gencode.simple.tsv"
	output: "output/statsFiles/" + "tmp/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.splicing_status-{splicedStatus}.endSupport-{endSupport}.vs.gencode.stats.tsv"
	shell:
		'''
uuidTmpOut=$(uuidgen)
cat {input} |cut -f2 | sort -T {TMPDIR} |uniq -c | awk -v s={wildcards.techname} -v c={wildcards.capDesign} -v si={wildcards.sizeFrac} -v b={wildcards.sampleRep} '{{print s"\t"c"\t"si"\t"b"\t"$2"\t"$1}}' > {TMPDIR}/$uuidTmpOut
mv {TMPDIR}/$uuidTmpOut {output}
		'''

rule getGffCompareGencodeSnPrStats:
	input:"output/mappings/mergedReads/gffcompare/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.splicing_status-{splicedStatus}.endSupport-{endSupport}.vs.gencode.adj.simple.tsv"
	output: "output/statsFiles/" + "tmp/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.splicing_status-{splicedStatus}.endSupport-{endSupport}.vs.gencode.SnPr.stats.tsv"
	shell:
		'''
uuidTmpOut=$(uuidgen)
uuid=$(uuidgen)
file=$(dirname {input})/$(basename {input} .simple.tsv)
for level in `echo Baselevel Exonlevel Intronchainlevel Intronlevel Locuslevel Transcriptlevel`; do
cat $file |grep "level:" > {TMPDIR}/$uuid || true
SnDEFAULT=0
SpDEFAULT='NA'


Sn=`cat {TMPDIR}/$uuid |sed 's/ //g'| sed 's/:/\\t/'|sed 's/|$//'|sed 's/|/\\t/g' | awk -v l=$level '$1==l' |cut -f2` 
Sn=${{Sn:-$SnDEFAULT}} #see https://vaneyckt.io/posts/safer_bash_scripts_with_set_euxo_pipefail/ default variable assignments

Sp=`cat {TMPDIR}/$uuid |sed 's/ //g'| sed 's/:/\\t/'|sed 's/|$//'|sed 's/|/\\t/g' | awk -v l=$level '$1==l' |cut -f3` 
Sp=${{Sp:-$SpDEFAULT}} #see https://vaneyckt.io/posts/safer_bash_scripts_with_set_euxo_pipefail/ default variable assignments

echo -e "{wildcards.techname}\t{wildcards.capDesign}\t{wildcards.sizeFrac}\t{wildcards.sampleRep}\t$level\t$Sn\t$Sp";
done |sed 's/level//g' > {TMPDIR}/$uuidTmpOut
mv {TMPDIR}/$uuidTmpOut {output}

		'''

rule aggGffCompareGencodeSnPrStats:
	input: lambda wildcards:expand("output/statsFiles/" + "tmp/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.splicing_status-{splicedStatus}.endSupport-{endSupport}.vs.gencode.SnPr.stats.tsv", filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, sampleRep=SAMPLEREPS, endSupport=wildcards.endSupport, minReadSupport=wildcards.minReadSupport, splicedStatus=wildcards.splicedStatus)
	output: "output/statsFiles/" + "all.tmerge.min{minReadSupport}reads.splicing_status-{splicedStatus}.endSupport-{endSupport}.vs.gencode.SnPr.stats.tsv"
	shell:
		'''
uuidTmpOut=$(uuidgen)
echo -e "seqTech\tcapDesign\tsizeFrac\tsampleRep\tlevel\tmetric\tvalue" > {TMPDIR}/$uuidTmpOut
cat {input} | awk '{{print $1"\\t"$2"\\t"$3"\\t"$4"\\t"$5"\\tSn\\t"$6"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\t"$5"\\tPr\\t"$7}}' | sort -T {TMPDIR}  >> {TMPDIR}/$uuidTmpOut
mv {TMPDIR}/$uuidTmpOut {output}

		'''

rule plotGffCompareGencodeSnPrStats:
	input:"output/statsFiles/" + "all.tmerge.min{minReadSupport}reads.splicing_status-{splicedStatus}.endSupport-{endSupport}.vs.gencode.SnPr.stats.tsv"
	output: returnPlotFilenames("output/plots/" + "tmerge.vs.gencode.SnPr.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.splicing_status-{splicedStatus}.endSupport-{endSupport}.vs.gencode.SnPr.stats")
	conda: "envs/R_env.yml"
	params:
		filterDat=lambda wildcards: multi_figures(wildcards.capDesign, wildcards.sizeFrac, wildcards.sampleRep, wildcards.techname)
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

cbPalette <- c('Sn'='#ffb366', 'Pr'='#2d8659')
dat <- read.table('{input}', header=T, as.is=T, sep='\\t')
{params.filterDat[technameFilterString]}
{params.filterDat[capDesignFilterString]}

{params.filterDat[sizeFracFilterString]}
{params.filterDat[sampleRepFilterString]}
{params.filterDat[substSeqTechString]}
{params.filterDat[substSampleRepString]}
{params.filterDat[graphDimensions]}

hXyPlot = hXyPlot +1
wXyPlot = wXyPlot +1

plotBase <- \\"p <- ggplot(dat, aes(x=level, y=value)) +
geom_point(aes(color=metric), shape=18, size=lineSize, alpha=0.7) +
scale_colour_manual (values=cbPalette, name='Metric', breaks=c('Sn', 'Pr'))+
ylab('Sn | Pr (%)') +
xlab('Evaluation level') +
scale_y_continuous() +
expand_limits(y=c(0,100))+
theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
{GGPLOT_PUB_QUALITY} + \\"

{params.filterDat[facetPlotSetup]}

save_plot('{output[0]}', legendOnly, base_width=wLegendOnly, base_height=hLegendOnly)
save_plot('{output[1]}', pXy, base_width=wXyPlot, base_height=hXyPlot)

save_plot('{output[2]}', pXyNoLegend, base_width=wXyNoLegendPlot, base_height=hXyNoLegendPlot)
save_plot('{output[3]}', pYx, base_width=wYxPlot, base_height=hYxPlot)

save_plot('{output[4]}', pYxNoLegend, base_width=wYxNoLegendPlot, base_height=hYxNoLegendPlot)


" > $(dirname {output[0]})/$(basename {output[0]} .legendOnly.png).r

cat $(dirname {output[0]})/$(basename {output[0]} .legendOnly.png).r | R --slave

		'''



rule aggGffCompareStats:
	input: lambda wildcards: expand("output/statsFiles/" + "tmp/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.splicing_status-{splicedStatus}.endSupport-{endSupport}.vs.gencode.stats.tsv",filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, sampleRep=SAMPLEREPS, endSupport=wildcards.endSupport, minReadSupport=wildcards.minReadSupport, splicedStatus=wildcards.splicedStatus)
	output: "output/statsFiles/" + "all.tmerge.min{minReadSupport}reads.splicing_status-{splicedStatus}.endSupport-{endSupport}.vs.gencode.stats.tsv"
	shell:
		'''
uuidTmpOut=$(uuidgen)
echo -e "seqTech\tcapDesign\tsizeFrac\tsampleRep\tcategory\tcount" > {TMPDIR}/$uuidTmpOut
cat {input} >> {TMPDIR}/$uuidTmpOut
mv {TMPDIR}/$uuidTmpOut {output}
		'''

rule plotGffCompareStats:
	input: "output/statsFiles/" + "all.tmerge.min{minReadSupport}reads.splicing_status-{splicedStatus}.endSupport-{endSupport}.vs.gencode.stats.tsv"
	output: returnPlotFilenames("output/plots/" + "tmerge.vs.gencode.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.splicing_status-{splicedStatus}.endSupport-{endSupport}.vs.gencode.stats")
	conda: "envs/R_env.yml"
	params:
		filterDat=lambda wildcards: multi_figures(wildcards.capDesign, wildcards.sizeFrac, wildcards.sampleRep, wildcards.techname)
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
{params.filterDat[sampleRepFilterString]}
{params.filterDat[substSeqTechString]}
{params.filterDat[substSampleRepString]}
{params.filterDat[graphDimensions]}


dat\$category<-factor(dat\$category, ordered=TRUE, levels=rev(c('Intergenic', 'Extends', 'Intronic', 'Overlaps', 'Antisense', 'Equal', 'Included')))
palette <- c('Intergenic' = '#0099cc', 'Extends' ='#00bfff', 'Intronic' = '#4dd2ff', 'Overlaps' = '#80dfff', 'Antisense' = '#ccf2ff', 'Equal' = '#c65353', 'Included' ='#d98c8c')

plotBase <- \\"p <- ggplot(dat[order(dat\$category), ], aes(x=1, y=count, fill=category)) +
geom_bar(stat='identity') +
scale_fill_manual(values=palette) +
ylab('# TMs') +
xlab('') +
guides(fill = guide_legend(title='Category'))+
scale_y_continuous(labels=scientific)+
{params.filterDat[hideXaxisLabels]}
{GGPLOT_PUB_QUALITY} + 
theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +

\\"

{params.filterDat[facetPlotSetup]}

save_plot('{output[0]}', legendOnly, base_width=wLegendOnly, base_height=hLegendOnly)
save_plot('{output[1]}', pXy, base_width=wXyPlot, base_height=hXyPlot)

save_plot('{output[2]}', pXyNoLegend, base_width=wXyNoLegendPlot, base_height=hXyNoLegendPlot)
save_plot('{output[3]}', pYx, base_width=wYxPlot, base_height=hYxPlot)

save_plot('{output[4]}', pYxNoLegend, base_width=wYxNoLegendPlot, base_height=hYxNoLegendPlot)


" > $(dirname {output[0]})/$(basename {output[0]} .legendOnly.png).r

cat $(dirname {output[0]})/$(basename {output[0]} .legendOnly.png).r | R --slave



		'''


rule getTmVsGencodeLengthStats:
	input: "output/mappings/mergedReads/gffcompare/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.splicing_status-{splicedStatus}.endSupport-{endSupport}.vs.gencode.simple.tsv"
	output: "output/statsFiles/" + "tmp/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.splicing_status-{splicedStatus}.endSupport-{endSupport}.vs.gencode.length.tsv"
	shell:
		'''
uuid=$(uuidgen)
file="$(dirname {input})/$(basename {input} .simple.tsv).tmap"
tail -n+2 $file | awk '$3=="c" || $3=="="' | cut -f10,12 > {TMPDIR}/$uuid
cat {TMPDIR}/$uuid | awk -v s={wildcards.techname} -v c={wildcards.capDesign} -v si={wildcards.sizeFrac} -v b={wildcards.sampleRep} -v sp={wildcards.splicedStatus} '{{print s"\t"c"\t"si"\t"b"\t"sp"\t"$0}}'  > {TMPDIR}/$uuid.TmpOut
mv {TMPDIR}/$uuid.TmpOut {output}

		'''


rule aggTmVsGencodeLengthStats:
	input: lambda wildcards:expand("output/statsFiles/" + "tmp/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.splicing_status-{splicedStatus}.endSupport-{endSupport}.vs.gencode.length.tsv", filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, sampleRep=SAMPLEREPS, endSupport=wildcards.endSupport, minReadSupport=wildcards.minReadSupport, splicedStatus=TMSPLICEDSTATUScategories)
	output: "output/statsFiles/" + "all.tmerge.min{minReadSupport}reads.endSupport-{endSupport}.vs.gencode.length.stats.tsv"
	shell:
		'''
uuidTmpOut=$(uuidgen)
echo -e "seqTech\tcapDesign\tsizeFrac\tsampleRep\tsplicingStatus\tlen\tref_match_len" > {TMPDIR}/$uuidTmpOut
cat {input} >> {TMPDIR}/$uuidTmpOut
mv {TMPDIR}/$uuidTmpOut {output}

		'''

rule plotTmVsGencodeLengthStats:
	input: "output/statsFiles/" + "all.tmerge.min{minReadSupport}reads.endSupport-{endSupport}.vs.gencode.length.stats.tsv"
	output: 
		bySS=returnPlotFilenames("output/plots/" + "tmerge.vs.gencode.length.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.splicing_status-bySplicingStatus.endSupport-{endSupport}.vs.gencode.length.stats"),
		all=returnPlotFilenames("output/plots/" + "tmerge.vs.gencode.length.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.splicing_status-all.endSupport-{endSupport}.vs.gencode.length.stats"),
	conda: "envs/R_env.yml"

	params:
		filterDat=lambda wildcards: multi_figures(wildcards.capDesign, wildcards.sizeFrac, wildcards.sampleRep, wildcards.techname)
	shell:
		'''
echo "
library(ggplot2)

library(cowplot)
library(dplyr)
library(scales)
library(gridExtra)
library(grid)
library(ggplotify)
library(ggExtra)
dat <- read.table('{input}', header=T, as.is=T, sep='\\t')
palette <- c('unspliced' = '#cc3300', 'spliced' = '#0099cc', 'all' = '#666666')
dat\$len <- as.numeric(dat\$len)
dat\$ref_match_len <- as.numeric(dat\$ref_match_len)

{params.filterDat[technameFilterString]}
{params.filterDat[capDesignFilterString]}

{params.filterDat[sizeFracFilterString]}
{params.filterDat[sampleRepFilterString]}
{params.filterDat[substSeqTechString]}
{params.filterDat[substSampleRepString]}
{params.filterDat[graphDimensions]}

datAll <- subset(dat, splicingStatus=='all')

wXyPlot = horizCats * 2.5
hXyPlot = vertCats * 2.5

# geom_textSize=geom_textSize-1
dat %>%
  group_by(splicingStatus) %>%
  summarise(n=n()) -> datSumm

summaryStats = transform(datSumm, Label = paste0('N= ', comma(n)) )

plotBase <- \\"p <- ggplot(dat, aes(x=ref_match_len, y=len, color=splicingStatus)) + 
geom_abline(intercept=0, alpha=0.09, size=lineSize) +
geom_point(alpha=0.1,size=0.5, stroke = 0) + 
#geom_density_2d(size=0.1, alpha=0.3) +
scale_y_log10(limits=c(100,10000)) +  
scale_x_log10(limits=c(100, 20000)) + 
geom_smooth() + 
annotate(x=100, y=10000,  label=paste('Pearson: ', round(cor(datAll\$len, datAll\$ref_match_len, method='pearson'),2)), geom='text', size=geom_textSize, color='#666666', hjust=0, vjust=1) +
annotate(x=100, y=7000,  label=paste('Spearman: ', round(cor(datAll\$len, datAll\$ref_match_len, method='spearman'),2)), geom='text', size=geom_textSize, color='#666666', hjust=0, vjust=1) +

geom_text(data = summaryStats[summaryStats\$splicingStatus=='all',], aes(label = Label, x = 100, y = 2000), hjust=0, vjust=-1,  size=geom_textSize, show.legend=FALSE) +
geom_text(data = summaryStats[summaryStats\$splicingStatus=='unspliced',], aes(label = Label, x = 100, y = 2000), hjust=0, vjust=0,  size=geom_textSize, show.legend=FALSE) +
geom_text(data = summaryStats[summaryStats\$splicingStatus=='spliced',], aes(label = Label, x = 100, y = 2000), hjust=0, vjust=1,  size=geom_textSize, show.legend=FALSE) +


guides(color = guide_legend(title='TM splicing\\nstatus'))+
xlab('Annotated length\\n(mature RNA, nts)') +
ylab('TM length\\n(mature RNA, nts)') +
scale_color_manual(values=palette) +
{GGPLOT_PUB_QUALITY} + theme(legend.position='left') \\"

plotFull <- parse(text =plotBase)
pXy <- eval(plotFull)
legend <- get_legend(pXy)
pXyNoLegend <- pXy + theme(legend.position='none')
pXyMar <- ggMarginal(pXy, groupColour = TRUE, groupFill = TRUE, xparams = list(size=0.1), yparams = list(size=0.1))
pXyMarNoLegend <- ggMarginal(pXyNoLegend, groupColour = TRUE, groupFill = TRUE, xparams = list(size=0.1), yparams = list(size=0.1))

legendOnly <- grid.arrange(legend)
pXyGrob <- as.grob(pXyMar)
pXyNoLegendGrob <- as.grob(pXyMarNoLegend)


hLegendOnly <- convertUnit(sum(legend\$heights), 'in', valueOnly=TRUE)
wLegendOnly <- convertUnit(sum(legend\$widths), 'in', valueOnly=TRUE)

hXyPlot <- hXyPlot
wXyPlot <- wXyPlot +2


hXyNoLegendPlot<- hXyPlot 
wXyNoLegendPlot<- wXyPlot - wLegendOnly




save_plot('{output.bySS[0]}', legendOnly, base_width=wLegendOnly, base_height=hLegendOnly)
save_plot('{output.bySS[1]}', pXy, base_width=wXyPlot, base_height=hXyPlot)

save_plot('{output.bySS[2]}', pXyMar, base_width=wXyPlot, base_height=hXyPlot)
save_plot('{output.bySS[3]}', pXyMar, base_width=wXyPlot, base_height=hXyPlot)

save_plot('{output.bySS[4]}', pXyMarNoLegend, base_width=wXyNoLegendPlot, base_height=hXyNoLegendPlot)

plotBase <- \\"p <- ggplot(datAll, aes(x=ref_match_len, y=len, color=splicingStatus)) + 
geom_abline(intercept=0, alpha=0.09, size=lineSize) +
geom_point(alpha=0.1,size=0.5, stroke = 0) + 
#geom_density_2d(size=0.1, alpha=0.3) +
scale_y_log10(limits=c(100,10000)) +  
scale_x_log10(limits=c(100, 20000)) + 
geom_smooth() + 
annotate(x=100, y=10000,  label=paste('Pearson: ', round(cor(datAll\$len, datAll\$ref_match_len, method='pearson'),2)), geom='text', size=geom_textSize, color='#666666', hjust=0, vjust=1) +
annotate(x=100, y=7000,  label=paste('Spearman: ', round(cor(datAll\$len, datAll\$ref_match_len, method='spearman'),2)), geom='text', size=geom_textSize, color='#666666', hjust=0, vjust=1) +
geom_text(data = summaryStats[summaryStats\$splicingStatus=='all',], aes(label = Label, x = 100, y = 2000), hjust=0, vjust=-1,  size=geom_textSize, show.legend=FALSE) +

guides(color = guide_legend(title='TM splicing\\nstatus'))+
xlab('Annotated length\\n(mature RNA, nts)') +
ylab('TM length\\n(mature RNA, nts)') +
scale_color_manual(values=palette) +
{GGPLOT_PUB_QUALITY} + theme(legend.position='left') \\"

plotFull <- parse(text =plotBase)
pXy <- eval(plotFull)
legend <- get_legend(pXy)
pXyNoLegend <- pXy + theme(legend.position='none')
pXyMar <- ggMarginal(pXy, groupColour = TRUE, groupFill = TRUE, xparams = list(size=0.1), yparams = list(size=0.1))
pXyMarNoLegend <- ggMarginal(pXyNoLegend, groupColour = TRUE, groupFill = TRUE, xparams = list(size=0.1), yparams = list(size=0.1))

legendOnly <- grid.arrange(legend)
pXyGrob <- as.grob(pXyMar)
pXyNoLegendGrob <- as.grob(pXyMarNoLegend)


hLegendOnly <- convertUnit(sum(legend\$heights), 'in', valueOnly=TRUE)
wLegendOnly <- convertUnit(sum(legend\$widths), 'in', valueOnly=TRUE)



hXyNoLegendPlot<- hXyPlot 
wXyNoLegendPlot<- wXyPlot - wLegendOnly




save_plot('{output.all[0]}', legendOnly, base_width=wLegendOnly, base_height=hLegendOnly)
save_plot('{output.all[1]}', pXy, base_width=wXyPlot, base_height=hXyPlot)

save_plot('{output.all[2]}', pXyMar, base_width=wXyPlot, base_height=hXyPlot)
save_plot('{output.all[3]}', pXyMar, base_width=wXyPlot, base_height=hXyPlot)

save_plot('{output.all[4]}', pXyMarNoLegend, base_width=wXyNoLegendPlot, base_height=hXyNoLegendPlot)

" > $(dirname {output.all[0]})/$(basename {output[0]} .legendOnly.png).r

cat $(dirname {output.all[0]})/$(basename {output[0]} .legendOnly.png).r | R --slave

		'''





rule mergeTmsWithGencode:
	input:
		annot="output/annotations/simplified/{capDesign}.gencode.simplified_biotypes.gtf",
		tm="output/mappings/mergedReads/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.HiSS.tmerge.min{minReadSupport}reads.splicing_status-all.endSupport-{endSupport}.gff.gz"
	output: "output/mappings/mergedReads/gencodeMerge/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.endSupport-{endSupport}.+gencode.gff.gz"
	threads:1
	shell:
		'''
uuidTmpOut=$(uuidgen)
zcat -f {input.annot} {input.tm}  | skipcomments | sort -T {TMPDIR}  -k1,1 -k4,4n -k5,5n  | tmerge --exonOverhangTolerance {ExonOverhangTolerance} - |sort -T {TMPDIR}  -k1,1 -k4,4n -k5,5n  |gzip > {TMPDIR}/$uuidTmpOut
mv {TMPDIR}/$uuidTmpOut {output}
		'''

rule makeClsGencodeLoci:
	input: "output/mappings/mergedReads/gencodeMerge/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.endSupport-{endSupport}.+gencode.gff.gz"
	params: locusPrefix=config["PROJECT_NAME"]
	output: temp("output/mappings/mergedReads/gencodeLociMerge/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.endSupport-{endSupport}.+gencode.loci.gff.gz")
	conda: "envs/xtools_env.yml"
	shell:
		'''
uuid=$(uuidgen)
uuidTmpOut=$(uuidgen)
zcat {input} > {TMPDIR}/$uuid


bedtools intersect -s -wao -a {TMPDIR}/$uuid -b {TMPDIR}/$uuid |awk '$1 !~ /ERCC/'| buildLoci.pl --locPrefix {params.locusPrefix}: - |sort -T {TMPDIR}  -k1,1 -k4,4n -k5,5n  | gzip> {TMPDIR}/$uuidTmpOut
mv {TMPDIR}/$uuidTmpOut {output}

		'''


rule mergeWithRef:
	input:
		clsGencode="output/mappings/mergedReads/gencodeLociMerge/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.endSupport-{endSupport}.+gencode.loci.gff.gz",
		gencode="output/annotations/simplified/{capDesign}.gencode.simplified_biotypes.gtf"
	output: "output/mappings/mergedReads/mergeToRef/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.endSupport-{endSupport}.+gencode.loci.refmerged.gff.gz"
	shell:
		'''
uuid=$(uuidgen)
uuidTmpOut=$(uuidgen)
zcat  {input.clsGencode} > {TMPDIR}/$uuid
mergeToRef.pl {input.gencode} {TMPDIR}/$uuid | sort -T {TMPDIR}  -k1,1 -k4,4n -k5,5n  |gzip > {TMPDIR}/$uuidTmpOut
mv {TMPDIR}/$uuidTmpOut {output}
		'''

rule getNovelIntergenicLoci:
	input:
		gencode="output/annotations/simplified/{capDesign}.gencode.simplified_biotypes.gtf",
		tmergeGencode="output/mappings/mergedReads/mergeToRef/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.endSupport-{endSupport}.+gencode.loci.refmerged.gff.gz"
	output:
		gtf="output/mappings/mergedReads/mergeToRef/novelLoci/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.endSupport-{endSupport}.novelIntergenicLoci.gff.gz",
		locusBed="output/mappings/mergedReads/mergeToRef/novelLoci/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.endSupport-{endSupport}.novelLoci.geneCoords.bed"
	conda: "envs/xtools_env.yml"
	shell:
		'''
uuid1=$(uuidgen)
uuid2=$(uuidgen)
uuid3=$(uuidgen)
uuid4=$(uuidgen)
uuid5=$(uuidgen)
uuidTmpOut=$(uuidgen)
cat {input.gencode} |awk '$3=="exon"' | extract_locus_coords.pl -| sort -T {TMPDIR}  -k1,1 -k2,2n -k3,3n  > {TMPDIR}/$uuid1
zcat {input.tmergeGencode} | tgrep -F 'gene_ref_status "novel";' > {TMPDIR}/$uuid4
cat {TMPDIR}/$uuid4 | extract_locus_coords.pl - | sort -T {TMPDIR}  -k1,1 -k2,2n -k3,3n  > {TMPDIR}/$uuid2



bedtools intersect -v -a {TMPDIR}/$uuid2 -b {TMPDIR}/$uuid1 > {TMPDIR}/$uuid5
cat {TMPDIR}/$uuid5 |awk '$1 !~ /ERCC/' |cut -f4 | sort -T {TMPDIR} |uniq > {TMPDIR}/$uuid3
zcat {input.tmergeGencode}| tgrep -F -w -f {TMPDIR}/$uuid3 - |gzip > {TMPDIR}/$uuidTmpOut
zcat {TMPDIR}/$uuidTmpOut |  extract_locus_coords.pl - | sort -T {TMPDIR}  -k1,1 -k2,2n -k3,3n  > {TMPDIR}/$uuid2.bed
mv {TMPDIR}/$uuidTmpOut {output.gtf}
mv {TMPDIR}/$uuid2.bed {output.locusBed}
		'''

rule getNovelIntergenicLociQcStats:
	input:
		tmergeGencode="output/mappings/mergedReads/mergeToRef/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.endSupport-{endSupport}.+gencode.loci.refmerged.gff.gz",
		repeats= lambda wildcards: "annotations/repeatMasker/" + CAPDESIGNTOGENOME[wildcards.capDesign] + ".repeatMasker.bed"
	conda: "envs/Redtools_env.yml"
	output: "output/statsFiles/" + "tmp/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.endSupport-{endSupport}.novelLoci.qc.stats.tsv"
	shell:
		'''
uuid=$(uuidgen)

# extract transcript models from novel genes
zcat {input.tmergeGencode} | tgrep -F 'gene_ref_status "novel";' > {TMPDIR}/$uuid
# extract gene_ids from novel genes and make a list of those
cat {TMPDIR}/$uuid | extractGffAttributeValue.pl gene_id |sort|uniq > {TMPDIR}/$uuid.genes.list

# extract bed12 file of transcript models (TMs) from novel genes
cat {TMPDIR}/$uuid | gff2bed_full.pl - |sort -T {TMPDIR}  -k1,1 -k2,2n -k3,3n > {TMPDIR}/$uuid.bed

############## spliced vs monoexonic novel genes:

# extract spliced TMs from bed12 file
cat {TMPDIR}/$uuid.bed | awk '$10>1' | cut -f4 | sort -T {TMPDIR} |uniq >  {TMPDIR}/$uuid.spliced.transcripts.list

# extract gene_ids of spliced TMs
tgrep -F -w -f {TMPDIR}/$uuid.spliced.transcripts.list {TMPDIR}/$uuid | extractGffAttributeValue.pl gene_id |sort|uniq > {TMPDIR}/$uuid.spliced.genes.list

total=$(cat {TMPDIR}/$uuid.genes.list |wc -l)
spliced=$(cat {TMPDIR}/$uuid.spliced.genes.list  |wc -l)
let mono=$total-$spliced || true

############## length statistics

cat {TMPDIR}/$uuid.bed | bed12ToTranscriptLength.pl - > {TMPDIR}/$uuid.stats.alltranscripts.tsv.tmp
tgrep -F -w -f {TMPDIR}/$uuid.spliced.transcripts.list {TMPDIR}/$uuid.stats.alltranscripts.tsv.tmp |cut -f2 > {TMPDIR}/$uuid.stats.splicedtranscripts.tsv
tgrep -F -w -v -f {TMPDIR}/$uuid.spliced.transcripts.list {TMPDIR}/$uuid.stats.alltranscripts.tsv.tmp | cut -f2 > {TMPDIR}/$uuid.stats.monotranscripts.tsv
cat {TMPDIR}/$uuid.stats.alltranscripts.tsv.tmp | cut -f2 >{TMPDIR}/$uuid.stats.alltranscripts.tsv



#use of readarray: see https://www.baeldung.com/linux/reading-output-into-array

readarray -t lengthStatsAll < <(minMedMaxStats.r {TMPDIR}/$uuid.stats.alltranscripts.tsv)
readarray -t lengthStatsMono < <(minMedMaxStats.r {TMPDIR}/$uuid.stats.monotranscripts.tsv)
readarray -t lenthStatsSpliced < <(minMedMaxStats.r {TMPDIR}/$uuid.stats.splicedtranscripts.tsv)


########### repeatmasked regions:

bedtools intersect -split -u -a {TMPDIR}/$uuid.bed -b {input.repeats} | cut -f4 | sort|uniq > {TMPDIR}/$uuid.repeats.transcripts.list
tgrep -F -w -f {TMPDIR}/$uuid.repeats.transcripts.list {TMPDIR}/$uuid | extractGffAttributeValue.pl gene_id |sort|uniq > {TMPDIR}/$uuid.repeats.genes.list

repeats=$(cat {TMPDIR}/$uuid.repeats.genes.list| wc -l)

############ write output
echo -e "{wildcards.techname}\\t{wildcards.capDesign}\\t{wildcards.sizeFrac}\\t{wildcards.sampleRep}\\t$total\\t$mono\\t$repeats\\t${{lengthStatsAll[*]}}\\t${{lengthStatsMono[*]}}\\t${{lenthStatsSpliced[*]}}" | ssv2tsv| awk '{{if ($5!=0) print $0"\t"$6/$5"\t"$7/$5; else print $0"\tNA\tNA"}}' > {TMPDIR}/$uuid.1
mv {TMPDIR}/$uuid.1 {output}

		'''

rule aggNovelIntergenicLociQcStats:
	input: lambda wildcards: expand("output/statsFiles/" + "tmp/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.endSupport-{endSupport}.novelLoci.qc.stats.tsv",filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, sampleRep=SAMPLEREPS, endSupport=wildcards.endSupport, minReadSupport=wildcards.minReadSupport)
	output: "output/statsFiles/" + "all.tmerge.min{minReadSupport}reads.endSupport-{endSupport}.novelLoci.qc.stats.tsv"
	shell:
		'''
uuid=$(uuidgen)
echo -e "seqTech\\tcapDesign\\tsizeFrac\\tsampleRep\\ttotal\\tcountMono\\tcountRepeats\\tminLengthAllTms\\tmedianLengthAllTms\\tmaxLengthAllTms\\tminLengthMonoTms\\tmedianLengthMonoTms\\tmaxLengthMonoTms\\tminLengthSplicedTms\\tmedianLengthSplicedTms\\tmaxLengthSplicedTms\\tpercentMono\\tpercentRepeats" > {TMPDIR}/$uuid
cat {input}  | sort -T {TMPDIR}  >> {TMPDIR}/$uuid
mv {TMPDIR}/$uuid {output}

		'''

rule getNovelIntergenicLociStats:
	input:
		tmergeGencode="output/mappings/mergedReads/mergeToRef/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.endSupport-{endSupport}.+gencode.loci.refmerged.gff.gz",
		intergenic="output/mappings/mergedReads/mergeToRef/novelLoci/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.endSupport-{endSupport}.novelIntergenicLoci.gff.gz"
	output: "output/statsFiles/" + "tmp/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.endSupport-{endSupport}.novelLoci.stats.tsv"
	shell:
		'''
uuidTmpOut=$(uuidgen)
totalNovel=$(zcat {input.tmergeGencode} | tgrep -F 'gene_ref_status "novel";' |extractGffAttributeValue.pl gene_id | sort -T {TMPDIR} | uniq | wc -l)
interg=$(zcat {input.intergenic} | extractGffAttributeValue.pl gene_id | sort -T {TMPDIR} | uniq | wc -l)
echo -e "{wildcards.techname}\t{wildcards.capDesign}\t{wildcards.sizeFrac}\t{wildcards.sampleRep}\t$totalNovel\t$interg"  > {TMPDIR}/$uuidTmpOut
mv {TMPDIR}/$uuidTmpOut {output}

		'''

rule aggNovelIntergenicLociStats:
	input: lambda wildcards: expand("output/statsFiles/" + "tmp/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.endSupport-{endSupport}.novelLoci.stats.tsv",filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, sampleRep=SAMPLEREPS, endSupport=wildcards.endSupport, minReadSupport=wildcards.minReadSupport)
	output: "output/statsFiles/" + "all.tmerge.min{minReadSupport}reads.endSupport-{endSupport}.novelLoci.stats.tsv"
	shell:
		'''
uuidTmpOut=$(uuidgen)
echo -e "seqTech\tcapDesign\tsizeFrac\tsampleRep\tcategory\tcount\tpercent" > {TMPDIR}/$uuidTmpOut
cat {input} | awk '{{if ($5!=0) print $1"\\t"$2"\\t"$3"\\t"$4"\\tintergenic\\t"$6"\\t"$6/$5"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\tintronic\\t"$5-$6"\\t"($5-$6)/$5; else print $1"\\t"$2"\\t"$3"\\t"$4"\\tintergenic\\t"$6"\\t0\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\tintronic\\t"$5-$6"\\t0"}}' | sort -T {TMPDIR}  >> {TMPDIR}/$uuidTmpOut
mv {TMPDIR}/$uuidTmpOut {output}
		'''


rule plotNovelIntergenicLociStats:
	input: "output/statsFiles/" + "all.tmerge.min{minReadSupport}reads.endSupport-{endSupport}.novelLoci.stats.tsv"
	output: returnPlotFilenames("output/plots/" + "tmerge.novelLoci.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.endSupport-{endSupport}.novelLoci.stats")
	conda: "envs/R_env.yml"
	params:
		filterDat=lambda wildcards: multi_figures(wildcards.capDesign, wildcards.sizeFrac, wildcards.sampleRep, wildcards.techname)
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
{params.filterDat[sampleRepFilterString]}
{params.filterDat[substSeqTechString]}
{params.filterDat[substSampleRepString]}
{params.filterDat[graphDimensions]}

dat\$category<-factor(dat\$category, ordered=TRUE, levels=rev(c('intronic', 'intergenic')))
plotBase <- \\"p <- ggplot(dat[order(dat\$category), ], aes(x=1, y=count, fill=category)) +
geom_bar(stat='identity') +
ylab('# Novel loci') +
scale_y_continuous(labels=comma)+
scale_fill_manual (values=c(intronic='#8cd9b3', intergenic='#00b386'))+
xlab('') +
guides(fill = guide_legend(title='Category\\n(w.r.t. GENCODE)'))+
geom_text(position = 'stack', size=geom_textSize, aes( y = count, label = comma(count), hjust = 0.5, vjust = 1))+
{params.filterDat[hideXaxisLabels]}
{GGPLOT_PUB_QUALITY}  + 
theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
\\"

{params.filterDat[facetPlotSetup]}

save_plot('{output[0]}', legendOnly, base_width=wLegendOnly, base_height=hLegendOnly)
save_plot('{output[1]}', pXy, base_width=wXyPlot, base_height=hXyPlot)

save_plot('{output[2]}', pXyNoLegend, base_width=wXyNoLegendPlot, base_height=hXyNoLegendPlot)
save_plot('{output[3]}', pYx, base_width=wYxPlot, base_height=hYxPlot)

save_plot('{output[4]}', pYxNoLegend, base_width=wYxNoLegendPlot, base_height=hYxNoLegendPlot)


" > $(dirname {output[0]})/$(basename {output[0]} .legendOnly.png).r

cat $(dirname {output[0]})/$(basename {output[0]} .legendOnly.png).r | R --slave



		'''



rule tmergeAll:
	input:
		tm=lambda wildcards: expand("output/mappings/mergedReads/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.HiSS.tmerge.min{minReadSupport}reads.splicing_status-{splicedStatus}.endSupport-{endSupport}.gff.gz", filtered_product, techname=TECHNAMES, capDesign=wildcards.capDesign, sizeFrac=SIZEFRACS, sampleRep=SAMPLEREPS, endSupport=wildcards.endSupport,  minReadSupport=wildcards.minReadSupport, splicedStatus=wildcards.splicedStatus)
		#gencode="annotations/simplified/{capDesign}.gencode.simplified_biotypes.gtf",
	output: 
		tm="output/mappings/mergedReads/tmergeAll/{capDesign}_min{minReadSupport}reads.splicing_status-{splicedStatus}.endSupport-{endSupport}.tmerge.gff.gz",
		quant="output/mappings/mergedReads/tmergeAll/{capDesign}_min{minReadSupport}reads.splicing_status-{splicedStatus}.endSupport-{endSupport}.tmerge.expQuant.tsv"
	shell:
		'''
uuid=$(uuidgen)
uuidTmpOutT=$(uuidgen)
uuidTmpOutQ=$(uuidgen)
for file in `ls {input.tm} | grep -v AlzhBrain`; do
bn=$(basename $file .HiSS.tmerge.min{wildcards.minReadSupport}reads.splicing_status-{wildcards.splicedStatus}.endSupport-{wildcards.endSupport}.gff.gz )

zcat $file | perl -sne '$_=~s/transcript_id \"(\S+)\"/transcript_id \"=$var=$1\"/g; print' -- -var=$bn
done > {TMPDIR}/$uuid
echo -e "transcript_id\tspliced\tflrpm\trpm" > {TMPDIR}/$uuidTmpOutQ
cat {TMPDIR}/$uuid | extractGffAttributeValue.pl transcript_id spliced flrpm rpm | sort|uniq >> {TMPDIR}/$uuidTmpOutQ
countDups=$(cat {TMPDIR}/$uuidTmpOutQ |cut -f1 |sort|uniq -d |wc -l)
if [ $countDups -gt 0 ]; then echoerr "$countDups duplicates found"; exit 1; fi;
cat {TMPDIR}/$uuid  | skipcomments | sort -T {TMPDIR} -k1,1 -k4,4n -k5,5n | tmerge --exonOverhangTolerance {ExonOverhangTolerance} - |sort -T {TMPDIR}  -k1,1 -k4,4n -k5,5n  > {TMPDIR}/$uuidTmpOutT
cat {TMPDIR}/$uuidTmpOutT | gzip > {output.tm}
mv {TMPDIR}/$uuidTmpOutQ {output.quant}

		'''




rule getAnnotatedGeneDetectionStats:
	input: "output/mappings/mergedReads/gffcompare/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.splicing_status-{splicedStatus}.endSupport-{endSupport}.vs.gencode.simple.tsv"
	output: "output/statsFiles/" + "tmp/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.splicing_status-{splicedStatus}.endSupport-{endSupport}.gencode.geneDetection.stats.tsv"
	shell:
		'''
uuid=$(uuidgen)
file="$(dirname {input})/$(basename {input} .simple.tsv).tmap"
totalGenes=$(tail -n+2 $file |cut -f1 | awk '$1!="-"'|sort|uniq|wc -l)
flGenes=$(tail -n+2 $file | awk '$3=="="' |cut -f1 |sort|uniq|wc -l)
let nonFlGenes=$totalGenes-$flGenes || true
echo -e "{wildcards.techname}\\t{wildcards.capDesign}\\t{wildcards.sizeFrac}\\t{wildcards.sampleRep}\\t$totalGenes\\t$flGenes\\t$nonFlGenes" | awk '{{print $0"\\t"$6/$5"\\t"$7/$5}}'> {TMPDIR}/$uuid.TmpOut
mv {TMPDIR}/$uuid.TmpOut {output}

		'''
rule aggAnnotatedGeneDetectionStats:
	input: lambda wildcards: expand("output/statsFiles/" + "tmp/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.splicing_status-{splicedStatus}.endSupport-{endSupport}.gencode.geneDetection.stats.tsv",filtered_product,  techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, sampleRep=SAMPLEREPS, endSupport=wildcards.endSupport, minReadSupport=wildcards.minReadSupport, splicedStatus=wildcards.splicedStatus)
	output: "output/statsFiles/" + "all.min{minReadSupport}reads.splicing_status-{splicedStatus}.endSupport-{endSupport}.gencode.geneDetection.stats.tsv"
	shell:
		'''
uuid=$(uuidgen)
echo -e "seqTech\tcapDesign\tsizeFrac\tsampleRep\tcategory\tcount\tpercent" > {TMPDIR}/$uuid
cat {input} | awk '{{print $1"\\t"$2"\\t"$3"\\t"$4"\\tflGenes\\t"$6"\\t"$8"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\tnonFlGenes\\t"$7"\\t"$9}}' | sort -T {TMPDIR}  >> {TMPDIR}/$uuid
mv {TMPDIR}/$uuid {output}

		'''

rule plotAnnotatedGeneDetectionStats:
	input: "output/statsFiles/" + "all.min{minReadSupport}reads.splicing_status-{splicedStatus}.endSupport-{endSupport}.gencode.geneDetection.stats.tsv"
	output: returnPlotFilenames("output/plots/" + "gencode.geneDetection.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.splicing_status-{splicedStatus}.endSupport-{endSupport}.gencode.geneDetection.stats")
	conda: "envs/R_env.yml"
	params:
		filterDat=lambda wildcards: multi_figures(wildcards.capDesign, wildcards.sizeFrac, wildcards.sampleRep, wildcards.techname)
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
dat\$category<-factor(dat\$category, ordered=TRUE, levels=rev(c('flGenes', 'nonFlGenes')))
{params.filterDat[technameFilterString]}
{params.filterDat[capDesignFilterString]}

{params.filterDat[sizeFracFilterString]}
{params.filterDat[sampleRepFilterString]}
{params.filterDat[substSeqTechString]}
{params.filterDat[substSampleRepString]}
{params.filterDat[graphDimensions]}

plotBase <- \\"p <- ggplot(dat[order(dat\$category), ], aes(x=1, y=count, fill=category)) +
geom_bar(stat='identity') +
scale_fill_manual(values=c('flGenes' = '#4da6ff', 'nonFlGenes' = '#cc9046ff'), labels = c(flGenes = 'Full-length', nonFlGenes = 'Partial')) +
ylab('# GENCODE genes detected') +
xlab('') +
guides(fill = guide_legend(title='Category'))+
scale_y_continuous(labels=scientific)+
{GGPLOT_PUB_QUALITY} + 
theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
\\"

{params.filterDat[facetPlotSetup]}

save_plot('{output[0]}', legendOnly, base_width=wLegendOnly, base_height=hLegendOnly)
save_plot('{output[1]}', pXy, base_width=wXyPlot, base_height=hXyPlot)

save_plot('{output[2]}', pXyNoLegend, base_width=wXyNoLegendPlot, base_height=hXyNoLegendPlot)
save_plot('{output[3]}', pYx, base_width=wYxPlot, base_height=hYxPlot)

save_plot('{output[4]}', pYxNoLegend, base_width=wYxNoLegendPlot, base_height=hYxNoLegendPlot)


" > $(dirname {output[0]})/$(basename {output[0]} .legendOnly.png).r

cat $(dirname {output[0]})/$(basename {output[0]} .legendOnly.png).r | R --slave

		'''


rule getNtCoverageByGenomePartition:
	input:
		tm="output/mappings/mergedReads/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.HiSS.tmerge.min{minReadSupport}reads.splicing_status-{splicedStatus}.endSupport-{endSupport}.gff.gz",
		gencodePart="output/annotations/{capDesign}.partition.gff"
	output: "output/statsFiles/" + "tmp/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.HiSS.tmerge.min{minReadSupport}reads.splicing_status-{splicedStatus}.endSupport-{endSupport}.ntCoverageByGenomePartition.stats.tsv"
	conda: "envs/xtools_env.yml"
	shell:
		'''


uuid=$(uuidgen)
zcat {input.tm} | bedtools merge -i stdin > {TMPDIR}/$uuid.tmp.tm.bed
bedtools coverage -a {input.gencodePart} -b {TMPDIR}/$uuid.tmp.tm.bed  > {TMPDIR}/$uuid.tmp.tm.cov.bedtsv

totalNts=$(cat {TMPDIR}/$uuid.tmp.tm.cov.bedtsv | awk '{{print $(NF-2)}}' | sum.sh)

for flag in `cat {TMPDIR}/$uuid.tmp.tm.cov.bedtsv | extractGffAttributeValue.pl region_flag|sort|uniq`; do 
nts=$(cat {TMPDIR}/$uuid.tmp.tm.cov.bedtsv |fgrep "region_flag \\"$flag\\";" | awk '{{print $(NF-2)}}'|sum.sh)
echo -e "{wildcards.techname}\\t{wildcards.capDesign}\\t{wildcards.sizeFrac}\\t{wildcards.sampleRep}\\t{wildcards.splicedStatus}\\t$flag\\t$nts" | awk -v t=$totalNts '{{print $0"\\t"$7/t}}'; done  > {TMPDIR}/$uuid.tmp.tm.cov.stats.tsv

mv {TMPDIR}/$uuid.tmp.tm.cov.stats.tsv {output}
		'''

rule aggNtCoverageByGenomePartition:
	input: lambda wildcards:expand("output/statsFiles/" + "tmp/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.HiSS.tmerge.min{minReadSupport}reads.splicing_status-{splicedStatus}.endSupport-{endSupport}.ntCoverageByGenomePartition.stats.tsv", filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, sampleRep=SAMPLEREPS, endSupport=wildcards.endSupport, minReadSupport=wildcards.minReadSupport, splicedStatus=TMSPLICEDSTATUScategories)
	output: "output/statsFiles/" + "all.tmerge.min{minReadSupport}reads.endSupport-{endSupport}.vs.ntCoverageByGenomePartition.stats.tsv"
	shell:
		'''
uuid=$(uuidgen)
echo -e "seqTech\\tcapDesign\\tsizeFrac\\tsampleRep\\tsplicingStatus\\tpartition\\tnt_coverage_count\\tnt_coverage_percent" > {TMPDIR}/$uuid
cat {input} | sort >> {TMPDIR}/$uuid
mv {TMPDIR}/$uuid {output}

		'''

rule plotNtCoverageByGenomePartition:
	input: "output/statsFiles/" + "all.tmerge.min{minReadSupport}reads.endSupport-{endSupport}.vs.ntCoverageByGenomePartition.stats.tsv"
	output: returnPlotFilenames("output/plots/" + "tmerge.ntCoverageByGenomePartition.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.endSupport-{endSupport}.ntCoverageByGenomePartition.stats")
	conda: "envs/R_env.yml"
	params:
		filterDat=lambda wildcards: multi_figures(wildcards.capDesign, wildcards.sizeFrac, wildcards.sampleRep, wildcards.techname)
	shell:
		'''
echo "
library(ggplot2)

library(cowplot)
library(dplyr)
library(scales)
library(gridExtra)
library(grid)
library(ggplotify)

dat <- read.table('{input}', header=T, as.is=T, sep='\\t')
{params.filterDat[technameFilterString]}
{params.filterDat[capDesignFilterString]}

{params.filterDat[sizeFracFilterString]}
{params.filterDat[sampleRepFilterString]}
{params.filterDat[substSeqTechString]}
{params.filterDat[substSampleRepString]}
{params.filterDat[graphDimensions]}

dat <- filter(dat, splicingStatus=='all')
dat <- select(dat, -splicingStatus)


dat\$partition<-factor(dat\$partition, ordered=TRUE)
plotBase <- \\"p <- ggplot(dat[order(dat\$partition), ], aes(x=1, y=nt_coverage_count, fill=partition)) +
geom_bar(stat='identity') +
ylab('# genomic nts covered') +
xlab('') +
#scale_y_continuous(labels=comma)+
scale_fill_manual (values=c(intron='#d98c8c', intergenic='#33ccff', CDS='#00e64d', exonOfNCT='#6666ff', exonOfPseudo='#e066ff', UTR='#999966'))+
guides(fill = guide_legend(title='Category\\n(w.r.t. GENCODE)'))+
#geom_text(position = 'stack', size=geom_textSize, aes( y = count, label = paste(sep='',percent(round(percent, digits=2)),' / ','(',comma(count),')'), hjust = 0.5, vjust = 1))+
{params.filterDat[hideXaxisLabels]}
{GGPLOT_PUB_QUALITY}  + 
theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
\\"

{params.filterDat[facetPlotSetup]}

save_plot('{output[0]}', legendOnly, base_width=wLegendOnly, base_height=hLegendOnly)
save_plot('{output[1]}', pXy, base_width=wXyPlot, base_height=hXyPlot)

save_plot('{output[2]}', pXyNoLegend, base_width=wXyNoLegendPlot, base_height=hXyNoLegendPlot)
save_plot('{output[3]}', pYx, base_width=wYxPlot, base_height=hYxPlot)

save_plot('{output[4]}', pYxNoLegend, base_width=wYxNoLegendPlot, base_height=hYxNoLegendPlot)


" > $(dirname {output[0]})/$(basename {output[0]} .legendOnly.png).r

cat $(dirname {output[0]})/$(basename {output[0]} .legendOnly.png).r | R --slave


		'''
