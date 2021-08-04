


rule getGencodeSupportedEnds:
	input:
		tm="annotations/simplified/{capDesign}.gencode.collapsed.simplified_biotypes.gtf",
		cagePeaks=lambda wildcards: GENOMETOCAGEPEAKS[CAPDESIGNTOGENOME[wildcards.capDesign]],
		genome = lambda wildcards: config["GENOMESDIR"] + CAPDESIGNTOGENOME[wildcards.capDesign] + ".sorted.genome",
		PAS=lambda wildcards: GENOMETOPAS[CAPDESIGNTOGENOME[wildcards.capDesign]],

	output:"mappings/nonAnchoredMergeReads/mergeWithPrevious/gencode/{capDesign}.gencode.cage+PASsupported.gff.gz"
	conda: "envs/xtools_env.yml"
	shell:
		'''
uuid5pEnds=$(uuidgen)
uuidTmpOut=$(uuidgen)
uuid2=$(uuidgen)
cat {input.genome} | cut -f1 | sort -T {config[TMPDIR]}  |uniq > {config[TMPDIR]}/$uuid2

cat {input.tm} | fgrep -w -f {config[TMPDIR]}/$uuid2 | gff2bed_full.pl - |extractTranscriptEndsFromBed12.pl 5 |sort -T {config[TMPDIR]}  -k1,1 -k2,2n -k3,3n > {config[TMPDIR]}/$uuid5pEnds
uuid3pEnds=$(uuidgen)
cat {input.tm} | fgrep -w -f {config[TMPDIR]}/$uuid2|gff2bed_full.pl - |extractTranscriptEndsFromBed12.pl 3 |sort -T {config[TMPDIR]}  -k1,1 -k2,2n -k3,3n > {config[TMPDIR]}/$uuid3pEnds

uuidCageSupported=$(uuidgen)
cat {config[TMPDIR]}/$uuid5pEnds | sort -T {config[TMPDIR]}  -k1,1 -k2,2n -k3,3n  | bedtools slop -s -l 50 -r 50 -i stdin -g {input.genome} | bedtools intersect -u -s -a stdin -b {input.cagePeaks} | cut -f4  |sort -T {config[TMPDIR]} |uniq > {config[TMPDIR]}/$uuidCageSupported

uuidPASsupported=$(uuidgen)
cat {config[TMPDIR]}/$uuid3pEnds | sort -T {config[TMPDIR]}  -k1,1 -k2,2n -k3,3n  | bedtools slop -s -l 50 -r -10 -i stdin -g {input.genome} | bedtools intersect -u -s -a stdin -b {input.PAS} |cut -f4 |sort -T {config[TMPDIR]} |uniq > {config[TMPDIR]}/$uuidPASsupported

uuidcagePASsupported=$(uuidgen)
comm -1 -2 {config[TMPDIR]}/$uuidCageSupported {config[TMPDIR]}/$uuidPASsupported |sort -T {config[TMPDIR]} |uniq | awk '{{print "transcript_id \\""$1"\\";"}}'> {config[TMPDIR]}/$uuidcagePASsupported

fgrep -w -f {config[TMPDIR]}/$uuidcagePASsupported {input.tm} |sort -T {config[TMPDIR]}  -k1,1 -k4,4n -k5,5n  | gzip> {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}

		'''

rule getCurrentPreviousPhaseTmsWithGencodeSupportedEnds:
	input:
		tm="mappings/nonAnchoredMergeReads/mergeWithPrevious+biotypes/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.tmerge.min{minReadSupport}reads.endSupport:all.gff.gz",
		cagePeaks=lambda wildcards: GENOMETOCAGEPEAKS[CAPDESIGNTOGENOME[wildcards.capDesign]],
		genome = lambda wildcards: config["GENOMESDIR"] + CAPDESIGNTOGENOME[wildcards.capDesign] + ".sorted.genome",
		PAS=lambda wildcards: GENOMETOPAS[CAPDESIGNTOGENOME[wildcards.capDesign]],

	output:"mappings/nonAnchoredMergeReads/mergeWithPrevious/cls/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.tmerge.min{minReadSupport}reads.cage+PASsupported.gff.gz"
	conda: "envs/xtools_env.yml"
	shell:
		'''
uuid5pEnds=$(uuidgen)
uuidChrList=$(uuidgen)
uuidTmpOut=$(uuidgen)

cut -f1 {input.genome} |sort|uniq > {config[TMPDIR]}/$uuidChrList
zcat {input.tm} | fgrep -w -f {config[TMPDIR]}/$uuidChrList |gff2bed_full.pl - |extractTranscriptEndsFromBed12.pl 5 |sort -T {config[TMPDIR]}  -k1,1 -k2,2n -k3,3n > {config[TMPDIR]}/$uuid5pEnds
uuid3pEnds=$(uuidgen)
zcat {input.tm} | fgrep -w -f {config[TMPDIR]}/$uuidChrList |gff2bed_full.pl - |extractTranscriptEndsFromBed12.pl 3 |sort -T {config[TMPDIR]}  -k1,1 -k2,2n -k3,3n > {config[TMPDIR]}/$uuid3pEnds

uuidCageSupported=$(uuidgen)
cat {config[TMPDIR]}/$uuid5pEnds | sort -T {config[TMPDIR]}  -k1,1 -k2,2n -k3,3n  | bedtools slop -s -l 50 -r 50 -i stdin -g {input.genome} | bedtools intersect -u -s -a stdin -b {input.cagePeaks} | cut -f4  |sort -T {config[TMPDIR]} |uniq > {config[TMPDIR]}/$uuidCageSupported

uuidPASsupported=$(uuidgen)
cat {config[TMPDIR]}/$uuid3pEnds | sort -T {config[TMPDIR]}  -k1,1 -k2,2n -k3,3n  | bedtools slop -s -l 50 -r -10 -i stdin -g {input.genome} | bedtools intersect -u -s -a stdin -b {input.PAS} |cut -f4 |sort -T {config[TMPDIR]} |uniq > {config[TMPDIR]}/$uuidPASsupported

uuidcagePASsupported=$(uuidgen)
comm -1 -2 {config[TMPDIR]}/$uuidCageSupported {config[TMPDIR]}/$uuidPASsupported |sort -T {config[TMPDIR]} |uniq | awk '{{print "transcript_id \\""$1"\\";"}}'> {config[TMPDIR]}/$uuidcagePASsupported

zcat {input.tm} | fgrep -w -f {config[TMPDIR]}/$uuidcagePASsupported - |sort -T {config[TMPDIR]}  -k1,1 -k4,4n -k5,5n  | gzip> {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}

		'''


rule getPreviousPhaseTmsWithGencodeSupportedEnds:
	input:
		tm="mappings/nonAnchoredMergeReads/previous+biotypes/{capDesign}.previous.tmerge.all.gff.gz",
		cagePeaks=lambda wildcards: GENOMETOCAGEPEAKS[CAPDESIGNTOGENOME[wildcards.capDesign]],
		genome = lambda wildcards: config["GENOMESDIR"] + CAPDESIGNTOGENOME[wildcards.capDesign] + ".sorted.genome",
		PAS=lambda wildcards: GENOMETOPAS[CAPDESIGNTOGENOME[wildcards.capDesign]],

	output:"mappings/nonAnchoredMergeReads/previous/{capDesign}.previous.tmerge.cage+PASsupported.gff.gz"
	conda: "envs/xtools_env.yml"
	shell:
		'''
uuid5pEnds=$(uuidgen)
uuidTmpOut=$(uuidgen)
uuid2=$(uuidgen)
cat {input.genome} | cut -f1 | sort -T {config[TMPDIR]}  |uniq > {config[TMPDIR]}/$uuid2

zcat {input.tm} | fgrep -w -f {config[TMPDIR]}/$uuid2|gff2bed_full.pl - |extractTranscriptEndsFromBed12.pl 5 |sort -T {config[TMPDIR]}  -k1,1 -k2,2n -k3,3n > {config[TMPDIR]}/$uuid5pEnds
uuid3pEnds=$(uuidgen)
zcat {input.tm} | fgrep -w -f {config[TMPDIR]}/$uuid2|gff2bed_full.pl - |extractTranscriptEndsFromBed12.pl 3 |sort -T {config[TMPDIR]}  -k1,1 -k2,2n -k3,3n > {config[TMPDIR]}/$uuid3pEnds


uuidCageSupported=$(uuidgen)
cat {config[TMPDIR]}/$uuid5pEnds | sort -T {config[TMPDIR]}  -k1,1 -k2,2n -k3,3n  | bedtools slop -s -l 50 -r 50 -i stdin -g {input.genome} | bedtools intersect -u -s -a stdin -b {input.cagePeaks} | cut -f4  |sort -T {config[TMPDIR]} |uniq > {config[TMPDIR]}/$uuidCageSupported

uuidPASsupported=$(uuidgen)
cat {config[TMPDIR]}/$uuid3pEnds | sort -T {config[TMPDIR]}  -k1,1 -k2,2n -k3,3n  | bedtools slop -s -l 50 -r -10 -i stdin -g {input.genome} | bedtools intersect -u -s -a stdin -b {input.PAS} |cut -f4 |sort -T {config[TMPDIR]} |uniq > {config[TMPDIR]}/$uuidPASsupported

uuidcagePASsupported=$(uuidgen)
comm -1 -2 {config[TMPDIR]}/$uuidCageSupported {config[TMPDIR]}/$uuidPASsupported |sort -T {config[TMPDIR]} |uniq | awk '{{print "transcript_id \\""$1"\\";"}}'> {config[TMPDIR]}/$uuidcagePASsupported

zcat {input.tm} | fgrep -w -f {config[TMPDIR]}/$uuidcagePASsupported - |sort -T {config[TMPDIR]}  -k1,1 -k4,4n -k5,5n  | gzip> {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}

		'''



rule getFlLocusStats:
	input:
		gencode="annotations/simplified/{capDesign}.gencode.collapsed.simplified_biotypes.gtf",
		clsGencode="mappings/nonAnchoredMergeReads/mergeWithPrevious+biotypes/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.tmerge.min{minReadSupport}reads.endSupport:all.gff.gz",
		gencodeFL="mappings/nonAnchoredMergeReads/mergeWithPrevious/gencode/{capDesign}.gencode.cage+PASsupported.gff.gz",
		previous="mappings/nonAnchoredMergeReads/previous+biotypes/{capDesign}.previous.tmerge.all.gff.gz",
		previousFL="mappings/nonAnchoredMergeReads/previous/{capDesign}.previous.tmerge.cage+PASsupported.gff.gz",
		clsGencodeFL="mappings/nonAnchoredMergeReads/mergeWithPrevious/cls/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.tmerge.min{minReadSupport}reads.cage+PASsupported.gff.gz"
	output: config["STATSDATADIR"] + "tmp/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.FLloci.stats.tsv"
	shell:
		'''
#PCG stats
gencodePcg=$(cat {input.gencode} | fgrep "gene_type \\"protein_coding\\";" | extractGffAttributeValue.pl gene_id|sort -T {config[TMPDIR]} |uniq|wc -l)
gencodePcgFL=$(zcat {input.gencodeFL} | fgrep "gene_type \\"protein_coding\\";" | extractGffAttributeValue.pl gene_id|sort -T {config[TMPDIR]} |uniq|wc -l)
clsGencodePcg=$(zcat {input.clsGencode} | fgrep "gene_type \\"protein_coding\\";" | extractGffAttributeValue.pl gene_id|sort -T {config[TMPDIR]} |uniq|wc -l)
clsGencodePcgFL=$(zcat {input.clsGencodeFL} | fgrep "gene_type \\"protein_coding\\";" | extractGffAttributeValue.pl gene_id|sort -T {config[TMPDIR]} |uniq|wc -l)
previousPcg=$(zcat {input.previous} | fgrep "gene_type \\"protein_coding\\";" | extractGffAttributeValue.pl gene_id|sort -T {config[TMPDIR]} |uniq|wc -l)
previousPcgFL=$(zcat {input.previousFL} | fgrep "gene_type \\"protein_coding\\";" | extractGffAttributeValue.pl gene_id|sort -T {config[TMPDIR]} |uniq|wc -l)

#lncRNA stats
uuidgencodeLnc=$(uuidgen)
uuidgencodeLncFL=$(uuidgen)
uuidclsGencodeLnc=$(uuidgen)
uuidclsGencodeLncFL=$(uuidgen)
uuidPreviousLnc=$(uuidgen)
uuidPreviousLncFL=$(uuidgen)
uuidTmpOut=$(uuidgen)

cat {input.gencode} | fgrep "gene_type \\"lncRNA\\";" |awk '$3=="exon"' > {config[TMPDIR]}/$uuidgencodeLnc
zcat {input.gencodeFL} | fgrep "gene_type \\"lncRNA\\";" |awk '$3=="exon"'> {config[TMPDIR]}/$uuidgencodeLncFL
zcat {input.clsGencode} | fgrep "gene_type \\"lncRNA\\";" |awk '$3=="exon"'> {config[TMPDIR]}/$uuidclsGencodeLnc
zcat {input.clsGencodeFL} | fgrep "gene_type \\"lncRNA\\";" |awk '$3=="exon"'> {config[TMPDIR]}/$uuidclsGencodeLncFL
zcat {input.previous} | fgrep "gene_type \\"lncRNA\\";" |awk '$3=="exon"'> {config[TMPDIR]}/$uuidPreviousLnc
zcat {input.previousFL} | fgrep "gene_type \\"lncRNA\\";" |awk '$3=="exon"'> {config[TMPDIR]}/$uuidPreviousLncFL


gencodeLncSpliced=$(cat {config[TMPDIR]}/$uuidgencodeLnc | extractGffAttributeValue.pl transcript_id |sort -T {config[TMPDIR]} |uniq -d | fgrep -w -f - {config[TMPDIR]}/$uuidgencodeLnc | extractGffAttributeValue.pl gene_id|sort -T {config[TMPDIR]} |uniq|wc -l)
gencodeLncFLSpliced=$(cat {config[TMPDIR]}/$uuidgencodeLncFL| extractGffAttributeValue.pl transcript_id |sort -T {config[TMPDIR]} |uniq -d | fgrep -w -f - {config[TMPDIR]}/$uuidgencodeLncFL | extractGffAttributeValue.pl gene_id|sort -T {config[TMPDIR]} |uniq|wc -l)
clsGencodeLncSpliced=$(cat {config[TMPDIR]}/$uuidclsGencodeLnc | extractGffAttributeValue.pl transcript_id |sort -T {config[TMPDIR]} |uniq -d | fgrep -w -f - {config[TMPDIR]}/$uuidclsGencodeLnc | extractGffAttributeValue.pl gene_id|sort -T {config[TMPDIR]} |uniq|wc -l)
clsGencodeLncFLSpliced=$(cat {config[TMPDIR]}/$uuidclsGencodeLncFL | extractGffAttributeValue.pl transcript_id |sort -T {config[TMPDIR]} |uniq -d | fgrep -w -f - {config[TMPDIR]}/$uuidclsGencodeLncFL | extractGffAttributeValue.pl gene_id|sort -T {config[TMPDIR]} |uniq|wc -l)
previousLncSpliced=$(cat {config[TMPDIR]}/$uuidPreviousLnc | extractGffAttributeValue.pl transcript_id |sort -T {config[TMPDIR]} |uniq -d | fgrep -w -f - {config[TMPDIR]}/$uuidPreviousLnc | extractGffAttributeValue.pl gene_id|sort -T {config[TMPDIR]} |uniq|wc -l)
previousLncFLSpliced=$(cat {config[TMPDIR]}/$uuidPreviousLncFL | extractGffAttributeValue.pl transcript_id |sort -T {config[TMPDIR]} |uniq -d | fgrep -w -f - {config[TMPDIR]}/$uuidPreviousLncFL | extractGffAttributeValue.pl gene_id|sort -T {config[TMPDIR]} |uniq|wc -l)

#clsGencodeLncKnownSpliced=$(cat {config[TMPDIR]}/$uuidclsGencodeLnc | fgrep "gene_ref_status \\"known\\";" |extractGffAttributeValue.pl transcript_id |sort -T {config[TMPDIR]} |uniq -d | fgrep -w -f - {config[TMPDIR]}/$uuidclsGencodeLnc | extractGffAttributeValue.pl gene_id|sort -T {config[TMPDIR]} |uniq|wc -l)
#clsGencodeLncFLKnownSpliced=$(cat {config[TMPDIR]}/$uuidclsGencodeLncFL | fgrep "gene_ref_status \\"known\\";" | extractGffAttributeValue.pl transcript_id |sort -T {config[TMPDIR]} |uniq -d | fgrep -w -f - {config[TMPDIR]}/$uuidclsGencodeLncFL | extractGffAttributeValue.pl gene_id|sort -T {config[TMPDIR]} |uniq|wc -l)



echo -e "{wildcards.techname}Corr{wildcards.corrLevel}\t{wildcards.capDesign}\t{wildcards.sizeFrac}\t{wildcards.barcodes}\t$gencodePcg\t$gencodePcgFL\t$clsGencodePcg\t$clsGencodePcgFL\t$gencodeLncSpliced\t$gencodeLncFLSpliced\t$clsGencodeLncSpliced\t$clsGencodeLncFLSpliced\t$previousLncSpliced\t$previousLncFLSpliced\t$previousPcg\t$previousPcgFL"  > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}


		'''


rule aggFlLocusStats:
	input: lambda wildcards: expand(config["STATSDATADIR"] + "tmp/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.FLloci.stats.tsv", techname='allSeqTechs', corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNSplusMERGED, sizeFrac=SIZEFRACS, barcodes='allTissues', minReadSupport=wildcards.minReadSupport)
	output: config["STATSDATADIR"] + "all.min{minReadSupport}reads.FLloci.stats.tsv"
	shell:
		'''
uuidTmpOut=$(uuidgen)
echo -e "seqTech\tcorrectionLevel\tcapDesign\tsizeFrac\ttissue\tannotation_set\tbiotype\tcategory\tcount\tpercent" > {config[TMPDIR]}/$uuidTmpOut

cat {input} | awk '{{print $1"\\t"$2"\\t"$3"\\t"$4"\\tGENCODE\\tprotein-coding\\tnon-FL\\t"$5-$6"\\t"($5-$6)/$5"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\tGENCODE\\tprotein-coding\\tFL\\t"$6"\\t"$6/$5"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\tGENCODE+previous+current\\tprotein-coding\\tnon-FL\\t"$7-$8"\\t"($7-$8)/$7"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\tGENCODE+previous+current\\tprotein-coding\\tFL\\t"$8"\\t"$8/$7"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\tGENCODE\\tlncRNA\\tnon-FL\\t"$9-$10"\\t"($9-$10)/$9"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\tGENCODE\\tlncRNA\\tFL\\t"$10"\\t"$10/$9"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\tGENCODE+previous+current\\tlncRNA\\tnon-FL\\t"$11-$12"\\t"($11-$12)/$11"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\tGENCODE+previous+current\\tlncRNA\\tFL\\t"$12"\\t"$12/$11"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\tGENCODE+previous\\tlncRNA\\tnon-FL\\t"$13-$14"\\t"($13-$14)/$13"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\tGENCODE+previous\\tlncRNA\\tFL\\t"$14"\\t"$14/$13"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\tGENCODE+previous\\tprotein-coding\\tnon-FL\\t"$15-$16"\\t"($15-$16)/$15"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\tGENCODE+previous\\tprotein-coding\\tFL\\t"$16"\\t"$16/$15}}' | sed 's/Corr0/\tNo/' | sed 's/Corr{lastK}/\tYes/' | sort -T {config[TMPDIR]}  >> {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
		'''

rule plotFlLocusGencodeOnlyStats:
	input: config["STATSDATADIR"] + "all.min{minReadSupport}reads.FLloci.stats.tsv"
	output: returnPlotFilenames(config["PLOTSDIR"] + "FLloci.gencodeOnly.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.FLloci.gencodeOnly.stats")
	conda: "envs/R_env.yml"
	params:
		filterDat=lambda wildcards: merge_figures_params(wildcards.capDesign, wildcards.sizeFrac, wildcards.barcodes, wildcards.corrLevel, wildcards.techname)
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
{params.filterDat[10]}
{params.filterDat[0]}
{params.filterDat[1]}
{params.filterDat[2]}
{params.filterDat[3]}
{params.filterDat[4]}
{params.filterDat[5]}
{params.filterDat[8]}


dat <- subset(dat, annotation_set=='GENCODE')

dat\$category<-factor(dat\$category, ordered=TRUE, levels=rev(c('FL', 'non-FL')))

dat\$annotbiotype <- paste(sep='', dat\$biotype, '\\n(', dat\$annotation_set, ')')
hXyPlot = hXyPlot +1
wXyPlot = wXyPlot +0.5

plotBase <- \\"p <- ggplot(dat[order(dat\$category), ], aes(x=factor(annotbiotype), y=count, fill=category)) +
geom_bar(stat='identity') +
ylab('# gene loci') +
scale_y_continuous(labels=comma)+
scale_fill_manual (values=c('FL'='#C453C4', 'non-FL'='#a6a6a6'))+
theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
xlab('Biotype (Annotation set)') +
guides(fill = guide_legend(title='Category'))+
geom_text(position = 'stack', size=geom_textSize, aes(x = factor(annotbiotype), y = count, label = paste(sep='',percent(round(percent, digits=2)),' / ','(',comma(count),')'), hjust = 0.5, vjust = 1))+
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


rule plotFlLocusStats:
	input: config["STATSDATADIR"] + "all.min{minReadSupport}reads.FLloci.stats.tsv"
	output: returnPlotFilenames(config["PLOTSDIR"] + "FLloci.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.FLloci.stats")
	conda: "envs/R_env.yml"
	params:
		filterDat=lambda wildcards: merge_figures_params(wildcards.capDesign, wildcards.sizeFrac, wildcards.barcodes, wildcards.corrLevel, wildcards.techname)
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
{params.filterDat[10]}
{params.filterDat[0]}
{params.filterDat[1]}
{params.filterDat[2]}
{params.filterDat[3]}
{params.filterDat[4]}
{params.filterDat[5]}
{params.filterDat[8]}



dat\$category<-factor(dat\$category, ordered=TRUE, levels=rev(c('FL', 'non-FL')))

dat\$annotbiotype <- paste(sep='', dat\$biotype, '\\n(', dat\$annotation_set, ')')

hXyPlot = hXyPlot +1
wXyPlot = wXyPlot +2.5

plotBase <- \\"p <- ggplot(dat[order(dat\$category), ], aes(x=factor(annotbiotype), y=count, fill=category)) +
geom_bar(stat='identity') +
ylab('# gene loci') +
scale_y_continuous(labels=comma)+
scale_fill_manual (values=c('FL'='#C453C4', 'non-FL'='#a6a6a6'))+
theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
xlab('Biotype (Annotation set)') +
guides(fill = guide_legend(title='Category'))+
geom_text(position = 'stack', size=geom_textSize, aes(x = factor(annotbiotype), y = count, label = paste(sep='',percent(round(percent, digits=2)),' ','(',comma(count),')'), hjust = 0.5, vjust = 1))+
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
