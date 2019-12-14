

# mapping of long reads:
rule readMapping:
	input:
		reads = lambda wildcards: expand(DEMULTIPLEXED_FASTQS + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.{barcodes}.fastq.gz", filtered_product, techname=wildcards.techname, corrLevel=wildcards.corrLevel, capDesign=wildcards.capDesign, sizeFrac=wildcards.sizeFrac,barcodes=wildcards.barcodes),
		genome = lambda wildcards: config["GENOMESDIR"] + CAPDESIGNTOGENOME[wildcards.capDesign] + ".fa"
	threads: 12
	params:
		minimap_preset = lambda wildcards: "splice" if wildcards.techname.find('ont') == 0 else "splice:hq" if wildcards.techname.find('pacBio') == 0 else None
	output:
		"mappings/readMapping/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bam"
	wildcard_constraints:
		barcodes='(?!allTissues).+',
		sizeFrac='[0-9-+\.]+',
		techname='(?!allSeqTechs).+'
	shell:
		'''
uuid=$(uuidgen)
uuidTmpOut=$(uuidgen)
echoerr "Mapping"
minimap2 -x {params.minimap_preset} -t {threads} --secondary=no -L -a {input.genome} {input.reads} > {config[TMPDIR]}/$uuid
echoerr "Mapping done"
echoerr "Creating/sorting BAM"
samtools view -H {config[TMPDIR]}/$uuid > {config[TMPDIR]}/$uuid.2
samtools view -F 256 -F4 -F 2048 {config[TMPDIR]}/$uuid >> {config[TMPDIR]}/$uuid.2
cat {config[TMPDIR]}/$uuid.2 | samtools sort -T {config[TMPDIR]}  --threads {threads}  -m 5G - > {config[TMPDIR]}/$uuidTmpOut
echoerr "Done creating/sorting BAM"
sleep 70s
samtools index {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
mv {config[TMPDIR]}/$uuidTmpOut.bai {output}.bai

		'''

rule getMappingStats:
	input:
		bams = "mappings/readMapping/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bam",
		fastqs = DEMULTIPLEXED_FASTQS + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.{barcodes}.fastq.gz"
	output: 
		basic=temp(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.{barcodes}.mapping.stats.tsv"),
		spikeIns=temp(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.{barcodes}.mapping.spikeIns.stats.tsv")
	shell:
		'''
uuidTmpOutB=$(uuidgen)
uuidTmpOutS=$(uuidgen)
totalReads=$(zcat {input.fastqs} | fastq2tsv.pl | wc -l)
mappedReads=$(samtools view  -F 4 {input.bams}|cut -f1|sort -T {config[TMPDIR]} |uniq|wc -l)
erccMappedReads=$(samtools view -F 4 {input.bams}|cut -f3| tgrep ERCC- | wc -l)
sirvMappedReads=$(samtools view -F 4 {input.bams}|cut -f3| tgrep SIRVome_isoforms | wc -l)
echo -e "{wildcards.techname}Corr{wildcards.corrLevel}\t{wildcards.capDesign}\t{wildcards.sizeFrac}\t{wildcards.barcodes}\t$totalReads\t$mappedReads" | awk '{{print $0"\t"$6/$5}}' > {config[TMPDIR]}/$uuidTmpOutB
mv {config[TMPDIR]}/$uuidTmpOutB {output.basic}
echo -e "{wildcards.techname}Corr{wildcards.corrLevel}\t{wildcards.capDesign}\t{wildcards.sizeFrac}\t{wildcards.barcodes}\t$totalReads\t$erccMappedReads\t$sirvMappedReads" | awk '{{print $0"\t"$6/$5"\t"$7/$5}}' > {config[TMPDIR]}/$uuidTmpOutS
mv {config[TMPDIR]}/$uuidTmpOutS {output.spikeIns}

		'''

rule aggMappingStats:
	input: lambda wildcards: expand(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.{barcodes}.mapping.stats.tsv",filtered_product, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES)
	output: config["STATSDATADIR"] + "all.basic.mapping.stats.tsv"
	shell:
		'''
uuidTmpOut=$(uuidgen)
echo -e "seqTech\tcorrectionLevel\tcapDesign\tsizeFrac\ttissue\ttotalReads\tmappedReads\tpercentMappedReads" > {config[TMPDIR]}/$uuidTmpOut
cat {input} | sed 's/Corr0/\tNo/' | sed 's/Corr{lastK}/\tYes/' | sort -T {config[TMPDIR]}  >> {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}

		'''

rule plotMappingStats:
	input: config["STATSDATADIR"] + "all.basic.mapping.stats.tsv"
	output: returnPlotFilenames(config["PLOTSDIR"] + "lrMapping.basic.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.lrMapping.basic.stats")
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

plotBase <- \\"ggplot(data=dat, aes(x=tissue, y=percentMappedReads, fill=sizeFrac)) +
geom_bar(width=0.75,stat='identity', position=position_dodge(width=0.9)) +
scale_fill_manual(values={sizeFrac_Rpalette}) +
geom_hline(aes(yintercept=1), linetype='dashed', alpha=0.7)+
geom_text(aes(group=sizeFrac, y=0.01, label = paste(sep='',percent(percentMappedReads),'\\n','(',comma(mappedReads),')')), angle=90, size=geom_textSize, hjust=0, vjust=0.5, position = position_dodge(width=0.9)) +
scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
xlab ('Sample') +
{GGPLOT_PUB_QUALITY} + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + \\"

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

rule aggSpikeInsMappingStats:
	input: lambda wildcards: expand(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.{barcodes}.mapping.spikeIns.stats.tsv",filtered_product, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES)
	output: config["STATSDATADIR"] + "all.spikeIns.mapping.stats.tsv"
	shell:
		'''
uuidTmpOut=$(uuidgen)
echo -e "seqTech\tcorrectionLevel\tcapDesign\tsizeFrac\ttissue\tcategory\tcount\tpercent" > {config[TMPDIR]}/$uuidTmpOut
cat {input} | awk '{{print $1"\\t"$2"\\t"$3"\\t"$4"\\tSIRVs\\t"$7"\\t"$9"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\tERCCs\\t"$6"\\t"$8}}' | sed 's/Corr0/\tNo/' | sed 's/Corr{lastK}/\tYes/' | sort -T {config[TMPDIR]}  >> {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}

		'''

rule plotSpikeInsMappingStats:
	input: config["STATSDATADIR"] + "all.spikeIns.mapping.stats.tsv"
	output: returnPlotFilenames(config["PLOTSDIR"] + "lrMapping.spikeIns.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.lrMapping.spikeIns.stats")
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

maxY <- max(dat\$percent)
plotBase <- \\"ggplot(data=dat, aes(x=factor(correctionLevel), y=percent, fill=category)) +
geom_bar(stat='identity', position=position_dodge()) +
geom_text(position = position_dodge(width = 0.9), size=geom_textSize, aes(x = factor(correctionLevel), y = 0, label = paste(sep='',percent(percent),'\\n','(',comma(count),')')), hjust = 0, vjust = 0.5, angle=90) +
scale_fill_manual(values=c('ERCCs' = '#e4b5ff', 'SIRVs' = '#5edba9')) +
ylab('% reads mapped on\\nspike-in sequences') +
xlab('{params.filterDat[6]}') +
guides(fill = guide_legend(title='Spike-in set'))+
scale_y_continuous(labels=percent)+
expand_limits(y=c(0,maxY))+
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


rule mergeBarcodeBams:
	input: lambda wildcards: expand("mappings/readMapping/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bam", filtered_product, techname=wildcards.techname, corrLevel=wildcards.corrLevel, capDesign=wildcards.capDesign, sizeFrac=wildcards.sizeFrac, barcodes=BARCODES)
	output: "mappings/readMapping/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_allTissues.bam"
	wildcard_constraints:
		sizeFrac='[0-9-+\.]+',
		techname='(?!allSeqTechs).+' #to avoid ambiguity with downstream merging rules
	shell:
		'''
uuidTmpOut=$(uuidgen)
samtools merge {config[TMPDIR]}/$uuidTmpOut {input}
sleep 120s
samtools index {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
mv {config[TMPDIR]}/$uuidTmpOut.bai {output}.bai

		'''


rule mergeAllSeqTechsFracsAndTissuesBams:
	input: lambda wildcards: expand("mappings/readMapping/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bam", filtered_product_merge, techname=TECHNAMES, corrLevel=wildcards.corrLevel, capDesign=wildcards.capDesign, sizeFrac=wildcards.sizeFrac, barcodes=wildcards.barcodes)
	output: "mappings/readMapping/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bam"
	wildcard_constraints:
#		techname='allSeqTechs',
#		sizeFrac='allFracs',
		barcodes='allTissues',
		capDesign='|'.join(CAPDESIGNS)
	shell:
		'''
uuidTmpOut=$(uuidgen)
samtools merge {config[TMPDIR]}/$uuidTmpOut {input}
sleep 120s
samtools index {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
mv {config[TMPDIR]}/$uuidTmpOut.bai {output}.bai

		'''


rule mergeAllCapDesignsAllSeqTechsFracsAndTissuesBams:
	input: lambda wildcards: expand("mappings/readMapping/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bam", filtered_capDesign_product_merge, techname=wildcards.techname, corrLevel=wildcards.corrLevel, capDesign=wildcards.capDesign, sizeFrac=wildcards.sizeFrac, barcodes=wildcards.barcodes)
	output: "mappings/readMapping/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bam"
	wildcard_constraints:
#		techname='allSeqTechs',
#		sizeFrac='allFracs',
		barcodes='allTissues',
		capDesign='|'.join(MERGEDCAPDESIGNS)
	shell:
		'''
uuidTmpOut=$(uuidgen)
samtools merge {config[TMPDIR]}/$uuidTmpOut {input}
sleep 120s
samtools index {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
mv {config[TMPDIR]}/$uuidTmpOut.bai {output}.bai

		'''


rule checkOnlyOneHit:
	input: "mappings/readMapping/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bam"
	output: "mappings/readMapping/qc/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bam.dupl.txt"
	shell:
		'''
uuidTmpOut=$(uuidgen)
samtools view {input} | cut -f1 | sort -T {config[TMPDIR]} | uniq -dc > {config[TMPDIR]}/$uuidTmpOut
count=$(cat {config[TMPDIR]}/$uuidTmpOut | wc -l)
if [ $count -gt 0 ]; then echo "$count duplicate read IDs found"; mv {input} {input}.dup.bkp; exit 1; fi
mv {config[TMPDIR]}/$uuidTmpOut {output}

		'''


rule readBamToBed:
	input: "mappings/readMapping/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bam"
	output: "mappings/readBamToBed/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bed.gz"
	shell:
		'''
uuidTmpOut=$(uuidgen)
#remove mappings with ultra-short exons after bamtobed
bedtools bamtobed -i {input} -bed12 | perl -ne '$line=$_; @line=split ("\\t", $line); @blockSizes=split(",", $line[10]); $allExonsOK=1; foreach $block (@blockSizes){{if ($block<2){{$allExonsOK=0; last;}}}}; if ($allExonsOK==1){{print $line}}'| sort -T {config[TMPDIR]}  -k1,1 -k2,2n -k3,3n  | gzip > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}

		'''

rule readBedToGff:
	input: "mappings/readBamToBed/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bed.gz"
	output: "mappings/readBedToGff/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.gff.gz"
	shell:
		'''
uuidTmpOut=$(uuidgen)
zcat {input} | awk -f ~jlagarde/julien_utils/bed12fields2gff.awk | sort -T {config[TMPDIR]}  -k1,1 -k4,4n -k5,5n  | gzip > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
		'''


# rule qualimap:
# 	input: "mappings/readMapping/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bam"
# 	output: "mappings/qualimap_reports/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}/genome_results.txt"
# 	shell:
# 		'''
# unset DISPLAY
# ~/bin/qualimap_v2.2.1/qualimap bamqc -bam {input} -outdir mappings/qualimap_reports/{wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}.merged2/ --java-mem-size=10G -outfile {wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}.merged2
# touch {output}
# 		'''
