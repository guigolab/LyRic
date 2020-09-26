

# mapping of long reads:
rule readMapping:
	input:
		reads = lambda wildcards: expand(DEMULTIPLEXED_FASTQS + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.fastq.gz", filtered_product, techname=wildcards.techname, corrLevel=wildcards.corrLevel, capDesign=wildcards.capDesign, sizeFrac=wildcards.sizeFrac,barcodes=wildcards.barcodes),
		genome = lambda wildcards: config["GENOMESDIR"] + CAPDESIGNTOGENOME[wildcards.capDesign] + ".sorted.fa",
		qc = FQ_CORR_PATH + "qc/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.dupl.txt" if config["DEMULTIPLEX"] else  FQ_CORR_PATH + "qc/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.{barcodes}.dupl.txt",

	threads: 12
	params:
		minimap_preset = lambda wildcards: "splice" if (sampleNameToSeqPlatform(wildcards.techname + "_" + wildcards.capDesign + "_" + wildcards.sizeFrac + "_" + wildcards.barcodes) == "ONT") else "splice:hq" if (sampleNameToSeqPlatform(wildcards.techname + "_" + wildcards.capDesign + "_" + wildcards.sizeFrac + "_" + wildcards.barcodes) == "pacBioSI" or sampleNameToSeqPlatform(wildcards.techname + "_" + wildcards.capDesign + "_" + wildcards.sizeFrac + "_" + wildcards.barcodes) == "pacBioSII") else None
	output:
		bam="mappings/readMapping/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bam",
		bai="mappings/readMapping/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bam.bai",
		bigwig="mappings/readMapping/bigWig/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bw"

	wildcard_constraints:
		barcodes='(?!allTissues).+',
		sizeFrac='[0-9-+\.]+',
		techname='(?!allSeqTechs).+'
	shell:
		'''
uuid=$(uuidgen)
uuidTmpOut=$(uuidgen)
set +eu
conda activate julenv
set -eu

echoerr "Mapping"

minimap2 --MD -x {params.minimap_preset} -t {threads} --secondary=no -L -a {input.genome} {input.reads} > {config[TMPDIR]}/$uuid
echoerr "Mapping done"
echoerr "Creating/sorting BAM"
samtools view -H {config[TMPDIR]}/$uuid > {config[TMPDIR]}/$uuid.2
samtools view -F 256 -F4 -F 2048 {config[TMPDIR]}/$uuid >> {config[TMPDIR]}/$uuid.2
cat {config[TMPDIR]}/$uuid.2 | samtools sort -T {config[TMPDIR]}  --threads {threads}  -m 5G - > {config[TMPDIR]}/$uuidTmpOut
echoerr "Done creating/sorting BAM"
sleep 200s
samtools index {config[TMPDIR]}/$uuidTmpOut

echoerr "Making BigWigs"
bamCoverage -b {config[TMPDIR]}/$uuidTmpOut -o {config[TMPDIR]}/$uuidTmpOut.bw --normalizeUsing CPM 

echoerr "Done making BigWigs"
mv {config[TMPDIR]}/$uuidTmpOut {output.bam}
mv {config[TMPDIR]}/$uuidTmpOut.bai {output.bai}
mv {config[TMPDIR]}/$uuidTmpOut.bw {output.bigwig}
		'''

rule bamqc:
	input: "mappings/readMapping/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bam"
	output: "mappings/readMapping/bamqc/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}/genome_results.txt"
	shell:
		'''
 unset DISPLAY
 ~/bin/qualimap_v2.2.1/qualimap bamqc -bam  {input} -outdir $(dirname {output}) --java-mem-size=25G 

		'''

rule getBamqcStats:
	input: "mappings/readMapping/bamqc/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}/genome_results.txt"
	output: config["STATSDATADIR"] + "tmp/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.sequencingError.stats.tsv"
	shell:
		'''
uuid=$(uuidgen)
qualimapReportToTsv.pl {input}  | cut -f2,3 |grep -v globalErrorRate| sed 's/PerMappedBase//' |awk -v t={wildcards.techname}Corr{wildcards.corrLevel} -v c={wildcards.capDesign} -v si={wildcards.sizeFrac} -v b={wildcards.barcodes} '{{print t"\t"c"\t"si"\t"b"\t"$1"\t"$2}}'| sed 's/Corr0/\tNo/' | sed 's/Corr{lastK}/\tYes/' > {config[TMPDIR]}/$uuid

mv {config[TMPDIR]}/$uuid {output}
		'''

rule aggBamqcStats:
	input: lambda wildcards: expand(config["STATSDATADIR"] + "tmp/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.sequencingError.stats.tsv",filtered_product_merge, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES)
	output: config["STATSDATADIR"] + "all.sequencingError.stats.tsv"
	shell:
		'''
uuid=$(uuidgen)
echo -e "seqTech\tcorrectionLevel\tcapDesign\tsizeFrac\ttissue\terrorCategory\terrorRate" > {config[TMPDIR]}/$uuid
cat {input} | sort >> {config[TMPDIR]}/$uuid
mv {config[TMPDIR]}/$uuid {output}

		'''

rule plotBamqcStats:
	input: config["STATSDATADIR"] + "all.sequencingError.stats.tsv"
	output: returnPlotFilenames(config["PLOTSDIR"] + "sequencingError.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.sequencingError.stats")
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
library(data.table)

dat <- read.table('{input}', header=T, as.is=T, sep='\\t')
{params.filterDat[10]}
{params.filterDat[0]}
{params.filterDat[1]}
{params.filterDat[2]}
{params.filterDat[3]}
{params.filterDat[4]}
{params.filterDat[5]}
{params.filterDat[8]}


plotBase <- \\"p <- ggplot(data=dat, aes(x=factor(correctionLevel), y=errorRate, fill=errorCategory)) +
geom_bar(stat='identity') + scale_fill_manual(values=c(deletions = '#999999', insertions = '#E69F00', mismatches = '#cc0099')) + ylab('# Errors per mapped base') + xlab('') + guides(fill = guide_legend(title='Error class')) +
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
 set +eu
conda activate R_env
set -eu
cat $(dirname {output[0]})/$(basename {output[0]} .legendOnly.png).r | R --slave



		'''


rule getReadProfileMatrix:
	input:
		annot="annotations/{capDesign}.bed",
		bw=lambda wildcards: expand("mappings/readMapping/bigWig/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bw", filtered_product, techname=TECHNAMES, corrLevel=wildcards.corrLevel, capDesign=wildcards.capDesign, sizeFrac=wildcards.sizeFrac, barcodes=wildcards.barcodes),
		sampleAnnot=config["SAMPLE_ANNOT"],

	output: 
		matrix=config["STATSDATADIR"] + "byTechCorr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.readProfileMatrix.tsv.gz",
		colorList=config["STATSDATADIR"] + "byTechCorr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.colors.txt",
		libraryPrepList=config["STATSDATADIR"] + "byTechCorr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.libraryPreps.txt"
	threads: 6
	shell:
		'''
uuid=$(uuidgen)

# extract libraryPrep names and matching colors, in the same order as input files

echo "
library(tidyverse); 
sampleAnnot <- read.table('{input.sampleAnnot}', header=T, as.is=T, sep='\\t') 
filesList <- read.table(file('stdin'), header=FALSE, sep=' '); 
colnames(filesList) <- c('sample_name');  
sampleAnnotation_colors = {config[SAMPLE_ANNOT_RPALETTE]}

#convert to dataframe
sampleAnnotation_colors_df <- setNames(stack(sampleAnnotation_colors\$libraryPrep)[2:1], c('libraryPrep','color'))
print(sampleAnnotation_colors_df)

join1 <- inner_join(filesList, sampleAnnot, by = 'sample_name')
join2 <- inner_join(join1, sampleAnnotation_colors_df, by='libraryPrep')

print(join1)
print(join2)

outTable <- select(join2, sample_name, libraryPrep, color)
write(outTable\$libraryPrep, '{output.libraryPrepList}', ncolumns=nrow(outTable))
write(outTable\$color, '{output.colorList}', ncolumns=nrow(outTable))


" > {output.matrix}.r

set +eu
conda activate R_env
set -eu

echo {input.bw} | xargs -n1 basename | sed 's/\.bw//' | sed 's/Corr0_/_/' | Rscript {output.matrix}.r

set +eu
conda activate julenv
set -eu

computeMatrix scale-regions -S {input.bw} -R {input.annot} -o {config[TMPDIR]}/$uuid.gz --upstream 1000 --downstream 1000 --sortRegions ascend  --missingDataAsZero --skipZeros --metagene -p {threads} --samplesLabel $(cat {output.libraryPrepList} | perl -ne 'chomp; print')

mv {config[TMPDIR]}/$uuid.gz {output.matrix}
		'''

rule plotReadProfileMatrix:
	input: 
		matrix=config["STATSDATADIR"] + "byTechCorr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.readProfileMatrix.tsv.gz",
		colorList=config["STATSDATADIR"] + "byTechCorr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.colors.txt"
	output: 
		profile=config["PLOTSDIR"] + "readProfile/byTechCorr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.readProfile.density.png",
		heatmap=config["PLOTSDIR"] + "readProfile/byTechCorr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.readProfile.heatmap.png"
	shell:
		'''
	
set +eu
conda activate julenv
set -eu
	
  plotProfile -m {input.matrix} -o {output.profile} --perGroup  --plotType se --yAxisLabel "mean CPM" --regionsLabel '' --colors $(cat {input.colorList} | perl -ne 'chomp; print')

  plotHeatmap -m {input.matrix} -o {output.heatmap} --perGroup  --plotType se --yAxisLabel "mean CPM" --regionsLabel '' --whatToShow 'heatmap and colorbar'

		'''





rule getMappingStats:
	input:
		bams = "mappings/readMapping/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bam",
		fastqs = DEMULTIPLEXED_FASTQS + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.fastq.gz"
	output: 
		basic=config["STATSDATADIR"] + "tmp/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.{barcodes}.mapping.stats.tsv",
		spikeIns=config["STATSDATADIR"] + "tmp/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.{barcodes}.mapping.spikeIns.stats.tsv"
	shell:
		'''
uuidTmpOutB=$(uuidgen)
uuidTmpOutS=$(uuidgen)
totalReads=$(zcat {input.fastqs} | fastq2tsv.pl | wc -l)
mappedReads=$(samtools view  -F 4 {input.bams}|cut -f1|sort -T {config[TMPDIR]} |uniq|wc -l)
erccMappedReads=$(samtools view -F 4 {input.bams}|cut -f3| tgrep ERCC- | wc -l)
sirvMappedReads=$(samtools view -F 4 {input.bams}|cut -f3| tgrep SIRV | wc -l)
echo -e "{wildcards.techname}Corr{wildcards.corrLevel}\t{wildcards.capDesign}\t{wildcards.sizeFrac}\t{wildcards.barcodes}\t$totalReads\t$mappedReads" | awk '{{print $0"\t"$6/$5}}' > {config[TMPDIR]}/$uuidTmpOutB
mv {config[TMPDIR]}/$uuidTmpOutB {output.basic}
echo -e "{wildcards.techname}Corr{wildcards.corrLevel}\t{wildcards.capDesign}\t{wildcards.sizeFrac}\t{wildcards.barcodes}\t$totalReads\t$erccMappedReads\t$sirvMappedReads" | awk '{{print $0"\t"$6/$5"\t"$7/$5}}' > {config[TMPDIR]}/$uuidTmpOutS
mv {config[TMPDIR]}/$uuidTmpOutS {output.spikeIns}

		'''

rule aggMappingStats:
	input: lambda wildcards: expand(config["STATSDATADIR"] + "tmp/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.{barcodes}.mapping.stats.tsv",filtered_product, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES)
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

plotBase <- \\"p <- ggplot(data=dat, aes(x=1, y=percentMappedReads, fill=sizeFrac)) +
geom_bar(width=0.75,stat='identity', position=position_dodge(width=0.9)) +
scale_fill_manual(values={sizeFrac_Rpalette}) +
geom_hline(aes(yintercept=1), linetype='dashed', alpha=0.7)+
geom_text(aes(group=sizeFrac, y=0.01, label = paste(sep='',percent(percentMappedReads),'\\n','(',comma(mappedReads),')')), angle=90, size=geom_textSize, hjust=0, vjust=0.5, position = position_dodge(width=0.9)) +
scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
xlab ('') +
{GGPLOT_PUB_QUALITY} + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) + \\"

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
 set +eu
conda activate R_env
set -eu
cat $(dirname {output[0]})/$(basename {output[0]} .legendOnly.png).r | R --slave

		'''

rule aggSpikeInsMappingStats:
	input: lambda wildcards: expand(config["STATSDATADIR"] + "tmp/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.{barcodes}.mapping.spikeIns.stats.tsv",filtered_product, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES)
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

maxY <- max(dat\$percent)
plotBase <- \\"p <- ggplot(data=dat, aes(x=factor(correctionLevel), y=percent, fill=category)) +
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
 set +eu
conda activate R_env
set -eu
cat $(dirname {output[0]})/$(basename {output[0]} .legendOnly.png).r | R --slave

		'''


rule mergeBarcodeBams:
	input: lambda wildcards: expand("mappings/readMapping/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bam", filtered_product, techname=wildcards.techname, corrLevel=wildcards.corrLevel, capDesign=wildcards.capDesign, sizeFrac=wildcards.sizeFrac, barcodes=BARCODES)
	output: 
		bam="mappings/readMapping/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_allTissues.bam",
		bai="mappings/readMapping/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_allTissues.bam.bai",
		bigwig="mappings/readMapping/bigWig/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_allTissues.bw"
	wildcard_constraints:
		sizeFrac='[0-9-+\.]+',
		techname='(?!allSeqTechs).+' #to avoid ambiguity with downstream merging rules
	shell:
		'''
uuidTmpOut=$(uuidgen)
samtools merge {config[TMPDIR]}/$uuidTmpOut {input}
sleep 200s
samtools index {config[TMPDIR]}/$uuidTmpOut
bamCoverage -b {config[TMPDIR]}/$uuidTmpOut -o {config[TMPDIR]}/$uuidTmpOut.bw --normalizeUsing CPM 

mv {config[TMPDIR]}/$uuidTmpOut {output.bam}
mv {config[TMPDIR]}/$uuidTmpOut.bai {output.bai}
mv {config[TMPDIR]}/$uuidTmpOut.bw {output.bigwig}

		'''


rule mergeAllSeqTechsFracsAndTissuesBams:
	input: lambda wildcards: expand("mappings/readMapping/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bam", filtered_product_merge, techname=TECHNAMES, corrLevel=wildcards.corrLevel, capDesign=wildcards.capDesign, sizeFrac=wildcards.sizeFrac, barcodes=wildcards.barcodes)
	output: 
		bam="mappings/readMapping/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bam",
		bai="mappings/readMapping/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bam.bai",
		bigwig="mappings/readMapping/bigWig/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bw"
	wildcard_constraints:
#		techname='allSeqTechs',
#		sizeFrac='allFracs',
		barcodes='allTissues',
		capDesign='|'.join(CAPDESIGNS)
	shell:
		'''
uuidTmpOut=$(uuidgen)
samtools merge {config[TMPDIR]}/$uuidTmpOut {input}
sleep 200s
samtools index {config[TMPDIR]}/$uuidTmpOut
bamCoverage -b {config[TMPDIR]}/$uuidTmpOut -o {config[TMPDIR]}/$uuidTmpOut.bw --normalizeUsing CPM 

mv {config[TMPDIR]}/$uuidTmpOut {output.bam}
mv {config[TMPDIR]}/$uuidTmpOut.bai {output.bai}
mv {config[TMPDIR]}/$uuidTmpOut.bw {output.bigwig}

		'''


rule mergeAllCapDesignsAllSeqTechsFracsAndTissuesBams:
	input: lambda wildcards: expand("mappings/readMapping/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bam", filtered_capDesign_product_merge, techname=wildcards.techname, corrLevel=wildcards.corrLevel, capDesign=wildcards.capDesign, sizeFrac=wildcards.sizeFrac, barcodes=wildcards.barcodes)
	output: 
		bam="mappings/readMapping/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bam",
		bai="mappings/readMapping/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bam.bai",
		bigwig="mappings/readMapping/bigWig/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bw"
	wildcard_constraints:
#		techname='allSeqTechs',
#		sizeFrac='allFracs',
		barcodes='allTissues',
		capDesign='|'.join(MERGEDCAPDESIGNS)
	shell:
		'''
uuidTmpOut=$(uuidgen)
samtools merge {config[TMPDIR]}/$uuidTmpOut {input}
sleep 200s
samtools index {config[TMPDIR]}/$uuidTmpOut
bamCoverage -b {config[TMPDIR]}/$uuidTmpOut -o {config[TMPDIR]}/$uuidTmpOut.bw --normalizeUsing CPM 

mv {config[TMPDIR]}/$uuidTmpOut {output.bam}
mv {config[TMPDIR]}/$uuidTmpOut.bai {output.bai}
mv {config[TMPDIR]}/$uuidTmpOut.bw {output.bigwig}

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
set +eu
conda activate julenv
set -eu
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
