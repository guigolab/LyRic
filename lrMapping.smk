

# mapping of long reads:
rule readMapping:
	input:
		reads = lambda wildcards: expand(LR_FASTQDIR + "{techname}_{capDesign}_{sizeFrac}_{barcodes}.fastq.gz", filtered_product, techname=wildcards.techname,  capDesign=wildcards.capDesign, sizeFrac=wildcards.sizeFrac,barcodes=wildcards.barcodes),
		genome = lambda wildcards: config["GENOMESDIR"] + CAPDESIGNTOGENOME[wildcards.capDesign] + ".sorted.fa",
		qc = LR_FASTQDIR + "qc/{techname}_{capDesign}_{sizeFrac}.{barcodes}.dupl.txt",

	threads: 12
	params:
		minimap_preset = lambda wildcards: "splice" if (sampleAnnotDict[wildcards.techname + "_" + wildcards.capDesign + "_" + wildcards.sizeFrac + "_" + wildcards.barcodes]['seqPlatform'] == "ONT") else "splice:hq" if ( sampleAnnotDict[wildcards.techname + "_" + wildcards.capDesign + "_" + wildcards.sizeFrac + "_" + wildcards.barcodes]['seqPlatform'] == "pacBioSI" or sampleAnnotDict[wildcards.techname + "_" + wildcards.capDesign + "_" + wildcards.sizeFrac + "_" + wildcards.barcodes]['seqPlatform'] == "pacBioSII") else None
	output:
		bam="mappings/readMapping/{techname}_{capDesign}_{sizeFrac}_{barcodes}.bam",
		bai="mappings/readMapping/{techname}_{capDesign}_{sizeFrac}_{barcodes}.bam.bai",
	conda: "envs/minimap2_env.yml"
	wildcard_constraints:
		sizeFrac='[0-9-+\.]+',
	shell:
		'''
uuid=$(uuidgen)
uuidTmpOut=$(uuidgen)

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

mv {config[TMPDIR]}/$uuidTmpOut {output.bam}
mv {config[TMPDIR]}/$uuidTmpOut.bai {output.bai}
		'''

rule makeBigWigs:
	input: "mappings/readMapping/{techname}_{capDesign}_{sizeFrac}_{barcodes}.bam"
	output: "mappings/readMapping/bigWig/{techname}_{capDesign}_{sizeFrac}_{barcodes}.bw"
	conda: "envs/xtools_env.yml"
	shell:
		'''
uuid=$(uuidgen)
bamCoverage --normalizeUsing CPM  -b {input} -o {config[TMPDIR]}/$uuid.bw 
mv {config[TMPDIR]}/$uuid.bw {output}
		'''



rule bamqc:
	input: "mappings/readMapping/{techname}_{capDesign}_{sizeFrac}_{barcodes}.bam"
	output: "mappings/readMapping/bamqc/{techname}_{capDesign}_{sizeFrac}_{barcodes}/genome_results.txt"
	shell:
		'''
 unset DISPLAY
 ~/bin/qualimap_v2.2.1/qualimap bamqc -bam  {input} -outdir $(dirname {output}) --java-mem-size=25G 

		'''

rule getBamqcStats:
	input: "mappings/readMapping/bamqc/{techname}_{capDesign}_{sizeFrac}_{barcodes}/genome_results.txt"
	output: config["STATSDATADIR"] + "tmp/{techname}_{capDesign}_{sizeFrac}_{barcodes}.sequencingError.stats.tsv"
	shell:
		'''
uuid=$(uuidgen)
qualimapReportToTsv.pl {input}  | cut -f2,3 |grep -v globalErrorRate| sed 's/PerMappedBase//' |awk -v t={wildcards.techname} -v c={wildcards.capDesign} -v si={wildcards.sizeFrac} -v b={wildcards.barcodes} '{{print t"\t"c"\t"si"\t"b"\t"$1"\t"$2}}' > {config[TMPDIR]}/$uuid

mv {config[TMPDIR]}/$uuid {output}
		'''

rule aggBamqcStats:
	input: lambda wildcards: expand(config["STATSDATADIR"] + "tmp/{techname}_{capDesign}_{sizeFrac}_{barcodes}.sequencingError.stats.tsv",filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES)
	output: config["STATSDATADIR"] + "all.sequencingError.stats.tsv"
	shell:
		'''
uuid=$(uuidgen)
echo -e "seqTech\tcapDesign\tsizeFrac\ttissue\terrorCategory\terrorRate" > {config[TMPDIR]}/$uuid
cat {input} | sort >> {config[TMPDIR]}/$uuid
mv {config[TMPDIR]}/$uuid {output}

		'''

rule plotBamqcStats:
	input: config["STATSDATADIR"] + "all.sequencingError.stats.tsv"
	output: 
		allErrors=returnPlotFilenames(config["PLOTSDIR"] + "sequencingError.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{barcodes}.sequencingError.allErrors.stats"),
		deletionsOnly=returnPlotFilenames(config["PLOTSDIR"] + "sequencingError.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{barcodes}.sequencingError.deletionsOnly.stats")
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
library(data.table)

dat <- read.table('{input}', header=T, as.is=T, sep='\\t')
{params.filterDat[technameFilterString]}
{params.filterDat[capDesignFilterString]}

{params.filterDat[sizeFracFilterString]}
{params.filterDat[tissueFilterString]}
{params.filterDat[substSeqTechString]}
{params.filterDat[substTissueString]}
{params.filterDat[graphDimensions]}


plotBase <- \\"p <- ggplot(data=dat, aes(x=1, y=errorRate, fill=errorCategory)) +
geom_bar(stat='identity') + scale_fill_manual(values=c(deletions = '#bfbfbf', insertions = '#ffa64d', mismatches = '#1a1aff')) + ylab('# Errors per mapped base') + xlab('') + guides(fill = guide_legend(title='Error class')) +
scale_y_continuous(labels = label_scientific(digits = 1)) +

{params.filterDat[hideXaxisLabels]}
{GGPLOT_PUB_QUALITY} + 
theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
\\"


{params.filterDat[facetPlotSetup]}

save_plot('{output.allErrors[0]}', legendOnly, base_width=wLegendOnly, base_height=hLegendOnly)
save_plot('{output.allErrors[1]}', legendOnly, base_width=wLegendOnly, base_height=hLegendOnly)

save_plot('{output.allErrors[2]}', pXy, base_width=wXyPlot, base_height=hXyPlot)
save_plot('{output.allErrors[3]}', pXy, base_width=wXyPlot, base_height=hXyPlot)

save_plot('{output.allErrors[4]}', pXyNoLegend, base_width=wXyNoLegendPlot, base_height=hXyNoLegendPlot)
save_plot('{output.allErrors[5]}', pXyNoLegend, base_width=wXyNoLegendPlot, base_height=hXyNoLegendPlot)

save_plot('{output.allErrors[6]}', pYx, base_width=wYxPlot, base_height=hYxPlot)
save_plot('{output.allErrors[7]}', pYx, base_width=wYxPlot, base_height=hYxPlot)

save_plot('{output.allErrors[8]}', pYxNoLegend, base_width=wYxNoLegendPlot, base_height=hYxNoLegendPlot)
save_plot('{output.allErrors[9]}', pYxNoLegend, base_width=wYxNoLegendPlot, base_height=hYxNoLegendPlot)



datDeletionsOnly <- subset(dat, errorCategory=='deletions')


plotBase <- \\"p <- ggplot(data=datDeletionsOnly, aes(x=1, y=errorRate, fill=errorCategory)) +
geom_bar(stat='identity') + scale_fill_manual(values=c(deletions = '#bfbfbf')) + ylab('# Errors per mapped base') + xlab('') + guides(fill = guide_legend(title='Error class')) +
scale_y_continuous(labels = label_scientific(digits = 1)) +

{params.filterDat[hideXaxisLabels]}
{GGPLOT_PUB_QUALITY} + 
theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
\\"


{params.filterDat[facetPlotSetup]}

save_plot('{output.deletionsOnly[0]}', legendOnly, base_width=wLegendOnly, base_height=hLegendOnly)
save_plot('{output.deletionsOnly[1]}', legendOnly, base_width=wLegendOnly, base_height=hLegendOnly)

save_plot('{output.deletionsOnly[2]}', pXy, base_width=wXyPlot, base_height=hXyPlot)
save_plot('{output.deletionsOnly[3]}', pXy, base_width=wXyPlot, base_height=hXyPlot)

save_plot('{output.deletionsOnly[4]}', pXyNoLegend, base_width=wXyNoLegendPlot, base_height=hXyNoLegendPlot)
save_plot('{output.deletionsOnly[5]}', pXyNoLegend, base_width=wXyNoLegendPlot, base_height=hXyNoLegendPlot)

save_plot('{output.deletionsOnly[6]}', pYx, base_width=wYxPlot, base_height=hYxPlot)
save_plot('{output.deletionsOnly[7]}', pYx, base_width=wYxPlot, base_height=hYxPlot)

save_plot('{output.deletionsOnly[8]}', pYxNoLegend, base_width=wYxNoLegendPlot, base_height=hYxNoLegendPlot)
save_plot('{output.deletionsOnly[9]}', pYxNoLegend, base_width=wYxNoLegendPlot, base_height=hYxNoLegendPlot)



" > $(dirname {output.allErrors[0]})/$(basename {output.allErrors[0]} .legendOnly.png).r

cat $(dirname {output.allErrors[0]})/$(basename {output.allErrors[0]} .legendOnly.png).r | R --slave



		'''


rule makeBigWigExonicRegions:
	input:
		bam= lambda wildcards: expand("mappings/readMapping/{techname}_{capDesign}_{sizeFrac}_{barcodes}.bam", filtered_product, techname=wildcards.techname,  capDesign=wildcards.capDesign, sizeFrac=wildcards.sizeFrac,barcodes=wildcards.barcodes),
		annotGff=lambda wildcards: CAPDESIGNTOANNOTGTF[wildcards.capDesign]
	conda: "envs/xtools_env.yml"
	output:
		"mappings/readMapping/bigWig_exonic/{techname}_{capDesign}_{sizeFrac}_{barcodes}.bw"
	shell:
		'''
uuid=$(uuidgen)
cat {input.annotGff} | awk '$3=="exon"' > {config[TMPDIR]}/$uuid.gff
bedtools intersect -split -u -a {input.bam} -b {config[TMPDIR]}/$uuid.gff > {config[TMPDIR]}/$uuid.bam
samtools index {config[TMPDIR]}/$uuid.bam

bamCoverage --normalizeUsing CPM  -b {config[TMPDIR]}/$uuid.bam -o {config[TMPDIR]}/$uuid.bw

mv {config[TMPDIR]}/$uuid.bw {output}

		'''


rule setupReadProfileMatrix:
	input: 
		bw=lambda wildcards: expand("mappings/readMapping/bigWig_exonic/{techname}_{capDesign}_{sizeFrac}_{barcodes}.bw", filtered_product, techname=TECHNAMES,  capDesign=wildcards.capDesign, sizeFrac=wildcards.sizeFrac, barcodes=wildcards.barcodes),
		sampleAnnot=config["SAMPLE_ANNOT"]
	output:
		colorList=config["STATSDATADIR"] + "byTech_{capDesign}_{sizeFrac}_{barcodes}.colors.txt",
		libraryPrepList=config["STATSDATADIR"] + "byTech_{capDesign}_{sizeFrac}_{barcodes}.libraryPreps.txt"
	conda: "envs/R_env.yml"
	shell:
		'''
# extract libraryPrep names and matching colors, in the same order as input files

echo "
library(tidyverse); 
sampleAnnot <- read.table('{input.sampleAnnot}', header=T, as.is=T, sep='\\t') 
filesList <- read.table(file('stdin'), header=FALSE, sep=' '); 
colnames(filesList) <- c('sample_name');  
sampleAnnotation_colors = {sampleAnnot_Rpalette}

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


" > {output.colorList}.r

echo {input.bw} | xargs -n1 basename | sed 's/\.bw//' | Rscript {output.colorList}.r

		'''


rule getReadProfileMatrix:
	input:
		annot="annotations/{capDesign}.bed",
		bw=lambda wildcards: expand("mappings/readMapping/bigWig_exonic/{techname}_{capDesign}_{sizeFrac}_{barcodes}.bw", filtered_product, techname=TECHNAMES,  capDesign=wildcards.capDesign, sizeFrac=wildcards.sizeFrac, barcodes=wildcards.barcodes),
		libraryPrepList=config["STATSDATADIR"] + "byTech_{capDesign}_{sizeFrac}_{barcodes}.libraryPreps.txt"

	output: 
		matrix=config["STATSDATADIR"] + "byTech_{capDesign}_{sizeFrac}_{barcodes}.readProfileMatrix.tsv.gz",
	conda: "envs/xtools_env.yml"
	threads: 6
	shell:
		'''
uuid=$(uuidgen)

computeMatrix scale-regions -S {input.bw} -R {input.annot} -o {config[TMPDIR]}/$uuid.gz --upstream 1000 --downstream 1000 --sortRegions ascend  --missingDataAsZero --skipZeros --metagene -p {threads} --samplesLabel $(cat {input.libraryPrepList} | perl -ne 'chomp; print')

mv {config[TMPDIR]}/$uuid.gz {output.matrix}
		'''

rule plotReadProfileMatrix:
	input: 
		matrix=config["STATSDATADIR"] + "byTech_{capDesign}_{sizeFrac}_{barcodes}.readProfileMatrix.tsv.gz",
		colorList=config["STATSDATADIR"] + "byTech_{capDesign}_{sizeFrac}_{barcodes}.colors.txt"
	output: 
		profile=config["PLOTSDIR"] + "readProfile/byTech_{capDesign}_{sizeFrac}_{barcodes}.readProfile.density.png",
		heatmap=config["PLOTSDIR"] + "readProfile/byTech_{capDesign}_{sizeFrac}_{barcodes}.readProfile.heatmap.png"
	conda: "envs/xtools_env.yml"
	shell:
		'''
	

  plotProfile -m {input.matrix} -o {output.profile} --perGroup  --plotType se --yAxisLabel "mean CPM" --regionsLabel '' --colors $(cat {input.colorList} | perl -ne 'chomp; print')

  plotHeatmap -m {input.matrix} -o {output.heatmap} --perGroup  --plotType se --yAxisLabel "mean CPM" --regionsLabel '' --whatToShow 'heatmap and colorbar'

		'''





rule getMappingStats:
	input:
		bams = "mappings/readMapping/{techname}_{capDesign}_{sizeFrac}_{barcodes}.bam",
		fastqs = LR_FASTQDIR + "{techname}_{capDesign}_{sizeFrac}_{barcodes}.fastq.gz"
	output: 
		basic=config["STATSDATADIR"] + "tmp/{techname}_{capDesign}_{sizeFrac}.{barcodes}.mapping.stats.tsv",
		spikeIns=config["STATSDATADIR"] + "tmp/{techname}_{capDesign}_{sizeFrac}.{barcodes}.mapping.spikeIns.stats.tsv"
	shell:
		'''
uuidTmpOutB=$(uuidgen)
uuidTmpOutS=$(uuidgen)
totalReads=$(zcat {input.fastqs} | fastq2tsv.pl | wc -l)
mappedReads=$(samtools view  -F 4 {input.bams}|cut -f1|sort -T {config[TMPDIR]} |uniq|wc -l)
erccMappedReads=$(samtools view -F 4 {input.bams}|cut -f3| tgrep ERCC- | wc -l)
sirvMappedReads=$(samtools view -F 4 {input.bams}|cut -f3| tgrep SIRV | wc -l)
echo -e "{wildcards.techname}\t{wildcards.capDesign}\t{wildcards.sizeFrac}\t{wildcards.barcodes}\t$totalReads\t$mappedReads" | awk '{{print $0"\t"$6/$5}}' > {config[TMPDIR]}/$uuidTmpOutB
mv {config[TMPDIR]}/$uuidTmpOutB {output.basic}
echo -e "{wildcards.techname}\t{wildcards.capDesign}\t{wildcards.sizeFrac}\t{wildcards.barcodes}\t$totalReads\t$erccMappedReads\t$sirvMappedReads" | awk '{{print $0"\t"$6/$5"\t"$7/$5}}' > {config[TMPDIR]}/$uuidTmpOutS
mv {config[TMPDIR]}/$uuidTmpOutS {output.spikeIns}

		'''

rule aggMappingStats:
	input: lambda wildcards: expand(config["STATSDATADIR"] + "tmp/{techname}_{capDesign}_{sizeFrac}.{barcodes}.mapping.stats.tsv",filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES)
	output: config["STATSDATADIR"] + "all.basic.mapping.stats.tsv"
	shell:
		'''
uuidTmpOut=$(uuidgen)
echo -e "seqTech\tcapDesign\tsizeFrac\ttissue\ttotalReads\tmappedReads\tpercentMappedReads" > {config[TMPDIR]}/$uuidTmpOut
cat {input} | sort -T {config[TMPDIR]}  >> {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}

		'''

rule plotMappingStats:
	input: config["STATSDATADIR"] + "all.basic.mapping.stats.tsv"
	output: returnPlotFilenames(config["PLOTSDIR"] + "lrMapping.basic.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{barcodes}.lrMapping.basic.stats")
	params:
		filterDat=lambda wildcards: multi_figures(wildcards.capDesign, wildcards.sizeFrac, wildcards.barcodes, wildcards.techname)
	conda: "envs/R_env.yml"
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

plotBase <- \\"p <- ggplot(data=dat, aes(x=1, y=percentMappedReads, fill=sizeFrac)) +
geom_bar(width=0.75,stat='identity', position=position_dodge(width=0.9)) +
scale_fill_manual(values={sizeFrac_Rpalette}) +
geom_hline(aes(yintercept=1), linetype='dashed', alpha=0.7, size=lineSize)+
geom_text(aes(group=sizeFrac, y=0.01, label = paste(sep='',percent(percentMappedReads),'\\n','(',comma(mappedReads),')')), angle=90, size=geom_textSize, hjust=0, vjust=0.5, position = position_dodge(width=0.9)) +
scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
xlab ('') +
{GGPLOT_PUB_QUALITY} + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) + \\"

{params.filterDat[facetPlotSetup]}


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
	input: lambda wildcards: expand(config["STATSDATADIR"] + "tmp/{techname}_{capDesign}_{sizeFrac}.{barcodes}.mapping.spikeIns.stats.tsv",filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES)
	output: config["STATSDATADIR"] + "all.spikeIns.mapping.stats.tsv"
	shell:
		'''
uuidTmpOut=$(uuidgen)
echo -e "seqTech\tcapDesign\tsizeFrac\ttissue\tcategory\tcount\tpercent" > {config[TMPDIR]}/$uuidTmpOut
cat {input} | awk '{{print $1"\\t"$2"\\t"$3"\\t"$4"\\tSIRVs\\t"$7"\\t"$9"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\tERCCs\\t"$6"\\t"$8}}' | sort -T {config[TMPDIR]}  >> {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}

		'''

rule plotSpikeInsMappingStats:
	input: config["STATSDATADIR"] + "all.spikeIns.mapping.stats.tsv"
	output: returnPlotFilenames(config["PLOTSDIR"] + "lrMapping.spikeIns.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{barcodes}.lrMapping.spikeIns.stats")
	params:
		filterDat=lambda wildcards: multi_figures(wildcards.capDesign, wildcards.sizeFrac, wildcards.barcodes, wildcards.techname)
	conda: "envs/R_env.yml"
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

maxY <- max(dat\$percent)
plotBase <- \\"p <- ggplot(data=dat, aes(x=1, y=percent, fill=category)) +
geom_bar(stat='identity', position=position_dodge()) +
geom_text(position = position_dodge(width = 0.9), size=geom_textSize, aes(y = 0, label = paste(sep='',percent(percent),'\\n','(',comma(count),')')), hjust = 0, vjust = 0.5, angle=90) +
scale_fill_manual(values=c('ERCCs' = '#e4b5ff', 'SIRVs' = '#5edba9')) +
ylab('% reads mapped on\\nspike-in sequences') +
xlab('') +
guides(fill = guide_legend(title='Spike-in set'))+
scale_y_continuous(labels=percent)+
expand_limits(y=c(0,maxY))+
{params.filterDat[hideXaxisLabels]}
{GGPLOT_PUB_QUALITY} + 
theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
\\"


{params.filterDat[facetPlotSetup]}

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




rule checkOnlyOneHit:
	input: "mappings/readMapping/{techname}_{capDesign}_{sizeFrac}_{barcodes}.bam"
	output: "mappings/readMapping/qc/{techname}_{capDesign}_{sizeFrac}_{barcodes}.bam.dupl.txt"
	shell:
		'''
uuidTmpOut=$(uuidgen)
samtools view {input} | cut -f1 | sort -T {config[TMPDIR]} | uniq -dc > {config[TMPDIR]}/$uuidTmpOut
count=$(cat {config[TMPDIR]}/$uuidTmpOut | wc -l)
if [ $count -gt 0 ]; then echo "$count duplicate read IDs found"; mv {input} {input}.dup.bkp; exit 1; fi
mv {config[TMPDIR]}/$uuidTmpOut {output}

		'''


rule readBamToBed:
	input: "mappings/readMapping/{techname}_{capDesign}_{sizeFrac}_{barcodes}.bam"
	output: "mappings/readBamToBed/{techname}_{capDesign}_{sizeFrac}_{barcodes}.bed.gz"
	conda: "envs/xtools_env.yml"
	shell:
		'''
uuidTmpOut=$(uuidgen)
#remove mappings with ultra-short exons after bamtobed

bedtools bamtobed -i {input} -bed12 | perl -ne '$line=$_; @line=split ("\\t", $line); @blockSizes=split(",", $line[10]); $allExonsOK=1; foreach $block (@blockSizes){{if ($block<2){{$allExonsOK=0; last;}}}}; if ($allExonsOK==1){{print $line}}'| sort -T {config[TMPDIR]}  -k1,1 -k2,2n -k3,3n  | gzip > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}

		'''

rule readBedToGff:
	input: "mappings/readBamToBed/{techname}_{capDesign}_{sizeFrac}_{barcodes}.bed.gz"
	output: "mappings/readBedToGff/{techname}_{capDesign}_{sizeFrac}_{barcodes}.gff.gz"
	shell:
		'''
uuidTmpOut=$(uuidgen)
zcat {input} | bed12togff | sort -T {config[TMPDIR]}  -k1,1 -k4,4n -k5,5n  | gzip > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
		'''

rule getReadBiotypeClassification:
	input: 
		reads="mappings/readMapping/{techname}_{capDesign}_{sizeFrac}_{barcodes}.bam",
		ann="annotations/simplified/{capDesign}.gencode.collapsed.simplified_biotypes.gtf"
	output: "mappings/readMapping/reads2biotypes/{techname}_{capDesign}_{sizeFrac}_{barcodes}.reads2biotypes.tsv.gz"
	conda: "envs/xtools_env.yml"
	shell:
		'''
uuidTmpOut=$(uuidgen)


bedtools intersect -split -wao -bed -a {input.reads} -b {input.ann} |perl -lane '$line=$_; $gid="NA"; $gt="nonExonic"; if($line=~/gene_id \"(\S+)\";/){{$gid=$1}}; if ($line=~/gene_type \"(\S+)\";/){{$gt=$1}}; print "$F[3]\\t$gid\\t$gt\\t$F[-1]"'|cut -f1,3|sort -T {config[TMPDIR]} |uniq | gzip > {config[TMPDIR]}/$uuidTmpOut.2
mv  {config[TMPDIR]}/$uuidTmpOut.2 {output}
		'''

rule getReadToBiotypeBreakdownStats:
	input: "mappings/readMapping/reads2biotypes/{techname}_{capDesign}_{sizeFrac}_{barcodes}.reads2biotypes.tsv.gz"
	output: config["STATSDATADIR"] + "tmp/{techname}_{capDesign}_{sizeFrac}_{barcodes}.readToBiotypeBreakdown.stats.tsv"
	shell:
		'''
totalPairs=$(zcat {input} | wc -l)
zcat {input} | cut -f2| sort -T {config[TMPDIR]} |uniq -c |ssv2tsv | awk -v t={wildcards.techname} -v c={wildcards.capDesign} -v si={wildcards.sizeFrac} -v b={wildcards.barcodes} -v tp=$totalPairs '{{print t\"\\t\"c\"\\t\"si\"\\t\"b\"\\t\"$2"\\t"$1"\\t"$1/tp}}' > {output}

		'''


rule aggReadToBiotypeBreakdownStats:
	input: expand(config["STATSDATADIR"] + "tmp/{techname}_{capDesign}_{sizeFrac}_{barcodes}.readToBiotypeBreakdown.stats.tsv", filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES)
	output: config["STATSDATADIR"] + "all.readToBiotypeBreakdown.stats.tsv"
	shell:
		'''
uuidTmpOut=$(uuidgen)
echo -e "seqTech\tcapDesign\tsizeFrac\ttissue\tbiotype\treadOverlapsCount\treadOverlapsPercent" > {config[TMPDIR]}/$uuidTmpOut
cat {input} | sort -T {config[TMPDIR]}  >> {config[TMPDIR]}/$uuidTmpOut
mv  {config[TMPDIR]}/$uuidTmpOut {output}

		'''

rule plotReadToBiotypeBreakdownStats:
	input: config["STATSDATADIR"] + "all.readToBiotypeBreakdown.stats.tsv"
	output: returnPlotFilenames(config["PLOTSDIR"] + "readToBiotypeBreakdown.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{barcodes}.readToBiotypeBreakdown.stats")
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
library(data.table)

dat <- read.table('{input}', header=T, as.is=T, sep='\\t')
{params.filterDat[technameFilterString]}
{params.filterDat[capDesignFilterString]}

{params.filterDat[sizeFracFilterString]}
{params.filterDat[tissueFilterString]}
{params.filterDat[substSeqTechString]}
{params.filterDat[substTissueString]}
{params.filterDat[graphDimensions]}

dat\$biotype=factor(dat\$biotype, levels=names({simpleBiotypes_Rpalette}), ordered=TRUE)  #otherwise the manual scale is not ordered correctly and "drop=FALSE" (include categories in scale that are absent from data frame) is ignored 

plotBase <- \\"p <- ggplot(data=dat, aes(x=1, y=readOverlapsPercent, fill=biotype)) +
geom_bar(stat='identity') + scale_fill_manual(values={simpleBiotypes_Rpalette}, drop = FALSE) + ylab('% read overlaps') + xlab('') + guides(fill = guide_legend(title='Region/biotype')) +
scale_y_continuous(labels = scales::percent) + coord_cartesian(ylim=c(0,1)) +

{params.filterDat[hideXaxisLabels]}
{GGPLOT_PUB_QUALITY} + 
theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
\\"


{params.filterDat[facetPlotSetup]}

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





