


rule polyAmapping:
	input:
		reads = "output/mappings/longReadMapping/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.bam",
		genome = lambda wildcards: config["GENOMESDIR"] + CAPDESIGNTOGENOME[wildcards.capDesign] + ".sorted.fa"
	params:
		minAcontent=0.8
	output: "output/mappings/polyAmapping/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.polyAsites.bed.gz"
	conda: "envs/perl_env.yml"
	shell:
		'''
uuidTmpOut=$(uuidgen)


samtools view {input.reads} | samToPolyA.pl --minClipped=10 --minAcontent={params.minAcontent}  --discardInternallyPrimed --minUpMisPrimeAlength=10 --genomeFasta={input.genome} - |sort -T {config[TMPDIR]}  -k1,1 -k2,2n -k3,3n  |gzip > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
		'''

rule removePolyAERCCs:
	input: "output/mappings/polyAmapping/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.polyAsites.bed.gz"
	output: temp("output/mappings/removePolyAERCCs/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.polyAsitesNoErcc.tmp.bed")
	shell:
		'''
uuidTmpOut=$(uuidgen)
zcat {input} | fgrep -v ERCC > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
		'''

rule getPolyAreadsList:
	input: "output/mappings/removePolyAERCCs/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.polyAsitesNoErcc.bed"
	output: "output/mappings/getPolyAreadsList/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.polyAreads.list"
	shell:
		'''
uuidTmpOut=$(uuidgen)
cat {input} | cut -f4 | sort -T {config[TMPDIR]} |uniq > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
		'''

rule getPolyAreadsStats:
	input:
		mappedReads= "output/mappings/longReadMapping/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.bam",
		polyAreads = "output/mappings/getPolyAreadsList/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.polyAreads.list"
	output: "output/statsFiles/" + "tmp/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.polyAreads.stats.tsv"
	shell:
		'''
uuidTmpOut=$(uuidgen)
mapped=$(samtools view -F4 {input.mappedReads} |cut -f1|sort -T {config[TMPDIR]} |uniq|wc -l)
polyA=$(cat {input.polyAreads} | wc -l)
echo -e "{wildcards.techname}\t{wildcards.capDesign}\t{wildcards.sizeFrac}\t{wildcards.sampleRep}\t$mapped\t$polyA" | awk '{{print $0"\t"$6/$5}}' > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}

		'''

rule aggPolyAreadsStats:
	input: expand("output/statsFiles/" + "tmp/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.polyAreads.stats.tsv",filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, sampleRep=SAMPLEREPS)
	output: "output/statsFiles/" + "all.polyAreads.stats.tsv"
	shell:
		'''
uuidTmpOut=$(uuidgen)
echo -e "seqTech\tcapDesign\tsizeFrac\tsampleRep\tcategory\tcount\tpercent" > {config[TMPDIR]}/$uuidTmpOut
cat {input} | awk '{{print $1"\\t"$2"\\t"$3"\\t"$4"\\tnonPolyA\\t"$5-$6"\\t"($5-$6)/$5"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\tpolyA\\t"$6"\\t"$6/$5}}' | sort -T {config[TMPDIR]}  >> {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
		'''


rule plotAllPolyAreadsStats:
	input: "output/statsFiles/" + "all.polyAreads.stats.tsv"
	output: returnPlotFilenames("output/plots/" + "polyAreads.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.polyAreads.stats")
	params:
		filterDat=lambda wildcards: multi_figures(wildcards.capDesign, wildcards.sizeFrac, wildcards.sampleRep, wildcards.techname)
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
library(data.table)

dat <- read.table('{input}', header=T, as.is=T, sep='\\t')
{params.filterDat[technameFilterString]}
{params.filterDat[capDesignFilterString]}

{params.filterDat[sizeFracFilterString]}
{params.filterDat[sampleRepFilterString]}
{params.filterDat[substSeqTechString]}
{params.filterDat[substSampleRepString]}
{params.filterDat[graphDimensions]}

plotBase <- \\"p <- ggplot(data=dat, aes(x=1, y=count, fill=category)) +
geom_bar(stat='identity') + 
scale_fill_manual(values=c('polyA' = '#c8e09e', 'nonPolyA' = '#e7a198')) +
ylab('# mapped reads') +
xlab('') +
guides(fill = guide_legend(title='Category'))+ scale_y_continuous(labels=scientific)+
geom_text(position = 'stack', size=geom_textSize, aes( y = count, label = paste(sep='',percent(round(percent, digits=2)),'\\n','(',comma(count),')'), hjust = 0.5, vjust = 1))+
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


rule clusterPolyAsites:
	input: "output/mappings/removePolyAERCCs/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.polyAsitesNoErcc.bed"
	output: "output/mappings/clusterPolyAsites/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.polyAsites.clusters.bed"
	conda: "envs/xtools_env.yml"
	shell:
		'''
uuidTmpOut=$(uuidgen)
cat {input} | bedtools merge -s -d 5 -c 4,6 -o distinct -i stdin | awk '{{print $1"\t"$2"\t"$3"\t"$4"\t0\t"$5}}'| perl -F"\\t" -lane 'if($F[5] eq "+"){{$F[1]=$F[2]-1}}elsif($F[5] eq "-"){{$F[2]=$F[1]+1}}else{{die}} print join("\t",@F);'|sort -T {config[TMPDIR]}  -k1,1 -k2,2n -k3,3n  > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
		'''
rule makePolyABigWigs:
	input:
		sites = "output/mappings/removePolyAERCCs/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.polyAsitesNoErcc.bed",
		genome = lambda wildcards: config["GENOMESDIR"] + CAPDESIGNTOGENOME[wildcards.capDesign] + ".sorted.genome"
	output: "output/mappings/makePolyABigWigs/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.polyAsitesNoErcc.{strand}.bw"
	conda: "envs/xtools_env.yml"
	shell:
		'''
tmpIn=$(uuidgen)
uuidTmpOut=$(uuidgen)
cat {input.sites} | grep -P "^chr" | grep -v "chrIS" > {config[TMPDIR]}/$tmpIn


bedtools genomecov -strand {wildcards.strand} -split -bg -i {config[TMPDIR]}/$tmpIn -g {input.genome} > {config[TMPDIR]}/$uuidTmpOut.bedgraph
bedGraphToBigWig {config[TMPDIR]}/$uuidTmpOut.bedgraph {input.genome} {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}

		'''
