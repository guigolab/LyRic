
rule getPolyAsitesTestPolyAmapping:
	input:
		reads="mappings/readMapping/{techname}Corr{corrLevel}_{capDesign}_0+_allTissues.bam",
		genome=lambda wildcards: config["GENOMESDIR"] + CAPDESIGNTOGENOME[wildcards.capDesign] + ".fa"
	params:
		minAcontent=0.8
	output: temp("mappings/polyAmapping/calibration/{techname}Corr{corrLevel}_{capDesign}.polyA.minClipped.{minA}.bed")
	shell:
		'''
uuidTmpOut=$(uuidgen)
samtools view {input.reads} | fgrep -v ERCC | fgrep -v SIRVome_isoforms |samToPolyA.pl --minClipped={wildcards.minA} --minAcontent={params.minAcontent} --discardInternallyPrimed --minUpMisPrimeAlength=10 --genomeFasta={input.genome} - | sort -T {config[TMPDIR]}  -k1,1 -k2,2n -k3,3n  | bedtools merge -s -d 5 -c 4,6 -o distinct -i stdin | awk '{{print $1"\\t"$2"\\t"$3"\\t"$4"\\t0\\t"$5}}'| perl -F"\\t" -lane 'if($F[5] eq "+"){{$F[1]=$F[2]-1}}elsif($F[5] eq "-"){{$F[2]=$F[1]+1}}else{{die}} print join("\\t",@F);'|sort -T {config[TMPDIR]}  -k1,1 -k2,2n -k3,3n  > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}

		'''

rule getPolyANonPolyAsitesTestPolyAmapping:
	input:
		reads=lambda wildcards: expand("mappings/strandGffs/{techname}Corr{corrLevel}_{capDesign}_{barcodes}.stranded.gff", filtered_product, techname=wildcards.techname, corrLevel=wildcards.corrLevel, capDesign=wildcards.capDesign, barcodes=BARCODES),
		inPolyA="mappings/polyAmapping/calibration/{techname}Corr{corrLevel}_{capDesign}.polyA.minClipped.{minA}.bed"
	params:
		minAcontent=0.8
	output:
		outPolyA=temp("mappings/polyAmapping/calibration/3pends/{techname}Corr{corrLevel}_{capDesign}.polyA.minClipped.{minA}.polyA.3pEnds.bed"),
		outNonPolyA=temp("mappings/polyAmapping/calibration/3pends/{techname}Corr{corrLevel}_{capDesign}.polyA.minClipped.{minA}.nonPolyA.3pEnds.bed")
	shell:
		'''
uuid=$(uuidgen)
uuidTmpOutN=$(uuidgen)
uuidTmpOutP=$(uuidgen)

cat {input.inPolyA} | cut -f4 | sed 's/,/\\n/g' | sort -T {config[TMPDIR]}  | uniq > {config[TMPDIR]}/$uuid.list

#non polyA:
cat {input.reads}| fgrep -v ERCC | fgrep -v SIRVome_isoforms| fgrep -v -w -f {config[TMPDIR]}/$uuid.list | gff2bed_full.pl - | extractTranscriptEndsFromBed12.pl 3 |sort -T {config[TMPDIR]}  -k1,1 -k2,2n -k3,3n | bedtools merge -s -d 5 -c 4,6 -o distinct -i stdin | awk '{{print $1"\\t"$2"\\t"$3"\\t"$4"\\t0\\t"$5}}'| perl -F"\\t" -lane 'if($F[5] eq "+"){{$F[1]=$F[2]-1}}elsif($F[5] eq "-"){{$F[2]=$F[1]+1}}else{{die}} print join("\\t",@F);' > {config[TMPDIR]}/$uuidTmpOutN
mv {config[TMPDIR]}/$uuidTmpOutN {output.outNonPolyA}

#polyA:
cat {input.reads}| fgrep -v ERCC | fgrep -v SIRVome_isoforms| fgrep -w -f {config[TMPDIR]}/$uuid.list | gff2bed_full.pl - | extractTranscriptEndsFromBed12.pl 3 |sort -T {config[TMPDIR]}  -k1,1 -k2,2n -k3,3n | bedtools merge -s -d 5 -c 4,6 -o distinct -i stdin | awk '{{print $1"\\t"$2"\\t"$3"\\t"$4"\\t0\\t"$5}}'| perl -F"\\t" -lane 'if($F[5] eq "+"){{$F[1]=$F[2]-1}}elsif($F[5] eq "-"){{$F[2]=$F[1]+1}}else{{die}} print join("\\t",@F);' > {config[TMPDIR]}/$uuidTmpOutP
mv {config[TMPDIR]}/$uuidTmpOutP {output.outPolyA}


		'''

rule compareToPASTestPolyAmapping:
	input:
		tpend="mappings/polyAmapping/calibration/3pends/{techname}Corr{corrLevel}_{capDesign}.polyA.minClipped.{minA}.{tpend}.3pEnds.bed",
		PAS=lambda wildcards: CAPDESIGNTOPAS[wildcards.capDesign],
		genome = lambda wildcards: config["GENOMESDIR"] + CAPDESIGNTOGENOME[wildcards.capDesign] + ".genome"
	output: temp("mappings/polyAmapping/calibration/3pends/vsPAS/{techname}Corr{corrLevel}_{capDesign}.polyA.minClipped.{minA}.{tpend}.3pEnds.bedtsv")

	shell:
		'''
uuidTmpOut=$(uuidgen)
cat {input.tpend} | sort -T {config[TMPDIR]}  -k1,1 -k2,2n -k3,3n  | bedtools slop -s -l 50 -r -10 -i stdin -g {input.genome} | bedtools intersect -u -s -a stdin -b {input.PAS} > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}

		'''

rule getPrecisionRecallTestPolyAmapping:
	input:
		pAvsPAS="mappings/polyAmapping/calibration/3pends/vsPAS/{techname}Corr{corrLevel}_{capDesign}.polyA.minClipped.{minA}.polyA.3pEnds.bedtsv",
		nonpAvsPAS="mappings/polyAmapping/calibration/3pends/vsPAS/{techname}Corr{corrLevel}_{capDesign}.polyA.minClipped.{minA}.nonPolyA.3pEnds.bedtsv",
		pA="mappings/polyAmapping/calibration/3pends/{techname}Corr{corrLevel}_{capDesign}.polyA.minClipped.{minA}.polyA.3pEnds.bed",
		nonpA="mappings/polyAmapping/calibration/3pends/{techname}Corr{corrLevel}_{capDesign}.polyA.minClipped.{minA}.nonPolyA.3pEnds.bed"

	output:temp(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{minA}.polyA.vs.PAS.precisionRecall.stats.tsv")
	shell:
		'''
closePolyA=$(cat {input.pAvsPAS} | cut -f4 | sort -T {config[TMPDIR]} |uniq | wc -l)
totalPolyA=$(cat {input.pA} | cut -f4 | sort -T {config[TMPDIR]} |uniq|wc -l)
let farPolyA=$totalPolyA-$closePolyA || true # "true" serves to avoid "let" exiting with status > 0 when its output is = 0
closeNonPolyA=$(cat {input.nonpAvsPAS} | cut -f4 | sort -T {config[TMPDIR]} |uniq | wc -l)
totalNonPolyA=$(cat {input.nonpA} | cut -f4 | sort -T {config[TMPDIR]} |uniq|wc -l)
let farNonPolyA=$totalNonPolyA-$closeNonPolyA || true
let denomPr=$closePolyA+$farPolyA || true
pr=$(echo $closePolyA $denomPr | awk '{{print $1/$2}}')
let denomSn=$closePolyA+$closeNonPolyA || true
sn=$(echo $closePolyA $denomSn | awk '{{print $1/$2}}')
snPr=$(echo $sn $pr | awk '{{print ($1+$2)/2}}')
echo -e "{wildcards.techname}Corr{wildcards.corrLevel}\t{wildcards.capDesign}\t{wildcards.minA}\tPr\t$pr
{wildcards.techname}Corr{wildcards.corrLevel}\t{wildcards.capDesign}\t{wildcards.minA}\tSn\t$sn
{wildcards.techname}Corr{wildcards.corrLevel}\t{wildcards.capDesign}\t{wildcards.minA}\tSnPr\t$snPr" | sed 's/Corr0/\tNo/' | sed 's/Corr{lastK}/\tYes/' > {output}
		'''

rule aggPrecisionRecallTestPolyAmapping:
	input: lambda wildcards: expand(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{minA}.polyA.vs.PAS.precisionRecall.stats.tsv", techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, minA=minPolyAlength)
	output: config["STATSDATADIR"] + "all.polyA.vs.PAS.precisionRecall.stats.tsv"
	shell:
		'''
cat {input} | sort -T {config[TMPDIR]}  > {output}
		'''

rule plotPrecisionRecallTestPolyAmapping:
	input: config["STATSDATADIR"] + "all.polyA.vs.PAS.precisionRecall.stats.tsv"
	output: config["PLOTSDIR"] + "all.polyA.vs.PAS.precisionRecall.stats.{ext}"
	shell:
		'''
echo "library(ggplot2)
dat <- read.table('{input}', header=F, as.is=T, sep='\\t')
colnames(dat)<-c('techname', 'correction', 'capDesign', 'minClipped','measure','value')
palette <- c('Sn' = '#b3ccff', 'Pr' = '#ffb399', 'SnPr' = '#009933')
ggplot(data=dat, aes(x=minClipped, y=value, group=measure, color=measure)) +
geom_point() +
scale_color_manual(values=palette) +
geom_line() + facet_grid( techname + capDesign ~ correction) + ylab('Pr/Sn') + xlab('minimum required A(n) length') + ylim(0, 1)
ggsave('{output}', width=7, height=8)
"  > {output}.r
cat {output}.r| R --slave

		'''



rule polyAmapping:
	input:
		reads = "mappings/readMapping/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bam",
		genome = lambda wildcards: config["GENOMESDIR"] + CAPDESIGNTOGENOME[wildcards.capDesign] + ".fa"
	params:
		minAcontent=0.8
	output: "mappings/polyAmapping/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.polyAsites.bed.gz"
	shell:
		'''
uuidTmpOut=$(uuidgen)
samtools view {input.reads} | samToPolyA.pl --minClipped=10 --minAcontent={params.minAcontent}  --discardInternallyPrimed --minUpMisPrimeAlength=10 --genomeFasta={input.genome} - |sort -T {config[TMPDIR]}  -k1,1 -k2,2n -k3,3n  |gzip > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
		'''

rule removePolyAERCCs:
	input: "mappings/polyAmapping/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.polyAsites.bed.gz"
	output: temp("mappings/removePolyAERCCs/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.polyAsitesNoErcc.tmp.bed")
	shell:
		'''
uuidTmpOut=$(uuidgen)
zcat {input} | fgrep -v ERCC > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
		'''

rule getPolyAreadsList:
	input: "mappings/removePolyAERCCs/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.polyAsitesNoErcc.bed"
	output: "mappings/getPolyAreadsList/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.polyAreads.list"
	shell:
		'''
uuidTmpOut=$(uuidgen)
cat {input} | cut -f4 | sort -T {config[TMPDIR]} |uniq > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
		'''

rule getPolyAreadsStats:
	input:
		mappedReads= "mappings/readMapping/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bam",
		polyAreads = "mappings/getPolyAreadsList/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.polyAreads.list"
	output: config["STATSDATADIR"] + "tmp/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.polyAreads.stats.tsv"
	shell:
		'''
uuidTmpOut=$(uuidgen)
mapped=$(samtools view -F4 {input.mappedReads} |cut -f1|sort -T {config[TMPDIR]} |uniq|wc -l)
polyA=$(cat {input.polyAreads} | wc -l)
echo -e "{wildcards.techname}Corr{wildcards.corrLevel}\t{wildcards.capDesign}\t{wildcards.sizeFrac}\t{wildcards.barcodes}\t$mapped\t$polyA" | awk '{{print $0"\t"$6/$5}}' > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}

		'''

rule aggPolyAreadsStats:
	input: expand(config["STATSDATADIR"] + "tmp/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.polyAreads.stats.tsv",filtered_product_merge, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODESpluSMERGED)
	output: config["STATSDATADIR"] + "all.polyAreads.stats.tsv"
	shell:
		'''
uuidTmpOut=$(uuidgen)
echo -e "seqTech\tcorrectionLevel\tcapDesign\tsizeFrac\ttissue\tcategory\tcount\tpercent" > {config[TMPDIR]}/$uuidTmpOut
cat {input} | awk '{{print $1"\\t"$2"\\t"$3"\\t"$4"\\tnonPolyA\\t"$5-$6"\\t"($5-$6)/$5"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\tpolyA\\t"$6"\\t"$6/$5}}'| sed 's/Corr0/\tNo/' | sed 's/Corr{lastK}/\tYes/' | sort -T {config[TMPDIR]}  >> {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
		'''


rule plotAllPolyAreadsStats:
	input: config["STATSDATADIR"] + "all.polyAreads.stats.tsv"
	output: returnPlotFilenames(config["PLOTSDIR"] + "polyAreads.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.polyAreads.stats")
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

plotBase <- \\"ggplot(data=dat, aes(x=factor(correctionLevel), y=count, fill=category)) +
geom_bar(stat='identity') + 
scale_fill_manual(values=c('polyA' = '#c8e09e', 'nonPolyA' = '#e7a198')) +
ylab('# mapped reads') +
xlab('{params.filterDat[6]}') +
guides(fill = guide_legend(title='Category'))+ scale_y_continuous(labels=scientific)+
geom_text(position = 'stack', size=geom_textSize, aes(x = factor(correctionLevel), y = count, label = paste(sep='',percent(round(percent, digits=2)),' / ','(',comma(count),')'), hjust = 0.5, vjust = 1))+
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


rule clusterPolyAsites:
	input: "mappings/removePolyAERCCs/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.polyAsitesNoErcc.bed"
	output: "mappings/clusterPolyAsites/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.polyAsites.clusters.bed"
	shell:
		'''
uuidTmpOut=$(uuidgen)
cat {input} | bedtools merge -s -d 5 -c 4,6 -o distinct -i stdin | awk '{{print $1"\t"$2"\t"$3"\t"$4"\t0\t"$5}}'| perl -F"\\t" -lane 'if($F[5] eq "+"){{$F[1]=$F[2]-1}}elsif($F[5] eq "-"){{$F[2]=$F[1]+1}}else{{die}} print join("\t",@F);'|sort -T {config[TMPDIR]}  -k1,1 -k2,2n -k3,3n  > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
		'''
rule makePolyABigWigs:
	input:
		sites = "mappings/removePolyAERCCs/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.polyAsitesNoErcc.bed",
		genome = lambda wildcards: config["GENOMESDIR"] + CAPDESIGNTOGENOME[wildcards.capDesign] + ".genome"
	output: "mappings/makePolyABigWigs/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.polyAsitesNoErcc.{strand}.bw"
	shell:
		'''
tmpIn=$(uuidgen)
uuidTmpOut=$(uuidgen)
cat {input.sites} | grep -P "^chr" | grep -v "chrIS" > {config[TMPDIR]}/$tmpIn

bedtools genomecov -strand {wildcards.strand} -split -bg -i {config[TMPDIR]}/$tmpIn -g {input.genome} > {config[TMPDIR]}/$uuidTmpOut.bedgraph
bedGraphToBigWig {config[TMPDIR]}/$uuidTmpOut.bedgraph {input.genome} {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}

		'''
