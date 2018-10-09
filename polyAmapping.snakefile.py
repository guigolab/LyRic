
rule getPolyAsitesTestPolyAmapping:
	input:
		reads="mappings/readMapping/{techname}Corr{corrLevel}_{capDesign}_allFracs_allTissues.bam",
		genome=lambda wildcards: config["GENOMESDIR"] + CAPDESIGNTOGENOME[wildcards.capDesign] + ".fa"
	params:
		minAcontent=0.8
	output: temp("mappings/polyAmapping/calibration/" + "{techname}Corr{corrLevel}_{capDesign}.polyA.minClipped.{minA}.bed")
	shell:
		'''
samtools view {input.reads} | fgrep -v ERCC | fgrep -v SIRVome_isoforms |samToPolyA.pl --minClipped={wildcards.minA} --minAcontent={params.minAcontent} --discardInternallyPrimed --minUpMisPrimeAlength=10 --genomeFasta={input.genome} - | sortbed | bedtools merge -s -d 5 -c 4,6 -o distinct -i stdin | awk '{{print $1"\\t"$2"\\t"$3"\\t"$4"\\t0\\t"$5}}'| perl -F"\\t" -lane 'if($F[5] eq "+"){{$F[1]=$F[2]-1}}elsif($F[5] eq "-"){{$F[2]=$F[1]+1}}else{{die}} print join("\\t",@F);'|sortbed > {output}

		'''

rule getPolyANonPolyAsitesTestPolyAmapping:
	input:
		reads=lambda wildcards: expand("mappings/" + "strandGffs/{techname}Corr{corrLevel}_{capDesign}_{barcodes}.stranded.gff", filtered_product, techname=wildcards.techname, corrLevel=wildcards.corrLevel, capDesign=wildcards.capDesign, barcodes=BARCODES),
		inPolyA="mappings/polyAmapping/calibration/" + "{techname}Corr{corrLevel}_{capDesign}.polyA.minClipped.{minA}.bed"
	params:
		minAcontent=0.8
	output:
		outPolyA=temp("mappings/polyAmapping/calibration/3pends/" + "{techname}Corr{corrLevel}_{capDesign}.polyA.minClipped.{minA}.polyA.3pEnds.bed"),
		outNonPolyA=temp("mappings/polyAmapping/calibration/3pends/" + "{techname}Corr{corrLevel}_{capDesign}.polyA.minClipped.{minA}.nonPolyA.3pEnds.bed")
	shell:
		'''
cat {input.inPolyA} | cut -f4 | sed 's/,/\\n/g' | sort | uniq > $TMPDIR/{wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}.polyA.minClipped.{wildcards.minA}.list

#non polyA:
cat {input.reads}| fgrep -v ERCC | fgrep -v SIRVome_isoforms| fgrep -v -w -f $TMPDIR/{wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}.polyA.minClipped.{wildcards.minA}.list | gff2bed_full.pl - | extractTranscriptEndsFromBed12.pl 3 |sortbed| bedtools merge -s -d 5 -c 4,6 -o distinct -i stdin | awk '{{print $1"\\t"$2"\\t"$3"\\t"$4"\\t0\\t"$5}}'| perl -F"\\t" -lane 'if($F[5] eq "+"){{$F[1]=$F[2]-1}}elsif($F[5] eq "-"){{$F[2]=$F[1]+1}}else{{die}} print join("\\t",@F);' > {output.outNonPolyA}

#polyA:
cat {input.reads}| fgrep -v ERCC | fgrep -v SIRVome_isoforms| fgrep -w -f $TMPDIR/{wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}.polyA.minClipped.{wildcards.minA}.list | gff2bed_full.pl - | extractTranscriptEndsFromBed12.pl 3 |sortbed| bedtools merge -s -d 5 -c 4,6 -o distinct -i stdin | awk '{{print $1"\\t"$2"\\t"$3"\\t"$4"\\t0\\t"$5}}'| perl -F"\\t" -lane 'if($F[5] eq "+"){{$F[1]=$F[2]-1}}elsif($F[5] eq "-"){{$F[2]=$F[1]+1}}else{{die}} print join("\\t",@F);' > {output.outPolyA}


		'''

rule compareToPASTestPolyAmapping:
	input:
		tpend="mappings/polyAmapping/calibration/3pends/" + "{techname}Corr{corrLevel}_{capDesign}.polyA.minClipped.{minA}.{tpend}.3pEnds.bed",
		PAS=lambda wildcards: CAPDESIGNTOPAS[wildcards.capDesign]
	params: genome = lambda wildcards: config["GENOMESDIR"] + CAPDESIGNTOGENOME[wildcards.capDesign] + ".genome"
	output: temp("mappings/polyAmapping/calibration/3pends/vsPAS/" + "{techname}Corr{corrLevel}_{capDesign}.polyA.minClipped.{minA}.{tpend}.3pEnds.bedtsv")

	shell:
		'''
cat {input.tpend} | sortbed | bedtools slop -s -l 50 -r -10 -i stdin -g {params.genome} | bedtools intersect -u -s -a stdin -b {input.PAS} > {output}

		'''

rule getPrecisionRecallTestPolyAmapping:
	input:
		pAvsPAS="mappings/" + "polyAmapping/calibration/3pends/vsPAS/" + "{techname}Corr{corrLevel}_{capDesign}.polyA.minClipped.{minA}.polyA.3pEnds.bedtsv",
		nonpAvsPAS="mappings/" + "polyAmapping/calibration/3pends/vsPAS/" + "{techname}Corr{corrLevel}_{capDesign}.polyA.minClipped.{minA}.nonPolyA.3pEnds.bedtsv",
		pA="mappings/" + "polyAmapping/calibration/3pends/" + "{techname}Corr{corrLevel}_{capDesign}.polyA.minClipped.{minA}.polyA.3pEnds.bed",
		nonpA="mappings/" + "polyAmapping/calibration/3pends/" + "{techname}Corr{corrLevel}_{capDesign}.polyA.minClipped.{minA}.nonPolyA.3pEnds.bed"

	output:temp(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{minA}.polyA.vs.PAS.precisionRecall.stats.tsv")
	shell:
		'''
closePolyA=$(cat {input.pAvsPAS} | cut -f4 | sort|uniq | wc -l)
totalPolyA=$(cat {input.pA} | cut -f4 | sort|uniq|wc -l)
let farPolyA=$totalPolyA-$closePolyA || true # "true" serves to avoid "let" exiting with status > 0 when its output is = 0
closeNonPolyA=$(cat {input.nonpAvsPAS} | cut -f4 | sort|uniq | wc -l)
totalNonPolyA=$(cat {input.nonpA} | cut -f4 | sort|uniq|wc -l)
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
cat {input} | sort > {output}
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
	output: "mappings/" + "polyAmapping/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.polyAsites.bed.gz"
	shell:
		'''
samtools view {input.reads} | samToPolyA.pl --minClipped=10 --minAcontent={params.minAcontent}  --discardInternallyPrimed --minUpMisPrimeAlength=10 --genomeFasta={input.genome} - |sortbed |gzip > {output}
		'''

rule removePolyAERCCs:
	input: "mappings/" + "polyAmapping/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.polyAsites.bed.gz"
	output: temp("mappings/" + "removePolyAERCCs/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.polyAsitesNoErcc.tmp.bed")
	shell:
		'''
zcat {input} | fgrep -v ERCC > {output}
		'''

rule getPolyAreadsList:
	input: "mappings/" + "removePolyAERCCs/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.polyAsitesNoErcc.bed"
	output: "mappings/" + "getPolyAreadsList/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.polyAreads.list"
	shell:
		'''
cat {input} | cut -f4 | sort|uniq > {output}
		'''

rule getPolyAreadsStats:
	input:
		mappedReads= "mappings/" + "readMapping/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bam",
		polyAreads = "mappings/" + "getPolyAreadsList/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.polyAreads.list"
	output: temp(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.polyAreads.stats.tsv")
	# wildcard_constraints:
	# 	sizeFrac='[^(allFracs)][\S]+', #to avoid ambiguity with downstream merging rules
	# 	barcodes='[^(allTissues)][\S]+',

	# wildcard_constraints:
	# 	techname='[^(allTechs)][\S]+', #to avoid ambiguity with downstream merging rules
	# 	corrLevel='[^(allCors)][\S]+',
	# 	capDesign='[^(allCapDesigns)][\S]+',
	shell:
		'''
mapped=$(samtools view -F4 {input.mappedReads} |cut -f1|sort|uniq|wc -l)
polyA=$(cat {input.polyAreads} | wc -l)
echo -e "{wildcards.techname}Corr{wildcards.corrLevel}\t{wildcards.capDesign}\t{wildcards.sizeFrac}\t{wildcards.barcodes}\t$mapped\t$polyA" | awk '{{print $0"\t"$6/$5}}'|sed 's/{wildcards.capDesign}_//' > {output}

		'''



rule aggPolyAreadsAllFracsAllTissuesStats:
	input: lambda wildcards: expand(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.polyAreads.stats.tsv", filtered_product, techname=wildcards.techname, corrLevel=wildcards.corrLevel, capDesign=wildcards.capDesign, sizeFrac=SIZEFRACS, barcodes=BARCODES)
	output: temp(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_allFracs_allTissues.polyAreads.agg.stats.tsv")
	shell:
		'''
totalMapped=$(cat {input} | cut -f5 | sum.sh)
totalPolyA=$(cat {input} | cut -f6 | sum.sh)
echo -e "{wildcards.techname}Corr{wildcards.corrLevel}\t{wildcards.capDesign}\tallFracs\tallTissues\t$totalMapped\t$totalPolyA" | awk '{{print $0"\t"$6/$5}}' > {output}
		'''

rule aggPolyAreadsAllTissuesStats:
	input: lambda wildcards: expand(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.polyAreads.stats.tsv", filtered_product, techname=wildcards.techname, corrLevel=wildcards.corrLevel, capDesign=wildcards.capDesign, sizeFrac=wildcards.sizeFrac, barcodes=BARCODES)
	output: temp(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_allTissues.polyAreads.agg.stats.tsv")
	wildcard_constraints:
		sizeFrac='[^(allFracs)][^_][\S]+', #to avoid ambiguity with downstream merging rules
	shell:
		'''
totalMapped=$(cat {input} | cut -f5 | sum.sh)
totalPolyA=$(cat {input} | cut -f6 | sum.sh)
echo -e "{wildcards.techname}Corr{wildcards.corrLevel}\t{wildcards.capDesign}\t{wildcards.sizeFrac}\tallTissues\t$totalMapped\t$totalPolyA" | awk '{{print $0"\t"$6/$5}}' > {output}
		'''


rule aggPolyAreadsAllFracsStats:
	input: lambda wildcards: expand(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.polyAreads.stats.tsv", filtered_product, techname=wildcards.techname, corrLevel=wildcards.corrLevel, capDesign=wildcards.capDesign, sizeFrac=SIZEFRACS, barcodes=wildcards.barcodes)
	output: temp(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_allFracs_{barcodes}.polyAreads.agg.stats.tsv")
	wildcard_constraints:
		barcodes='[^(allTissues)][\S]+'
	shell:
		'''
totalMapped=$(cat {input} | cut -f5 | sum.sh)
totalPolyA=$(cat {input} | cut -f6 | sum.sh)
echo -e "{wildcards.techname}Corr{wildcards.corrLevel}\t{wildcards.capDesign}\tallFracs\t{wildcards.barcodes}\t$totalMapped\t$totalPolyA" | awk '{{print $0"\t"$6/$5}}' |sed 's/{wildcards.capDesign}_//' > {output}
		'''

rule aggPolyAreadsStats:
	#input: config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.polyAreads.stats.tsv"
#	input: lambda wildcards: expand(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.polyAreads.agg.stats.tsv",filtered_product_merge, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACSpluSMERGED, barcodes=BARCODESpluSMERGED)
	input: lambda wildcards: expand(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_allFracs_{barcodes}.polyAreads.agg.stats.tsv", filtered_product, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, barcodes=BARCODES),
		lambda wildcards: expand(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_allTissues.polyAreads.agg.stats.tsv", filtered_product, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS),
		lambda wildcards: expand(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_allFracs_allTissues.polyAreads.agg.stats.tsv", techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS)
	output: config["STATSDATADIR"] + "all.polyAreads.stats.tsv"
	shell:
		'''
echo -e "seqTech\tcorrectionLevel\tcapDesign\tsizeFrac\ttissue\tcategory\tcount\tpercent" > {output}
cat {input} | awk '{{print $1"\\t"$2"\\t"$3"\\t"$4"\\tnonPolyA\\t"$5-$6"\\t"($5-$6)/$5"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\tpolyA\\t"$6"\\t"$6/$5}}'| sed 's/Corr0/\tNo/' | sed 's/Corr{lastK}/\tYes/' | sort >> {output}
		'''



# rule aggPolyAreadsStats:
# 	input: lambda wildcards: expand(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.polyAreads.stats.tsv",filtered_product_merge, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACSpluSMERGED, barcodes=BARCODESpluSMERGED)
# 	output: config["STATSDATADIR"] + "all.polyAreads.stats.tsv"
# 	shell:
# 		'''
# echo -e "seqTech\tcorrectionLevel\tcapDesign\tsizeFrac\ttissue\tcategory\tcount" > {output}
# cat {input} | awk '{{print $1"\\t"$2"\\t"$3"\\t"$4"\\tnonPolyA\\t"$5-$6"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\tpolyA\\t"$6}}' | sed 's/Corr0/\tNo/' | sed 's/Corr{lastK}/\tYes/' | sort -r >> {output}
# 		'''

rule plotAllPolyAreadsStats:
	input: config["STATSDATADIR"] + "all.polyAreads.stats.tsv"
	output: config["PLOTSDIR"] + "{capDesign}_{byFrac}_{byTissue}.polyAreads.stats.{ext}"
	params:
		filterDat=lambda wildcards: merge_figures_params(wildcards.capDesign, wildcards.byFrac, wildcards.byTissue)
	shell:
		'''
echo "library(ggplot2)
library(plyr)
library(scales)
dat <- read.table('{input}', header=T, as.is=T, sep='\\t')
dat\$seqTech <- gsub(':', '\\n', dat\$seqTech)
{params.filterDat}
horizCats <- length(unique(dat\$correctionLevel)) * length(unique(dat\$capDesign)) * length(unique(dat\$tissue))
vertCats <- length(unique(dat\$seqTech)) * length(unique(dat\$sizeFrac))
plotWidth = horizCats + 1
plotHeight = vertCats + 1
geom_textSize=0.1 * plotWidth * plotHeight
plotWidth
plotHeight
geom_textSize
ggplot(data=dat, aes(x=factor(correctionLevel), y=count, fill=category)) +
geom_bar(stat='identity') + scale_fill_manual(values=c('polyA' = '#c8e09e', 'nonPolyA' = '#e7a198')) + facet_grid( seqTech + sizeFrac ~ capDesign + tissue)+ ylab('# mapped reads') + xlab('Error correction') + guides(fill = guide_legend(title='Category'))+ scale_y_continuous(labels=scientific)+
geom_text(position = 'stack', size=geom_textSize, aes(x = factor(correctionLevel), y = count, ymax=count, label = paste(sep='',percent(round(percent, digits=2)),' / ','(',comma(count),')'), hjust = 0.5, vjust = 1))+
{GGPLOT_PUB_QUALITY}
ggsave('{output}', width=plotWidth, height=plotHeight)
" > {output}.r
cat {output}.r | R --slave

		'''


rule clusterPolyAsites:
	input: "mappings/" + "removePolyAERCCs/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.polyAsitesNoErcc.bed"
	output: "mappings/" + "clusterPolyAsites/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.polyAsites.clusters.bed"
	shell:
		'''
cat {input} | bedtools merge -s -d 5 -c 4,6 -o distinct -i stdin | awk '{{print $1"\t"$2"\t"$3"\t"$4"\t0\t"$5}}'| perl -F"\\t" -lane 'if($F[5] eq "+"){{$F[1]=$F[2]-1}}elsif($F[5] eq "-"){{$F[2]=$F[1]+1}}else{{die}} print join("\t",@F);'|sortbed > {output}
		'''
rule makePolyABigWigs:
	input:
		sites = "mappings/" + "removePolyAERCCs/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.polyAsitesNoErcc.bed",
		genome = lambda wildcards: config["GENOMESDIR"] + CAPDESIGNTOGENOME[wildcards.capDesign] + ".genome"
	output: "mappings/" + "makePolyABigWigs/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.polyAsitesNoErcc.{strand}.bw"
	shell:
		'''
bedtools genomecov -strand {wildcards.strand} -split -bg -i {input.sites} -g {input.genome} > {output}.bedgraph
bedGraphToBigWig {output}.bedgraph {input.genome} {output}
rm -f {output}.bedgraph
		'''
