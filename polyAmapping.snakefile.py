rule polyAmapping:
	input:
		reads = "mappings/" + "mergeSizeFracBams/{techname}_{capDesign}_{barcodes}.merged.bam",
		genome = lambda wildcards: config["GENOMESDIR"] + CAPDESIGNTOGENOME[wildcards.capDesign] + ".fa"

	output: "mappings/" + "polyAmapping/{techname}_{capDesign}_{barcodes}.polyAsites.bed.gz"
	shell:
		'''
samtools view {input.reads} | samToPolyA.pl --minClipped=20 --minAcontent=0.9 --discardInternallyPrimed --minUpMisPrimeAlength=10 --genomeFasta={input.genome} - |sortbed |gzip > {output}
		'''

rule removePolyAERCCs:
	input: "mappings/" + "polyAmapping/{techname}_{capDesign}_{barcodes}.polyAsites.bed.gz"
	output: "mappings/" + "removePolyAERCCs/{techname}_{capDesign}_{barcodes}.polyAsitesNoErcc.bed"
	shell:
		'''
zcat {input} | fgrep -v ERCC > {output}
		'''

rule getPolyAreadsList:
	input: "mappings/" + "removePolyAERCCs/{techname}_{capDesign}_{barcodes}.polyAsitesNoErcc.bed"
	output: "mappings/" + "getPolyAreadsList/{techname}_{capDesign}_{barcodes}.polyAreads.list"
	shell:
		'''
cat {input} | cut -f4 | sort|uniq > {output}
		'''

rule getPolyAreadsStats:
	input:
		mappedReads= "mappings/" + "mergeSizeFracBams/{techname}_{capDesign}_{barcodes}.merged.bam",
		polyAreads = "mappings/" + "getPolyAreadsList/{techname}_{capDesign}_{barcodes}.polyAreads.list"
	output: config["STATSDATADIR"] + "{techname}_{capDesign}_{barcodes}.polyAreads.stats.tsv"
	shell:
		'''
mapped=$(samtools view -F4 {input.mappedReads} |cut -f1|sort|uniq|wc -l)
polyA=$(cat {input.polyAreads} | wc -l)
echo -e "{wildcards.techname}\t{wildcards.capDesign}\t{wildcards.barcodes}\t$mapped\t$polyA" | awk '{{print $0"\t"$5/$4}}' > {output}

		'''

rule aggPolyAreadsCapDesignStats:
	input: lambda wildcards: expand(config["STATSDATADIR"] + "{techname}_{capDesign}_{barcodes}.polyAreads.stats.tsv", filtered_product, techname=wildcards.techname, capDesign=wildcards.capDesign, barcodes=BARCODES)
	output: config["STATSDATADIR"] + "{techname}_{capDesign}.polyAreads.stats.tsv"
	shell:
		'''
totalMapped=$(cat {input} | cut -f4 | sum.sh)
totalPolyA=$(cat {input} | cut -f5 | sum.sh)
echo -e "{wildcards.techname}\t{wildcards.capDesign}\t$totalMapped\t$totalPolyA" | awk '{{print $0"\t"$4/$3}}' > {output}
		'''

rule aggPolyAreadsStats:
	input: lambda wildcards: expand(config["STATSDATADIR"] + "{techname}_{capDesign}.polyAreads.stats.tsv", techname=TECHNAMES, capDesign=CAPDESIGNS)
	output: config["STATSDATADIR"] + "all.polyAreads.stats.tsv"
	shell:
		'''
echo -e "seqTech\tcorrectionLevel\tcapDesign\tcategory\tcount" > {output}
cat {input} | awk '{{print $1"\\t"$2"\\tnonPolyA\\t"$3-$4"\\n"$1"\\t"$2"\\tpolyA\\t"$4}}' | sed 's/Corr0/\tNo/' | sed 's/Corr90/\tYes/' | sort -r >> {output}
		'''

rule plotPolyAreadsStats:
	input: config["STATSDATADIR"] + "all.polyAreads.stats.tsv"
	output: config["PLOTSDIR"] + "all.polyAreads.stats.{ext}"
	shell:
		'''
echo "library(ggplot2)
library(plyr)
library(scales)
dat <- read.table('{input}', header=T, as.is=T, sep='\\t')
ggplot(data=dat, aes(x=factor(correctionLevel), y=count, fill=category)) +
geom_bar(stat='identity') + scale_fill_manual(values=c('polyA' = '#25804C', 'nonPolyA' = '#FB3B24')) + facet_grid( seqTech ~ capDesign)+ ylab('# mapped reads') + xlab('Error correction') + guides(fill = guide_legend(title='Category'))+ scale_y_continuous(labels=scientific)+
{GGPLOT_PUB_QUALITY}
ggsave('{output}', width=7, height=3)
" > {output}.r
cat {output}.r | R --slave

		'''


rule clusterPolyAsites:
	input: "mappings/" + "removePolyAERCCs/{techname}_{capDesign}_{barcodes}.polyAsitesNoErcc.bed"
	output: "mappings/" + "clusterPolyAsites/{techname}_{capDesign}_{barcodes}.polyAsites.clusters.bed"
	shell:
		'''
cat {input} | bedtools merge -s -d 5 -c 4,6 -o distinct -i stdin | awk '{{print $1"\t"$2"\t"$3"\t"$4"\t0\t"$5}}'| perl -F"\\t" -lane 'if($F[5] eq "+"){{$F[1]=$F[2]-1}}elsif($F[5] eq "-"){{$F[2]=$F[1]+1}}else{{die}} print join("\t",@F);'|sortbed > {output}
		'''
rule makePolyABigWigs:
	input:
		sites = "mappings/" + "removePolyAERCCs/{techname}_{capDesign}_{barcodes}.polyAsitesNoErcc.bed",
		genome = lambda wildcards: config["GENOMESDIR"] + CAPDESIGNTOGENOME[wildcards.capDesign] + ".genome"
	output: "mappings/" + "makePolyABigWigs/{techname}_{capDesign}_{barcodes}.polyAsitesNoErcc.{strand}.bw"
	shell:
		'''
bedtools genomecov -strand {wildcards.strand} -split -bg -i {input.sites} -g {input.genome} > {output}.bedgraph
bedGraphToBigWig {output}.bedgraph {input.genome} {output}
rm -f {output}.bedgraph
		'''
