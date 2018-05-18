

# mapping of long reads:
rule readMapping:
#	wildcard_constraints:
#		 barcodesU = lambda wildcards: {wildcards.capDesign} + "_.+"
	input:
#		reads = returnCapDesignBarcodesFastqs,
#		reads = config["DEMULTIPLEX_DIR"] + "demultiplexFastqs/{techname}_{capDesign}_{sizeFrac}.{barcodes}.fastq.gz",
		reads = lambda wildcards: expand(config["DEMULTIPLEX_DIR"] + "demultiplexFastqs/{techname}_{capDesign}_{sizeFrac}.{barcodes}.fastq.gz", filtered_product, techname=wildcards.techname, capDesign=wildcards.capDesign, sizeFrac=wildcards.sizeFrac,barcodes=wildcards.barcodes),
		genome = lambda wildcards: config["GENOMESDIR"] + CAPDESIGNTOGENOME[wildcards.capDesign] + ".fa"
#	params:
#		reference=  lambda wildcards: CAPDESIGNTOGENOME[wildcards.capDesign]
	threads: 12
	output:
		temp("mappings/" + "readMapping/{techname}_{capDesign}_{sizeFrac}_{barcodes}.bam")
	shell:
		'''
echoerr "Mapping"
minimap2 --cs -t {threads} --secondary=no -L -ax splice {input.genome} {input.reads} > {output}.tmp
echoerr "Mapping done"
echoerr "Creating/sorting BAM"
cat {output}.tmp | samtools view -F 256 -F4 -F 2048 -b -u -S - | samtools sort --threads {threads} -T $TMPDIR -m 5G - >{output}
echoerr "Done creating/sorting BAM"
rm {output}.tmp
		'''

rule getMappingStats:
	input:
		bams = "mappings/" + "readMapping/{techname}_{capDesign}_{sizeFrac}_{barcodes}.bam",
		fastqs = config["DEMULTIPLEX_DIR"] + "demultiplexFastqs/{techname}_{capDesign}_{sizeFrac}.{barcodes}.fastq.gz"
	output: temp(config["STATSDATADIR"] + "{techname}_{capDesign}_{sizeFrac}.{barcodes}.mapping.perSample.perFraction.stats.tsv")
	shell:
		'''
totalReads=$(zcat {input.fastqs} | fastq2tsv.pl | wc -l)
mappedReads=$(samtools view  -F 4 {input.bams}|cut -f1|sort|uniq|wc -l)
echo -e "{wildcards.capDesign}\t{wildcards.sizeFrac}\t{wildcards.barcodes}\t$totalReads\t$mappedReads" | awk '{{print $0"\t"$5/$4}}' > {output}
		'''
rule aggMappingStats:
	input: lambda wildcards: expand(config["STATSDATADIR"] + "{techname}_{capDesign}_{sizeFrac}.{barcodes}.mapping.perSample.perFraction.stats.tsv",filtered_product, techname=wildcards.techname, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES)
	output: config["STATSDATADIR"] + "{techname}.mapping.perSample.perFraction.stats.tsv"
	shell:
		'''
cat {input} | sort > {output}
		'''

rule plotMappingStats:
	input: config["STATSDATADIR"] + "{techname}.mapping.perSample.perFraction.stats.tsv"
	output: config["PLOTSDIR"] + "{techname}.mapping.perSample.perFraction.stats.{ext}"
	shell:
		'''
echo "library(ggplot2)
library(plyr)
library(scales)
dat <- read.table('{input}', header=F, as.is=T, sep='\\t')
colnames(dat)<-c('capDesign', 'sizeFraction','barcode','totalReads', 'mappedReads', 'percentMappedReads')
ggplot(dat, aes(x=barcode, y=percentMappedReads, fill=sizeFraction)) +
geom_bar(width=0.75,stat='identity', position=position_dodge(width=0.9)) +
scale_fill_manual(values={sizeFrac_Rpalette}) +
geom_hline(aes(yintercept=1), linetype='dashed', alpha=0.7)+
facet_grid(sizeFraction ~ capDesign, scales='free') +
geom_text(aes(group=sizeFraction, y=0.01, label = paste(sep='',percent(percentMappedReads),' / ','(',comma(mappedReads),')')), angle=90, size=5, hjust=0, vjust=0.5, position = position_dodge(width=0.9)) +
scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
xlab ('Sample (barcode)') +
theme_bw(base_size=17) +
{GGPLOT_PUB_QUALITY} + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('{output}', width=13, height=9)
" > {output}.r
cat {output}.r | R --slave
#dropbox_uploader.sh upload {output} {DROPBOX_PLOTS};

		'''


rule mergeSizeFracBams:
	input: lambda wildcards: expand("mappings/" + "readMapping/{techname}_{capDesign}_{sizeFrac}_{barcodes}.bam", filtered_product, techname=wildcards.techname, capDesign=wildcards.capDesign, sizeFrac=SIZEFRACS, barcodes=wildcards.barcodes)
	output: "mappings/" + "mergeSizeFracBams/{techname}_{capDesign}_{barcodes}.merged.bam"
	shell:
		'''
samtools merge {output} {input}
samtools index {output}
		'''

rule checkOnlyOneHit:
	input: "mappings/" + "mergeSizeFracBams/{techname}_{capDesign}_{barcodes}.merged.bam"
	output: "mappings/" + "qc/{techname}_{capDesign}_{barcodes}.merged.bam.dupl.txt"
	shell:
		'''
samtools view {input} | cut -f1 | sort| uniq -dc > {output}
count=$(cat {output} | wc -l)
if [ $count -gt 0 ]; then echo "$count duplicate read IDs found"; mv {output} {output}.tmp; exit 1; fi
		'''


rule readBamToBed:
	input: lambda wildcards: expand("mappings/" + "mergeSizeFracBams/{techname}_{capDesign}_{barcodes}.merged.bam", filtered_product, techname=wildcards.techname, capDesign=wildcards.capDesign, barcodes=wildcards.barcodes)
	output: "mappings/" + "readBamToBed/{techname}_{capDesign}_{barcodes}.merged.bed"
	shell:
		'''
bamToBed -i {input} -bed12 > {output}

		'''

rule readBedToGff:
	input: lambda wildcards: expand("mappings/" + "readBamToBed/{techname}_{capDesign}_{barcodes}.merged.bed", filtered_product, techname=wildcards.techname, capDesign=wildcards.capDesign, barcodes=wildcards.barcodes)
	output: "mappings/" + "readBedToGff/{techname}_{capDesign}_{barcodes}.merged.gff"
	shell:
		'''
cat {input} | awk -f ~jlagarde/julien_utils/bed12fields2gff.awk | sortgff> {output}
		'''


rule mergeCapDesignBams:
	input: lambda wildcards: expand("mappings/" + "mergeSizeFracBams/{techname}_{capDesign}_{barcodes}.merged.bam", filtered_product, techname=wildcards.techname, capDesign=wildcards.capDesign, barcodes=BARCODES)
	output: "mappings/" + "mergeCapDesignBams/{techname}_{capDesign}.merged2.bam"
	shell:
		'''
samtools merge {output} {input}
#sleep so latency doesn't make bai younger than bam
sleep 60s
samtools index {output}

		'''

rule qualimap:
	input: "mappings/" + "mergeCapDesignBams/{techname}_{capDesign}.merged2.bam"
	output: "mappings/qualimap_reports/" + "{techname}_{capDesign}.merged2/genome_results.txt"
	shell:
		'''
unset DISPLAY
~/bin/qualimap_v2.2.1/qualimap bamqc -bam {input} -outdir mappings/qualimap_reports/{wildcards.techname}_{wildcards.capDesign}.merged2/ --java-mem-size=10G -outfile {wildcards.techname}_{wildcards.capDesign}.merged2
touch {output}
		'''

