

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
echoerr "Mapping"
minimap2 -x {params.minimap_preset} -t {threads} --secondary=no -L -a {input.genome} {input.reads} > {config[TMPDIR]}/$uuid
echoerr "Mapping done"
echoerr "Creating/sorting BAM"
samtools view -H {config[TMPDIR]}/$uuid > {config[TMPDIR]}/$uuid.2
samtools view -F 256 -F4 -F 2048 {config[TMPDIR]}/$uuid >> {config[TMPDIR]}/$uuid.2
cat {config[TMPDIR]}/$uuid.2 | samtools sort -T {config[TMPDIR]}  --threads {threads}  -m 5G - > {output}
echoerr "Done creating/sorting BAM"
sleep 70s
samtools index {output}
		'''

rule getMappingStats:
	input:
		bams = "mappings/readMapping/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bam",
		fastqs = DEMULTIPLEXED_FASTQS + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.{barcodes}.fastq.gz"
	output: temp(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.{barcodes}.mapping.perSample.perFraction.stats.tsv")
	shell:
		'''
totalReads=$(zcat {input.fastqs} | fastq2tsv.pl | wc -l)
mappedReads=$(samtools view  -F 4 {input.bams}|cut -f1|sort -T {config[TMPDIR]} |uniq|wc -l)
echo -e "{wildcards.capDesign}\t{wildcards.sizeFrac}\t{wildcards.barcodes}\t$totalReads\t$mappedReads" | awk '{{print $0"\t"$5/$4}}' > {output}
		'''
rule aggMappingStats:
	input: lambda wildcards: expand(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.{barcodes}.mapping.perSample.perFraction.stats.tsv",filtered_product, techname=wildcards.techname, corrLevel=wildcards.corrLevel, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES)
	output: config["STATSDATADIR"] + "{techname}Corr{corrLevel}.mapping.perSample.perFraction.stats.tsv"
	shell:
		'''
cat {input} | sort -T {config[TMPDIR]}  > {output}
		'''

rule plotMappingStats:
	input: config["STATSDATADIR"] + "{techname}Corr{corrLevel}.mapping.perSample.perFraction.stats.tsv"
	output: config["PLOTSDIR"] + "lrMapping.stats/{techname}Corr{corrLevel}.mapping.perSample.perFraction.stats.{ext}"
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

		'''


rule mergeSizeFracBams:
	input: lambda wildcards: expand("mappings/readMapping/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bam", filtered_product, techname=wildcards.techname, corrLevel=wildcards.corrLevel, capDesign=wildcards.capDesign, sizeFrac=SIZEFRACS, barcodes=wildcards.barcodes)
	output: "mappings/readMapping/{techname}Corr{corrLevel}_{capDesign}_allFracs_{barcodes}.bam"
	wildcard_constraints:
		barcodes='(?!allTissues).+',
		techname='(?!allSeqTechs).+'
	shell:
		'''

samtools merge {output} {input}
sleep 120s
samtools index {output}
		'''

rule mergeBarcodeBams:
	input: lambda wildcards: expand("mappings/readMapping/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bam", filtered_product, techname=wildcards.techname, corrLevel=wildcards.corrLevel, capDesign=wildcards.capDesign, sizeFrac=wildcards.sizeFrac, barcodes=BARCODES)
	output: "mappings/readMapping/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_allTissues.bam"
	wildcard_constraints:
		sizeFrac='[0-9-+\.]+',
		techname='(?!allSeqTechs).+' #to avoid ambiguity with downstream merging rules
	shell:
		'''

samtools merge {output} {input}
sleep 120s
samtools index {output}
		'''


rule mergeSizeFracAndBarcodeBams:
	input: lambda wildcards: expand("mappings/readMapping/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bam", filtered_product, techname=wildcards.techname, corrLevel=wildcards.corrLevel, capDesign=wildcards.capDesign, sizeFrac=SIZEFRACS, barcodes=BARCODES)
	output: "mappings/readMapping/{techname}Corr{corrLevel}_{capDesign}_allFracs_allTissues.bam"
	wildcard_constraints:
		techname='(?!allSeqTechs).+',
		capDesign='|'.join(CAPDESIGNS)
	shell:
		'''
samtools merge {output} {input}
sleep 120s
samtools index {output}

		'''

rule mergeAllSeqTechsFracsAndTissuesBams:
	input: lambda wildcards: expand("mappings/readMapping/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bam", filtered_product_merge, techname=TECHNAMES, corrLevel=wildcards.corrLevel, capDesign=wildcards.capDesign, sizeFrac=wildcards.sizeFrac, barcodes=wildcards.barcodes)
	output: "mappings/readMapping/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bam"
	wildcard_constraints:
#		techname='allSeqTechs',
		sizeFrac='allFracs',
		barcodes='allTissues',
		capDesign='|'.join(CAPDESIGNS)
	shell:
		'''
samtools merge {output} {input}
sleep 120s
samtools index {output}

		'''


rule mergeAllCapDesignsAllSeqTechsFracsAndTissuesBams:
	input: lambda wildcards: expand("mappings/readMapping/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bam", filtered_capDesign_product_merge, techname=wildcards.techname, corrLevel=wildcards.corrLevel, capDesign=wildcards.capDesign, sizeFrac=wildcards.sizeFrac, barcodes=wildcards.barcodes)
	output: "mappings/readMapping/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bam"
	wildcard_constraints:
#		techname='allSeqTechs',
		sizeFrac='allFracs',
		barcodes='allTissues',
		capDesign='|'.join(MERGEDCAPDESIGNS)
	shell:
		'''
samtools merge {output} {input}
sleep 120s
samtools index {output}

		'''


rule checkOnlyOneHit:
	input: "mappings/readMapping/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bam"
	output: "mappings/readMapping/qc/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bam.dupl.txt"
	shell:
		'''
samtools view {input} | cut -f1 | sort -T {config[TMPDIR]} | uniq -dc > {output}
count=$(cat {output} | wc -l)
if [ $count -gt 0 ]; then echo "$count duplicate read IDs found"; mv {input} {input}.dup.bkp; exit 1; fi
		'''


rule readBamToBed:
	input: "mappings/readMapping/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bam"
	output: "mappings/readBamToBed/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bed.gz"
	#input: "mappings/readMapping/{basename}.bam"
	#output: "mappings/readBamToBed/{basename}.bed"
	shell:
		'''
bedtools bamtobed -i {input} -bed12 | sort -T {config[TMPDIR]}  -k1,1 -k2,2n -k3,3n  | gzip > {output}

		'''

rule readBedToGff:
	input: "mappings/readBamToBed/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bed.gz"
	output: "mappings/readBedToGff/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.gff.gz"
	shell:
		'''
zcat {input} | awk -f ~jlagarde/julien_utils/bed12fields2gff.awk | sort -T {config[TMPDIR]}  -k1,1 -k4,4n -k5,5n  | gzip > {output}
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
