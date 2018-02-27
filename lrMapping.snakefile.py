
# # STAR mapping of reads:
# rule readMapping:
# #	wildcard_constraints:
# #		 barcodesU = lambda wildcards: {wildcards.capDesign} + "_.+"
# 	input:
# #		reads = returnCapDesignBarcodesFastqs,
# 		reads = config["DEMULTIPLEX_DIR"] + "{capDesign}_{sizeFrac}.{barcodesU}.fastq",
# 		reference=  lambda wildcards: "/users/rg/buszczynska/Projects/RNA.capture/pacbio.postCapture.20150317/mapping/reference/" + CAPDESIGNTOGENOME[wildcards.capDesign] + "/star.index/"
# 	output:
# 		config["PB_MAPPINGS"] + "{capDesign}_{sizeFrac}.{barcodesU}.bam"
# 	shell:
# 		'''
# # (Recommended STAR parameters from https://github.com/PacificBiosciences/cDNA_primer/wiki/Bioinfx-study:-Optimizing-STAR-aligner-for-Iso-Seq-data )

# sameBarcodeCapDesign=$(echo "{wildcards.capDesign} {wildcards.barcodesU}" | perl -slane 'if($F[1]=~/$F[0]/ || $F[1] eq "Undeter"){{print "1"}} else {{print "0"}}')

# if [ {TECHNAME} == "ont" && {wildcards.barcodesU} == "Undeter" ]; then
# run=0
# else
# run=1
# fi


# if [ $sameBarcodeCapDesign == 1 && $run == 1 ]; then
# /users/rg/jlagarde/bin/STAR-2.5.3a/bin/Linux_x86_64/STARlong \
# --runMode alignReads \
# --outSAMattributes NH HI NM MD \
# --readNameSeparator space \
# --outFilterMultimapScoreRange 1 \
# --outFilterMismatchNmax 2000 \
# --scoreGapNoncan -20 \
# --scoreGapGCAG -4 \
# --scoreGapATAC -8 \
# --scoreDelOpen -1 \
# --scoreDelBase -1 \
# --scoreInsOpen -1 \
# --scoreInsBase -1 \
# --alignEndsType Local \
# --seedSearchStartLmax 50 \
# --seedPerReadNmax 100000 \
# --seedPerWindowNmax 1000 \
# --alignTranscriptsPerReadNmax 100000 \
# --alignTranscriptsPerWindowNmax 10000 \
# --readFilesIn {input.reads} \
# --genomeDir {input.reference} \
# --outFileNamePrefix {pacBioMappingDir}/{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.barcodesU} \
# --outStd SAM \
# | samtools view -b -u -S - \
# | samtools sort -m 20000000000 - >{output}
# samtools index {output}
# else
# touch {output}
# fi
# 		'''

# mapping of long reads:
rule readMapping:
#	wildcard_constraints:
#		 barcodesU = lambda wildcards: {wildcards.capDesign} + "_.+"
	input:
#		reads = returnCapDesignBarcodesFastqs,
		reads = config["DEMULTIPLEX_DIR"] + "{techname}_{capDesign}_{sizeFrac}.{barcodesU}.fastq.gz",
		genome = lambda wildcards: "/users/rg/jlagarde/genomes/" + CAPDESIGNTOGENOME[wildcards.capDesign] + ".fa"
#	params:
#		reference=  lambda wildcards: CAPDESIGNTOGENOME[wildcards.capDesign]
	threads: 12
	output:
		temp("mappings/" + "{techname}_{capDesign}_{sizeFrac}.{barcodesU}.bam")
	shell:
		'''
sameBarcodeCapDesign=$(echo "{wildcards.capDesign} {wildcards.barcodesU}" | perl -slane 'if($F[1]=~/$F[0]/ || $F[1] eq "Undeter"){{print "1"}} else {{print "0"}}')
run=1
if [[ $sameBarcodeCapDesign == 1 && $run == 1 ]]; then
echoerr "Mapping"
minimap2 --cs -t {threads} --secondary=no -L -ax splice {input.genome} {input.reads} > {output}.tmp
echoerr "Mapping done"
echoerr "Creating/sorting BAM"
cat {output}.tmp | samtools view -F 256 -F4 -F 2048 -b -u -S - | samtools sort --threads {threads} -T $TMPDIR -m 5G - >{output}
echoerr "Done creating/sorting BAM"
#echoerr "Indexing BAM"
#samtools index {output}
rm {output}.tmp
else
touch {output}
fi
		'''

# gmap -D /users/rg/jlagarde/genomes/gmapdb/ -d {params.reference} -f samse -n 0 -t {threads} {input.reads} > {output}.tmp


rule getMappingStats:
	input:
		bams = "mappings/" + "{techname}_{capDesign}_{sizeFrac}.{barcodesU}.bam",
		fastqs = config["DEMULTIPLEX_DIR"] + "{techname}_{capDesign}_{sizeFrac}.{barcodesU}.fastq.gz"
	output: temp(config["STATSDATADIR"] + "{techname}_{capDesign}_{sizeFrac}.{barcodesU}.mapping.perSample.perFraction.stats.tsv")
	shell:
		'''
sameBarcodeCapDesign=$(echo "{wildcards.capDesign} {wildcards.barcodesU}" | perl -slane 'if($F[1]=~/$F[0]/ || $F[1] eq "Undeter"){{print "1"}} else {{print "0"}}')
if [ $sameBarcodeCapDesign == 1 ]; then
totalReads=$(zcat {input.fastqs} | fastq2tsv.pl | wc -l)
mappedReads=$(samtools view  -F 4 {input.bams}|cut -f1|sort|uniq|wc -l)
echo -e "{wildcards.capDesign}\t{wildcards.sizeFrac}\t{wildcards.barcodesU}\t$totalReads\t$mappedReads" | awk '{{print $0"\t"$5/$4}}' > {output}
else
touch {output}
fi
		'''
rule aggMappingStats:
	input: lambda wildcards: expand(config["STATSDATADIR"] + "{techname}_{capDesign}_{sizeFrac}.{barcodesU}.mapping.perSample.perFraction.stats.tsv", filtered_product, techname=wildcards.techname, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodesU=BARCODESUNDETER)
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
dropbox_uploader.sh upload {output} {DROPBOX_PLOTS};

		'''


rule mergeSizeFracBams:
	input: lambda wildcards: expand("mappings/" + "{techname}_{capDesign}_{sizeFrac}.{barcodesU}.bam", filtered_product, techname=wildcards.techname, capDesign=wildcards.capDesign, sizeFrac=SIZEFRACS, barcodesU=wildcards.barcodesU)
	output: "mappings/" + "{techname}_{capDesign}_{barcodesU}.merged.bam"
	shell:
		'''
sameBarcodeCapDesign=$(echo "{wildcards.capDesign} {wildcards.barcodesU}" | perl -slane 'if($F[1]=~/$F[0]/ || $F[1] eq "Undeter"){{print "1"}} else {{print "0"}}')
if [ $sameBarcodeCapDesign == 1 ]; then
samtools merge {output} {input}
samtools index {output}
else
touch {output}
fi
		'''

rule checkOnlyOneHit:
	input: "mappings/" + "{techname}_{capDesign}_{barcodesU}.merged.bam"
	output: "mappings/" + "qc/{techname}_{capDesign}_{barcodesU}.merged.bam.dupl.txt"
	shell:
		'''
samtools view {input} | cut -f1 | sort| uniq -dc > {output}
count=$(cat {output} | wc -l)
if [ $count -gt 0 ]; then echo "$count duplicate read IDs found"; mv {output} {output}.tmp; exit 1; fi
		'''


rule readBamToBed:
	input: lambda wildcards: expand("mappings/" + "{techname}_{capDesign}_{barcodesU}.merged.bam", filtered_product, techname=wildcards.techname, capDesign=wildcards.capDesign, barcodesU=wildcards.barcodesU)
	output: "mappings/bed/" + "{techname}_{capDesign}_{barcodesU}.merged.bed"
	shell:
		'''
sameBarcodeCapDesign=$(echo "{wildcards.capDesign} {wildcards.barcodesU}" | perl -slane 'if($F[1]=~/$F[0]/ || $F[1] eq "Undeter"){{print "1"}} else {{print "0"}}')
if [ $sameBarcodeCapDesign == 1 ]; then
bamToBed -i {input} -bed12 > {output}
else
touch {output}
fi

		'''

rule readBedToGff:
	input: lambda wildcards: expand("mappings/bed/" + "{techname}_{capDesign}_{barcodesU}.merged.bed", filtered_product, techname=wildcards.techname, capDesign=wildcards.capDesign, barcodesU=wildcards.barcodesU)
	output: "mappings/gff/" + "{techname}_{capDesign}_{barcodesU}.merged.gff"
	shell:
		'''
sameBarcodeCapDesign=$(echo "{wildcards.capDesign} {wildcards.barcodesU}" | perl -slane 'if($F[1]=~/$F[0]/ || $F[1] eq "Undeter"){{print "1"}} else {{print "0"}}')
if [ $sameBarcodeCapDesign == 1 ]; then
cat {input} | awk -f ~jlagarde/julien_utils/bed12fields2gff.awk | sortgff> {output}
else
touch {output}
fi
		'''


rule mergeCapDesignBams:
	input: lambda wildcards: expand("mappings/" + "{techname}_{capDesign}_{barcodes}.merged.bam", techname=wildcards.techname, capDesign=wildcards.capDesign, barcodes=BARCODES)
	output: temp("mappings/" + "{techname}_{capDesign}.merged2.bam")
	shell:
		'''
samtools merge {output} {input}
samtools index {output}

		'''

rule qualimap:
	input: "mappings/" + "{techname}_{capDesign}.merged2.bam"
	output: "mappings/qualimap_reports/" + "{techname}_{capDesign}.merged2/genome_results.txt"
	shell:
		'''
~/bin/qualimap_v2.2.1/qualimap bamqc -bam {input} -outdir mappings/qualimap_reports/{wildcards.techname}_{wildcards.capDesign}.merged2/ --java-mem-size=10G -outfile {wildcards.techname}_{wildcards.capDesign}.merged2
touch {output}
		'''

