rule makeStarIndex:
	input: genome = config["GENOMESDIR"] +"{genome}.sorted.fa"
	output: config["GENOMESDIR"] + "STARshort_indices/" + "{genome}/SA"
	shell:
		'''
uuid=$(uuidgen)
mkdir -p {config[TMPDIR]}/$uuid ; 
mkdir -p $(dirname {output}); 
conda activate star
STAR --runMode genomeGenerate --runThreadN 3 --genomeDir {config[TMPDIR]}/$uuid --genomeFastaFiles {input}
mv -f {config[TMPDIR]}/$uuid/* $(dirname {output})
		'''


rule hiSeqReadMapping:
	input:
		reads1 = config["MATCHED_HISEQ_PATH"] + "hiSeq_{capDesign}_1.fastq.gz",
		reads2 = config["MATCHED_HISEQ_PATH"] + "hiSeq_{capDesign}_2.fastq.gz",
		genome = lambda wildcards: config["GENOMESDIR"] + "STARshort_indices/" + CAPDESIGNTOGENOME[wildcards.capDesign] + "/SA",
#		referenceAnnot = lambda wildcards: CAPDESIGNTOANNOTGTF[wildcards.capDesign]
	threads: 12
	output:
		"mappings/hiSeq_{capDesign}.bam"
	wildcard_constraints:
		barcodes='[^(allTissues)][\S]+',
	shell:
		'''
uuidTmpOut=$(uuidgen)
 conda activate star

echoerr "Mapping"
mkdir -p mappings/STAR/`basename {output}`/
STAR \
--runThreadN {threads} \
--readFilesIn {input.reads1} {input.reads2} \
--genomeDir $(dirname {input.genome}) \
--readFilesCommand zcat \
--sjdbOverhang 124 \
--outFileNamePrefix "mappings/STAR/`basename {output}`/" \
--outStd SAM \
--genomeLoad NoSharedMemory \
--outSAMunmapped Within \
--outFilterType BySJout \
--outFilterMultimapNmax 20 \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverLmax 0.04 \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--outFilterIntronMotifs RemoveNoncanonical \
--outSAMstrandField intronMotif \
--outSAMattributes NH HI NM MD AS nM XS \
| samtools view -b -u -S - | samtools sort -T {config[TMPDIR]}  -@ 2   -m 15000000000 - > {config[TMPDIR]}/$uuidTmpOut
sleep 200s
samtools index {config[TMPDIR]}/$uuidTmpOut
echoerr "Mapping done"
mv {config[TMPDIR]}/$uuidTmpOut {output}
mv {config[TMPDIR]}/$uuidTmpOut.bai {output}.bai

		'''


rule getHiSeqMappingStats:
	input:
		bams = "mappings/hiSeq_{capDesign}.bam",
		reads1 = config["MATCHED_HISEQ_PATH"] + "hiSeq_{capDesign}_1.fastq.gz",
		reads2 = config["MATCHED_HISEQ_PATH"] + "hiSeq_{capDesign}_2.fastq.gz",
	output: config["STATSDATADIR"] + "tmp/{capDesign}_tmp.hiSeq.mapping.stats.tsv"
	shell:
		'''
uuidTmpOut=$(uuidgen)
totalReads=$(zcat {input.reads1} {input.reads2}  | fastq2tsv.pl | wc -l)
mappedReads=$(samtools view  -F 4 {input.bams}|cut -f1|sort -T {config[TMPDIR]} |uniq|wc -l)
echo -e "{wildcards.capDesign}\t$totalReads\t$mappedReads" | awk '{{print $0"\t"$3/$2}}' > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
		'''

rule aggHiSeqMappingStats:
	input: lambda wildcards: expand(config["STATSDATADIR"] + "tmp/{capDesign}_tmp.hiSeq.mapping.stats.tsv",capDesign=CAPDESIGNS)
	output: config["STATSDATADIR"] + "all.hiSeq.mapping.stats.tsv"
	shell:
		'''
uuidTmpOut=$(uuidgen)
cat {input} | sort -T {config[TMPDIR]}  > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
		'''


rule plotHiSeqMappingStats:
	input: config["STATSDATADIR"] + "all.hiSeq.mapping.stats.tsv"
	output: config["PLOTSDIR"] + "hiSeq.mapping.stats/all.hiSeq.mapping.stats.{ext}"
	shell:
		'''
echo "library(ggplot2)
library(plyr)
library(scales)
themeSize = 9
lineSize=0.175
minorLineSize=lineSize/2

dat <- read.table('{input}', header=F, as.is=T, sep='\\t')
colnames(dat)<-c('capDesign', 'totalReads', 'mappedReads', 'percentMappedReads')
ggplot(dat, aes(x=capDesign, y=percentMappedReads), fill='#b0a19d') +
geom_bar(width=0.75,stat='identity', position=position_dodge(width=0.9)) +
geom_hline(aes(yintercept=1), linetype='dashed', alpha=0.7)+
geom_text(aes(group=capDesign, y=0.01, label = paste(sep='',percent(percentMappedReads),' / ','(',comma(mappedReads),')')), angle=90, size=5, hjust=0, vjust=0.5, position = position_dodge(width=0.9)) +
scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
xlab ('capDesign') +
theme_bw(base_size=17) +
{GGPLOT_PUB_QUALITY} + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('{output}', width=7, height=9)
" > {output}.r
 set +eu
conda activate R_env
set -eu
cat {output}.r | R --slave

		'''


rule getHiSeqCanonicalIntronsList:
	input:
		bam="mappings/hiSeq_{capDesign}.bam",
		genome = lambda wildcards: config["GENOMESDIR"] + CAPDESIGNTOGENOME[wildcards.capDesign] + ".sorted.fa"
	threads: 6
	output:
		list="mappings/hiSeqIntrons/hiSeq_{capDesign}.canonicalIntrons.list",
		stats=config["STATSDATADIR"] + "tmp/{capDesign}_tmp.hiSeq.SJs.stats.tsv"
	shell:
		'''
uuid=$(uuidgen)
uuidTmpOutL=$(uuidgen)
uuidTmpOutS=$(uuidgen)

echo $PATH

echo 

which perl

echoerr "making bed"
set +eu
# conda activate julenv
conda env list
set -eu

samtools view -b -F 256 -F4 -F 2048 {input.bam}  |bedtools bamtobed -i stdin -bed12 | fgrep -v ERCC- > {config[TMPDIR]}/$uuid.hiSeq_{wildcards.capDesign}.bed
set +eu

conda deactivate

conda env list
set -eu

echo $PATH

echo 

which perl




echoerr "splitting"
split -a 3 -d -e -n l/24 {config[TMPDIR]}/$uuid.hiSeq_{wildcards.capDesign}.bed {config[TMPDIR]}/$uuid.hiSeq_{wildcards.capDesign}.bed.split
rm {config[TMPDIR]}/$uuid.hiSeq_{wildcards.capDesign}.bed
for file in `ls {config[TMPDIR]}/$uuid.hiSeq_{wildcards.capDesign}.bed.split*`; do
echo "cat $file | awk -f ~/julien_utils/bed12fields2gff.awk > $file.gff; sort -T {config[TMPDIR]} -k1,1 -k4,4n -k5,5n $file.gff | makeIntrons.pl - | extract_intron_strand_motif.pl - {input.genome} {config[TMPDIR]}/$(basename $file); rm $file $file.gff $file.transcripts.tsv"
done > {config[TMPDIR]}/$uuid.parallelIntrons.sh
echoerr "extracting introns on split files"

parallel -j {threads} < {config[TMPDIR]}/$uuid.parallelIntrons.sh
echoerr "getting SJs and merging into output..."
cat {config[TMPDIR]}/$uuid.hiSeq_{wildcards.capDesign}.bed.split*.introns.gff | perl -lane '$start=$F[3]-1; $end=$F[4]+1; print $F[0]."_".$start."_".$end."_".$F[6]' | sort  -T {config[TMPDIR]} |uniq> {config[TMPDIR]}/$uuidTmpOutL
echo -e "{wildcards.capDesign}\t$(cat {config[TMPDIR]}/$uuidTmpOutL | wc -l )" > {config[TMPDIR]}/$uuidTmpOutS
mv {config[TMPDIR]}/$uuidTmpOutS {output.stats}
mv {config[TMPDIR]}/$uuidTmpOutL {output.list}

		'''

rule aggHiSeqSjStats:
	input: lambda wildcards: expand(config["STATSDATADIR"] + "tmp/{capDesign}_tmp.hiSeq.SJs.stats.tsv",capDesign=CAPDESIGNS)
	output: config["STATSDATADIR"] + "all.hiSeq.SJs.stats.tsv"
	shell:
		'''
uuidTmpOut=$(uuidgen)
cat {input} | sort -T {config[TMPDIR]}  > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
		'''

rule plotHiSeqSjStats:
	input: config["STATSDATADIR"] + "all.hiSeq.SJs.stats.tsv"
	output: config["PLOTSDIR"] + "hiSeq.SJs.stats/all.hiSeq.SJs.stats.{ext}"
	shell:
		'''
echo "library(ggplot2)
library(plyr)
library(scales)
themeSize=9
lineSize=0.175
minorLineSize=lineSize/2

dat <- read.table('{input}', header=F, as.is=T, sep='\\t')
colnames(dat)<-c('capDesign', 'totalSJs')
ggplot(dat, aes(x=capDesign, y=totalSJs), fill='#BCB5B5') +
geom_bar(width=0.75,stat='identity', position=position_dodge(width=0.9)) +
xlab ('capDesign') +
{GGPLOT_PUB_QUALITY} + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('{output}', width=7, height=9)
" > {output}.r
 set +eu
conda activate R_env
set -eu
cat {output}.r | R --slave

		'''
