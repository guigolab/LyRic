rule makeStarIndex:
	input: genome = config["GENOMESDIR"] +"{genome}.sorted.fa"
	output: config["GENOMESDIR"] + "STARshort_indices/" + "{genome}/SA"
	conda: "envs/star_env.yml"
	shell:
		'''
uuid=$(uuidgen)
mkdir -p {TMPDIR}/$uuid ; 
mkdir -p $(dirname {output}); 
STAR --runMode genomeGenerate --runThreadN 3 --genomeDir {TMPDIR}/$uuid --genomeFastaFiles {input}
mv -f {TMPDIR}/$uuid/* $(dirname {output})
		'''


rule hiSeqReadMapping:
	input:
		reads1 = "fastqs/hiSeq/" + "hiSeq_{capDesign}_1.fastq.gz",
		reads2 = "fastqs/hiSeq/" + "hiSeq_{capDesign}_2.fastq.gz",
		genome = lambda wildcards: config["GENOMESDIR"] + "STARshort_indices/" + CAPDESIGNTOGENOME[wildcards.capDesign] + "/SA",
	threads: 12
	conda: "envs/star_env.yml"
	output:
		"output/mappings/shortReadMappings/hiSeq_{capDesign}.bam"
	shell:
		'''
uuidTmpOut=$(uuidgen)

echoerr "Mapping"
mkdir -p output/mappings/STAR/`basename {output}`/
STAR \
--runThreadN {threads} \
--readFilesIn {input.reads1} {input.reads2} \
--genomeDir $(dirname {input.genome}) \
--readFilesCommand zcat \
--sjdbOverhang 124 \
--outFileNamePrefix "output/mappings/STAR/`basename {output}`/" \
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
| samtools view -b -u -S - | samtools sort -T {TMPDIR}  -@ 2   -m 15000000000 - > {TMPDIR}/$uuidTmpOut
sleep 200s
samtools index {TMPDIR}/$uuidTmpOut
echoerr "Mapping done"
mv {TMPDIR}/$uuidTmpOut {output}
mv {TMPDIR}/$uuidTmpOut.bai {output}.bai

		'''


rule getHiSeqMappingStats:
	input:
		bams = "output/mappings/shortReadMappings/hiSeq_{capDesign}.bam",
		reads1 = "fastqs/hiSeq/" + "hiSeq_{capDesign}_1.fastq.gz",
		reads2 = "fastqs/hiSeq/" + "hiSeq_{capDesign}_2.fastq.gz",
	output: "output/statsFiles/" + "tmp/{capDesign}_tmp.hiSeq.mapping.stats.tsv"
	conda: "envs/xtools_env.yml"
	shell:
		'''
uuidTmpOut=$(uuidgen)
totalReads=$(zcat {input.reads1} {input.reads2}  | fastq2tsv.pl | wc -l)
mappedReads=$(samtools view  -F 4 {input.bams}|cut -f1|sort -T {TMPDIR} |uniq|wc -l)
echo -e "{wildcards.capDesign}\t$totalReads\t$mappedReads" | awk '{{print $0"\t"$3/$2}}' > {TMPDIR}/$uuidTmpOut
mv {TMPDIR}/$uuidTmpOut {output}
		'''

rule aggHiSeqMappingStats:
	input: lambda wildcards: expand("output/statsFiles/" + "tmp/{capDesign}_tmp.hiSeq.mapping.stats.tsv",capDesign=CAPDESIGNS)
	output: "output/statsFiles/" + "all.hiSeq.mapping.stats.tsv"
	shell:
		'''
uuidTmpOut=$(uuidgen)
cat {input} | sort -T {TMPDIR}  > {TMPDIR}/$uuidTmpOut
mv {TMPDIR}/$uuidTmpOut {output}
		'''


rule plotHiSeqMappingStats:
	input: "output/statsFiles/" + "all.hiSeq.mapping.stats.tsv"
	output: "output/plots/" + "hiSeq.mapping.stats/all.hiSeq.mapping.stats.{ext}"
	conda: "envs/R_env.yml"
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

cat {output}.r | R --slave

		'''


rule getHiSeqCanonicalIntronsList:
	input:
		bam="output/mappings/shortReadMappings/hiSeq_{capDesign}.bam",
	threads: 6
	output: temp("output/mappings/hiSeqIntrons/hiSeq_{capDesign}.canonicalIntrons.bed")
	conda: "envs/xtools_env.yml"
	shell:
		'''
uuid=$(uuidgen)

echoerr "making bed"

samtools view -b -F 256 -F4 -F 2048 {input.bam}  |bedtools bamtobed -i stdin -bed12 | fgrep -v ERCC- > {TMPDIR}/$uuid.hiSeq_{wildcards.capDesign}.bed

mv {TMPDIR}/$uuid.hiSeq_{wildcards.capDesign}.bed {output}
		'''

rule getHiSeqCanonicalIntronsList2:
	input:
		bed="output/mappings/hiSeqIntrons/hiSeq_{capDesign}.canonicalIntrons.bed",
		genome = lambda wildcards: config["GENOMESDIR"] + CAPDESIGNTOGENOME[wildcards.capDesign] + ".sorted.fa"
	output:	
		list="output/mappings/hiSeqIntrons/hiSeq_{capDesign}.canonicalIntrons.list",
		stats="output/statsFiles/" + "tmp/{capDesign}_tmp.hiSeq.SJs.stats.tsv"
	conda: "envs/perl_env.yml"
	
	shell:
		'''

uuid=$(uuidgen)
uuidTmpOutL=$(uuidgen)
uuidTmpOutS=$(uuidgen)

echoerr "splitting"
split -a 3 -d -e -n l/24 {input.bed} {TMPDIR}/$uuid.hiSeq_{wildcards.capDesign}.bed.split
for file in `ls {TMPDIR}/$uuid.hiSeq_{wildcards.capDesign}.bed.split*`; do
echo "cat $file | bed12togff > $file.gff; sort -T {TMPDIR} -k1,1 -k4,4n -k5,5n $file.gff | makeIntrons.pl - | extract_intron_strand_motif.pl - {input.genome} {TMPDIR}/$(basename $file); rm $file $file.gff $file.transcripts.tsv"
done > {TMPDIR}/$uuid.parallelIntrons.sh
echoerr "extracting introns on split files"

parallel -j {threads} < {TMPDIR}/$uuid.parallelIntrons.sh
echoerr "getting SJs and merging into output..."
cat {TMPDIR}/$uuid.hiSeq_{wildcards.capDesign}.bed.split*.introns.gff | perl -lane '$start=$F[3]-1; $end=$F[4]+1; print $F[0]."_".$start."_".$end."_".$F[6]' | sort  -T {TMPDIR} |uniq> {TMPDIR}/$uuidTmpOutL
echo -e "{wildcards.capDesign}\t$(cat {TMPDIR}/$uuidTmpOutL | wc -l )" > {TMPDIR}/$uuidTmpOutS
mv {TMPDIR}/$uuidTmpOutS {output.stats}
mv {TMPDIR}/$uuidTmpOutL {output.list}

		'''

rule aggHiSeqSjStats:
	input: lambda wildcards: expand("output/statsFiles/" + "tmp/{capDesign}_tmp.hiSeq.SJs.stats.tsv",capDesign=CAPDESIGNS)
	output: "output/statsFiles/" + "all.hiSeq.SJs.stats.tsv"
	shell:
		'''
uuidTmpOut=$(uuidgen)
cat {input} | sort -T {TMPDIR}  > {TMPDIR}/$uuidTmpOut
mv {TMPDIR}/$uuidTmpOut {output}
		'''

rule plotHiSeqSjStats:
	input: "output/statsFiles/" + "all.hiSeq.SJs.stats.tsv"
	output: "output/plots/" + "hiSeq.SJs.stats/all.hiSeq.SJs.stats.{ext}"
	conda: "envs/R_env.yml"
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

cat {output}.r | R --slave

		'''
