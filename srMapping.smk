rule hiSeqReadMapping:
	input:
		reads1 = config["MATCHED_HISEQ_PATH"] + "hiSeq_{capDesign}_1.fastq.gz",
		reads2 = config["MATCHED_HISEQ_PATH"] + "hiSeq_{capDesign}_2.fastq.gz",
		genome = lambda wildcards: config["GENOMESDIR"] + "STARshort_indices/" + CAPDESIGNTOGENOME[wildcards.capDesign] + "/",
		referenceAnnot = lambda wildcards: CAPDESIGNTOANNOTGTF[wildcards.capDesign]
	threads: 12
	output:
		"mappings/hiSeq_{capDesign}.bam"
	wildcard_constraints:
		barcodes='[^(allTissues)][\S]+',
	shell:
		'''
echoerr "Mapping"
mkdir -p mappings/STAR/`basename {output}`/
STAR \
--runThreadN {threads} \
--readFilesIn {input.reads1} {input.reads2} \
--genomeDir {input.genome} \
--readFilesCommand zcat \
--sjdbOverhang 124 \
--sjdbGTFfile {input.referenceAnnot} \
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
| samtools view -b -u -S - | samtools sort -T {config[TMPDIR]}  -@ 2   -m 15000000000 - > {output}
sleep 120s
samtools index {output}
echoerr "Mapping done"

		'''


rule getHiSeqCanonicalIntronsList:
	input:
		bam="mappings/hiSeq_{capDesign}.bam",
		genome = lambda wildcards: config["GENOMESDIR"] + CAPDESIGNTOGENOME[wildcards.capDesign] + ".fa"
	threads: 6
	output:"mappings/hiSeqIntrons/hiSeq_{capDesign}.canonicalIntrons.list"
	shell:
		'''
uuid=$(uuidgen)
echoerr "making bed"
samtools view -b -F 256 -F4 -F 2048 {input.bam}  |bedtools bamtobed -i stdin -bed12 | fgrep -v ERCC- > {config[TMPDIR]}/$uuid.hiSeq_{wildcards.capDesign}.bed
echoerr "splitting"
split -a 3 -d -e -n l/24 {config[TMPDIR]}/$uuid.hiSeq_{wildcards.capDesign}.bed {config[TMPDIR]}/$uuid.hiSeq_{wildcards.capDesign}.bed.split
rm {config[TMPDIR]}/$uuid.hiSeq_{wildcards.capDesign}.bed
for file in `ls {config[TMPDIR]}/$uuid.hiSeq_{wildcards.capDesign}.bed.split*`; do
echo "cat $file | awk -f ~/julien_utils/bed12fields2gff.awk > $file.gff; sort -T {config[TMPDIR]}  -T {config[TMPDIR]} -k1,1 -k4,4n -k5,5n $file.gff | makeIntrons.pl - | extract_intron_strand_motif.pl - {input.genome} {config[TMPDIR]}/$(basename $file); rm $file $file.gff $file.transcripts.tsv"
done > {config[TMPDIR]}/$uuid.parallelIntrons.sh
echoerr "extracting introns on split files"

parallel -j {threads} < {config[TMPDIR]}/$uuid.parallelIntrons.sh
echoerr "getting SJs and merging into output..."
cat {config[TMPDIR]}/$uuid.hiSeq_{wildcards.capDesign}.bed.split*.introns.gff | perl -lane '$start=$F[3]-1; $end=$F[4]+1; print $F[0]."_".$start."_".$end."_".$F[6]' | sort -T {config[TMPDIR]}  -T {config[TMPDIR]} |uniq> {output}

		'''
