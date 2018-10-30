rule hiSeqReadMapping:
	input:
		reads1 = config["MATCHED_HISEQ_PATH"] + "hiSeq_{capDesign}_1.fastq.gz" if config["DEMULTIPLEX"] else lambda wildcards: expand(config["MATCHED_HISEQ_PATH"] + "{techname}_{capDesign}.{barcodes}_1.fastq.gz", filtered_product, techname=wildcards.techname, capDesign=wildcards.capDesign, barcodes=wildcards.barcodes),
		reads2 = config["MATCHED_HISEQ_PATH"] + "hiSeq_{capDesign}_2.fastq.gz" if config["DEMULTIPLEX"] else lambda wildcards: expand(config["MATCHED_HISEQ_PATH"] + "{techname}_{capDesign}.{barcodes}_2.fastq.gz", filtered_product, techname=wildcards.techname, capDesign=wildcards.capDesign, barcodes=wildcards.barcodes),
		genome = lambda wildcards: config["GENOMESDIR"] + "STARshort_indices/" + CAPDESIGNTOGENOME[wildcards.capDesign] + "/",
		referenceAnnot = lambda wildcards: CAPDESIGNTOANNOTGTF[wildcards.capDesign]
	threads: 12
	output:
		"mappings/" + "hiSeq_{capDesign}.bam" if config["DEMULTIPLEX"] else "mappings/" + "hiSeq_{techname}_{capDesign}.{barcodes}.bam"
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
| samtools view -b -u -S - | samtools sort -@ 2 -T $TMPDIR  -m 15000000000 - > {output}
echoerr "Mapping done"

		'''

if not config["DEMULTIPLEX"]:
	rule mergeHiSeqBarcodeBams:
		input: lambda wildcards: expand("mappings/" + "hiSeq_{techname}_{capDesign}.{barcodes}.bam", filtered_product, techname=wildcards.techname, capDesign=wildcards.capDesign, barcodes=BARCODES)
		output: "mappings/" + "hiSeq_{techname}_{capDesign}.allTissues.bam"
		shell:
			'''
samtools merge {output} {input}
sleep 120s
samtools index {output}
			'''


rule getHiSeqCanonicalIntronsList:
	input:
		bam="mappings/" + "hiSeq_{capDesign}.bam" if config["DEMULTIPLEX"] else "mappings/" + "hiSeq_{techname}_{capDesign}.{barcodes}.bam",
		genome = lambda wildcards: config["GENOMESDIR"] + CAPDESIGNTOGENOME[wildcards.capDesign] + ".fa"
	threads: 6
	output:"mappings/hiSeqIntrons/" + "hiSeq_{capDesign}.canonicalIntrons.list" if config["DEMULTIPLEX"] else "mappings/hiSeqIntrons/hiSeq_{techname}_{capDesign}.{barcodes}.canonicalIntrons.list"
	shell:
		'''
echoerr "making bed"
samtools view -b -F 256 -F4 -F 2048 {input.bam}  |bedtools bamtobed -i stdin -bed12 | fgrep -v ERCC- > $TMPDIR/hiSeq_{wildcards.capDesign}.bed
echoerr "splitting"
split -a 3 -d -e -n l/24 $TMPDIR/hiSeq_{wildcards.capDesign}.bed $TMPDIR/hiSeq_{wildcards.capDesign}.bed.split
rm $TMPDIR/hiSeq_{wildcards.capDesign}.bed
for file in `ls $TMPDIR/hiSeq_{wildcards.capDesign}.bed.split*`; do
echo "cat $file | awk -f ~/julien_utils/bed12fields2gff.awk > $file.gff; sort -T $TMPDIR -k12,12 -k4,4n -k5,5n $file.gff | awk -v fldgn=10 -v fldtr=12 -f ~/julien_utils/make_introns.awk | extract_intron_strand_motif.pl - {input.genome} $TMPDIR/$(basename $file); rm $file $file.gff $file.transcripts.tsv"
done > $TMPDIR/parallelIntrons.sh
echoerr "extracting introns on split files"

parallel -j {threads} < $TMPDIR/parallelIntrons.sh
echoerr "getting SJs and merging into output..."
cat $TMPDIR/hiSeq_{wildcards.capDesign}.bed.split*.introns.gff | perl -lane '$start=$F[3]-1; $end=$F[4]+1; print $F[0]."_".$start."_".$end."_".$F[6]' | sort -T $TMPDIR |uniq> {output}

		'''

