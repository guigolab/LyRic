rule hiSeqReadMapping:
	input:
		reads = "fastqs/" + "hiSeq_{capDesign}.fastq.gz",
		genome = lambda wildcards: "/users/rg/jlagarde/genomes/STARshort_indices/" + CAPDESIGNTOGENOME[wildcards.capDesign] + "/"
	params:
		referenceAnnot = lambda wildcards: CAPDESIGNTOANNOTGTF[wildcards.capDesign]
	threads: 12
	output:
		"mappings/" + "hiSeq_{capDesign}.bam"
	shell:
		'''
echoerr "Mapping"
mkdir -p mappings/STAR/`basename {output}`/
STAR \
--runThreadN {threads} \
--readFilesIn {input.reads} \
--genomeDir {input.genome} \
--readFilesCommand zcat \
--sjdbOverhang 124 \
--sjdbGTFfile {params.referenceAnnot} \
--outFileNamePrefix "mappings/STAR/`basename {output}`/" \
--outStd SAM \
--genomeLoad NoSharedMemory \
--outSAMunmapped Within \
--outFilterType BySJout \
-- outFilterMultimapNmax 20 \
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


