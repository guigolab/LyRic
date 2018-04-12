rule makeIntrons:
	input: "mappings/" + "readBedToGff/{techname}_{capDesign}_{barcodes}.merged.gff"
	output: "mappings/" + "makeIntrons/{techname}_{capDesign}_{barcodes}.introns.gff.gz"
	shell:
		'''
cat {input} | sort -k12,12 -k4,4n -k5,5n | awk -v fldgn=10 -v fldtr=12 -f ~/julien_utils/make_introns.awk | sortgff |gzip> {output}
		'''

rule getIntronMotif:
	input:
		introns = "mappings/" + "makeIntrons/{techname}_{capDesign}_{barcodes}.introns.gff.gz",
		genome = lambda wildcards: GENOMESDIR + CAPDESIGNTOGENOME[wildcards.capDesign] + ".fa"
	output:
		gff = "mappings/" + "getIntronMotif/{techname}_{capDesign}_{barcodes}.introns.gff",
		tsv = "mappings/" + "getIntronMotif/{techname}_{capDesign}_{barcodes}.transcripts.tsv"
	shell:
		'''
zcat {input.introns} | grep -vP "^ERCC"| extract_intron_strand_motif.pl - {input.genome} $(dirname {output.gff})/$(basename {output.gff} .introns.gff)

		'''
