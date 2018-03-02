rule polyAmapping:
	input:
		reads = "mappings/" + "mergeSizeFracBams/{techname}_{capDesign}_{barcodes}.merged.bam",
		genome = lambda wildcards: "/users/rg/jlagarde/genomes/" + CAPDESIGNTOGENOME[wildcards.capDesign] + ".fa"

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


rule clusterPolyAsites:
	input: "mappings/" + "removePolyAERCCs/{techname}_{capDesign}_{barcodes}.polyAsitesNoErcc.bed"
	output: "mappings/" + "clusterPolyAsites/{techname}_{capDesign}_{barcodes}.polyAsites.clusters.bed"
	shell:
		'''
cat {input} | $HOME/bin/bedtools2/bin/bedtools merge -s -d 5 -c 4 -o collapse -i stdin | awk '{{print $1"\t"$2"\t"$3"\t"$5"\t0\t"$4}}'| perl -F"\\t" -lane 'if($F[5] eq "+"){{$F[1]=$F[2]-1}}elsif($F[5] eq "-"){{$F[2]=$F[1]+1}}else{{die}} print join("\t",@F);'|sortbed > {output}
		'''
rule makePolyABigWigs:
	input:
		sites = "mappings/" + "removePolyAERCCs/{techname}_{capDesign}_{barcodes}.polyAsitesNoErcc.bed",
		genome = lambda wildcards: "/users/rg/jlagarde/genomes/" + CAPDESIGNTOGENOME[wildcards.capDesign] + ".genome"
	output: "mappings/" + "makePolyABigWigs/{techname}_{capDesign}_{barcodes}.polyAsitesNoErcc.{strand}.bw"
	shell:
		'''
genomeCoverageBed -strand {wildcards.strand} -split -bg -i {input.sites} -g {input.genome} > {output}.bedgraph
bedGraphToBigWig {output}.bedgraph {input.genome} {output}
rm -f {output}.bedgraph
		'''
