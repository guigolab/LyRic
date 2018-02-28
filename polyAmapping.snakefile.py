rule polyAmapping:
	input:
		reads = "mappings/" + "{techname}_{capDesign}_{barcodesU}.merged.bam",
		genome = lambda wildcards: "/users/rg/jlagarde/genomes/" + CAPDESIGNTOGENOME[wildcards.capDesign] + ".fa"

	output: "mappings/polyA/" + "{techname}_{capDesign}_{barcodesU}.polyAsites.bed.gz"
	shell:
		'''
samtools view {input.reads} | samToPolyA.pl --minClipped=20 --minAcontent=0.9 --discardInternallyPrimed --minUpMisPrimeAlength=10 --genomeFasta={input.genome} - |sortbed |gzip > {output}
		'''

rule removePolyAERCCs:
	input: "mappings/polyA/" + "{techname}_{capDesign}_{barcodesU}.polyAsites.bed.gz"
	output: "mappings/polyA/" + "{techname}_{capDesign}_{barcodesU}.polyAsitesNoErcc.bed"
	shell:
		'''
zcat {input} | fgrep -v ERCC > {output}
		'''

rule getPolyAreadsList:
	input: "mappings/polyA/" + "{techname}_{capDesign}_{barcodesU}.polyAsitesNoErcc.bed"
	output: "mappings/polyA/" + "{techname}_{capDesign}_{barcodesU}.polyAreads.list"
	shell:
		'''
cat {input} | cut -f4 | sort|uniq > {output}
		'''


rule clusterPolyAsites:
	input: "mappings/polyA/" + "{techname}_{capDesign}_{barcodesU}.polyAsitesNoErcc.bed"
	output: "mappings/polyA/" + "{techname}_{capDesign}_{barcodesU}.polyAsites.clusters.bed"
	shell:
		'''
cat {input} | bedtools merge -s -c 4 -d 5 -o collapse -i stdin | awk '{{print $1"\t"$2"\t"$3"\t"$5"\t0\t"$4}}'| perl -F"\t" -lane 'if($F[5] eq "+"){{$F[1]=$F[2]-1}}elsif($F[5] eq "-"){{$F[2]=$F[1]+1}}else{{die}} print join("\t",@F);'|sortbed > {output}
		'''

rule makePolyABigWigs:
	input:
		sites = "mappings/polyA/" + "{techname}_{capDesign}_{barcodesU}.polyAsitesNoErcc.bed",
		genome = lambda wildcards: "/users/rg/jlagarde/genomes/" + CAPDESIGNTOGENOME[wildcards.capDesign] + ".genome"
	output: "mappings/polyA/" + "{techname}_{capDesign}_{barcodesU}.polyAsitesNoErcc.{strand}.bw"
	shell:
		'''
genomeCoverageBed -strand {wildcards.strand} -split -bg -i {input.sites} -g {input.genome} > {output}.bedgraph
bedGraphToBigWig {output}.bedgraph {input.genome} {output}
rm -f {output}.bedgraph
		'''
