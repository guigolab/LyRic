rule extractPooledTmsFivepEnds:
	input: "mappings/" + "nonAnchoredMergeReads/pooled/{techname}_{capDesign}_pooled.tmerge.bed"
	output: "mappings/" + "nonAnchoredMergeReads/pooled/5pEnds/{techname}_{capDesign}_pooled.tmerge.5pEnds.bed"
	shell:
		'''
cat {input} |extractTranscriptEndsFromBed12.pl 5 |sortbed> {output}
		'''

rule cageSupportedfivepEnds:
	input:
		fivePends="mappings/" + "nonAnchoredMergeReads/pooled/5pEnds/{techname}_{capDesign}_pooled.tmerge.5pEnds.bed",
		tms="mappings/" + "nonAnchoredMergeReads/pooled/{techname}_{capDesign}_pooled.tmerge.bed",
		cagePeaks=lambda wildcards: CAPDESIGNTOCAGEPEAKS[wildcards.capDesign]
	params: genome = lambda wildcards: config["GENOMESDIR"] + CAPDESIGNTOGENOME[wildcards.capDesign] + ".genome"
	output: "mappings/" + "nonAnchoredMergeReads/pooled/cageSupported/{techname}_{capDesign}_pooled.tmerge.cageSupported.bed"
	shell:
		'''
cat {input.fivePends} | sortbed | bedtools slop -s -l 50 -r 50 -i stdin -g {params.genome} | bedtools intersect -u -s -a stdin -b {input.cagePeaks} | cut -f4 | fgrep -w -f - {input.tms} > {output}
		'''