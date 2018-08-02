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


rule extractPooledTmsThreepEnds:
	input: "mappings/" + "nonAnchoredMergeReads/pooled/{techname}_{capDesign}_pooled.tmerge.bed"
	output: "mappings/" + "nonAnchoredMergeReads/pooled/3pEnds/{techname}_{capDesign}_pooled.tmerge.3pEnds.bed"
	shell:
		'''
cat {input} |extractTranscriptEndsFromBed12.pl 3 |sortbed> {output}
		'''

rule polyASupportedthreepEnds:
	input:
		threePends="mappings/" + "nonAnchoredMergeReads/pooled/3pEnds/{techname}_{capDesign}_pooled.tmerge.3pEnds.bed",
		tms="mappings/" + "nonAnchoredMergeReads/pooled/{techname}_{capDesign}_pooled.tmerge.bed",
		polyAsites=lambda wildcards: expand("mappings/" + "removePolyAERCCs/{techname}_{capDesign}_{barcodes}.polyAsitesNoErcc.bed", filtered_product, techname=wildcards.techname, capDesign=wildcards.capDesign, barcodes=BARCODES)
	params: genome = lambda wildcards: config["GENOMESDIR"] + CAPDESIGNTOGENOME[wildcards.capDesign] + ".genome"
	output: "mappings/nonAnchoredMergeReads/pooled/polyASupported/{techname}_{capDesign}_pooled.tmerge.polyASupported.bed"
	shell:
		'''
cat {input.polyAsites} |sortbed > $TMPDIR/polyAsites.bed
cat {input.threePends} | sortbed | bedtools slop -s -l 5 -r 5 -i stdin -g {params.genome} | bedtools intersect -u -s -a stdin -b $TMPDIR/polyAsites.bed | cut -f4 | fgrep -w -f - {input.tms} > {output}
		'''

