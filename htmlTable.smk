
rule filterSampleAnnot:
	input: config['SAMPLE_ANNOT']
	output: temp("output/tmp/{subProject}_samples.tsv")
	run:
		if wildcards.subProject == 'ALL':
			sampleAnnot.to_csv(output[0], sep="\t")
		else:
			sampleAnnotFiltered=sampleAnnot[sampleAnnot.subProject == wildcards.subProject]
			sampleAnnotFiltered.to_csv(output[0], sep="\t")

rule makeHtmlSummaryDashboard:
	input:
		sampleAnnotationFile="output/tmp/{subProject}_samples.tsv",
		allFastqTimeStamps="output/statsFiles/all.fastq.timestamps.tsv",
		allReadLengths="output/statsFiles/all.readlength.summary.tsv",
		allBasicMappingStats="output/statsFiles/all.basic.mapping.stats.tsv",
		allHissStats="output/statsFiles/all.HiSS.stats.tsv",
		allMergedStats="output/statsFiles/all.min{minReadSupport}reads.merged.stats.tsv",
		allMatureRnaLengthStats="output/statsFiles/all.min{minReadSupport}reads.matureRNALengthSummary.stats.tsv",
		allTmergeVsSirvStats="output/statsFiles/all.HiSS.tmerge.min{minReadSupport}reads.vs.SIRVs.stats.tsv",
		allCagePolyASupportStats="output/statsFiles/all.min{minReadSupport}reads.splicing_status-all.cagePolyASupport.stats.tsv",
		allNovelLociStats="output/statsFiles/all.tmerge.min{minReadSupport}reads.endSupport-all.novelLoci.stats.tsv",
		allNovelFlLociStats="output/statsFiles/all.tmerge.min{minReadSupport}reads.endSupport-cagePolyASupported.novelLoci.stats.tsv",
		allNovelLociQcStats="output/statsFiles/all.tmerge.min{minReadSupport}reads.endSupport-all.novelLoci.qc.stats.tsv",
		allNovelFlLociQcStats="output/statsFiles/all.tmerge.min{minReadSupport}reads.endSupport-cagePolyASupported.novelLoci.qc.stats.tsv",
		allNtCoverageStats="output/statsFiles/all.tmerge.min{minReadSupport}reads.endSupport-all.vs.ntCoverageByGenomePartition.stats.tsv",
	conda: "envs/R_env.yml"
	params: indexEntryPart= lambda wildcards: "(no filter, all samples)" if wildcards.subProject == 'ALL' else "sub-project"
	output:
		html="output/html/summary_table_min{minReadSupport}reads_{subProject}.html",
		tsv="output/html/summary_table_min{minReadSupport}reads_{subProject}.tsv",
		index=temp("output/html/summary_table_min{minReadSupport}reads_{subProject}.index.tmp.html")
	shell:
		'''
./LyRic/makeHtmlDashboard.r {output.html} {input.sampleAnnotationFile} {input.allFastqTimeStamps} {input.allReadLengths} {input.allBasicMappingStats} {input.allHissStats} {input.allMergedStats} {input.allMatureRnaLengthStats} {input.allTmergeVsSirvStats} {input.allCagePolyASupportStats} {input.allNovelLociStats} {input.allNovelFlLociStats} {input.allNovelLociQcStats} {input.allNovelFlLociQcStats} {input.allNtCoverageStats}

htmlBn=$(basename {output.html})
tsvBn=$(basename {output.tsv})

echo "<li> <b>{wildcards.subProject}</b> {params.indexEntryPart}: <a href='$htmlBn'>HTML</a> / <a href='$tsvBn'>TSV</a></li>" > {output.index}
		'''

rule makeHtmlSummaryDashboardIndex:
	input:  lambda wildcards: expand("output/html/summary_table_min{minReadSupport}reads_{subProject}.index.tmp.html", subProject=subProjects, minReadSupport='|'.join(MINIMUM_TMERGE_READ_SUPPORT))
	output: "output/html/index.html"
	shell:
		'''
uuid=$(uuidgen)
# print start of html:
printf '<!DOCTYPE html\>
' > {TMPDIR}/$uuid.html
printf "<html>
<head>
<title>Summary statistics tables for {config[PROJECT_NAME]} project</title>
</head><body>
<h1>Summary statistics tables for {config[PROJECT_NAME]} project</h1>
<ul> " >> {TMPDIR}/$uuid.html

# print list:
cat {input} | sort >> {TMPDIR}/$uuid.html 

# print end of html:
printf '</ul>
'  >> {TMPDIR}/$uuid.html 

printf "<br>Produced with <a href="https://github.com/julienlag/LyRic/">LyRic</a>.<br>" >> {TMPDIR}/$uuid.html 
date=$(date)
printf "<br>(Last updated $date)"  >> {TMPDIR}/$uuid.html 
printf '</body>
</html> ' >> {TMPDIR}/$uuid.html 

mv {TMPDIR}/$uuid.html {output}
		'''