rule make_hubtxt:
	output: config["TRACK_HUB_DIR"] + "hub.txt"
	shell:
		'''
echo "hub {config[PROJECT_NAME]}
shortLabel {config[PROJECT_LONG_NAME]}
longLabel {config[PROJECT_LONG_NAME]}
genomesFile genomes.txt
email {config[PROJECT_CONTACT_EMAIL]}
" > {output}

		'''

rule make_genomestxt:
	output:temp(config["TRACK_HUB_DIR"] + "tmp/{genome}.txt")
	shell:
		'''

echo "genome {wildcards.genome}
trackDb {wildcards.genome}/trackDb.txt
" > {output}
		'''

rule agg_genomestxt:
	input: expand(config["TRACK_HUB_DIR"] + "tmp/{genome}.txt", genome=GENOMES)
	output: config["TRACK_HUB_DIR"] + "genomes.txt"
	shell:
		'''
cat {input} > {output}
		'''

rule makeCapDesignSupertrackHeader:
	output: temp(config["TRACK_HUB_DIR"] + "{capDesign}_SupertrackHeader.txt")
	shell:
		'''
echo "track {wildcards.capDesign}
superTrack on show
shortLabel {config[PROJECT_NAME]} {wildcards.capDesign} design
longLabel {config[PROJECT_NAME]} {wildcards.capDesign} capture design
visibility full
" > {output}

		'''

rule makeTrackDbInputHeader:
	output: temp(config["TRACK_HUB_DIR"] + "{capDesign}_trackDbInputHeader.txt")
	shell:
		'''
echo "track {wildcards.capDesign}_input
compositeTrack on
type bigBed 12
shortLabel {config[PROJECT_NAME]} {wildcards.capDesign} Targets
longLabel {config[PROJECT_NAME]} Targets (Captured regions in {wildcards.capDesign} design)
priority 9
parent {wildcards.capDesign}
visibility dense
color 110,30,0
" > {output}
		'''






rule makeTrackDbInputTracks:
	input:
		gff=config["TARGETSDIR"] + "{capDesign}_primary_targets.exons.reduced.gene_type.segments.gtf",
		genome=lambda wildcards: config["GENOMESDIR"] + CAPDESIGNTOGENOME[wildcards.capDesign] + ".genome"
	output:
		bigBed=config["TRACK_HUB_DIR"] + "dataFiles/{capDesign}_primary_targets.bb",
		trackDb=temp(config["TRACK_HUB_DIR"] + "{capDesign}_primary_targets.trackDb.txt")
	shell:
		'''
uuid=$(uuidgen)
uuid2=$(uuidgen)
uuid3=$(uuidgen)
cat {input.genome} | cut -f1 | sort |uniq > $TMPDIR/$uuid2

cat {input.gff} | grep -P "^chr" |gff2bed_full.pl -|sortbed > $TMPDIR/$uuid3
join  -j1 $TMPDIR/$uuid2 $TMPDIR/$uuid3 |ssv2tsv > $TMPDIR/$uuid
bedToBigBed -as={config[BED12_AUTOSQL]} -type=bed12 -extraIndex=name $TMPDIR/$uuid {input.genome} {output.bigBed}

echo -e "
\ttrack {wildcards.capDesign}_targets
\tbigDataUrl {TRACK_HUB_DATA_URL}$(basename {output.bigBed})
\tshortLabel {config[PROJECT_NAME]} {wildcards.capDesign} targets
\tlongLabel {config[PROJECT_NAME]} Targeted regions in {wildcards.capDesign} capture design
\ttype bigBed 12
\titemRgb on
\tsearchIndex name
\tparent {wildcards.capDesign}_input
\tvisibility dense

" > {output.trackDb}

		'''

if config["USE_MATCHED_ILLUMINA"]:
	rule makeTrackDbHiSeqOutputHeader:
		output: temp(config["TRACK_HUB_DIR"] + "{capDesign}_trackDbHiSeqOutputHeader.txt")
		shell:
			'''
echo "track {wildcards.capDesign}_hiSeq_output
compositeTrack on
type bam
shortLabel {config[PROJECT_NAME]} {wildcards.capDesign} Illumina HiSeq sequencing output
longLabel {config[PROJECT_NAME]} {wildcards.capDesign} Illumina HiSeq sequencing output
priority 9
parent {wildcards.capDesign}
visibility hide
" > {output}
			'''

if config["USE_MATCHED_ILLUMINA"]:
	rule makeHiSeqBamOutputTracks:
		input: "mappings/hiSeq_{capDesign}.bam" if config["DEMULTIPLEX"] else "mappings/hiSeq_{techname}_{capDesign}.{barcodes}.bam"
		output:
			bam=config["TRACK_HUB_DIR"] + "dataFiles/hiSeq_{capDesign}.bam" if config["DEMULTIPLEX"] else config["TRACK_HUB_DIR"] + "dataFiles/hiSeq_{techname}_{capDesign}.{barcodes}.bam",
			bai=config["TRACK_HUB_DIR"] + "dataFiles/hiSeq_{capDesign}.bam.bai" if config["DEMULTIPLEX"] else config["TRACK_HUB_DIR"] + "dataFiles/hiSeq_{techname}_{capDesign}.{barcodes}.bam.bai"
		shell:
			'''
uuid=$(uuidgen)
samtools view -H {input} > $TMPDIR/$uuid
samtools view {input} | grep -P "\tchr" >> $TMPDIR/$uuid
samtools view -b $TMPDIR/$uuid > {output.bam}
sleep 60s
samtools index {output.bam}

			'''

if config["USE_MATCHED_ILLUMINA"]:
	rule makeHiSeqBamOutputTrackDb:
		input: config["TRACK_HUB_DIR"] + "dataFiles/hiSeq_{capDesign}.bam" if config["DEMULTIPLEX"] else config["TRACK_HUB_DIR"] + "dataFiles/hiSeq_{techname}_{capDesign}.{barcodes}.bam",
		output:
			trackDb=temp(config["TRACK_HUB_DIR"] + "hiSeq_{capDesign}_BamOutput.trackDb.txt" if config["DEMULTIPLEX"] else config["TRACK_HUB_DIR"] +"hiSeq_{techname}_{capDesign}.{barcodes}_BamOutput.trackDb.txt")
		params:
#			track= lambda wildcards: "hiSeq_{wildcards.capDesign}_BAM" if config["DEMULTIPLEX"] else "hiSeq_{wildcards.techname}_{wildcards.capDesign}.{wildcards.barcodes}_BAM",
#			shortLabel=lambda wildcards: config["PROJECT_NAME"] + " Illumina HiSeq ({wildcards.capDesign}) BAMs" if config["DEMULTIPLEX"] else config["PROJECT_NAME"] + "Illumina HiSeq ({wildcards.techname} {wildcards.capDesign} {wildcards.barcodes}) BAMs"
			track="hiSeq_{capDesign}_BAM" if config["DEMULTIPLEX"] else "hiSeq_{techname}_{capDesign}.{barcodes}_BAM",
			shortLabel=config["PROJECT_NAME"] + " Illumina HiSeq ({capDesign}) BAMs" if config["DEMULTIPLEX"] else config["PROJECT_NAME"] + " Illumina HiSeq ({techname} {capDesign} {barcodes}) BAMs"
		shell:
			'''
echo -e "
\ttrack {params.track}
\tbigDataUrl {TRACK_HUB_DATA_URL}$(basename {input})
\tshortLabel {params.shortLabel}
\tlongLabel {params.shortLabel}
\ttype bam
\tbamColorMode gray
\tbamGrayMode aliQual
\tdoWiggle on
\tindelPolyA on
\talwaysZero on
\tmaxHeightPixels 100:35:11
\tparent {wildcards.capDesign}_hiSeq_output

" > {output.trackDb}

			'''


rule makeTrackDbTechnameHeader:
	params:
		subGroupString=lambda wildcards: trackHubSubGroupString(wildcards.techname, wildcards.capDesign, SIZEFRACSpluSMERGED,BARCODESpluSMERGED, FINALCORRECTIONLEVELS, config["MINIMUM_TMERGE_READ_SUPPORT"]),
		viewsString="allTMs=All_TMs flTMs=Full-length_TMs BAMs=BAMs"
	output: temp(config["TRACK_HUB_DIR"] + "{techname}_{capDesign}_trackDbOutputHeader.txt")
	shell:
		'''
echo "track {wildcards.techname}_{wildcards.capDesign}_output
compositeTrack on
type bigBed 12
shortLabel {config[PROJECT_NAME]} {wildcards.capDesign} {wildcards.techname} sequencing output
longLabel {config[PROJECT_NAME]} {wildcards.techname} sequencing output in {wildcards.capDesign} design
itemRgb on
priority 10
parent {wildcards.capDesign}
{params.subGroupString}
subGroup4 filter Filter {params.viewsString}
dimensions dimX=sample dimY=filter dimA=sizeFraction dimB=corrLevel dimC=minReadSupport
filterComposite dimA dimB dimC
visibility squish
" > {output}
		'''

rule makeBamOutputTracks:
	input: "mappings/readMapping/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bam"
	output:
		bam=config["TRACK_HUB_DIR"] + "dataFiles/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bam",
		bai=config["TRACK_HUB_DIR"] + "dataFiles/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bam.bai"
	shell:
		'''
uuid=$(uuidgen)
samtools view -H {input} > $TMPDIR/$uuid
samtools view {input} | grep -P "\tchr" >> $TMPDIR/$uuid
samtools view -b $TMPDIR/$uuid > {output.bam}
sleep 60s
samtools index {output.bam}

		'''


rule makeBamOutputTrackDb:
	input: config["TRACK_HUB_DIR"] + "dataFiles/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bam"
	output:
		trackDb=temp(config["TRACK_HUB_DIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}_BamOutput.trackDb.txt")
	params:
		trackOnOff=lambda wildcards: trackChecked(wildcards.techname, wildcards.corrLevel, wildcards.capDesign, wildcards.sizeFrac, wildcards.barcodes, None)
	shell:
		'''
echo -e "
\ttrack {wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.barcodes}_BAM
\tbigDataUrl {TRACK_HUB_DATA_URL}$(basename {input})
\tshortLabel {config[PROJECT_NAME]} {wildcards.techname} Corr{wildcards.corrLevel} {wildcards.capDesign} {wildcards.sizeFrac} {wildcards.barcodes} BAMs
\tlongLabel {config[PROJECT_NAME]} BAM (read alignments) (Protocol: {wildcards.techname}, Correction: {wildcards.corrLevel}, Capture design: {wildcards.capDesign}, Size fraction: {wildcards.sizeFrac}, Sample: {wildcards.barcodes})
\ttype bam
\tbamColorMode gray
\tbamGrayMode aliQual
\tdoWiggle on
\tindelPolyA on
\talwaysZero on
\tmaxHeightPixels 100:35:11
\tparent {wildcards.techname}_{wildcards.capDesign}_output {params.trackOnOff}
\tsubGroups filter=BAMs sample={wildcards.barcodes} protocol={wildcards.techname} sizeFraction={wildcards.sizeFrac} corrLevel={wildcards.corrLevel}

" > {output.trackDb}

		'''





rule makeTmergeOutputTracks:
	input:
		bed="mappings/nonAnchoredMergeReads/colored/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.all.bed",
		genome=lambda wildcards: config["GENOMESDIR"] + CAPDESIGNTOGENOME[wildcards.capDesign] + ".genome"
	output:
		bigBed=config["TRACK_HUB_DIR"] + "dataFiles/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.bb",
		trackDb=temp(config["TRACK_HUB_DIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}_TmergeOutput.min{minReadSupport}reads.trackDb.txt")
	params:
		trackOnOff=lambda wildcards: trackChecked(wildcards.techname, wildcards.corrLevel, wildcards.capDesign, wildcards.sizeFrac, wildcards.barcodes, wildcards.minReadSupport)
	shell:
		'''
uuid=$(uuidgen)

cat {input.bed} | grep -P "^chr"  > $TMPDIR/$uuid
bedToBigBed -as={config[BED12_AUTOSQL]} -type=bed12 -extraIndex=name $TMPDIR/$uuid {input.genome} {output.bigBed}

echo -e "
\ttrack {wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.barcodes}_{wildcards.minReadSupport}_TMs
\tbigDataUrl {TRACK_HUB_DATA_URL}$(basename {output.bigBed})
\tshortLabel {config[PROJECT_NAME]} {wildcards.techname} Corr{wildcards.corrLevel} {wildcards.capDesign} {wildcards.sizeFrac} {wildcards.barcodes} min{wildcards.minReadSupport}reads TMs
\tlongLabel {config[PROJECT_NAME]} Transcript models (Protocol: {wildcards.techname}, Correction: {wildcards.corrLevel}, Capture design: {wildcards.capDesign}, Size fraction: {wildcards.sizeFrac}, Sample: {wildcards.barcodes}, Min. read support per TM: {wildcards.minReadSupport})
\ttype bigBed 12
\titemRgb on
\tsearchIndex name
\tparent {wildcards.techname}_{wildcards.capDesign}_output {params.trackOnOff}
\tsubGroups filter=allTMs sample={wildcards.barcodes} protocol={wildcards.techname} sizeFraction={wildcards.sizeFrac} corrLevel={wildcards.corrLevel} minReadSupport={wildcards.minReadSupport}

" > {output.trackDb}

		'''

rule makeFLTmergeOutputTracks:
	input:
		bed="mappings/nonAnchoredMergeReads/colored/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.cage+polyASupported.bed",
		genome=lambda wildcards: config["GENOMESDIR"] + CAPDESIGNTOGENOME[wildcards.capDesign] + ".genome"
	output:
		bigBed=config["TRACK_HUB_DIR"] + "dataFiles/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}_FL.tmerge.min{minReadSupport}reads.bb",
		trackDb=temp(config["TRACK_HUB_DIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}_FLTmergeOutput.min{minReadSupport}reads.trackDb.txt")
	params:
		trackOnOff=lambda wildcards: trackChecked(wildcards.techname, wildcards.corrLevel, wildcards.capDesign, wildcards.sizeFrac, wildcards.barcodes, wildcards.minReadSupport)
	shell:
		'''
uuid=$(uuidgen)

cat {input.bed} | grep -P "^chr"  > $TMPDIR/$uuid
bedToBigBed -as={config[BED12_AUTOSQL]} -type=bed12 -extraIndex=name $TMPDIR/$uuid {input.genome} {output.bigBed}

echo -e "
\ttrack {wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.barcodes}_{wildcards.minReadSupport}_FL_TMs
\tbigDataUrl {TRACK_HUB_DATA_URL}$(basename {output.bigBed})
\tshortLabel {config[PROJECT_NAME]} {wildcards.techname} Corr{wildcards.corrLevel} {wildcards.capDesign} {wildcards.sizeFrac} {wildcards.barcodes} min{wildcards.minReadSupport}reads FL TMs
\tlongLabel {config[PROJECT_NAME]} Full-length transcript models (Protocol: {wildcards.techname}, Correction: {wildcards.corrLevel}, Capture design: {wildcards.capDesign}, Size fraction: {wildcards.sizeFrac}, Sample: {wildcards.barcodes}, Min. read support per TM: {wildcards.minReadSupport})
\ttype bigBed 12
\titemRgb on
\tsearchIndex name
\tparent {wildcards.techname}_{wildcards.capDesign}_output {params.trackOnOff}
\tsubGroups filter=flTMs sample={wildcards.barcodes} protocol={wildcards.techname} sizeFraction={wildcards.sizeFrac} corrLevel={wildcards.corrLevel} minReadSupport={wildcards.minReadSupport}

" > {output.trackDb}

		'''




rule catTrackDb:
	input:
		lambda wildcards: expand(config["TRACK_HUB_DIR"] + "{capDesign}_SupertrackHeader.txt", capDesign=GENOMETOCAPDESIGNS[wildcards.genome]),
		lambda wildcards: expand(config["TRACK_HUB_DIR"] + "{capDesign}_trackDbInputHeader.txt", capDesign=GENOMETOCAPDESIGNS[wildcards.genome]),
		lambda wildcards: expand(config["TRACK_HUB_DIR"] + "{capDesign}_trackDbHiSeqOutputHeader.txt", capDesign=GENOMETOCAPDESIGNS[wildcards.genome]),
		lambda wildcards: expand(config["TRACK_HUB_DIR"] + "{capDesign}_primary_targets.trackDb.txt", capDesign=GENOMETOCAPDESIGNS[wildcards.genome]),
		lambda wildcards: expand(config["TRACK_HUB_DIR"] + "{techname}_{capDesign}_trackDbOutputHeader.txt", filtered_product, techname=TECHNAMES, capDesign=GENOMETOCAPDESIGNS[wildcards.genome]),
		lambda wildcards: expand(config["TRACK_HUB_DIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}_TmergeOutput.min{minReadSupport}reads.trackDb.txt", filtered_product_merge, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=GENOMETOCAPDESIGNS[wildcards.genome], sizeFrac=SIZEFRACSpluSMERGED, barcodes=BARCODESpluSMERGED, minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"]),
		lambda wildcards: expand(config["TRACK_HUB_DIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}_FLTmergeOutput.min{minReadSupport}reads.trackDb.txt", filtered_product_merge, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=GENOMETOCAPDESIGNS[wildcards.genome], sizeFrac=SIZEFRACSpluSMERGED, barcodes=BARCODESpluSMERGED, minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"]),
		lambda wildcards: expand(config["TRACK_HUB_DIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}_BamOutput.trackDb.txt", filtered_product_merge, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=GENOMETOCAPDESIGNS[wildcards.genome], sizeFrac=SIZEFRACSpluSMERGED, barcodes=BARCODESpluSMERGED),
		lambda wildcards: expand(config["TRACK_HUB_DIR"] + "hiSeq_{capDesign}_BamOutput.trackDb.txt", capDesign=GENOMETOCAPDESIGNS[wildcards.genome]) if config["USE_MATCHED_ILLUMINA"] and config["DEMULTIPLEX"] else expand(config["TRACK_HUB_DIR"] + "hiSeq_{techname}_{capDesign}.{barcodes}_BamOutput.trackDb.txt", filtered_product, techname=TECHNAMES, capDesign=GENOMETOCAPDESIGNS[wildcards.genome], barcodes=BARCODES) if config["USE_MATCHED_ILLUMINA"] and not config["DEMULTIPLEX"] else None
	output:
		config["TRACK_HUB_DIR"] + "{genome}/trackDb.txt"
	shell:
		'''
cat {input} > {output}
		'''

