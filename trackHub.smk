rule make_hubtxt:
	output: config["TRACK_HUB_DIR"] + "hub.txt"
	shell:
		'''
uuidTmpOut=$(uuidgen)
echo "hub {config[PROJECT_NAME]}
shortLabel {config[PROJECT_LONG_NAME]}
longLabel {config[PROJECT_LONG_NAME]}
genomesFile genomes.txt
email {config[PROJECT_CONTACT_EMAIL]}
" > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}

		'''

rule make_genomestxt:
	output:temp(config["TRACK_HUB_DIR"] + "tmp/{genome}.txt")
	shell:
		'''
uuidTmpOut=$(uuidgen)
echo "genome {wildcards.genome}
trackDb {wildcards.genome}/trackDb.txt
" > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
		'''

rule agg_genomestxt:
	input: expand(config["TRACK_HUB_DIR"] + "tmp/{genome}.txt", genome=GENOMES)
	output: config["TRACK_HUB_DIR"] + "genomes.txt"
	shell:
		'''
uuidTmpOut=$(uuidgen)
cat {input} > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
		'''

rule makeCapDesignSupertrackHeader:
	output: temp(config["TRACK_HUB_DIR"] + "{capDesign}_SupertrackHeader.txt")
	shell:
		'''
uuidTmpOut=$(uuidgen)
echo "track {wildcards.capDesign}
superTrack on show
shortLabel {config[PROJECT_NAME]} {wildcards.capDesign} design
longLabel {config[PROJECT_NAME]} {wildcards.capDesign} capture design
visibility full
" > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}

		'''

rule makeTrackDbInputHeader:
	output: temp(config["TRACK_HUB_DIR"] + "{capDesign}_trackDbInputHeader.txt")
	shell:
		'''
uuidTmpOut=$(uuidgen)
echo "track {wildcards.capDesign}_input
compositeTrack on
type bigBed 12
shortLabel {config[PROJECT_NAME]} {wildcards.capDesign} Targets
longLabel {config[PROJECT_NAME]} Targets (Captured regions in {wildcards.capDesign} design)
priority 9
parent {wildcards.capDesign}
visibility dense
color 110,30,0
" > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
		'''






rule makeTrackDbInputTracks:
	input:
		gff=config["TARGETSDIR"] + "{capDesign}_primary_targets.exons.reduced.gene_type.segments.gtf",
		genome=lambda wildcards: config["GENOMESDIR"] + CAPDESIGNTOGENOME[wildcards.capDesign] + ".sorted.genome"
	output:
		bigBed=config["TRACK_HUB_DIR"] + "dataFiles/{capDesign}_primary_targets.bb",
		trackDb=temp(config["TRACK_HUB_DIR"] + "{capDesign}_primary_targets.trackDb.txt")
	shell:
		'''
uuid=$(uuidgen)
uuid2=$(uuidgen)
uuid3=$(uuidgen)
uuidTmpOutB=$(uuidgen)
uuidTmpOutT=$(uuidgen)
cat {input.genome} | cut -f1 | sort -T {config[TMPDIR]}  |uniq > {config[TMPDIR]}/$uuid2

cat {input.gff} | grep -P "^chr" | grep -v "chrIS" |gff2bed_full.pl -|sort -T {config[TMPDIR]}  -k1,1 -k2,2n -k3,3n  > {config[TMPDIR]}/$uuid3
join  -j1 {config[TMPDIR]}/$uuid2 {config[TMPDIR]}/$uuid3 |ssv2tsv > {config[TMPDIR]}/$uuid
bedToBigBed -as={config[BED12_AUTOSQL]} -type=bed12 -extraIndex=name {config[TMPDIR]}/$uuid {input.genome} {config[TMPDIR]}/$uuidTmpOutB
mv {config[TMPDIR]}/$uuidTmpOutB {output.bigBed}

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

" > {config[TMPDIR]}/$uuidTmpOutT
mv {config[TMPDIR]}/$uuidTmpOutT {output.trackDb}

		'''

rule makeTrackDbHiSeqOutputHeader:
	output: temp(config["TRACK_HUB_DIR"] + "{capDesign}_trackDbHiSeqOutputHeader.txt")
	shell:
		'''
uuidTmpOut=$(uuidgen)
echo "track {wildcards.capDesign}_hiSeq_output
compositeTrack on
type bam
shortLabel {config[PROJECT_NAME]} {wildcards.capDesign} Illumina HiSeq sequencing output
longLabel {config[PROJECT_NAME]} {wildcards.capDesign} Illumina HiSeq sequencing output
priority 9
parent {wildcards.capDesign}
visibility hide
" > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
		'''

rule makeHiSeqBamOutputTracks:
	input: "mappings/hiSeq_{capDesign}.bam"
	output:
		bam=config["TRACK_HUB_DIR"] + "dataFiles/hiSeq_{capDesign}.bam",
		bai=config["TRACK_HUB_DIR"] + "dataFiles/hiSeq_{capDesign}.bam.bai"
	shell:
		'''
uuidTmpOut=$(uuidgen)
uuid=$(uuidgen)
samtools view -H {input} > {config[TMPDIR]}/$uuid
samtools view {input} | grep -P "\tchr" | grep -v "chrIS" >> {config[TMPDIR]}/$uuid
samtools view -b {config[TMPDIR]}/$uuid > {config[TMPDIR]}/$uuidTmpOut
sleep 200s
samtools index {config[TMPDIR]}/$uuidTmpOut

mv {config[TMPDIR]}/$uuidTmpOut {output.bam}
mv {config[TMPDIR]}/$uuidTmpOut.bai  {output.bai}
		'''

rule makeHiSeqBamOutputTrackDb:
	input: config["TRACK_HUB_DIR"] + "dataFiles/hiSeq_{capDesign}.bam",
	output:
		trackDb=temp(config["TRACK_HUB_DIR"] + "hiSeq_{capDesign}_BamOutput.trackDb.txt")
	params:
		track="hiSeq_{capDesign}_BAM" ,
		shortLabel=config["PROJECT_NAME"] + " Illumina HiSeq ({capDesign}) BAMs"
	shell:
		'''
uuidTmpOut=$(uuidgen)
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

" > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output.trackDb}

		'''


rule makeTrackDbTechnameHeader:
	params:
		subGroupString=lambda wildcards: trackHubSubGroupString(wildcards.techname, wildcards.capDesign, SIZEFRACS,BARCODES, config["MINIMUM_TMERGE_READ_SUPPORT"]),
		viewsString="allTMs=All_TMs flTMs=Full-length_TMs BAMs=BAMs"
	output: temp(config["TRACK_HUB_DIR"] + "{techname}_{capDesign}_trackDbOutputHeader.txt")
	shell:
		'''
uuidTmpOut=$(uuidgen)
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
dimensions dimX=sample dimY=filter dimA=sizeFraction dimB=minReadSupport
filterComposite dimA dimB
visibility squish
" > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
		'''

rule makeBamOutputTracks:
	input: "mappings/readMapping/{techname}_{capDesign}_{sizeFrac}_{barcodes}.bam"
	output:
		bam=config["TRACK_HUB_DIR"] + "dataFiles/{techname}_{capDesign}_{sizeFrac}_{barcodes}.bam",
		bai=config["TRACK_HUB_DIR"] + "dataFiles/{techname}_{capDesign}_{sizeFrac}_{barcodes}.bam.bai"
	shell:
		'''
uuidTmpOut=$(uuidgen)
uuid=$(uuidgen)
samtools view -H {input} > {config[TMPDIR]}/$uuid
samtools view {input} | grep -P "\\tchr" | grep -v "chrIS" >> {config[TMPDIR]}/$uuid
samtools view -b {config[TMPDIR]}/$uuid > {config[TMPDIR]}/$uuidTmpOut
sleep 200s
samtools index {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output.bam}
mv {config[TMPDIR]}/$uuidTmpOut.bai {output.bai}

		'''


rule makeBamOutputTrackDb:
	input: config["TRACK_HUB_DIR"] + "dataFiles/{techname}_{capDesign}_{sizeFrac}_{barcodes}.bam"
	output:
		trackDb=temp(config["TRACK_HUB_DIR"] + "{techname}_{capDesign}_{sizeFrac}_{barcodes}_BamOutput.trackDb.txt")

	shell:
		'''
uuidTmpOut=$(uuidgen)
echo -e "
\ttrack {wildcards.techname}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.barcodes}_BAM
\tbigDataUrl {TRACK_HUB_DATA_URL}$(basename {input})
\tshortLabel {config[PROJECT_NAME]} {wildcards.techname} {wildcards.capDesign} {wildcards.sizeFrac} {wildcards.barcodes} BAMs
\tlongLabel {config[PROJECT_NAME]} BAM (read alignments) (Protocol: {wildcards.techname}, Capture design: {wildcards.capDesign}, Size fraction: {wildcards.sizeFrac}, Sample: {wildcards.barcodes})
\ttype bam
\tbamColorMode gray
\tbamGrayMode aliQual
\tdoWiggle on
\tindelPolyA on
\talwaysZero on
\tmaxHeightPixels 100:35:11
\tparent {wildcards.techname}_{wildcards.capDesign}_output
\tsubGroups filter=BAMs sample={wildcards.barcodes} protocol={wildcards.techname} sizeFraction={wildcards.sizeFrac} 

" > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output.trackDb}

		'''





rule makeTmergeOutputTracks:
	input:
		bed="mappings/nonAnchoredMergeReads/colored/{techname}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.splicing_status:all.endSupport:all.bed",
		gff="mappings/nonAnchoredMergeReads/{techname}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.tmerge.min{minReadSupport}reads.splicing_status:all.endSupport:all.gff.gz",
		genome=lambda wildcards: config["GENOMESDIR"] + CAPDESIGNTOGENOME[wildcards.capDesign] + ".sorted.genome"
	output:
		bigBed=config["TRACK_HUB_DIR"] + "dataFiles/{techname}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.bb",
		gff=config["TRACK_HUB_DIR"] + "dataFiles/{techname}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.gff.gz",
		trackDb=temp(config["TRACK_HUB_DIR"] + "{techname}_{capDesign}_{sizeFrac}_{barcodes}_TmergeOutput.min{minReadSupport}reads.trackDb.txt")
	shell:
		'''
uuid=$(uuidgen)
uuidTmpOutB=$(uuidgen)
uuidTmpOutT=$(uuidgen)

cat {input.bed} | grep -P "^chr" | grep -v "chrIS"  > {config[TMPDIR]}/$uuid
bedToBigBed -as={config[BED12_AUTOSQL]} -type=bed12 -extraIndex=name {config[TMPDIR]}/$uuid {input.genome} {config[TMPDIR]}/$uuidTmpOutB
mv {config[TMPDIR]}/$uuidTmpOutB {output.bigBed}
cp -f {input.gff} {output.gff}
echo -e "
\ttrack {wildcards.techname}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.barcodes}_{wildcards.minReadSupport}_TMs
\tbigDataUrl {TRACK_HUB_DATA_URL}$(basename {output.bigBed})
\tshortLabel {config[PROJECT_NAME]} {wildcards.techname} {wildcards.capDesign} {wildcards.sizeFrac} {wildcards.barcodes} min{wildcards.minReadSupport}reads TMs
\tlongLabel {config[PROJECT_NAME]} Transcript models (Protocol: {wildcards.techname}, Capture design: {wildcards.capDesign}, Size fraction: {wildcards.sizeFrac}, Sample: {wildcards.barcodes}, Min. read support per TM: {wildcards.minReadSupport})
\ttype bigBed 12
\titemRgb on
\tsearchIndex name
\tparent {wildcards.techname}_{wildcards.capDesign}_output
\tsubGroups filter=allTMs sample={wildcards.barcodes} protocol={wildcards.techname} sizeFraction={wildcards.sizeFrac} minReadSupport={wildcards.minReadSupport}

" > {config[TMPDIR]}/$uuidTmpOutT
mv {config[TMPDIR]}/$uuidTmpOutT {output.trackDb}

		'''

rule makeFLTmergeOutputTracks:
	input:
		bed="mappings/nonAnchoredMergeReads/colored/{techname}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.splicing_status:all.endSupport:cagePolyASupported.bed",
		gff="mappings/nonAnchoredMergeReads/{techname}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.tmerge.min{minReadSupport}reads.splicing_status:all.endSupport:cagePolyASupported.gff.gz",
		genome=lambda wildcards: config["GENOMESDIR"] + CAPDESIGNTOGENOME[wildcards.capDesign] + ".sorted.genome"
	output:
		bigBed=config["TRACK_HUB_DIR"] + "dataFiles/{techname}_{capDesign}_{sizeFrac}_{barcodes}_FL.tmerge.min{minReadSupport}reads.bb",
		gff=config["TRACK_HUB_DIR"] + "dataFiles/{techname}_{capDesign}_{sizeFrac}_{barcodes}_FL.tmerge.min{minReadSupport}reads.gff.gz",
		trackDb=temp(config["TRACK_HUB_DIR"] + "{techname}_{capDesign}_{sizeFrac}_{barcodes}_FLTmergeOutput.min{minReadSupport}reads.trackDb.txt")
	shell:
		'''
uuid=$(uuidgen)
uuidTmpOutB=$(uuidgen)
uuidTmpOutT=$(uuidgen)

cat {input.bed} | grep -P "^chr" | grep -v "chrIS"  > {config[TMPDIR]}/$uuid
bedToBigBed -as={config[BED12_AUTOSQL]} -type=bed12 -extraIndex=name {config[TMPDIR]}/$uuid {input.genome} {config[TMPDIR]}/$uuidTmpOutB
mv {config[TMPDIR]}/$uuidTmpOutB {output.bigBed}
cp -f {input.gff} {output.gff}

echo -e "
\ttrack {wildcards.techname}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.barcodes}_{wildcards.minReadSupport}_FL_TMs
\tbigDataUrl {TRACK_HUB_DATA_URL}$(basename {output.bigBed})
\tshortLabel {config[PROJECT_NAME]} {wildcards.techname} {wildcards.capDesign} {wildcards.sizeFrac} {wildcards.barcodes} min{wildcards.minReadSupport}reads FL TMs
\tlongLabel {config[PROJECT_NAME]} Full-length transcript models (Protocol: {wildcards.techname}, Capture design: {wildcards.capDesign}, Size fraction: {wildcards.sizeFrac}, Sample: {wildcards.barcodes}, Min. read support per TM: {wildcards.minReadSupport})
\ttype bigBed 12
\titemRgb on
\tsearchIndex name
\tparent {wildcards.techname}_{wildcards.capDesign}_output
\tsubGroups filter=flTMs sample={wildcards.barcodes} protocol={wildcards.techname} sizeFraction={wildcards.sizeFrac} minReadSupport={wildcards.minReadSupport}

" > {config[TMPDIR]}/$uuidTmpOutT
mv {config[TMPDIR]}/$uuidTmpOutT {output.trackDb}

		'''




rule catTrackDb:
	input:
		lambda wildcards: expand(config["TRACK_HUB_DIR"] + "{capDesign}_SupertrackHeader.txt", capDesign=GENOMETOCAPDESIGNS[wildcards.genome]),
		lambda wildcards: expand(config["TRACK_HUB_DIR"] + "{capDesign}_trackDbInputHeader.txt", nonPreCapOnly, capDesign=GENOMETOCAPDESIGNS[wildcards.genome]),
		lambda wildcards: expand(config["TRACK_HUB_DIR"] + "{capDesign}_trackDbHiSeqOutputHeader.txt", capDesign=GENOMETOCAPDESIGNS[wildcards.genome]) if config["USE_MATCHED_ILLUMINA"] else None,
		lambda wildcards: expand(config["TRACK_HUB_DIR"] + "{capDesign}_primary_targets.trackDb.txt", nonPreCapOnly, capDesign=GENOMETOCAPDESIGNS[wildcards.genome]),
		lambda wildcards: expand(config["TRACK_HUB_DIR"] + "{techname}_{capDesign}_trackDbOutputHeader.txt", filtered_product, techname=TECHNAMES, capDesign=GENOMETOCAPDESIGNS[wildcards.genome]),
		lambda wildcards: expand(config["TRACK_HUB_DIR"] + "{techname}_{capDesign}_{sizeFrac}_{barcodes}_TmergeOutput.min{minReadSupport}reads.trackDb.txt", filtered_product, techname=TECHNAMES, capDesign=GENOMETOCAPDESIGNS[wildcards.genome], sizeFrac=SIZEFRACS, barcodes=BARCODES, minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"]),
		lambda wildcards: expand(config["TRACK_HUB_DIR"] + "{techname}_{capDesign}_{sizeFrac}_{barcodes}_FLTmergeOutput.min{minReadSupport}reads.trackDb.txt", filtered_product, techname=TECHNAMES, capDesign=GENOMETOCAPDESIGNS[wildcards.genome], sizeFrac=SIZEFRACS, barcodes=BARCODES, minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"]),
		lambda wildcards: expand(config["TRACK_HUB_DIR"] + "{techname}_{capDesign}_{sizeFrac}_{barcodes}_BamOutput.trackDb.txt", filtered_product, techname=TECHNAMES, capDesign=GENOMETOCAPDESIGNS[wildcards.genome], sizeFrac=SIZEFRACS, barcodes=BARCODES),
		lambda wildcards: expand(config["TRACK_HUB_DIR"] + "hiSeq_{capDesign}_BamOutput.trackDb.txt", capDesign=GENOMETOCAPDESIGNS[wildcards.genome]) if config["USE_MATCHED_ILLUMINA"] else None
	output:
		config["TRACK_HUB_DIR"] + "{genome}/trackDb.txt"
	shell:
		'''
uuidTmpOut=$(uuidgen)
cat {input} > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
		'''
