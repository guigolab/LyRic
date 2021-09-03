rule make_hubtxt:
	output: "output/trackHub/" + "hub.txt"
	shell:
		'''
uuidTmpOut=$(uuidgen)
echo "hub {config[PROJECT_NAME]}
shortLabel {config[PROJECT_LONG_NAME]}
longLabel {config[PROJECT_LONG_NAME]}
genomesFile genomes.txt
email {config[PROJECT_CONTACT_EMAIL]}
" > {TMPDIR}/$uuidTmpOut
mv {TMPDIR}/$uuidTmpOut {output}

		'''

rule make_genomestxt:
	output:temp("output/trackHub/" + "tmp/{genome}.txt")
	shell:
		'''
uuidTmpOut=$(uuidgen)
echo "genome {wildcards.genome}
trackDb {wildcards.genome}/trackDb.txt
" > {TMPDIR}/$uuidTmpOut
mv {TMPDIR}/$uuidTmpOut {output}
		'''

rule agg_genomestxt:
	input: expand("output/trackHub/" + "tmp/{genome}.txt", genome=GENOMES)
	output: "output/trackHub/" + "genomes.txt"
	shell:
		'''
uuidTmpOut=$(uuidgen)
cat {input} > {TMPDIR}/$uuidTmpOut
mv {TMPDIR}/$uuidTmpOut {output}
		'''

rule makeCapDesignSupertrackHeader:
	output: temp("output/trackHub/" + "{capDesign}_SupertrackHeader.txt")
	shell:
		'''
uuidTmpOut=$(uuidgen)
echo "track {wildcards.capDesign}
superTrack on show
shortLabel {config[PROJECT_NAME]} {wildcards.capDesign} design
longLabel {config[PROJECT_NAME]} {wildcards.capDesign} capture design
visibility full
" > {TMPDIR}/$uuidTmpOut
mv {TMPDIR}/$uuidTmpOut {output}

		'''

rule makeTrackDbInputHeader:
	output: temp("output/trackHub/" + "{capDesign}_trackDbInputHeader.txt")
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
" > {TMPDIR}/$uuidTmpOut
mv {TMPDIR}/$uuidTmpOut {output}
		'''






rule makeTrackDbInputTracks:
	input:
		gff=lambda wildcards: CAPDESIGNTOTARGETSGFF[CAPDESIGNTOCAPDESIGN[wildcards.capDesign]],
		genome=lambda wildcards: config["GENOMESDIR"] + CAPDESIGNTOGENOME[wildcards.capDesign] + ".sorted.genome"
	output:
		bigBed="output/trackHub/" + "dataFiles/{capDesign}_primary_targets.bb",
		trackDb=temp("output/trackHub/" + "{capDesign}_primary_targets.trackDb.txt")
	shell:
		'''
uuid=$(uuidgen)
uuid2=$(uuidgen)
uuid3=$(uuidgen)
uuidTmpOutB=$(uuidgen)
uuidTmpOutT=$(uuidgen)
cat {input.genome} | cut -f1 | sort -T {TMPDIR}  |uniq > {TMPDIR}/$uuid2

cat {input.gff} | grep -P "^chr" | grep -v "chrIS" |gff2bed_full.pl -|sort -T {TMPDIR}  -k1,1 -k2,2n -k3,3n  > {TMPDIR}/$uuid3
join  -j1 {TMPDIR}/$uuid2 {TMPDIR}/$uuid3 |ssv2tsv > {TMPDIR}/$uuid
bedToBigBed -as=./LyRic/misc/bed12.as -type=bed12 -extraIndex=name {TMPDIR}/$uuid {input.genome} {TMPDIR}/$uuidTmpOutB
mv {TMPDIR}/$uuidTmpOutB {output.bigBed}

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

" > {TMPDIR}/$uuidTmpOutT
mv {TMPDIR}/$uuidTmpOutT {output.trackDb}

		'''

rule makeTrackDbHiSeqOutputHeader:
	output: temp("output/trackHub/" + "{capDesign}_trackDbHiSeqOutputHeader.txt")
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
" > {TMPDIR}/$uuidTmpOut
mv {TMPDIR}/$uuidTmpOut {output}
		'''

rule makeHiSeqBamOutputTracks:
	input: "output/mappings/shortReadMappings/hiSeq_{capDesign}.bam"
	output:
		bam="output/trackHub/" + "dataFiles/hiSeq_{capDesign}.bam",
		bai="output/trackHub/" + "dataFiles/hiSeq_{capDesign}.bam.bai"
	shell:
		'''
uuidTmpOut=$(uuidgen)
uuid=$(uuidgen)
samtools view -H {input} > {TMPDIR}/$uuid
samtools view {input} | grep -P "\tchr" | grep -v "chrIS" >> {TMPDIR}/$uuid
samtools view -b {TMPDIR}/$uuid > {TMPDIR}/$uuidTmpOut
sleep 200s
samtools index {TMPDIR}/$uuidTmpOut

mv {TMPDIR}/$uuidTmpOut {output.bam}
mv {TMPDIR}/$uuidTmpOut.bai  {output.bai}
		'''

rule makeHiSeqBamOutputTrackDb:
	input: "output/trackHub/" + "dataFiles/hiSeq_{capDesign}.bam",
	output:
		trackDb=temp("output/trackHub/" + "hiSeq_{capDesign}_BamOutput.trackDb.txt")
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

" > {TMPDIR}/$uuidTmpOut
mv {TMPDIR}/$uuidTmpOut {output.trackDb}

		'''


rule makeTrackDbTechnameHeader:
	params:
		subGroupString=lambda wildcards: trackHubSubGroupString(wildcards.techname, wildcards.capDesign, SIZEFRACS,SAMPLEREPS, MINIMUM_TMERGE_READ_SUPPORT),
		viewsString="allTMs=All_TMs flTMs=Full-length_TMs BAMs=BAMs"
	output: temp("output/trackHub/" + "{techname}_{capDesign}_trackDbOutputHeader.txt")
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
" > {TMPDIR}/$uuidTmpOut
mv {TMPDIR}/$uuidTmpOut {output}
		'''

rule makeBamOutputTracks:
	input: "output/mappings/longReadMapping/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.bam"
	output:
		bam="output/trackHub/" + "dataFiles/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.bam",
		bai="output/trackHub/" + "dataFiles/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.bam.bai"
	shell:
		'''
uuidTmpOut=$(uuidgen)
uuid=$(uuidgen)
samtools view -H {input} > {TMPDIR}/$uuid
samtools view {input} | grep -P "\\tchr" | grep -v "chrIS" >> {TMPDIR}/$uuid
samtools view -b {TMPDIR}/$uuid > {TMPDIR}/$uuidTmpOut
sleep 200s
samtools index {TMPDIR}/$uuidTmpOut
mv {TMPDIR}/$uuidTmpOut {output.bam}
mv {TMPDIR}/$uuidTmpOut.bai {output.bai}

		'''


rule makeBamOutputTrackDb:
	input: "output/trackHub/" + "dataFiles/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.bam"
	output:
		trackDb=temp("output/trackHub/" + "{techname}_{capDesign}_{sizeFrac}_{sampleRep}_BamOutput.trackDb.txt")

	shell:
		'''
uuidTmpOut=$(uuidgen)
echo -e "
\ttrack {wildcards.techname}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.sampleRep}_BAM
\tbigDataUrl {TRACK_HUB_DATA_URL}$(basename {input})
\tshortLabel {config[PROJECT_NAME]} {wildcards.techname} {wildcards.capDesign} {wildcards.sizeFrac} {wildcards.sampleRep} BAMs
\tlongLabel {config[PROJECT_NAME]} BAM (read alignments) (Protocol: {wildcards.techname}, Capture design: {wildcards.capDesign}, Size fraction: {wildcards.sizeFrac}, Sample: {wildcards.sampleRep})
\ttype bam
\tbamColorMode gray
\tbamGrayMode aliQual
\tdoWiggle on
\tindelPolyA on
\talwaysZero on
\tmaxHeightPixels 100:35:11
\tparent {wildcards.techname}_{wildcards.capDesign}_output
\tsubGroups filter=BAMs sample={wildcards.sampleRep} protocol={wildcards.techname} sizeFraction={wildcards.sizeFrac} 

" > {TMPDIR}/$uuidTmpOut
mv {TMPDIR}/$uuidTmpOut {output.trackDb}

		'''





rule makeTmergeOutputTracks:
	input:
		bed="output/mappings/mergedReads/colored/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.splicing_status-all.endSupport-all.bed",
		gff="output/mappings/mergedReads/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.HiSS.tmerge.min{minReadSupport}reads.splicing_status-all.endSupport-all.gff.gz",
		genome=lambda wildcards: config["GENOMESDIR"] + CAPDESIGNTOGENOME[wildcards.capDesign] + ".sorted.genome"
	output:
		bigBed="output/trackHub/" + "dataFiles/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.bb",
		gff="output/trackHub/" + "dataFiles/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.gff.gz",
		trackDb=temp("output/trackHub/" + "{techname}_{capDesign}_{sizeFrac}_{sampleRep}_TmergeOutput.min{minReadSupport}reads.trackDb.txt")
	shell:
		'''
uuid=$(uuidgen)
uuidTmpOutB=$(uuidgen)
uuidTmpOutT=$(uuidgen)

cat {input.bed} | grep -P "^chr" | grep -v "chrIS"  > {TMPDIR}/$uuid
bedToBigBed -as=./LyRic/misc/bed12.as -type=bed12 -extraIndex=name {TMPDIR}/$uuid {input.genome} {TMPDIR}/$uuidTmpOutB
mv {TMPDIR}/$uuidTmpOutB {output.bigBed}
cp -f {input.gff} {output.gff}
echo -e "
\ttrack {wildcards.techname}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.sampleRep}_{wildcards.minReadSupport}_TMs
\tbigDataUrl {TRACK_HUB_DATA_URL}$(basename {output.bigBed})
\tshortLabel {config[PROJECT_NAME]} {wildcards.techname} {wildcards.capDesign} {wildcards.sizeFrac} {wildcards.sampleRep} min{wildcards.minReadSupport}reads TMs
\tlongLabel {config[PROJECT_NAME]} Transcript models (Protocol: {wildcards.techname}, Capture design: {wildcards.capDesign}, Size fraction: {wildcards.sizeFrac}, Sample: {wildcards.sampleRep}, Min. read support per TM: {wildcards.minReadSupport})
\ttype bigBed 12
\titemRgb on
\tsearchIndex name
\tparent {wildcards.techname}_{wildcards.capDesign}_output
\tsubGroups filter=allTMs sample={wildcards.sampleRep} protocol={wildcards.techname} sizeFraction={wildcards.sizeFrac} minReadSupport={wildcards.minReadSupport}

" > {TMPDIR}/$uuidTmpOutT
mv {TMPDIR}/$uuidTmpOutT {output.trackDb}

		'''

rule makeFLTmergeOutputTracks:
	input:
		bed="output/mappings/mergedReads/colored/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.splicing_status-all.endSupport-cagePolyASupported.bed",
		gff="output/mappings/mergedReads/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.HiSS.tmerge.min{minReadSupport}reads.splicing_status-all.endSupport-cagePolyASupported.gff.gz",
		genome=lambda wildcards: config["GENOMESDIR"] + CAPDESIGNTOGENOME[wildcards.capDesign] + ".sorted.genome"
	output:
		bigBed="output/trackHub/" + "dataFiles/{techname}_{capDesign}_{sizeFrac}_{sampleRep}_FL.tmerge.min{minReadSupport}reads.bb",
		gff="output/trackHub/" + "dataFiles/{techname}_{capDesign}_{sizeFrac}_{sampleRep}_FL.tmerge.min{minReadSupport}reads.gff.gz",
		trackDb=temp("output/trackHub/" + "{techname}_{capDesign}_{sizeFrac}_{sampleRep}_FLTmergeOutput.min{minReadSupport}reads.trackDb.txt")
	shell:
		'''
uuid=$(uuidgen)
uuidTmpOutB=$(uuidgen)
uuidTmpOutT=$(uuidgen)

cat {input.bed} | grep -P "^chr" | grep -v "chrIS"  > {TMPDIR}/$uuid
bedToBigBed -as=./LyRic/misc/bed12.as -type=bed12 -extraIndex=name {TMPDIR}/$uuid {input.genome} {TMPDIR}/$uuidTmpOutB
mv {TMPDIR}/$uuidTmpOutB {output.bigBed}
cp -f {input.gff} {output.gff}

echo -e "
\ttrack {wildcards.techname}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.sampleRep}_{wildcards.minReadSupport}_FL_TMs
\tbigDataUrl {TRACK_HUB_DATA_URL}$(basename {output.bigBed})
\tshortLabel {config[PROJECT_NAME]} {wildcards.techname} {wildcards.capDesign} {wildcards.sizeFrac} {wildcards.sampleRep} min{wildcards.minReadSupport}reads FL TMs
\tlongLabel {config[PROJECT_NAME]} Full-length transcript models (Protocol: {wildcards.techname}, Capture design: {wildcards.capDesign}, Size fraction: {wildcards.sizeFrac}, Sample: {wildcards.sampleRep}, Min. read support per TM: {wildcards.minReadSupport})
\ttype bigBed 12
\titemRgb on
\tsearchIndex name
\tparent {wildcards.techname}_{wildcards.capDesign}_output
\tsubGroups filter=flTMs sample={wildcards.sampleRep} protocol={wildcards.techname} sizeFraction={wildcards.sizeFrac} minReadSupport={wildcards.minReadSupport}

" > {TMPDIR}/$uuidTmpOutT
mv {TMPDIR}/$uuidTmpOutT {output.trackDb}

		'''




rule catTrackDb:
	input:
		lambda wildcards: expand("output/trackHub/" + "{capDesign}_SupertrackHeader.txt", capDesign=GENOMETOCAPDESIGNS[wildcards.genome]),
		lambda wildcards: expand("output/trackHub/" + "{capDesign}_trackDbInputHeader.txt", nonPreCapOnly, capDesign=GENOMETOCAPDESIGNS[wildcards.genome]),
		lambda wildcards: expand("output/trackHub/" + "{capDesign}_trackDbHiSeqOutputHeader.txt", capDesign=GENOMETOCAPDESIGNS[wildcards.genome]) if config["USE_MATCHED_ILLUMINA"] else None,
		lambda wildcards: expand("output/trackHub/" + "{capDesign}_primary_targets.trackDb.txt", nonPreCapOnly, capDesign=GENOMETOCAPDESIGNS[wildcards.genome]),
		lambda wildcards: expand("output/trackHub/" + "{techname}_{capDesign}_trackDbOutputHeader.txt", filtered_product, techname=TECHNAMES, capDesign=GENOMETOCAPDESIGNS[wildcards.genome]),
		lambda wildcards: expand("output/trackHub/" + "{techname}_{capDesign}_{sizeFrac}_{sampleRep}_TmergeOutput.min{minReadSupport}reads.trackDb.txt", filtered_product, techname=TECHNAMES, capDesign=GENOMETOCAPDESIGNS[wildcards.genome], sizeFrac=SIZEFRACS, sampleRep=SAMPLEREPS, minReadSupport=MINIMUM_TMERGE_READ_SUPPORT),
		lambda wildcards: expand("output/trackHub/" + "{techname}_{capDesign}_{sizeFrac}_{sampleRep}_FLTmergeOutput.min{minReadSupport}reads.trackDb.txt", filtered_product, techname=TECHNAMES, capDesign=GENOMETOCAPDESIGNS[wildcards.genome], sizeFrac=SIZEFRACS, sampleRep=SAMPLEREPS, minReadSupport=MINIMUM_TMERGE_READ_SUPPORT),
		lambda wildcards: expand("output/trackHub/" + "{techname}_{capDesign}_{sizeFrac}_{sampleRep}_BamOutput.trackDb.txt", filtered_product, techname=TECHNAMES, capDesign=GENOMETOCAPDESIGNS[wildcards.genome], sizeFrac=SIZEFRACS, sampleRep=SAMPLEREPS),
		lambda wildcards: expand("output/trackHub/" + "hiSeq_{capDesign}_BamOutput.trackDb.txt", capDesign=GENOMETOCAPDESIGNS[wildcards.genome]) if config["USE_MATCHED_ILLUMINA"] else None
	output:
		"output/trackHub/" + "{genome}/trackDb.txt"
	shell:
		'''
uuidTmpOut=$(uuidgen)
cat {input} > {TMPDIR}/$uuidTmpOut
mv {TMPDIR}/$uuidTmpOut {output}
		'''
