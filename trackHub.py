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

rule makeTrackDbInputHeader:
	output: temp(config["TRACK_HUB_DIR"] + "{genome}/trackDbInputHeader.txt")
	shell:
		'''
echo "track input
superTrack on
shortLabel Targets
longLabel Targets (Captured regions)
priority 10
visibility dense
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
cat {input.genome} | cut -f1 | sort |uniq > $TMPDIR/$uuid2

cat {input.gff} | grep -P "^chr" | grep -w -f $TMPDIR/$uuid2 | gff2bed_full.pl - |sortbed > $TMPDIR/$uuid
bedToBigBed -as={config[BED12_AUTOSQL]} -type=bed12 -extraIndex=name $TMPDIR/$uuid {input.genome} {output.bigBed}

echo -e "
\ttrack {wildcards.capDesign}_targets
\tbigDataUrl {TRACK_HUB_DATA_URL}$(basename {output.bigBed})
\tshortLabel {wildcards.capDesign} targets
\tlongLabel Targeted regions in {wildcards.capDesign} capture design
\ttype bigBed 12
\titemRgb on
\tsearchIndex name
\tparent input
\tvisibility dense

" > {output.trackDb}

		'''

rule makeBam









rule catTrackDb:
	input:
		config["TRACK_HUB_DIR"] + "{genome}/trackDbInputHeader.txt",
		lambda wildcards: expand(config["TRACK_HUB_DIR"] + "{capDesign}" + "_primary_targets.trackDb.txt", capDesign=GENOMETOCAPDESIGNS[wildcards.genome])
	output:
		config["TRACK_HUB_DIR"] + "{genome}/trackDb.txt"
	shell:
		'''
cat {input} > {output}
		'''