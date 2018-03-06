rule integratePolyaAndSjInfo:
	input:
		polyA = "mappings/" + "removePolyAERCCs/{techname}_{capDesign}_{barcodes}.polyAsitesNoErcc.bed",
		SJs = "mappings/" + "getIntronMotif/{techname}_{capDesign}_{barcodes}.transcripts.tsv"
	output: "mappings/" + "integratePolyaAndSjInfo/{techname}_{capDesign}_{barcodes}.polyA+SJ.strandInfo.tsv"
	shell:
		'''
cat {input.polyA} | cut -f4,6 | awk '$2!="."'| sort > $TMPDIR/reads.polyA.strandInfo.tsv
cat {input.SJs} | skipcomments | cut -f 1,2 | awk '$2!="."' | sort> $TMPDIR/reads.SJ.strandInfo.tsv
cat $TMPDIR/reads.polyA.strandInfo.tsv $TMPDIR/reads.SJ.strandInfo.tsv | awk '$2!="."'| sort|uniq > {output}
		'''

rule strandGffs:
	input:
		gff = "mappings/" + "readBedToGff/{techname}_{capDesign}_{barcodes}.merged.gff",
		strandInfo = "mappings/" + "integratePolyaAndSjInfo/{techname}_{capDesign}_{barcodes}.polyA+SJ.strandInfo.tsv"
	output: "mappings/" + "strandGffs/{techname}_{capDesign}_{barcodes}.stranded.gff"
	shell:
		'''
get_right_transcript_strand.pl {input.gff} {input.strandInfo} | fgrep -v ERCC- | sortgff> {output}

		'''

rule highConfidenceReads:
	input:
		strandInfo = "mappings/" + "integratePolyaAndSjInfo/{techname}_{capDesign}_{barcodes}.polyA+SJ.strandInfo.tsv",
		strandedReads = "mappings/" + "strandGffs/{techname}_{capDesign}_{barcodes}.stranded.gff"
	output:
		"mappings/" + "highConfidenceReads/{techname}_{capDesign}_{barcodes}.strandedHCGMs.gff"
	shell:
		'''
#select read IDs with canonical GT|GC/AG and high-confidence SJs
cat {input.strandInfo} | skipcomments | awk '$6==1 && $7==1' | cut -f1 | sort|uniq > $TMPDIR/reads.hcSJs.list
echo $?
fgrep -w -f $TMPDIR/reads.hcSJs.list {input.strandedReads} || true > $TMPDIR/gtag.gff
 echo $?
cat {input.strandedReads} | extractGffAttributeValue.pl transcript_id | sort|uniq -u > $TMPDIR/tmp
fgrep -w -f $TMPDIR/tmp {input.strandedReads} || true > $TMPDIR/tmp2
cat $TMPDIR/tmp2 | fgrep -v ERCC- || true > $TMPDIR/monoPolyA.gff
 echo $?
cat $TMPDIR/gtag.gff $TMPDIR/monoPolyA.gff | sortgff> {output}
 echo $?
		'''