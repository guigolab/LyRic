import glob
from collections import defaultdict
import os
import itertools
import sys


GGPLOT_PUB_QUALITY=config["GGPLOT_PUB_QUALITY"]

TECHNAMES=["pacBio:Cshl:Smarter:", "ont:Cshl:Smarter:", "ont:Cshl:dRNA:"]
CAPDESIGNS=["HpreCap",]
HISEQ_LIB_PROTOCOLS=["RzDirSmarter", "RzKapaCapTrap", "RzKapaSmarter", "TotDirSmarter"]
SIZEFRACS=["0+",]
BARCODES=["Heart",]
FINALCORRECTIONLEVELS=["0",]
MINSEQQUAL=["0", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20"]
MINREADSUPPORT=["1","2","3","4","5","6","7","8","9","10"]
EXONOVERHANG=["0","2","4","6","8","10","12","14", "16"]
CAPDESIGNTOGENOME=config["capDesignToGenome"]
#TmergeEndFuzz="100"

optimalMINSEQQUAL="10"
optimalEXONOVERHANG="8"
optimalMINREADSUPPORT="2"


def filtered_product(*args): #to avoid useless combinations of wildcards
	found=False
	for wc_comb in itertools.product(*args):
		#yield(wc_comb)
		#print(wc_comb)

		if (wc_comb[5][0] in ('minReadSupport') and wc_comb[5][1] in ('1')) and (wc_comb[6][0] in ('minSeqQual') and wc_comb[6][1] in ('0')) and (wc_comb[7][0] in ('exonOverhang') and wc_comb[7][1] in ('0')):
			found=True
			#print ("AUTH")
			yield(wc_comb)
		elif (wc_comb[5][0] in ('minReadSupport') and wc_comb[5][1] not in ('1')) and (wc_comb[6][0] in ('minSeqQual') and wc_comb[6][1] in ('0')) and (wc_comb[7][0] in ('exonOverhang') and wc_comb[7][1] in ('0')):
			found=True
			#print ("AUTH")
			yield(wc_comb)
		elif (wc_comb[5][0] in ('minReadSupport') and wc_comb[5][1] in ('1')) and (wc_comb[6][0] in ('minSeqQual') and wc_comb[6][1] not in ('0')) and (wc_comb[7][0] in ('exonOverhang') and wc_comb[7][1] in ('0')):
			found=True
			#print ("AUTH")
			yield(wc_comb)
		elif (wc_comb[5][0] in ('minReadSupport') and wc_comb[5][1] in ('1')) and (wc_comb[6][0] in ('minSeqQual') and wc_comb[6][1] in ('0')) and (wc_comb[7][0] in ('exonOverhang') and wc_comb[7][1] not in ('0')):
			found=True
			#print ("AUTH")
			yield(wc_comb)
		elif (wc_comb[5][0] in ('minReadSupport') and wc_comb[5][1] == optimalMINREADSUPPORT) and (wc_comb[6][0] in ('minSeqQual') and wc_comb[6][1] == optimalMINSEQQUAL) and (wc_comb[7][0] in ('exonOverhang') and wc_comb[7][1] == optimalEXONOVERHANG):
			found=True
			#print ("AUTH")
			yield(wc_comb)
		#else:
			#print ("NO")
	if not found:
		print(" Error in function filtered_product. Args were:")
		print((args))
		quit(" Error. Could not yield any input file.")



rule all:
	input:
		expand("plots/non_snakemake/all.tmerge.minByreads.minQual0.exonOverhang0.vs.SIRVs_lib{libProt}.stats.png", libProt=HISEQ_LIB_PROTOCOLS),
		expand("plots/non_snakemake/all.tmerge.min1reads.minQualBy.exonOverhang0.vs.SIRVs_lib{libProt}.stats.png", libProt=HISEQ_LIB_PROTOCOLS),
		expand("plots/non_snakemake/all.tmerge.min1reads.minQual0.exonOverhangBy.vs.SIRVs_lib{libProt}.stats.png", libProt=HISEQ_LIB_PROTOCOLS),
		expand("plots/non_snakemake/all.tmerge.optimalParams.vs.SIRVs_lib{libProt}.stats.png", libProt=HISEQ_LIB_PROTOCOLS),
		expand("test/sirvMappings/nonAnchoredMergeReads/colored/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.SIRVome.minQual{minSeqQual}.exonOverhang{exonOverhang}_lib{libProt}.all.bed", filtered_product, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES, minReadSupport=MINREADSUPPORT, minSeqQual=MINSEQQUAL, exonOverhang=EXONOVERHANG, libProt=HISEQ_LIB_PROTOCOLS)



rule getHiSeqCanonicalIntronsList:
	input:
		bam="fastqs/hiSeq/mappings/{barcodes}_{libProt}.bam",
		genome = lambda wildcards: config["GENOMESDIR"] + CAPDESIGNTOGENOME["HpreCap"] + ".fa"
	threads: 6
	output:
		list="test/sirvMappings/hiSeqIntrons/hiSeq_{barcodes}_lib{libProt}.canonicalIntrons.list",
		#stats=temp(config["STATSDATADIR"] + "{capDesign}_tmp.hiSeq.SJs.stats.tsv")
	shell:
		'''
uuid=$(uuidgen)
echoerr "making bed"
samtools view -b -F 256 -F4 -F 2048 {input.bam}  |bedtools bamtobed -i stdin -bed12 | fgrep -v ERCC- > {config[TMPDIR]}/$uuid.hiSeq.bed
echoerr "splitting"
split -a 3 -d -e -n l/24 {config[TMPDIR]}/$uuid.hiSeq.bed {config[TMPDIR]}/$uuid.hiSeq.bed.split
rm {config[TMPDIR]}/$uuid.hiSeq.bed
for file in `ls {config[TMPDIR]}/$uuid.hiSeq.bed.split*`; do
echo "cat $file | awk -f ~/julien_utils/bed12fields2gff.awk > $file.gff; sort -T {config[TMPDIR]} -k1,1 -k4,4n -k5,5n $file.gff | makeIntrons.pl - | extract_intron_strand_motif.pl - {input.genome} {config[TMPDIR]}/$(basename $file); rm $file $file.gff $file.transcripts.tsv"
done > {config[TMPDIR]}/$uuid.parallelIntrons.sh
echoerr "extracting introns on split files"

parallel -j {threads} < {config[TMPDIR]}/$uuid.parallelIntrons.sh
echoerr "getting SJs and merging into output..."
cat {config[TMPDIR]}/$uuid.hiSeq.bed.split*.introns.gff | perl -lane '$start=$F[3]-1; $end=$F[4]+1; print $F[0]."_".$start."_".$end."_".$F[6]' | sort -T {config[TMPDIR]} |uniq> {output.list}

		'''




rule readMapping:
	input:
		reads="fastqs/{techname}_{capDesign}_{sizeFrac}.{barcodes}.fastq.gz",
		genome="/users/rg/jlagarde/genomes/lexogen_SIRVs/SIRV_Set1_Lot00141_Sequences_181206a/SIRVome_isoforms_170612a.unix.fasta"
	output: "test/sirvMappings/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bam"
	params:
		minimap_preset = lambda wildcards: "splice" if wildcards.techname.find('ont') == 0 else "splice:hq" if wildcards.techname.find('pacBio') == 0 else None
	threads: 4
	shell:
		'''
uuid=$(uuidgen)
minimap2 -x {params.minimap_preset} --cs -t {threads} --secondary=no --splice-flank=no -L -a {input.genome}  {input.reads} > {config[TMPDIR]}/$uuid
samtools view -H {config[TMPDIR]}/$uuid > {config[TMPDIR]}/$uuid.2
samtools view -F 256 -F4 -F 2048 {config[TMPDIR]}/$uuid >> {config[TMPDIR]}/$uuid.2
cat {config[TMPDIR]}/$uuid.2 | samtools sort --threads {threads} -m 5G - > {output}

		'''

rule readBamToBed:
	input: "test/sirvMappings/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bam"
	output: "test/sirvMappings/readBamToBed/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bed.gz"
	#input: "mappings/readMapping/{basename}.bam"
	#output: "mappings/readBamToBed/{basename}.bed"
	shell:
		'''
bedtools bamtobed -i {input} -bed12 | sortbed | gzip > {output}

		'''

rule readBedToGff:
	input: "test/sirvMappings/readBamToBed/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bed.gz"
	output: "test/sirvMappings/readBedToGff/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.gff.gz"
	shell:
		'''
zcat {input} | awk -f ~jlagarde/julien_utils/bed12fields2gff.awk | sortgff | gzip > {output}
		'''


rule polyAmapping:
	input:
		reads = "test/sirvMappings/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bam",
		genome = lambda wildcards: config["GENOMESDIR"] + CAPDESIGNTOGENOME[wildcards.capDesign] + ".fa"
	params:
		minAcontent=0.8
	output: temp("test/sirvMappings/polyAmapping/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.polyAsites.bed")
	shell:
		'''
samtools view {input.reads} | samToPolyA.pl --minClipped=10 --minAcontent={params.minAcontent}  --discardInternallyPrimed --minUpMisPrimeAlength=10 --genomeFasta={input.genome} - |sortbed > {output}
		'''


rule makeIntrons:
	input: "test/sirvMappings/readBedToGff/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.gff.gz"
	output: temp("test/sirvMappings/makeIntrons/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.introns.gff.gz")
	shell:
		'''
zcat {input} | makeIntrons.pl - | sortgff |gzip> {output}
		'''

rule getIntronMotif:
	input:
		introns = "test/sirvMappings/makeIntrons/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.introns.gff.gz",
		genome = lambda wildcards: config["GENOMESDIR"] + CAPDESIGNTOGENOME[wildcards.capDesign] + ".fa"
	output:
		gff = temp("test/sirvMappings/getIntronMotif/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.introns.gff"),
		tsv = temp("test/sirvMappings/getIntronMotif/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.transcripts.tsv")
	shell:
		'''
zcat {input.introns} | extract_intron_strand_motif.pl - {input.genome} $(dirname {output.gff})/$(basename {output.gff} .introns.gff)

		'''


rule integratePolyaAndSjInfo:
	input:
		polyA = "test/sirvMappings/polyAmapping/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.polyAsites.bed",
		SJs = "test/sirvMappings/getIntronMotif/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.transcripts.tsv"
	output:
		strandInfo=temp("test/sirvMappings/integratePolyaAndSjInfo/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.polyA+SJ.strandInfo.tsv"),
		wrongPolyAs=temp("test/sirvMappings/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.wrongPolyAs.list")
	shell:
		'''
uuid=$(uuidgen)
cat {input.SJs} | tgrep -v -P "^#" | cut -f 1,2 | awk '$2!="."' | sort> {config[TMPDIR]}/$uuid.reads.SJ.strandInfo.tsv
cat {input.polyA} | cut -f4,6 | awk '$2!="."'| sort > {config[TMPDIR]}/$uuid.reads.polyA.strandInfo.tsv
join -a1 -a2 -e '.' -o '0,1.2,2.2' {config[TMPDIR]}/$uuid.reads.SJ.strandInfo.tsv  {config[TMPDIR]}/$uuid.reads.polyA.strandInfo.tsv > {config[TMPDIR]}/$uuid.reads.SJ.polyA.strandInfo.tsv

#make list of reads with wrongly called polyA sites (i.e. their strand is different from the one inferred using SJs):
cat {config[TMPDIR]}/$uuid.reads.SJ.polyA.strandInfo.tsv | perl -slane 'if($F[2] ne "." && $F[1] ne "." && $F[1] ne $F[2]){{print join ("\\t", $F[0])}}' |sort|uniq > {output.wrongPolyAs}
#get strand info, prioritizing SJ inference
cat {config[TMPDIR]}/$uuid.reads.SJ.polyA.strandInfo.tsv | perl -slane 'if($F[1] ne "."){{print "$F[0]\\t$F[1]"}} else{{print "$F[0]\\t$F[2]"}}' |sort|uniq > {output.strandInfo}
		'''


rule removeWrongPolyAs:
	input:
		polyA="test/sirvMappings/polyAmapping/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.polyAsites.bed",
		wrongPolyAs="test/sirvMappings/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.wrongPolyAs.list"
	output: temp("test/sirvMappings/removePolyAERCCs/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.polyAsitesNoErcc.bed")
	shell:
		'''
cat {input.polyA} | fgrep -v -w -f {input.wrongPolyAs} > {output}
		'''

rule strandGffs:
	input:
		gff = "test/sirvMappings/readBedToGff/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.gff.gz",
		strandInfo = "test/sirvMappings/integratePolyaAndSjInfo/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.polyA+SJ.strandInfo.tsv"
	output: temp("test/sirvMappings/strandGffs/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.stranded.gff.gz")
	shell:
		'''
uuid=$(uuidgen)
zcat {input.gff} > {config[TMPDIR]}/$uuid.in.gff
get_right_transcript_strand.pl {config[TMPDIR]}/$uuid.in.gff {input.strandInfo} | sortgff |gzip> {output}

		'''



rule highConfidenceReads:
	input:
		transcriptStrandInfo = "test/sirvMappings/getIntronMotif/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.transcripts.tsv",
		strandedReads = "test/sirvMappings/strandGffs/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.stranded.gff.gz",
		bam = "test/sirvMappings/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bam"
	output:
		temp("test/sirvMappings/highConfidenceReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.strandedHCGMs.SIRVome.minQual{minSeqQual}.gff.gz")
	shell:
		'''
#select read IDs with canonical GT|GC/AG and high-confidence SJs
uuid=$(uuidgen)
cat {input.transcriptStrandInfo} | tgrep -v -P "^#" | awk '$6==1 && $7==1' | cut -f1 | sort|uniq > {config[TMPDIR]}/$uuid.reads.hcSJs.list.tmp
wc -l {config[TMPDIR]}/$uuid.reads.hcSJs.list.tmp
zcat {input.strandedReads} > {config[TMPDIR]}/$uuid.str.gff

##high-quality Phred SJs
samtools view -F 256 -F4 -F 2048 {input.bam} | samHQintrons.pl --minQual {wildcards.minSeqQual} - |cut -f1|sort|uniq > {config[TMPDIR]}/$uuid.HQintrons.reads.list
tgrep -F -w -f {config[TMPDIR]}/$uuid.HQintrons.reads.list {config[TMPDIR]}/$uuid.reads.hcSJs.list.tmp > {config[TMPDIR]}/$uuid.reads.hcSJs.list

tgrep -F -w -f {config[TMPDIR]}/$uuid.reads.hcSJs.list {config[TMPDIR]}/$uuid.str.gff > {config[TMPDIR]}/$uuid.gtag.gff
wc -l {config[TMPDIR]}/$uuid.gtag.gff
cat {config[TMPDIR]}/$uuid.str.gff | extractGffAttributeValue.pl transcript_id | sort|uniq -u > {config[TMPDIR]}/$uuid.tmp
tgrep -F -w -f {config[TMPDIR]}/$uuid.tmp {config[TMPDIR]}/$uuid.str.gff > {config[TMPDIR]}/$uuid.tmp2
cat {config[TMPDIR]}/$uuid.tmp2  > {config[TMPDIR]}/$uuid.monoPolyA.gff
 echo $?
cat {config[TMPDIR]}/$uuid.gtag.gff {config[TMPDIR]}/$uuid.monoPolyA.gff | sortgff |gzip> {output}
 echo $?
		'''


rule getHCGMintrons:
	input: "test/sirvMappings/highConfidenceReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.strandedHCGMs.SIRVome.minQual{minSeqQual}.gff.gz"
	output: temp("test/sirvMappings/highConfidenceReads/introns/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.strandedHCGMs.introns.SIRVome.minQual{minSeqQual}.tsv")
	shell:
		'''
zcat {input} | makeIntrons.pl -| perl -lane '$F[11]=~s/"|;//g; $start=$F[3]-1; $end=$F[4]+1; print $F[11]."\t".$F[0]."_".$start."_".$end."_".$F[6]' | sort -k2,2 > {output}

		'''

rule getHiSeqSupportedHCGMs:
	input:
		hiSeqIntrons="test/sirvMappings/hiSeqIntrons/hiSeq_{barcodes}_lib{libProt}.canonicalIntrons.list",
		#hiSeqIntrons="mappings/hiSeqIntrons/hiSeq_{techname}_{capDesign}.{barcodes}.canonicalIntrons.list",
		lrIntrons="test/sirvMappings/highConfidenceReads/introns/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.strandedHCGMs.introns.SIRVome.minQual{minSeqQual}.tsv",
		hcgmGTF= "test/sirvMappings/highConfidenceReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.strandedHCGMs.SIRVome.minQual{minSeqQual}.gff.gz"
	output:temp("test/sirvMappings/highConfidenceReads/HiSS/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.HiSS_lib{libProt}.SIRVome.minQual{minSeqQual}.gff.gz")
	shell:
		'''
uuid=$(uuidgen)
join -v1 -1 2 -2 1 {input.lrIntrons} {input.hiSeqIntrons} |awk '{{print $2"\t"$1}}' > {config[TMPDIR]}/$uuid.{wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.barcodes}.introns.noHiSeq.tsv
zcat {input.hcgmGTF} > {config[TMPDIR]}/$uuid.$(basename {input.hcgmGTF} .gz)
cut -f1 {config[TMPDIR]}/$uuid.{wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.barcodes}.introns.noHiSeq.tsv | sort|uniq | fgrep -wv -f - {config[TMPDIR]}/$uuid.$(basename {input.hcgmGTF} .gz) |sortgff |gzip > {output}

		'''



rule nonAnchoredMergeReads:
	input: "test/sirvMappings/highConfidenceReads/HiSS/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.HiSS_lib{libProt}.SIRVome.minQual{minSeqQual}.gff.gz"
	output: "test/sirvMappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.SIRVome.minQual{minSeqQual}.exonOverhang{exonOverhang}_lib{libProt}.all.gff",
	threads:1
	wildcard_constraints:
		barcodes='(?!allTissues).+',
		sizeFrac='[0-9-+\.]+',
		techname='(?!allSeqTechs).+'
	shell:
		'''
#--endFuzz {wildcards.exonOverhang}
zcat {input} | tmerge --exonOverhangTolerance {wildcards.exonOverhang} --endFuzz {wildcards.exonOverhang} --minReadSupport {wildcards.minReadSupport} --tmPrefix {wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.barcodes}.NAM_ - |sortgff > {output}
		'''

rule nonAnchoredMergeReadsToBed:
	input: "test/sirvMappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.SIRVome.minQual{minSeqQual}.exonOverhang{exonOverhang}_lib{libProt}.all.gff"
	output: temp("test/sirvMappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.SIRVome.minQual{minSeqQual}.exonOverhang{exonOverhang}_lib{libProt}.all.bed")
	shell:
		'''
cat {input} | gff2bed_full.pl - > {output}
		'''


rule colorBedAccordingToGffCompare:
	input:
		classes="test/sirvMappings/nonAnchoredMergeReads/gffcompare/SIRVs/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.vs.SIRVs.SIRVome.minQual{minSeqQual}.exonOverhang{exonOverhang}_lib{libProt}.simple.tsv",
		tm="test/sirvMappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.SIRVome.minQual{minSeqQual}.exonOverhang{exonOverhang}_lib{libProt}.all.bed"
	output: "test/sirvMappings/nonAnchoredMergeReads/colored/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.SIRVome.minQual{minSeqQual}.exonOverhang{exonOverhang}_lib{libProt}.all.bed"
	shell:
			'''
colorNovelTxBed.pl {input.classes} {input.tm} > {output}

			'''


rule gffcompareToSirvAnnotation:
	input:
		annot=config["SIRVgff"],
		tm="test/sirvMappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.SIRVome.minQual{minSeqQual}.exonOverhang{exonOverhang}_lib{libProt}.all.gff"
	output: "test/sirvMappings/nonAnchoredMergeReads/gffcompare/SIRVs/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.vs.SIRVs.SIRVome.minQual{minSeqQual}.exonOverhang{exonOverhang}_lib{libProt}.simple.tsv"
	shell:
		'''
pref=$(basename {output} .simple.tsv)
annotFullPath=$(fullpath {input.annot})
cat {input.tm} | awk '$1=="SIRVome_isoforms" ' > $(dirname {output})/$(basename {input.tm})
cd $(dirname {output})
gffcompare -o $pref -r $annotFullPath $(basename {input.tm})
cat $pref.tracking | simplifyGffCompareClasses.pl - > $(basename {output})

		'''

rule getGffCompareSirvStats:
	input:"test/sirvMappings/nonAnchoredMergeReads/gffcompare/SIRVs/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.vs.SIRVs.SIRVome.minQual{minSeqQual}.exonOverhang{exonOverhang}_lib{libProt}.simple.tsv"
	output: temp("test/Rinput/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.SIRVome.minQual{minSeqQual}.exonOverhang{exonOverhang}_lib{libProt}.vs.SIRVs.stats.tsv")
	shell:
		'''
file=$(dirname {input})/$(basename {input} .simple.tsv)
for level in `echo Baselevel Exonlevel Intronchainlevel Intronlevel Locuslevel Transcriptlevel`; do
Sn=`cat $file |tgrep "level:" |sed 's/ //g'| sed 's/:/\\t/'|sed 's/|$//'|sed 's/|/\\t/g' | awk -v l=$level '$1==l' |cut -f2` || Sn='NaN'
Sp=`cat $file |tgrep "level:" |sed 's/ //g'| sed 's/:/\\t/'|sed 's/|$//'|sed 's/|/\\t/g' | awk -v l=$level '$1==l' |cut -f3` || Sp='NaN'
echo -e "{wildcards.techname}Corr{wildcards.corrLevel}\t{wildcards.capDesign}\t{wildcards.sizeFrac}\t{wildcards.barcodes}\t{wildcards.minReadSupport}\t{wildcards.exonOverhang}\t{wildcards.minSeqQual}\t$level\t$Sn\t$Sp";
done |sed 's/level//g' > {output}

		'''

rule aggGffCompareSirvStats:
	input: lambda wildcards:expand("test/Rinput/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.SIRVome.minQual{minSeqQual}.exonOverhang{exonOverhang}_lib{libProt}.vs.SIRVs.stats.tsv", filtered_product, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES, minReadSupport=MINREADSUPPORT, minSeqQual=MINSEQQUAL, exonOverhang=EXONOVERHANG, libProt=wildcards.libProt)
	output: "test/Rinput/all.tmerge.vs.SIRVs_lib{libProt}.stats.tsv"
	shell:
		'''
echo -e "seqTech\tcorrectionLevel\tcapDesign\tsizeFrac\ttissue\tminReadSupport\texonOverhang\tminSeqQual\tlevel\tmetric\tvalue" > {output}
echo "{input}" > tmp.1
xargs -a tmp.1 -n 20000 -L 20000 cat | awk '{{print $1"\\t"$2"\\t"$3"\\t"$4"\\t"$5"\\t"$6"\\t"$7"\\t"$8"\\tSn\\t"$9"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\t"$5"\\t"$6"\\t"$7"\\t"$8"\\tPr\\t"$10"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\t"$5"\\t"$6"\\t"$7"\\t"$8"\\tSnPr\\t"($9+$10)/2}}'| sed 's/Corr0/\tNo/'  | sort >> {output}

		'''

rule plotGffCompareSirvByminQualStats:
	input: "test/Rinput/all.tmerge.vs.SIRVs_lib{libProt}.stats.tsv"
	output: "plots/non_snakemake/all.tmerge.min1reads.minQualBy.exonOverhang0.vs.SIRVs_lib{libProt}.stats.png"
	shell:
		'''
echo "library(ggplot2)
dat <- read.table('{input}', header=T, as.is=T, sep='\\t')
dat <- subset(dat, level=='Transcript')
dat <- subset(dat, minReadSupport==1)
dat <- subset(dat, exonOverhang==0)
dat <- subset(dat, metric!='SnPr')
dat\$seqM <- paste(dat\$seqTech, dat\$metric)
palette <- c('Sn' = '#cc6600', 'Pr' = '#2d8659', 'SnPr' = '#999999')
ggplot(data=dat, aes(x=minSeqQual, y=value, group=seqM, color=metric)) +
geom_point() +
#xlim(1, 51) +
scale_color_manual(values=palette) +
geom_line(aes(linetype=seqTech)) + ylab('Pr & Sn (%)') + xlab('Min avg. Phred score required +/- 3 bases around SJs') +
ylim(0, 100) +
{GGPLOT_PUB_QUALITY}
ggsave('{output}', width=7, height=5)
"  | R --slave
		'''


rule plotGffCompareSirvByminReadSupport:
	input: "test/Rinput/all.tmerge.vs.SIRVs_lib{libProt}.stats.tsv"
	output: "plots/non_snakemake/all.tmerge.minByreads.minQual0.exonOverhang0.vs.SIRVs_lib{libProt}.stats.png"
	shell:
		'''
echo "library(ggplot2)
dat <- read.table('{input}', header=T, as.is=T, sep='\\t')
dat <- subset(dat, level=='Transcript')
dat <- subset(dat, minSeqQual==0)
dat <- subset(dat, exonOverhang==0)
dat <- subset(dat, metric!='SnPr')
dat\$seqM <- paste(dat\$seqTech, dat\$metric)
palette <- c('Sn' = '#cc6600', 'Pr' = '#2d8659', 'SnPr' = '#999999')
ggplot(data=dat, aes(x=minReadSupport, y=value, group=seqM, color=metric)) +
geom_point() +
#xlim(1, 51) +
scale_color_manual(values=palette) +
geom_line(aes(linetype=seqTech)) + ylab('Pr & Sn (%)') + xlab('Min read support required per TM') +
ylim(0, 100) +
{GGPLOT_PUB_QUALITY}
ggsave('{output}', width=7, height=5)
"  | R --slave
		'''



rule plotGffCompareSirvByExonOverhang:
	input: "test/Rinput/all.tmerge.vs.SIRVs_lib{libProt}.stats.tsv"
	output: "plots/non_snakemake/all.tmerge.min1reads.minQual0.exonOverhangBy.vs.SIRVs_lib{libProt}.stats.png"
	shell:
		'''
echo "library(ggplot2)
dat <- read.table('{input}', header=T, as.is=T, sep='\\t')
dat <- subset(dat, level=='Transcript')
dat <- subset(dat, minSeqQual==0)
dat <- subset(dat, minReadSupport==1)
dat <- subset(dat, metric!='SnPr')
dat\$seqM <- paste(dat\$seqTech, dat\$metric)
palette <- c('Sn' = '#cc6600', 'Pr' = '#2d8659', 'SnPr' = '#999999')
ggplot(data=dat, aes(x=exonOverhang, y=value, group=seqM, color=metric)) +
geom_point() +
#xlim(1, 51) +
scale_color_manual(values=palette) +
geom_line(aes(linetype=seqTech)) + ylab('Pr & Sn (%)') + xlab('Exon overhang tolerance (nts)') +
ylim(0, 100) +
{GGPLOT_PUB_QUALITY}
ggsave('{output}', width=7, height=5)
"  | R --slave
		'''

rule plotGffCompareSirvStats:
	input:"test/Rinput/all.tmerge.vs.SIRVs_lib{libProt}.stats.tsv"
	output: "plots/non_snakemake/all.tmerge.optimalParams.vs.SIRVs_lib{libProt}.stats.png"
	shell:
		'''
echo "library(ggplot2)
library(plyr)
library(scales)
cbPalette <- c('Sn'='#cc6600', 'Pr'='#2d8659')
dat <- read.table('{input}', header=T, as.is=T, sep='\\t')
dat\$seqTech <- gsub(':$', '', dat\$seqTech)
dat\$seqTech <- gsub(':', '\\n', dat\$seqTech)

dat <- subset(dat, minReadSupport=={optimalMINREADSUPPORT})
dat <- subset(dat, exonOverhang=={optimalEXONOVERHANG})
dat <- subset(dat, minSeqQual=={optimalMINSEQQUAL})
dat <- subset(dat, metric!='SnPr')

#dat\$levelCorrlevel <- paste(sep='', dat\$level, ' (Corr: ', dat\$correctionLevel, ')')
ggplot(dat, aes(x=level, y=value)) +
geom_point(aes(color=metric, shape=correctionLevel), size=3, alpha=0.8) +
scale_colour_manual (values=cbPalette, name='Metric', breaks=c('Sn', 'Pr'))+
scale_shape_manual(values=c(16,21), name='Error correction') +
ylab('Sn | Pr (%)') +
xlab('Evaluation level') +
scale_y_continuous() +
expand_limits(y=c(0,100))+
theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
facet_grid(capDesign + tissue ~  seqTech + sizeFrac )+
{GGPLOT_PUB_QUALITY}
ggsave('{output}', width=5.5, height=3)
" > {output}.r
cat {output}.r | R --slave
		'''
