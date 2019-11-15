rule integratePolyaAndSjInfo:
	input:
		polyA = "mappings/removePolyAERCCs/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.polyAsitesNoErcc.tmp.bed",
		SJs = "mappings/getIntronMotif/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.transcripts.tsv"
	output:
		strandInfo="mappings/integratePolyaAndSjInfo/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.polyA+SJ.strandInfo.tsv",
		wrongPolyAs="mappings/removePolyAERCCs/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.wrongPolyAs.list"
	shell:
		'''
uuid=$(uuidgen)
cat {input.SJs} | skipcomments | cut -f 1,2 | awk '$2!="."' | sort -T {config[TMPDIR]} > {config[TMPDIR]}/$uuid.reads.SJ.strandInfo.tsv
cat {input.polyA} | cut -f4,6 | awk '$2!="."'| sort -T {config[TMPDIR]}  > {config[TMPDIR]}/$uuid.reads.polyA.strandInfo.tsv
join -a1 -a2 -e '.' -o '0,1.2,2.2' {config[TMPDIR]}/$uuid.reads.SJ.strandInfo.tsv  {config[TMPDIR]}/$uuid.reads.polyA.strandInfo.tsv > {config[TMPDIR]}/$uuid.reads.SJ.polyA.strandInfo.tsv

#make list of reads with wrongly called polyA sites (i.e. their strand is different from the one inferred using SJs):
cat {config[TMPDIR]}/$uuid.reads.SJ.polyA.strandInfo.tsv | perl -slane 'if($F[2] ne "." && $F[1] ne "." && $F[1] ne $F[2]){{print join ("\\t", $F[0])}}' |sort -T {config[TMPDIR]} |uniq > {output.wrongPolyAs}
#get strand info, prioritizing SJ inference
cat {config[TMPDIR]}/$uuid.reads.SJ.polyA.strandInfo.tsv | perl -slane 'if($F[1] ne "."){{print "$F[0]\\t$F[1]"}} else{{print "$F[0]\\t$F[2]"}}' |sort -T {config[TMPDIR]} |uniq > {output.strandInfo}
		'''

rule removeWrongPolyAs:
	input:
		polyA="mappings/removePolyAERCCs/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.polyAsitesNoErcc.tmp.bed",
		wrongPolyAs="mappings/removePolyAERCCs/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.wrongPolyAs.list"
	output: "mappings/removePolyAERCCs/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.polyAsitesNoErcc.bed"
	shell:
		'''
cat {input.polyA} | fgrep -v -w -f {input.wrongPolyAs} > {output}
		'''

rule strandGffs:
	input:
		gff = "mappings/readBedToGff/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.gff.gz",
		strandInfo = "mappings/integratePolyaAndSjInfo/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.polyA+SJ.strandInfo.tsv"
	output: "mappings/strandGffs/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.stranded.gff.gz"
	shell:
		'''
uuid=$(uuidgen)
zcat {input.gff} > {config[TMPDIR]}/$uuid.in.gff
get_right_transcript_strand.pl {config[TMPDIR]}/$uuid.in.gff {input.strandInfo} | fgrep -v ERCC- | sort -T {config[TMPDIR]}  -k1,1 -k4,4n -k5,5n  |gzip> {output}

		'''

rule highConfidenceReads:
	input:
		transcriptStrandInfo = "mappings/getIntronMotif/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.transcripts.tsv",
		strandedReads = "mappings/strandGffs/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.stranded.gff.gz",
		bam = "mappings/readMapping/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bam"
	output:
		gff="mappings/highConfidenceReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.strandedHCGMs.gff.gz",
		stats=temp(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.highConfSplicedReads.stats.tsv")
	shell:
		'''
#### reads QC stats
totalSplicedReads=$(cat {input.transcriptStrandInfo} | tgrep -v -P "^#" | wc -l)
# reads with 100% canonical SJs
canonSjReads=$(cat {input.transcriptStrandInfo} | awk '$6==1'| wc -l)
# reads with no fishy SJs (i.e. not surrounded by direct repeats)
noFishySjReads=$(cat {input.transcriptStrandInfo} | awk '$7==1'| wc -l)
# reads with no fishy SJs and canonical SJs
noFishyCanonSjReads=$(cat {input.transcriptStrandInfo} | awk '$6==1 && $7==1'| wc -l)
echo -e "{wildcards.techname}Corr{wildcards.corrLevel}\t{wildcards.capDesign}\t{wildcards.sizeFrac}\t{wildcards.barcodes}\t$totalSplicedReads\t$canonSjReads\t$noFishySjReads\t$noFishyCanonSjReads" | awk '{{print $0"\\t"$6/$5"\\t"$7/$5"\\t"$8/$5}}' > {output.stats}


#select read IDs with canonical GT|GC/AG
uuid=$(uuidgen)
cat {input.transcriptStrandInfo} | tgrep -v -P "^#" | awk '$6==1' | cut -f1 | sort -T {config[TMPDIR]} |uniq > {config[TMPDIR]}/$uuid.reads.hcSJs.list.tmp
wc -l {config[TMPDIR]}/$uuid.reads.hcSJs.list.tmp
zcat {input.strandedReads} > {config[TMPDIR]}/$uuid.str.gff

##high-quality Phred SJs
samtools view -F 256 -F4 -F 2048 {input.bam} | samHQintrons.pl --minQual {config[minPhredQualAroundSJ]} - |cut -f1|sort -T {config[TMPDIR]} |uniq > {config[TMPDIR]}/$uuid.HQintrons.reads.list
tgrep -F -w -f {config[TMPDIR]}/$uuid.HQintrons.reads.list {config[TMPDIR]}/$uuid.reads.hcSJs.list.tmp > {config[TMPDIR]}/$uuid.reads.hcSJs.list

tgrep -F -w -f {config[TMPDIR]}/$uuid.reads.hcSJs.list {config[TMPDIR]}/$uuid.str.gff > {config[TMPDIR]}/$uuid.gtag.gff
wc -l {config[TMPDIR]}/$uuid.gtag.gff
cat {config[TMPDIR]}/$uuid.str.gff | extractGffAttributeValue.pl transcript_id | sort -T {config[TMPDIR]} |uniq -u > {config[TMPDIR]}/$uuid.tmp
tgrep -F -w -f {config[TMPDIR]}/$uuid.tmp {config[TMPDIR]}/$uuid.str.gff > {config[TMPDIR]}/$uuid.tmp2
cat {config[TMPDIR]}/$uuid.tmp2  > {config[TMPDIR]}/$uuid.monoPolyA.gff
 echo $?
cat {config[TMPDIR]}/$uuid.gtag.gff {config[TMPDIR]}/$uuid.monoPolyA.gff | sort -T {config[TMPDIR]}  -k1,1 -k4,4n -k5,5n  |gzip> {output.gff}
 echo $?
		'''

rule aggHighConfSplicedReadsStats:
	input: expand(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.highConfSplicedReads.stats.tsv",filtered_product_merge, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACSpluSMERGED, barcodes=BARCODESpluSMERGED)
	output: config["STATSDATADIR"] + "all.highConfSplicedReads.stats.tsv"
	shell:
		'''
echo -e "seqTech\tcorrectionLevel\tcapDesign\tsizeFrac\ttissue\ttotalSplicedReads\tcanonSjReads\tnoFishySjReads\tnoFishyCanonSjReads" > {output}
cat {input} | sed 's/Corr0/\tNo/' | sed 's/Corr{lastK}/\tYes/' | sort -T {config[TMPDIR]}  >> {output}
		'''


rule getHCGMintrons:
	input: "mappings/highConfidenceReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.strandedHCGMs.gff.gz"
	output: temp("mappings/highConfidenceReads/introns/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.strandedHCGMs.introns.tsv")
	shell:
		'''
zcat {input} | makeIntrons.pl -| perl -lane '$F[11]=~s/"|;//g; $start=$F[3]-1; $end=$F[4]+1; print $F[11]."\t".$F[0]."_".$start."_".$end."_".$F[6]' | sort -T {config[TMPDIR]}  -k2,2 > {output}

		'''

rule getHiSeqSupportedHCGMs:
	input:
		hiSeqIntrons="mappings/hiSeqIntrons/hiSeq_{capDesign}.canonicalIntrons.list" if config["USE_MATCHED_ILLUMINA"] else config["SUPPORT_INTRONS_DB"],
		#hiSeqIntrons="mappings/hiSeqIntrons/hiSeq_{techname}_{capDesign}.{barcodes}.canonicalIntrons.list",
		lrIntrons="mappings/highConfidenceReads/introns/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.strandedHCGMs.introns.tsv",
		hcgmGTF= "mappings/highConfidenceReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.strandedHCGMs.gff.gz"
	output:"mappings/highConfidenceReads/HiSS/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.gff.gz"
	shell:
		'''
uuid=$(uuidgen)
join -v1 -1 2 -2 1 {input.lrIntrons} {input.hiSeqIntrons} |awk '{{print $2"\t"$1}}' > {config[TMPDIR]}/$uuid.{wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.barcodes}.introns.noHiSeq.tsv
zcat {input.hcgmGTF} > {config[TMPDIR]}/$uuid.$(basename {input.hcgmGTF} .gz)
cut -f1 {config[TMPDIR]}/$uuid.{wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.barcodes}.introns.noHiSeq.tsv | sort -T {config[TMPDIR]} |uniq | fgrep -wv -f - {config[TMPDIR]}/$uuid.$(basename {input.hcgmGTF} .gz) |sort -T {config[TMPDIR]}  -k1,1 -k4,4n -k5,5n  |gzip > {output}

		'''

rule getHiSSStats:
	input:
		reads = "mappings/readMapping/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bam",
		HiSSGTF="mappings/highConfidenceReads/HiSS/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.gff.gz"
	output: temp(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.stats.tsv")
	shell:
		'''
uuid=$(uuidgen)
bedtools bamtobed -i {input.reads} -bed12 > {config[TMPDIR]}/$uuid.{wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.barcodes}.merged.bed
zcat {input.HiSSGTF} | gff2bed_full.pl - > {config[TMPDIR]}/$uuid.{wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.barcodes}.HiSS.bed

mappedReadsMono=$(cat {config[TMPDIR]}/$uuid.{wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.barcodes}.merged.bed | awk '$10<=1'|cut -f4 |sort -T {config[TMPDIR]} |uniq|wc -l)
mappedReadsSpliced=$(cat {config[TMPDIR]}/$uuid.{wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.barcodes}.merged.bed | awk '$10>1'|cut -f4 |sort -T {config[TMPDIR]} |uniq|wc -l)

HiSSMono=$(cat {config[TMPDIR]}/$uuid.{wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.barcodes}.HiSS.bed| awk '$10<=1'|cut -f4  | sort -T {config[TMPDIR]} |uniq|wc -l)
HiSSSpliced=$(cat {config[TMPDIR]}/$uuid.{wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.barcodes}.HiSS.bed| awk '$10>1'|cut -f4  | sort -T {config[TMPDIR]} |uniq|wc -l)

#let totalMapped=$mappedReadsMono+$mappedReadsSpliced || true
let nonHiSSMono=$mappedReadsMono-$HiSSMono || true
let nonHiSSSPliced=$mappedReadsSpliced-$HiSSSpliced || true
echo -e "{wildcards.techname}Corr{wildcards.corrLevel}\t{wildcards.capDesign}\t{wildcards.sizeFrac}\t{wildcards.barcodes}\t$HiSSMono\t$HiSSSpliced\t$nonHiSSMono\t$nonHiSSSPliced" > {output}

		'''


rule aggHiSSStats:
	input: expand(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.stats.tsv",filtered_product_merge, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACSpluSMERGED, barcodes=BARCODESpluSMERGED)
	output: config["STATSDATADIR"] + "all.HiSS.stats.tsv"

	shell:
		'''
echo -e "seqTech\tcorrectionLevel\tcapDesign\tsizeFrac\ttissue\tcategory\tcount" > {output}
cat {input} | awk '{{print $1"\\t"$2"\\t"$3"\\t"$4"\\tHCGM-mono\\t"$5"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\tHCGM-spliced\\t"$6"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\tnonHCGM-mono\\t"$7"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\tnonHCGM-spliced\\t"$8}}'| sed 's/Corr0/\tNo/' | sed 's/Corr{lastK}/\tYes/' | sort -T {config[TMPDIR]}  >> {output}

		'''

rule plotAllHiSSStats:
	input: config["STATSDATADIR"] + "all.HiSS.stats.tsv"
	output:  config["PLOTSDIR"] + "HiSS.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.stats.{ext}"
	params:
		filterDat=lambda wildcards: merge_figures_params(wildcards.capDesign, wildcards.sizeFrac, wildcards.barcodes, wildcards.corrLevel, wildcards.techname)
	shell:
		'''
echo "library(ggplot2)
library(plyr)
library(scales)
dat <- read.table('{input}', header=T, as.is=T, sep='\\t')
dat\$category<-factor(dat\$category, ordered=TRUE, levels=rev(c('HCGM-mono', 'HCGM-spliced', 'nonHCGM-mono', 'nonHCGM-spliced')))
{params.filterDat[10]}
{params.filterDat[0]}
{params.filterDat[1]}
{params.filterDat[2]}
{params.filterDat[3]}
{params.filterDat[4]}
{params.filterDat[5]}
{params.filterDat[8]}

ggplot(dat[order(dat\$category), ], aes(x=factor(correctionLevel), y=count, fill=category)) +
geom_bar(stat='identity') +
scale_fill_manual(values=c('HCGM-mono' = '#9ce2bb', 'HCGM-spliced' = '#39c678', 'nonHCGM-mono' = '#fda59b', 'nonHCGM-spliced' = '#fa341e')) +
facet_grid( seqTech + sizeFrac ~ capDesign + tissue, scales = 'free_y')+ ylab('# mapped reads') +
xlab('{params.filterDat[6]}') +
guides(fill = guide_legend(title='Category'))+
geom_text(position = 'stack', size=geom_textSize, aes(x = factor(correctionLevel), y = count, label = comma(count), hjust = 0.5, vjust = 1))+
scale_y_continuous(labels=scientific)+
{params.filterDat[7]}
{GGPLOT_PUB_QUALITY}
ggsave('{output}', width=plotWidth, height=plotHeight)
" > {output}.r
cat {output}.r | R --slave

		'''



rule nonAnchoredMergeReads:
	input: "mappings/highConfidenceReads/HiSS/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.gff.gz"
	output: "mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.all.gff",
	threads:1
	wildcard_constraints:
		barcodes='(?!allTissues).+',
		sizeFrac='[0-9-+\.]+',
		techname='(?!allSeqTechs).+'
	shell:
		'''
zcat {input} | tmerge --exonOverhangTolerance {config[exonOverhangTolerance]} --minReadSupport {wildcards.minReadSupport} --endFuzz {config[exonOverhangTolerance]} --tmPrefix {wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.barcodes}.NAM_ - |sort -T {config[TMPDIR]}  -k1,1 -k4,4n -k5,5n  > {output}
		'''

rule mergeTissuesNonAnchoredMergeReads:
	input: lambda wildcards: expand("mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.all.gff", filtered_product, techname=wildcards.techname, corrLevel=wildcards.corrLevel, capDesign=wildcards.capDesign, sizeFrac=wildcards.sizeFrac, barcodes=BARCODES, minReadSupport=wildcards.minReadSupport)
	output: "mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.all.gff"
	threads:1
	wildcard_constraints:
		barcodes='allTissues',
		sizeFrac='[0-9-+\.]+',
		techname='(?!allSeqTechs).+'
	shell:
		'''
uuid=$(uuidgen)
cat {input} |sort -T {config[TMPDIR]}  -k1,1 -k4,4n -k5,5n  > {config[TMPDIR]}/$uuid
cat {config[TMPDIR]}/$uuid | tmerge --exonOverhangTolerance {config[exonOverhangTolerance]} --tmPrefix {wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.barcodes}.NAM_ - |sort -T {config[TMPDIR]}  -k1,1 -k4,4n -k5,5n  > {output}
#cp {output} {output}.bkp
#checkTmergeOutput.sh {config[TMPDIR]}/$uuid {output} &> {output}.qc.txt
#rm {output}.bkp

		'''

rule mergeFracsNonAnchoredMergeReads:
	input: lambda wildcards: expand("mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.all.gff", filtered_product, techname=wildcards.techname, corrLevel=wildcards.corrLevel, capDesign=wildcards.capDesign, sizeFrac=SIZEFRACS, barcodes=wildcards.barcodes, minReadSupport=wildcards.minReadSupport)
	output: "mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.all.gff"
	threads:1
	wildcard_constraints:
		sizeFrac='allFracs',
		barcodes='(?!allTissues).+',
		techname='(?!allSeqTechs).+'
	shell:
		'''
uuid=$(uuidgen)
cat {input} |sort -T {config[TMPDIR]}  -k1,1 -k4,4n -k5,5n  > {config[TMPDIR]}/$uuid
cat {config[TMPDIR]}/$uuid | tmerge --exonOverhangTolerance {config[exonOverhangTolerance]} --tmPrefix {wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.barcodes}.NAM_ - |sort -T {config[TMPDIR]}  -k1,1 -k4,4n -k5,5n  > {output}
#cp {output} {output}.bkp
#checkTmergeOutput.sh {config[TMPDIR]}/$uuid {output} &> {output}.qc.txt
#rm {output}.bkp
		'''


rule mergeFracsAndTissuesNonAnchoredMergeReads:
	input:
		lambda wildcards: expand("mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.all.gff", filtered_product, techname=wildcards.techname, corrLevel=wildcards.corrLevel, capDesign=wildcards.capDesign, sizeFrac=SIZEFRACS, barcodes=BARCODES, minReadSupport=wildcards.minReadSupport),

	output: "mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_allFracs_allTissues.tmerge.min{minReadSupport}reads.all.gff"
	threads:1
	wildcard_constraints:
		techname='(?!allSeqTechs).+',
		capDesign='|'.join(CAPDESIGNS)
	shell:
		'''
uuid=$(uuidgen)
cat {input} |sort -T {config[TMPDIR]}  -k1,1 -k4,4n -k5,5n  > {config[TMPDIR]}/$uuid
cat {config[TMPDIR]}/$uuid | tmerge --exonOverhangTolerance {config[exonOverhangTolerance]} --tmPrefix {wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}_allFracs_allTissues.NAM_ - |sort -T {config[TMPDIR]}  -k1,1 -k4,4n -k5,5n  > {output}
#cp {output} {output}.bkp
#checkTmergeOutput.sh {config[TMPDIR]}/$uuid {output} &> {output}.qc.txt
#rm {output}.bkp

		'''

rule mergeAllSeqTechsFracsAndTissuesNonAnchoredMergeReads:
	input:
		lambda wildcards: expand("mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.all.gff", filtered_product_merge, techname=TECHNAMES, corrLevel=wildcards.corrLevel, capDesign=wildcards.capDesign, sizeFrac=wildcards.sizeFrac, barcodes=wildcards.barcodes, minReadSupport=wildcards.minReadSupport),

	output: "mappings/nonAnchoredMergeReads/allSeqTechsCorr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.all.gff"
	threads:1
	wildcard_constraints:
#		techname='allSeqTechs',
		sizeFrac='allFracs',
		barcodes='allTissues',
		capDesign='|'.join(CAPDESIGNS)
	shell:
		'''
uuid=$(uuidgen)
cat {input} |sort -T {config[TMPDIR]}  -k1,1 -k4,4n -k5,5n  > {config[TMPDIR]}/$uuid
cat {config[TMPDIR]}/$uuid | tmerge --exonOverhangTolerance {config[exonOverhangTolerance]} --tmPrefix allSeqTechsCorr{wildcards.corrLevel}_{wildcards.capDesign}_allFracs_allTissues.NAM_ - |sort -T {config[TMPDIR]}  -k1,1 -k4,4n -k5,5n  > {output}
#cp {output} {output}.bkp
#checkTmergeOutput.sh {config[TMPDIR]}/$uuid {output} &> {output}.qc.txt
#rm {output}.bkp

		'''


rule mergeAllCapDesignsSeqTechsFracsAndTissuesNonAnchoredMergeReads:
	input:
		lambda wildcards: expand("mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.all.gff", filtered_capDesign_product_merge, techname=wildcards.techname, corrLevel=wildcards.corrLevel, capDesign=wildcards.capDesign, sizeFrac=wildcards.sizeFrac, barcodes=wildcards.barcodes, minReadSupport=wildcards.minReadSupport),

	output: "mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.all.gff"
	threads:1
	wildcard_constraints:
		sizeFrac='allFracs',
		barcodes='allTissues',
		capDesign='|'.join(MERGEDCAPDESIGNS)
	shell:
		'''
uuid=$(uuidgen)
cat {input} |sort -T {config[TMPDIR]}  -k1,1 -k4,4n -k5,5n  > {config[TMPDIR]}/$uuid
cat {config[TMPDIR]}/$uuid | tmerge --exonOverhangTolerance {config[exonOverhangTolerance]} --tmPrefix allSeqTechsCorr{wildcards.corrLevel}_{wildcards.capDesign}_allFracs_allTissues.NAM_ - |sort -T {config[TMPDIR]}  -k1,1 -k4,4n -k5,5n  > {output}
#cp {output} {output}.bkp
#checkTmergeOutput.sh {config[TMPDIR]}/$uuid {output} &> {output}.qc.txt
#rm {output}.bkp

		'''




rule nonAnchoredMergeReadsToBed:
	input: "mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.all.gff"
	output: temp("mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.all.bed")
	shell:
		'''
cat {input} | gff2bed_full.pl - > {output}
		'''



rule getTmStats:
	input: "mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.{endSupport}.gff"
	output: temp(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.{endSupport}.TmStats.stats.tsv")
	# wildcard_constraints:
	# 	barcodes='(?!allTissues).+',
	# 	sizeFrac='[0-9-+\.]+',
	# 	techname='(?!allSeqTechs).+'
	shell:
		'''
cat {input} | extractGffAttributeValue.pl transcript_id spliced mature_RNA_length contains_count 3p_dists_to_3p 5p_dists_to_5p meta_3p_dists_to_5p meta_5p_dists_to_5p |sort|uniq | awk -v t={wildcards.techname} -v c={wildcards.corrLevel} -v ca={wildcards.capDesign} -v s={wildcards.sizeFrac} -v b={wildcards.barcodes} '{{print t"Corr"c"\\t"ca"\\t"s"\\t"b"\t"$0}}'> {output}
		'''

rule aggTmStats:
	input: lambda wildcards: expand(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.{endSupport}.TmStats.stats.tsv", filtered_product, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES, minReadSupport=wildcards.minReadSupport, endSupport=wildcards.endSupport)
	output: config["STATSDATADIR"] + "all.min{minReadSupport}.{endSupport}.TmStats.stats.tsv"
	shell:
		'''
echo -e "seqTech\tcorrectionLevel\tcapDesign\tsizeFrac\ttissue\ttranscript_id\tspliced\tcontains_count\tend\tdistance\tnormDistance" > {output}

cat {input} | perl -F"\\t" -slane '@ara=split(",", $F[8]); @arb=split(",", $F[9]); @arc=split(",", $F[10]), @ard=split(",", $F[11]); for ($i=0; $i<=$#ara; $i++){{$threepMinusDist=-$ara[$i];print  "$F[0]\\t$F[1]\\t$F[2]\\t$F[3]\\t$F[4]\\t$F[5]\\t$F[7]\\t3p\\t$threepMinusDist\\t$arc[$i]\\n$F[0]\\t$F[1]\\t$F[2]\\t$F[3]\\t$F[4]\\t$F[5]\\t$F[7]\\t5p\\t$arb[$i]\\t$ard[$i]"}}'| sed 's/Corr0/\\tNo/' | sed 's/Corr{lastK}/\\tYes/' | sort -T {config[TMPDIR]}  >> {output}

		'''



rule getMergingStats:
	input:
		hcgms = "mappings/highConfidenceReads/HiSS/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.gff.gz",
		pooledMerged = "mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.all.gff"
	output: temp(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.merged.stats.tsv")
	shell:
		'''
hcgms=$(zcat {input.hcgms} | extractGffAttributeValue.pl transcript_id | sort -T {config[TMPDIR]} |uniq|wc -l)
merged=$(cat {input.pooledMerged} | extractGffAttributeValue.pl transcript_id | sort -T {config[TMPDIR]} |uniq|wc -l)
echo -e "{wildcards.techname}Corr{wildcards.corrLevel}\t{wildcards.capDesign}\t{wildcards.sizeFrac}\t{wildcards.barcodes}\t$hcgms\t$merged" | awk '{{print $0"\t"$6/$5}}' > {output}

		'''


rule aggMergingStats:
	input: lambda wildcards: expand(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.merged.stats.tsv",filtered_product_merge, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACSpluSMERGED, barcodes=BARCODESpluSMERGED, minReadSupport=wildcards.minReadSupport)

	output: config["STATSDATADIR"] + "all.min{minReadSupport}reads.merged.stats.tsv"
	shell:
		'''
echo -e "seqTech\tcorrectionLevel\tcapDesign\tsizeFrac\ttissue\tcategory\tcount" > {output}
cat {input} | awk '{{print $1"\\t"$2"\\t"$3"\\t"$4"\\tHCGMreads\\t"$5"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\tmergedTMs\\t"$6}}' | sed 's/Corr0/\tNo/' | sed 's/Corr{lastK}/\tYes/' | sort -T {config[TMPDIR]}  >> {output}

		'''
rule plotMergingStats:
	input:  config["STATSDATADIR"] + "all.min{minReadSupport}reads.merged.stats.tsv"
	output: config["PLOTSDIR"] + "merged.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.merged.stats.{ext}"
	params:
		filterDat=lambda wildcards: merge_figures_params(wildcards.capDesign, wildcards.sizeFrac, wildcards.barcodes, wildcards.corrLevel, wildcards.techname)
	shell:
		'''
echo "library(ggplot2)
library(plyr)
library(scales)
dat <- read.table('{input}', header=T, as.is=T, sep='\\t')
{params.filterDat[10]}
{params.filterDat[0]}
{params.filterDat[1]}
{params.filterDat[2]}
{params.filterDat[3]}
{params.filterDat[4]}
{params.filterDat[5]}
{params.filterDat[8]}


ggplot(data=dat, aes(x=factor(correctionLevel), y=count, fill=category)) +
geom_bar(stat='identity', position=position_dodge()) +
geom_text(position = position_dodge(width = 0.9), size=geom_textSize, aes(x = factor(correctionLevel), y = 1, label = comma(count), hjust = 0, vjust = 0.5), angle=90) +
scale_fill_manual(values=c('HCGMreads' = '#d98cb3', 'mergedTMs' = '#cc9966')) +
facet_grid( seqTech + sizeFrac ~ capDesign + tissue)+ ylab('# objects') +
xlab('{params.filterDat[6]}') +
guides(fill = guide_legend(title='Category'))+
scale_y_continuous(labels=comma)+
{params.filterDat[7]}
{GGPLOT_PUB_QUALITY}
ggsave('{output}', width=plotWidth, height=plotHeight)
" > {output}.r
cat {output}.r | R --slave

		'''

# rule getTmLengthStats:
# 	input:
# 		gencode=lambda wildcards: CAPDESIGNTOANNOTGTF[wildcards.capDesign],
# 		tms="mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.all.bed",
# 		flTms="mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.cage+polyASupported.bed"
# 	output:temp(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.splicedLength.stats.tsv")
# 	shell:
# 		'''
# uuid=$(uuidgen)
# cat {input.gencode} | awk '$3=="exon"' | fgrep "transcript_type \\"protein_coding\\";" |gff2bed_full.pl - | bed12ToTranscriptLength.pl - | awk -v t={wildcards.techname}Corr{wildcards.corrLevel} -v c={wildcards.capDesign} -v si={wildcards.sizeFrac} -v b={wildcards.barcodes} '{{print t"\t"c"\t"si"\t"b"\tGENCODE_protein_coding\t"$2}}'| sed 's/Corr0/\tNo/' | sed 's/Corr{lastK}/\tYes/' > {config[TMPDIR]}/$uuid.gencode.pcg.tsv
#
# cat {input.tms} | bed12ToTranscriptLength.pl - | awk -v t={wildcards.techname}Corr{wildcards.corrLevel} -v c={wildcards.capDesign} -v si={wildcards.sizeFrac} -v b={wildcards.barcodes} '{{print t"\t"c"\t"si"\t"b"\tCLS_TMs\t"$2}}'| sed 's/Corr0/\tNo/' | sed 's/Corr{lastK}/\tYes/'  > {config[TMPDIR]}/$uuid.tms.tsv
#
# cat {input.flTms} | bed12ToTranscriptLength.pl - | awk -v t={wildcards.techname}Corr{wildcards.corrLevel} -v c={wildcards.capDesign} -v si={wildcards.sizeFrac} -v b={wildcards.barcodes} '{{print t"\t"c"\t"si"\t"b"\tCLS_FL_TMs\t"$2}}'| sed 's/Corr0/\tNo/' | sed 's/Corr{lastK}/\tYes/'  > {config[TMPDIR]}/$uuid.flTms.tsv
#
# cat {config[TMPDIR]}/$uuid.gencode.pcg.tsv {config[TMPDIR]}/$uuid.tms.tsv {config[TMPDIR]}/$uuid.flTms.tsv | sed 's/Corr0/\tNo/' | sed 's/Corr{lastK}/\tYes/'  > {output}
# 		'''
#
#

rule aggTmLengthStats:
	input:
		all=lambda wildcards: expand(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.all.TmStats.stats.tsv", filtered_product, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES, minReadSupport=wildcards.minReadSupport),
		fl=lambda wildcards: expand(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.cage+polyASupported.TmStats.stats.tsv", filtered_product, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES, minReadSupport=wildcards.minReadSupport)
	output: config["STATSDATADIR"] + "all.min{minReadSupport}reads.splicedLength.stats.tsv"
	shell:
		'''
echo -e "seqTech\tcorrectionLevel\tcapDesign\tsizeFrac\ttissue\ttranscript_id\tspliced\tmature_RNA_length\tcategory" > {output}

cat {input.all} |cut -f1-7| awk '{{print $0"\\tCLS_TMs"}}' | sed 's/Corr0/\\tNo/' | sed 's/Corr{lastK}/\\tYes/' |sort -T {config[TMPDIR]}  >> {output}
cat {input.fl} |cut -f1-7| awk '{{print $0"\\tCLS_FL_TMs"}}' | sed 's/Corr0/\\tNo/' | sed 's/Corr{lastK}/\\tYes/' |sort -T {config[TMPDIR]}  >> {output}


		'''


rule plotTmLengthStats:
	input: config["STATSDATADIR"] + "all.min{minReadSupport}reads.splicedLength.stats.tsv"
	output: config["PLOTSDIR"] + "splicedLength.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.splicing_status:{splicedStatus}.splicedLength.stats.{ext}"
	params:
		filterDat=lambda wildcards: merge_figures_params(wildcards.capDesign, wildcards.sizeFrac, wildcards.barcodes, wildcards.corrLevel, wildcards.techname, wildcards.splicedStatus)
	shell:
		'''
echo "
library(data.table)
library(ggplot2)
library(scales)
library(plyr)
palette <- c('GENCODE_protein_coding' = '#009900', 'CLS_TMs' = '#cc9966', 'CLS_FL_TMs' = '#cc00cc')
dat<-fread('{input}', header=T, sep='\\t')
{params.filterDat[10]}
{params.filterDat[0]}
{params.filterDat[1]}
{params.filterDat[2]}
{params.filterDat[3]}
{params.filterDat[4]}
{params.filterDat[5]}
{params.filterDat[8]}
{params.filterDat[11]}


fun_length <- function(x){{
return(data.frame(y=-8.5,label= paste0('N=', comma(length(x)))))
}}
fun_median <- function(x){{
return(data.frame(y=-8.5,label= paste0('Median=', comma(median(x)))))
}}
ggplot(dat, aes(x=factor(correctionLevel), y=mature_RNA_length, color=category)) +
geom_boxplot(position=position_dodge(0.9), outlier.shape=NA) +
coord_cartesian(ylim=c(100, 5000)) +
scale_y_continuous(labels=comma)+
scale_color_manual(values=palette, name='Category', labels = c('GENCODE_protein_coding' = 'GENCODE\nprotein-coding', 'CLS_TMs'='CLS TMs', 'CLS_FL_TMs'='CLS FL TMs')) +
facet_grid( seqTech + sizeFrac ~ capDesign + tissue)+

stat_summary(aes(x=factor(correctionLevel), group=category), position=position_dodge(0.9), fun.data = fun_length, geom = 'text', vjust = +1, hjust=0, show.legend=FALSE, size=geom_textSize, color='#666666') +
stat_summary(aes(x=factor(correctionLevel), group=category), position=position_dodge(0.9), fun.data = fun_median, geom = 'text', vjust = 0, hjust=0, show.legend=FALSE, size=geom_textSize, color='#666666') +
ylab('Mature RNA length') +
xlab('{params.filterDat[6]}') +
coord_flip(ylim=c(100, 5000)) +
{params.filterDat[9]}
{GGPLOT_PUB_QUALITY} +
theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('{output}', width=plotWidth, height=plotHeight)

" > {output}.r
cat {output}.r | R --slave

		'''
