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
uuidTmpOutS=$(uuidgen)
uuidTmpOutW=$(uuidgen)
cat {input.SJs} | skipcomments | cut -f 1,2 | awk '$2!="."' | sort -T {config[TMPDIR]} > {config[TMPDIR]}/$uuid.reads.SJ.strandInfo.tsv
cat {input.polyA} | cut -f4,6 | awk '$2!="."'| sort -T {config[TMPDIR]}  > {config[TMPDIR]}/$uuid.reads.polyA.strandInfo.tsv
join -a1 -a2 -e '.' -o '0,1.2,2.2' {config[TMPDIR]}/$uuid.reads.SJ.strandInfo.tsv  {config[TMPDIR]}/$uuid.reads.polyA.strandInfo.tsv > {config[TMPDIR]}/$uuid.reads.SJ.polyA.strandInfo.tsv

#make list of reads with wrongly called polyA sites (i.e. their strand is different from the one inferred using SJs):
cat {config[TMPDIR]}/$uuid.reads.SJ.polyA.strandInfo.tsv | perl -slane 'if($F[2] ne "." && $F[1] ne "." && $F[1] ne $F[2]){{print join ("\\t", $F[0])}}' |sort -T {config[TMPDIR]} |uniq > {config[TMPDIR]}/$uuidTmpOutW
mv {config[TMPDIR]}/$uuidTmpOutW {output.wrongPolyAs}
#get strand info, prioritizing SJ inference
cat {config[TMPDIR]}/$uuid.reads.SJ.polyA.strandInfo.tsv | perl -slane 'if($F[1] ne "."){{print "$F[0]\\t$F[1]"}} else{{print "$F[0]\\t$F[2]"}}' |sort -T {config[TMPDIR]} |uniq > {config[TMPDIR]}/$uuidTmpOutS
mv {config[TMPDIR]}/$uuidTmpOutS {output.strandInfo}
		'''

rule removeWrongPolyAs:
	input:
		polyA="mappings/removePolyAERCCs/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.polyAsitesNoErcc.tmp.bed",
		wrongPolyAs="mappings/removePolyAERCCs/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.wrongPolyAs.list"
	output: "mappings/removePolyAERCCs/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.polyAsitesNoErcc.bed"
	shell:
		'''
uuidTmpOut=$(uuidgen)
cat {input.polyA} | fgrep -v -w -f {input.wrongPolyAs} > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
		'''

rule strandGffs:
	input:
		gff = "mappings/readBedToGff/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.gff.gz",
		strandInfo = "mappings/integratePolyaAndSjInfo/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.polyA+SJ.strandInfo.tsv"
	output: "mappings/strandGffs/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.stranded.gff.gz"
	shell:
		'''
uuid=$(uuidgen)
uuidTmpOut=$(uuidgen)
zcat {input.gff} > {config[TMPDIR]}/$uuid.in.gff
get_right_transcript_strand.pl {config[TMPDIR]}/$uuid.in.gff {input.strandInfo} | fgrep -v ERCC- | sort -T {config[TMPDIR]}  -k1,1 -k4,4n -k5,5n  |gzip> {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}

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
uuidTmpOutG=$(uuidgen)
uuidTmpOutS=$(uuidgen)
#### reads QC stats
totalSplicedReads=$(cat {input.transcriptStrandInfo} | tgrep -v -P "^#" | wc -l)
# reads with 100% canonical SJs
canonSjReads=$(cat {input.transcriptStrandInfo} | awk '$6==1'| wc -l)
# reads with no fishy SJs (i.e. not surrounded by direct repeats)
noFishySjReads=$(cat {input.transcriptStrandInfo} | awk '$7==1'| wc -l)
# reads with no fishy SJs and canonical SJs
noFishyCanonSjReads=$(cat {input.transcriptStrandInfo} | awk '$6==1 && $7==1'| wc -l)
echo -e "{wildcards.techname}Corr{wildcards.corrLevel}\t{wildcards.capDesign}\t{wildcards.sizeFrac}\t{wildcards.barcodes}\t$totalSplicedReads\t$canonSjReads\t$noFishySjReads\t$noFishyCanonSjReads" | awk '{{print $0"\\t"$6/$5"\\t"$7/$5"\\t"$8/$5}}' > {config[TMPDIR]}/$uuidTmpOutS
mv {config[TMPDIR]}/$uuidTmpOutS {output.stats}


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
cat {config[TMPDIR]}/$uuid.gtag.gff {config[TMPDIR]}/$uuid.monoPolyA.gff | sort -T {config[TMPDIR]}  -k1,1 -k4,4n -k5,5n  |gzip> {config[TMPDIR]}/$uuidTmpOutG
mv {config[TMPDIR]}/$uuidTmpOutG {output.gff}
		'''

rule aggHighConfSplicedReadsStats:
	input: expand(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.highConfSplicedReads.stats.tsv",filtered_product_merge, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODESpluSMERGED)
	output: config["STATSDATADIR"] + "all.highConfSplicedReads.stats.tsv"
	shell:
		'''
uuidTmpOut=$(uuidgen)
echo -e "seqTech\tcorrectionLevel\tcapDesign\tsizeFrac\ttissue\ttotalSplicedReads\tcanonSjReads\tnoFishySjReads\tnoFishyCanonSjReads" > {config[TMPDIR]}/$uuidTmpOut
cat {input} | sed 's/Corr0/\tNo/' | sed 's/Corr{lastK}/\tYes/' | sort -T {config[TMPDIR]}  >> {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
		'''


rule getHCGMintrons:
	input: "mappings/highConfidenceReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.strandedHCGMs.gff.gz"
	output: temp("mappings/highConfidenceReads/introns/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.strandedHCGMs.introns.tsv")
	shell:
		'''
uuidTmpOut=$(uuidgen)
zcat {input} | makeIntrons.pl -| perl -lane '$F[11]=~s/"|;//g; $start=$F[3]-1; $end=$F[4]+1; print $F[11]."\t".$F[0]."_".$start."_".$end."_".$F[6]' | sort -T {config[TMPDIR]}  -k2,2 > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}

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
uuidTmpOut=$(uuidgen)
join -v1 -1 2 -2 1 {input.lrIntrons} {input.hiSeqIntrons} |awk '{{print $2"\t"$1}}' > {config[TMPDIR]}/$uuid.{wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.barcodes}.introns.noHiSeq.tsv
zcat {input.hcgmGTF} > {config[TMPDIR]}/$uuid.$(basename {input.hcgmGTF} .gz)
cut -f1 {config[TMPDIR]}/$uuid.{wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.barcodes}.introns.noHiSeq.tsv | sort -T {config[TMPDIR]} |uniq | fgrep -wv -f - {config[TMPDIR]}/$uuid.$(basename {input.hcgmGTF} .gz) |sort -T {config[TMPDIR]}  -k1,1 -k4,4n -k5,5n  |gzip > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}

		'''

rule getHiSSStats:
	input:
		reads = "mappings/readMapping/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bam",
		HiSSGTF="mappings/highConfidenceReads/HiSS/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.gff.gz"
	output: temp(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.stats.tsv")
	shell:
		'''
uuid=$(uuidgen)
uuidTmpOut=$(uuidgen)
bedtools bamtobed -i {input.reads} -bed12 > {config[TMPDIR]}/$uuid.{wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.barcodes}.merged.bed
zcat {input.HiSSGTF} | gff2bed_full.pl - > {config[TMPDIR]}/$uuid.{wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.barcodes}.HiSS.bed

mappedReadsMono=$(cat {config[TMPDIR]}/$uuid.{wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.barcodes}.merged.bed | awk '$10<=1'|cut -f4 |sort -T {config[TMPDIR]} |uniq|wc -l)
mappedReadsSpliced=$(cat {config[TMPDIR]}/$uuid.{wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.barcodes}.merged.bed | awk '$10>1'|cut -f4 |sort -T {config[TMPDIR]} |uniq|wc -l)

HiSSMono=$(cat {config[TMPDIR]}/$uuid.{wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.barcodes}.HiSS.bed| awk '$10<=1'|cut -f4  | sort -T {config[TMPDIR]} |uniq|wc -l)
HiSSSpliced=$(cat {config[TMPDIR]}/$uuid.{wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.barcodes}.HiSS.bed| awk '$10>1'|cut -f4  | sort -T {config[TMPDIR]} |uniq|wc -l)

#let totalMapped=$mappedReadsMono+$mappedReadsSpliced || true
let nonHiSSMono=$mappedReadsMono-$HiSSMono || true
let nonHiSSSPliced=$mappedReadsSpliced-$HiSSSpliced || true
echo -e "{wildcards.techname}Corr{wildcards.corrLevel}\t{wildcards.capDesign}\t{wildcards.sizeFrac}\t{wildcards.barcodes}\t$HiSSMono\t$HiSSSpliced\t$nonHiSSMono\t$nonHiSSSPliced" > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}

		'''


rule aggHiSSStats:
	input: expand(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.stats.tsv",filtered_product_merge, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODESpluSMERGED)
	output: config["STATSDATADIR"] + "all.HiSS.stats.tsv"

	shell:
		'''
uuidTmpOut=$(uuidgen)
echo -e "seqTech\tcorrectionLevel\tcapDesign\tsizeFrac\ttissue\tcategory\tcount" > {config[TMPDIR]}/$uuidTmpOut
cat {input} | awk '{{print $1"\\t"$2"\\t"$3"\\t"$4"\\tHCGM-mono\\t"$5"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\tHCGM-spliced\\t"$6"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\tnonHCGM-mono\\t"$7"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\tnonHCGM-spliced\\t"$8}}'| sed 's/Corr0/\tNo/' | sed 's/Corr{lastK}/\tYes/' | sort -T {config[TMPDIR]}  >> {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}

		'''

rule plotAllHiSSStats:
	input: config["STATSDATADIR"] + "all.HiSS.stats.tsv"
	output:  returnPlotFilenames(config["PLOTSDIR"] + "HiSS.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.stats")
	params:
		filterDat=lambda wildcards: merge_figures_params(wildcards.capDesign, wildcards.sizeFrac, wildcards.barcodes, wildcards.corrLevel, wildcards.techname)
	shell:
		'''
echo "
library(cowplot)
library(plyr)
library(scales)
library(gridExtra)
library(grid)
library(ggplotify)
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

plotBase <- \\"ggplot(dat[order(dat\$category), ], aes(x=factor(correctionLevel), y=count, fill=category)) +
geom_bar(stat='identity') +
scale_fill_manual(values=c('HCGM-mono' = '#9ce2bb', 'HCGM-spliced' = '#39c678', 'nonHCGM-mono' = '#fda59b', 'nonHCGM-spliced' = '#fa341e')) +
ylab('# mapped reads') +
xlab('{params.filterDat[6]}') +
guides(fill = guide_legend(title='Category'))+
#geom_text(position = 'stack', size=geom_textSize, aes(x = factor(correctionLevel), y = count, label = comma(count), hjust = 0.5, vjust = 1))+
scale_y_continuous(labels=scientific)+
{params.filterDat[7]}
{GGPLOT_PUB_QUALITY} + \\"

{params.filterDat[12]}

save_plot('{output[0]}', legendOnly, base_width=wLegendOnly, base_height=hLegendOnly)
save_plot('{output[1]}', legendOnly, base_width=wLegendOnly, base_height=hLegendOnly)

save_plot('{output[2]}', pXy, base_width=wXyPlot, base_height=hXyPlot)
save_plot('{output[3]}', pXy, base_width=wXyPlot, base_height=hXyPlot)

save_plot('{output[4]}', pXyNoLegend, base_width=wXyNoLegendPlot, base_height=hXyNoLegendPlot)
save_plot('{output[5]}', pXyNoLegend, base_width=wXyNoLegendPlot, base_height=hXyNoLegendPlot)

save_plot('{output[6]}', pYx, base_width=wYxPlot, base_height=hYxPlot)
save_plot('{output[7]}', pYx, base_width=wYxPlot, base_height=hYxPlot)

save_plot('{output[8]}', pYxNoLegend, base_width=wYxNoLegendPlot, base_height=hYxNoLegendPlot)
save_plot('{output[9]}', pYxNoLegend, base_width=wYxNoLegendPlot, base_height=hYxNoLegendPlot)


" > $(dirname {output[0]})/$(basename {output[0]} .legendOnly.png).r
cat $(dirname {output[0]})/$(basename {output[0]} .legendOnly.png).r | R --slave


		'''



rule nonAnchoredMergeReads:
	input: "mappings/highConfidenceReads/HiSS/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.gff.gz"
	output: temp("mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.tmerge.min{minReadSupport}reads.endSupport:all.gff"),
	threads:1
	wildcard_constraints:
		barcodes='(?!allTissues).+',
		sizeFrac='[0-9-+\.]+',
		techname='(?!allSeqTechs).+'
	shell:
		'''
uuidTmpOut=$(uuidgen)
zcat {input} | tmerge --exonOverhangTolerance {config[exonOverhangTolerance]} --minReadSupport {wildcards.minReadSupport} --endFuzz {config[exonOverhangTolerance]} --tmPrefix {wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.barcodes}.NAM_ - |sort -T {config[TMPDIR]}  -k1,1 -k4,4n -k5,5n  > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
		'''

rule nonAnchoredMergeUnfilteredSirvReads:
	input: "mappings/strandGffs/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.stranded.gff.gz"
	output: "mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.noFilt.tmerge.min{minReadSupport}reads.splicing_status:all.endSupport:all.gff",
	threads:1
	wildcard_constraints:
		barcodes='(?!allTissues).+',
		sizeFrac='[0-9-+\.]+',
		techname='(?!allSeqTechs).+'
	shell:
		'''
uuidTmpOut=$(uuidgen)
uuid=$(uuidgen)
zcat {input} | awk '$1=="SIRVome_isoforms"' > {config[TMPDIR]}/$uuid

cat {config[TMPDIR]}/$uuid| tmerge --minReadSupport {wildcards.minReadSupport} --tmPrefix {wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.barcodes}.NAM_ - |sort -T {config[TMPDIR]}  -k1,1 -k4,4n -k5,5n  > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
		'''



rule splitTmsBySplicedStatus:
	input: "mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.tmerge.min{minReadSupport}reads.endSupport:all.gff"
	output: "mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.tmerge.min{minReadSupport}reads.splicing_status:{splicedStatus}.endSupport:all.gff"
	params:
		grepSpliced = lambda wildcards: '| fgrep \'spliced \"1\"\'' if wildcards.splicedStatus == "spliced" else '| fgrep \'spliced \"0\"\'' if wildcards.splicedStatus == "unspliced" else ''
	shell:
		'''
uuidTmpOut=$(uuidgen)
cat {input} {params.grepSpliced} > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
		'''

rule mergeTissuesNonAnchoredMergeReads:
	input: lambda wildcards: expand("mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.tmerge.min{minReadSupport}reads.splicing_status:all.endSupport:all.gff", filtered_product, techname=wildcards.techname, corrLevel=wildcards.corrLevel, capDesign=wildcards.capDesign, sizeFrac=wildcards.sizeFrac, barcodes=BARCODES, minReadSupport=wildcards.minReadSupport)
	output: "mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.tmerge.min{minReadSupport}reads.splicing_status:all.endSupport:all.gff"
	threads:1
	wildcard_constraints:
		barcodes='allTissues',
		sizeFrac='[0-9-+\.]+',
		techname='(?!allSeqTechs).+'
	shell:
		'''
uuidTmpOut=$(uuidgen)
uuid=$(uuidgen)
cat {input} |sort -T {config[TMPDIR]}  -k1,1 -k4,4n -k5,5n  > {config[TMPDIR]}/$uuid
cat {config[TMPDIR]}/$uuid | tmerge --exonOverhangTolerance {config[exonOverhangTolerance]} --tmPrefix {wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.barcodes}.NAM_ - |sort -T {config[TMPDIR]}  -k1,1 -k4,4n -k5,5n  > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}


		'''


rule mergeAllSeqTechsFracsAndTissuesNonAnchoredMergeReads:
	input:
		lambda wildcards: expand("mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.tmerge.min{minReadSupport}reads.splicing_status:all.endSupport:all.gff", filtered_product_merge, techname=TECHNAMES, corrLevel=wildcards.corrLevel, capDesign=wildcards.capDesign, sizeFrac=wildcards.sizeFrac, barcodes=wildcards.barcodes, minReadSupport=wildcards.minReadSupport),

	output: "mappings/nonAnchoredMergeReads/allSeqTechsCorr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.tmerge.min{minReadSupport}reads.endSupport:all.gff"
	threads:1
	wildcard_constraints:
#		techname='allSeqTechs',
#		sizeFrac='allFracs',
		barcodes='allTissues',
		capDesign='|'.join(CAPDESIGNS)
	shell:
		'''
uuidTmpOut=$(uuidgen)
uuid=$(uuidgen)
cat {input} |sort -T {config[TMPDIR]}  -k1,1 -k4,4n -k5,5n  > {config[TMPDIR]}/$uuid
cat {config[TMPDIR]}/$uuid | tmerge --exonOverhangTolerance {config[exonOverhangTolerance]} --tmPrefix allSeqTechsCorr{wildcards.corrLevel}_{wildcards.capDesign}_{wildcards.sizeFrac}_allTissues.NAM_ - |sort -T {config[TMPDIR]}  -k1,1 -k4,4n -k5,5n  > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}


		'''


rule mergeAllCapDesignsSeqTechsFracsAndTissuesNonAnchoredMergeReads:
	input:
		lambda wildcards: expand("mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.tmerge.min{minReadSupport}reads.splicing_status:all.endSupport:all.gff", filtered_capDesign_product_merge, techname=wildcards.techname, corrLevel=wildcards.corrLevel, capDesign=wildcards.capDesign, sizeFrac=wildcards.sizeFrac, barcodes=wildcards.barcodes, minReadSupport=wildcards.minReadSupport),

	output: "mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.tmerge.min{minReadSupport}reads.splicing_status:all.endSupport:all.gff"
	threads:1
	wildcard_constraints:
#		sizeFrac='0+',
		barcodes='allTissues',
		capDesign='|'.join(MERGEDCAPDESIGNS)
	shell:
		'''
uuidTmpOut=$(uuidgen)
uuid=$(uuidgen)
cat {input} |sort -T {config[TMPDIR]}  -k1,1 -k4,4n -k5,5n  > {config[TMPDIR]}/$uuid
cat {config[TMPDIR]}/$uuid | tmerge --exonOverhangTolerance {config[exonOverhangTolerance]} --tmPrefix allSeqTechsCorr{wildcards.corrLevel}_{wildcards.capDesign}_{wildcards.sizeFrac}_allTissues.NAM_ - |sort -T {config[TMPDIR]}  -k1,1 -k4,4n -k5,5n  > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}


		'''




rule nonAnchoredMergeReadsToBed:
	input: "mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.tmerge.min{minReadSupport}reads.splicing_status:all.endSupport:all.gff"
	output: temp("mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.splicing_status:all.endSupport:all.bed")
	shell:
		'''
uuidTmpOut=$(uuidgen)
cat {input} | gff2bed_full.pl - > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
		'''



rule getTmStats:
	input: "mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.tmerge.min{minReadSupport}reads.splicing_status:all.endSupport:{endSupport}.gff"
	output: temp(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.splicing_status:all.endSupport:{endSupport}.TmStats.stats.tsv")
	# wildcard_constraints:
	# 	barcodes='(?!allTissues).+',
	# 	sizeFrac='[0-9-+\.]+',
	# 	techname='(?!allSeqTechs).+'
	shell:
		'''
uuidTmpOut=$(uuidgen)
cat {input} | extractGffAttributeValue.pl transcript_id spliced mature_RNA_length contains_count 3p_dists_to_3p 5p_dists_to_5p meta_3p_dists_to_5p meta_5p_dists_to_5p |sort|uniq | awk -v t={wildcards.techname} -v c={wildcards.corrLevel} -v ca={wildcards.capDesign} -v s={wildcards.sizeFrac} -v b={wildcards.barcodes} '{{print t"Corr"c"\\t"ca"\\t"s"\\t"b"\t"$0}}'> {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
		'''

rule aggTmStats:
	input: lambda wildcards: expand(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.splicing_status:all.endSupport:{endSupport}.TmStats.stats.tsv", filtered_product, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES, minReadSupport=wildcards.minReadSupport, endSupport=wildcards.endSupport)
	output: config["STATSDATADIR"] + "all.min{minReadSupport}.endSupport:{endSupport}.TmStats.stats.tsv"
	shell:
		'''
uuidTmpOut=$(uuidgen)
echo -e "seqTech\tcorrectionLevel\tcapDesign\tsizeFrac\ttissue\ttranscript_id\tspliced\tcontains_count\tend\tdistance\tnormDistance" > {config[TMPDIR]}/$uuidTmpOut

cat {input} | perl -F"\\t" -slane '@ara=split(",", $F[8]); @arb=split(",", $F[9]); @arc=split(",", $F[10]), @ard=split(",", $F[11]); for ($i=0; $i<=$#ara; $i++){{$threepMinusDist=-$ara[$i];print  "$F[0]\\t$F[1]\\t$F[2]\\t$F[3]\\t$F[4]\\t$F[5]\\t$F[7]\\t3p\\t$threepMinusDist\\t$arc[$i]\\n$F[0]\\t$F[1]\\t$F[2]\\t$F[3]\\t$F[4]\\t$F[5]\\t$F[7]\\t5p\\t$arb[$i]\\t$ard[$i]"}}'| sed 's/Corr0/\\tNo/' | sed 's/Corr{lastK}/\\tYes/' | sort -T {config[TMPDIR]}  >> {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}

		'''



rule getMergingStats:
	input:
		hcgms = "mappings/highConfidenceReads/HiSS/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.gff.gz",
		pooledMerged = "mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.tmerge.min{minReadSupport}reads.splicing_status:all.endSupport:all.gff"
	output: temp(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.merged.stats.tsv")
	shell:
		'''
uuidTmpOut=$(uuidgen)
hcgms=$(zcat {input.hcgms} | extractGffAttributeValue.pl transcript_id | sort -T {config[TMPDIR]} |uniq|wc -l)
merged=$(cat {input.pooledMerged} | extractGffAttributeValue.pl transcript_id | sort -T {config[TMPDIR]} |uniq|wc -l)
echo -e "{wildcards.techname}Corr{wildcards.corrLevel}\t{wildcards.capDesign}\t{wildcards.sizeFrac}\t{wildcards.barcodes}\t$hcgms\t$merged" | awk '{{print $0"\t"$6/$5}}' > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}

		'''


rule aggMergingStats:
	input: lambda wildcards: expand(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.merged.stats.tsv",filtered_product_merge, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODESpluSMERGED, minReadSupport=wildcards.minReadSupport)

	output: config["STATSDATADIR"] + "all.min{minReadSupport}reads.merged.stats.tsv"
	shell:
		'''
uuidTmpOut=$(uuidgen)
echo -e "seqTech\tcorrectionLevel\tcapDesign\tsizeFrac\ttissue\tcategory\tcount" > {config[TMPDIR]}/$uuidTmpOut
cat {input} | awk '{{print $1"\\t"$2"\\t"$3"\\t"$4"\\tHCGMreads\\t"$5"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\tmergedTMs\\t"$6}}' | sed 's/Corr0/\tNo/' | sed 's/Corr{lastK}/\tYes/' | sort -T {config[TMPDIR]}  >> {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}

		'''
rule plotMergingStats:
	input:  config["STATSDATADIR"] + "all.min{minReadSupport}reads.merged.stats.tsv"
	output: returnPlotFilenames(config["PLOTSDIR"] + "merged.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.merged.stats")
	params:
		filterDat=lambda wildcards: merge_figures_params(wildcards.capDesign, wildcards.sizeFrac, wildcards.barcodes, wildcards.corrLevel, wildcards.techname)
	shell:
		'''
echo "
library(cowplot)
library(plyr)
library(scales)
library(gridExtra)
library(grid)
library(ggplotify)
dat <- read.table('{input}', header=T, as.is=T, sep='\\t')
{params.filterDat[10]}
{params.filterDat[0]}
{params.filterDat[1]}
{params.filterDat[2]}
{params.filterDat[3]}
{params.filterDat[4]}
{params.filterDat[5]}
{params.filterDat[8]}


plotBase <- \\"ggplot(data=dat, aes(x=factor(correctionLevel), y=count, fill=category)) +
geom_bar(stat='identity', position=position_dodge()) +
geom_text(position = position_dodge(width = 0.9), size=geom_textSize, aes(x = factor(correctionLevel), y = 1, label = comma(count), hjust = 0, vjust = 0.5), angle=90) +
scale_fill_manual(values=c('HCGMreads' = '#d98cb3', 'mergedTMs' = '#cc9966')) +
ylab('# objects') +
xlab('{params.filterDat[6]}') +
guides(fill = guide_legend(title='Category'))+
scale_y_continuous(labels=comma)+
{params.filterDat[7]}
{GGPLOT_PUB_QUALITY} + \\"

{params.filterDat[12]}

save_plot('{output[0]}', legendOnly, base_width=wLegendOnly, base_height=hLegendOnly)
save_plot('{output[1]}', legendOnly, base_width=wLegendOnly, base_height=hLegendOnly)

save_plot('{output[2]}', pXy, base_width=wXyPlot, base_height=hXyPlot)
save_plot('{output[3]}', pXy, base_width=wXyPlot, base_height=hXyPlot)

save_plot('{output[4]}', pXyNoLegend, base_width=wXyNoLegendPlot, base_height=hXyNoLegendPlot)
save_plot('{output[5]}', pXyNoLegend, base_width=wXyNoLegendPlot, base_height=hXyNoLegendPlot)

save_plot('{output[6]}', pYx, base_width=wYxPlot, base_height=hYxPlot)
save_plot('{output[7]}', pYx, base_width=wYxPlot, base_height=hYxPlot)

save_plot('{output[8]}', pYxNoLegend, base_width=wYxNoLegendPlot, base_height=hYxNoLegendPlot)
save_plot('{output[9]}', pYxNoLegend, base_width=wYxNoLegendPlot, base_height=hYxNoLegendPlot)


" > $(dirname {output[0]})/$(basename {output[0]} .legendOnly.png).r
cat $(dirname {output[0]})/$(basename {output[0]} .legendOnly.png).r | R --slave


		'''


rule aggTmLengthStats:
	input:
		all=lambda wildcards: expand(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.splicing_status:all.endSupport:all.TmStats.stats.tsv", filtered_product, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES, minReadSupport=wildcards.minReadSupport),
		fl=lambda wildcards: expand(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.splicing_status:all.endSupport:cagePolyASupported.TmStats.stats.tsv", filtered_product, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES, minReadSupport=wildcards.minReadSupport)
	output: config["STATSDATADIR"] + "all.min{minReadSupport}reads.matureRNALength.stats.tsv"
	shell:
		'''
uuidTmpOut=$(uuidgen)
echo -e "seqTech\tcorrectionLevel\tcapDesign\tsizeFrac\ttissue\ttranscript_id\tspliced\tmature_RNA_length\tcategory" > {config[TMPDIR]}/$uuidTmpOut

cat {input.all} |cut -f1-7| awk '{{print $0"\\tCLS_TMs"}}' | sed 's/Corr0/\\tNo/' | sed 's/Corr{lastK}/\\tYes/' |sort -T {config[TMPDIR]}  >> {config[TMPDIR]}/$uuidTmpOut
cat {input.fl} |cut -f1-7| awk '{{print $0"\\tCLS_FL_TMs"}}' | sed 's/Corr0/\\tNo/' | sed 's/Corr{lastK}/\\tYes/' |sort -T {config[TMPDIR]}  >> {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
		'''


rule plotTmLengthStats:
	input: config["STATSDATADIR"] + "all.min{minReadSupport}reads.matureRNALength.stats.tsv"
	output: returnPlotFilenames(config["PLOTSDIR"] + "matureRNALength.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.splicing_status:{splicedStatus}.matureRNALength.stats")
	params:
		filterDat=lambda wildcards: merge_figures_params(wildcards.capDesign, wildcards.sizeFrac, wildcards.barcodes, wildcards.corrLevel, wildcards.techname, wildcards.splicedStatus)
	shell:
		'''
echo "
library(cowplot)
library(plyr)
library(scales)
library(gridExtra)
library(grid)
library(ggplotify)
library(data.table)

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
plotBase <- \\"ggplot(dat, aes(x=factor(correctionLevel), y=mature_RNA_length, color=category)) +
geom_boxplot(position=position_dodge(0.9), outlier.shape=NA) +
coord_cartesian(ylim=c(100, 3000)) +
scale_y_continuous(labels=comma)+
scale_color_manual(values=palette, name='Category', labels = c('GENCODE_protein_coding' = 'GENCODE\nprotein-coding', 'CLS_TMs'='TMs', 'CLS_FL_TMs'='FL TMs')) +

stat_summary(aes(x=factor(correctionLevel), group=category), position=position_dodge(0.9), fun.data = fun_length, geom = 'text', vjust = +1, hjust=0, show.legend=FALSE, size=geom_textSize, colour='black') +
stat_summary(aes(x=factor(correctionLevel), group=category), position=position_dodge(0.9), fun.data = fun_median, geom = 'text', vjust = 0, hjust=0, show.legend=FALSE, size=geom_textSize, colour='black') +
ylab('Mature RNA length') +
xlab('{params.filterDat[6]}') +
coord_flip(ylim=c(100, 3000)) +
{params.filterDat[9]}
{GGPLOT_PUB_QUALITY} +
theme(axis.text.x = element_text(angle = 45, hjust = 1)) + \\"

{params.filterDat[12]}

save_plot('{output[0]}', legendOnly, base_width=wLegendOnly, base_height=hLegendOnly)
save_plot('{output[1]}', legendOnly, base_width=wLegendOnly, base_height=hLegendOnly)

save_plot('{output[2]}', pXy, base_width=wXyPlot, base_height=hXyPlot)
save_plot('{output[3]}', pXy, base_width=wXyPlot, base_height=hXyPlot)

save_plot('{output[4]}', pXyNoLegend, base_width=wXyNoLegendPlot, base_height=hXyNoLegendPlot)
save_plot('{output[5]}', pXyNoLegend, base_width=wXyNoLegendPlot, base_height=hXyNoLegendPlot)

save_plot('{output[6]}', pYx, base_width=wYxPlot, base_height=hYxPlot)
save_plot('{output[7]}', pYx, base_width=wYxPlot, base_height=hYxPlot)

save_plot('{output[8]}', pYxNoLegend, base_width=wYxNoLegendPlot, base_height=hYxNoLegendPlot)
save_plot('{output[9]}', pYxNoLegend, base_width=wYxNoLegendPlot, base_height=hYxNoLegendPlot)


" > $(dirname {output[0]})/$(basename {output[0]} .legendOnly.png).r
cat $(dirname {output[0]})/$(basename {output[0]} .legendOnly.png).r | R --slave

		'''
