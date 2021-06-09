rule integratePolyaAndSjInfo:
	input:
		polyA = "mappings/removePolyAERCCs/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.polyAsitesNoErcc.tmp.bed",
		SJs = "mappings/getIntronMotif/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.transcripts.tsv.gz"
	output:
		strandInfo="mappings/integratePolyaAndSjInfo/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.polyA+SJ.strandInfo.tsv",
		wrongPolyAs="mappings/removePolyAERCCs/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.wrongPolyAs.list"
	shell:
		'''
uuid=$(uuidgen)
uuidTmpOutS=$(uuidgen)
uuidTmpOutW=$(uuidgen)
zcat {input.SJs} | skipcomments | cut -f 1,2 | awk '$2!="."' | sort -T {config[TMPDIR]} > {config[TMPDIR]}/$uuid.reads.SJ.strandInfo.tsv
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
rule removeIntraPriming:
	input: 
		strandedGff="mappings/strandGffs/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.stranded.gff.gz",
		genome = lambda wildcards: config["GENOMESDIR"] + CAPDESIGNTOGENOME[wildcards.capDesign] + ".sorted.fa"
	output: 
		list="mappings/intraPriming/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.list",
		stats=config["STATSDATADIR"] + "tmp/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.intraPriming.stats.tsv"
	shell:
		'''
uuid=$(uuidgen)
echo $PATH

echo 

which perl
set +eu
# conda activate bedtools_env
conda env list

conda deactivate

conda env list
set -eu

echo $PATH

echo 

which perl

zcat {input.strandedGff} | awk '$3=="exon"'>  {config[TMPDIR]}/$uuid.gff
gtfToGenePred {config[TMPDIR]}/$uuid.gff {config[TMPDIR]}/$uuid.gp
genePredToBed {config[TMPDIR]}/$uuid.gp {config[TMPDIR]}/$uuid.bed
rm {config[TMPDIR]}/$uuid.gp
totalReads=$(cat {config[TMPDIR]}/$uuid.bed | wc -l)
set +eu
conda activate perl_env
set -eu

cat {config[TMPDIR]}/$uuid.bed | findIntraPriming --genomeFa {input.genome} --downSeqLength 20 - | gzip > $(dirname {output.list})/$(basename {output.list} .list).bed.gz
rm {config[TMPDIR]}/$uuid.bed
zcat $(dirname {output.list})/$(basename {output.list} .list).bed.gz | awk '$5>0.6'|cut -f4 | sort|uniq > {config[TMPDIR]}/$uuid
mv {config[TMPDIR]}/$uuid {output.list}
intraPrimed=$(cat {output.list} | wc -l)
let nonIntraPrimed=$totalReads-$intraPrimed || true
echo -e "{wildcards.techname}Corr{wildcards.corrLevel}\t{wildcards.capDesign}\t{wildcards.sizeFrac}\t{wildcards.barcodes}\t$totalReads\t$intraPrimed" | awk '{{print $0"\t"$6/$5}}' > {config[TMPDIR]}/$uuid.2
mv {config[TMPDIR]}/$uuid.2 {output.stats}
rm {config[TMPDIR]}/$uuid*
		'''

rule aggIntraPrimingStats:
	input: lambda wildcards: expand(config["STATSDATADIR"] + "tmp/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.intraPriming.stats.tsv",filtered_product_merge, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES)
	output: config["STATSDATADIR"] + "all.intraPriming.stats.tsv"
	shell:
		'''
uuid=$(uuidgen)
echo -e "seqTech\tcorrectionLevel\tcapDesign\tsizeFrac\ttissue\ttotalReads\tintraPrimed\tpercentIntraPrimed" > {config[TMPDIR]}/$uuid
cat {input} | sed 's/Corr0/\tNo/' | sed 's/Corr{lastK}/\tYes/' | sort -T {config[TMPDIR]}  >> {config[TMPDIR]}/$uuid
mv {config[TMPDIR]}/$uuid {output}
		'''

rule plotIntraPrimingStats:
	input: config["STATSDATADIR"] + "all.intraPriming.stats.tsv"
	output: returnPlotFilenames(config["PLOTSDIR"] + "intraPriming.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.intraPriming.stats")
	params:
		filterDat=lambda wildcards: merge_figures_params(wildcards.capDesign, wildcards.sizeFrac, wildcards.barcodes, wildcards.corrLevel, wildcards.techname)
	shell:
		'''
echo "
library(ggplot2)

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

plotBase <- \\"p <- ggplot(data=dat, aes(x=1, y=percentIntraPrimed, fill=sizeFrac)) +
geom_bar(width=0.75,stat='identity', position=position_dodge(width=0.9)) +
scale_fill_manual(values={sizeFrac_Rpalette}) +
#geom_hline(aes(yintercept=1), linetype='dashed', alpha=0.7, size=lineSize)+
geom_text(aes(group=sizeFrac, y=0, label = paste(sep='',percent(percentIntraPrimed),'\\n','(',comma(intraPrimed),')')), angle=90, size=geom_textSize, hjust=0, vjust=0.5, position = position_dodge(width=0.9)) +
scale_y_continuous(labels = scales::percent) +
xlab ('') +
ylab ('% intra-primed reads') +

{GGPLOT_PUB_QUALITY} + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) + \\"

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
 set +eu
conda activate R_env
set -eu
cat $(dirname {output[0]})/$(basename {output[0]} .legendOnly.png).r | R --slave

		'''

rule highConfidenceReads:
	input:
		transcriptStrandInfo = "mappings/getIntronMotif/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.transcripts.tsv.gz",
		strandedReads = "mappings/strandGffs/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.stranded.gff.gz",
		bam = "mappings/readMapping/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bam",
		intraPriming="mappings/intraPriming/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.list"
	params: 
		minSJsPhredScore = lambda wildcards: config['minPhredQualAroundSJwithHiSeq'] if sampleAnnotDict[wildcards.techname + "_" + wildcards.capDesign + "_" + wildcards.sizeFrac + "_" + wildcards.barcodes]['use_matched_HiSeq'] else config['minPhredQualAroundSJnoHiSeq']

	output:
		gff="mappings/highConfidenceReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.strandedHCGMs.gff.gz",
		stats=config["STATSDATADIR"] + "tmp/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.highConfSplicedReads.stats.tsv"
	shell:
		'''
uuidTmpOutG=$(uuidgen)
uuidTmpOutS=$(uuidgen)
#### reads QC stats
totalSplicedReads=$(zcat {input.transcriptStrandInfo} | tgrep -v -P "^#" | wc -l)
# reads with 100% canonical SJs
canonSjReads=$(zcat {input.transcriptStrandInfo} | awk '$6==1'| wc -l)
# reads with no fishy SJs (i.e. not surrounded by direct repeats)
noFishySjReads=$(zcat {input.transcriptStrandInfo} | awk '$7==1'| wc -l)
# reads with no fishy SJs and canonical SJs
noFishyCanonSjReads=$(zcat {input.transcriptStrandInfo} | awk '$6==1 && $7==1'| wc -l)
echo -e "{wildcards.techname}Corr{wildcards.corrLevel}\t{wildcards.capDesign}\t{wildcards.sizeFrac}\t{wildcards.barcodes}\t$totalSplicedReads\t$canonSjReads\t$noFishySjReads\t$noFishyCanonSjReads" | awk '{{print $0"\\t"$6/$5"\\t"$7/$5"\\t"$8/$5}}' > {config[TMPDIR]}/$uuidTmpOutS
mv {config[TMPDIR]}/$uuidTmpOutS {output.stats}


#select read IDs with canonical GT|GC/AG
uuid=$(uuidgen)
zcat {input.transcriptStrandInfo} | tgrep -v -P "^#" | awk '$6==1' | cut -f1 | sort -T {config[TMPDIR]} |uniq > {config[TMPDIR]}/$uuid.reads.hcSJs.list.tmp
wc -l {config[TMPDIR]}/$uuid.reads.hcSJs.list.tmp
zcat {input.strandedReads} > {config[TMPDIR]}/$uuid.str.gff

##high-quality Phred SJs
samtools view -F 256 -F4 -F 2048 {input.bam} | samHQintrons.pl --minQual {params.minSJsPhredScore} - |cut -f1|sort -T {config[TMPDIR]} |uniq > {config[TMPDIR]}/$uuid.HQintrons.reads.list
#tgrep -F -w -f {config[TMPDIR]}/$uuid.HQintrons.reads.list {config[TMPDIR]}/$uuid.reads.hcSJs.list.tmp > {config[TMPDIR]}/$uuid.reads.hcSJs.list
comm -1 -2 {config[TMPDIR]}/$uuid.HQintrons.reads.list {config[TMPDIR]}/$uuid.reads.hcSJs.list.tmp > {config[TMPDIR]}/$uuid.reads.hcSJs.list



tgrep -F -w -f {config[TMPDIR]}/$uuid.reads.hcSJs.list {config[TMPDIR]}/$uuid.str.gff > {config[TMPDIR]}/$uuid.gtag.gff
wc -l {config[TMPDIR]}/$uuid.gtag.gff
cat {config[TMPDIR]}/$uuid.str.gff | extractGffAttributeValue.pl transcript_id | sort -T {config[TMPDIR]} |uniq -u > {config[TMPDIR]}/$uuid.tmp
tgrep -F -w -f {config[TMPDIR]}/$uuid.tmp {config[TMPDIR]}/$uuid.str.gff > {config[TMPDIR]}/$uuid.tmp2
cat {config[TMPDIR]}/$uuid.tmp2  > {config[TMPDIR]}/$uuid.monoPolyA.gff
 echo $?
cat {config[TMPDIR]}/$uuid.gtag.gff {config[TMPDIR]}/$uuid.monoPolyA.gff | fgrep -vw -f {input.intraPriming} -| sort -T {config[TMPDIR]}  -k1,1 -k4,4n -k5,5n  |gzip> {config[TMPDIR]}/$uuidTmpOutG
mv {config[TMPDIR]}/$uuidTmpOutG {output.gff}
		'''

rule aggHighConfSplicedReadsStats:
	input: expand(config["STATSDATADIR"] + "tmp/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.highConfSplicedReads.stats.tsv",filtered_product_merge, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODESpluSMERGED)
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
		hiSeqIntrons=lambda wildcards: "mappings/hiSeqIntrons/hiSeq_{capDesign}.canonicalIntrons.list" if sampleAnnotDict[wildcards.techname + "_" + wildcards.capDesign + "_" + wildcards.sizeFrac + "_" + wildcards.barcodes]['use_matched_HiSeq'] else "/dev/null",
		lrIntrons="mappings/highConfidenceReads/introns/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.strandedHCGMs.introns.tsv",
		hcgmGTF= "mappings/highConfidenceReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.strandedHCGMs.gff.gz"
	params:
		useHiSeq = lambda wildcards: 'pleasedo' if sampleAnnotDict[wildcards.techname + "_" + wildcards.capDesign + "_" + wildcards.sizeFrac + "_" + wildcards.barcodes]['use_matched_HiSeq'] else 'donot'
	output:"mappings/highConfidenceReads/HiSS/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.gff.gz"
	shell:
		'''

echoerr "#########################"
echoerr "Use HiSeq SJs: {params.useHiSeq}"
echoerr "#########################"
if [ "{params.useHiSeq}" = "pleasedo" ]; then
uuid=$(uuidgen)
uuidTmpOut=$(uuidgen)
join -v1 -1 2 -2 1 {input.lrIntrons} {input.hiSeqIntrons} |awk '{{print $2"\t"$1}}' > {config[TMPDIR]}/$uuid.{wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.barcodes}.introns.noHiSeq.tsv
zcat {input.hcgmGTF} > {config[TMPDIR]}/$uuid.$(basename {input.hcgmGTF} .gz)
cut -f1 {config[TMPDIR]}/$uuid.{wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.barcodes}.introns.noHiSeq.tsv | sort -T {config[TMPDIR]} |uniq | fgrep -wv -f - {config[TMPDIR]}/$uuid.$(basename {input.hcgmGTF} .gz) |sort -T {config[TMPDIR]}  -k1,1 -k4,4n -k5,5n  |gzip > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
else
ln -srf {input.hcgmGTF} {output}
fi

		'''

rule getHiSSStats:
	input:
		reads = "mappings/readMapping/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bam",
		HiSSGTF="mappings/highConfidenceReads/HiSS/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.gff.gz"
	output: config["STATSDATADIR"] + "tmp/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.stats.tsv"
	shell:
		'''
uuid=$(uuidgen)
uuidTmpOut=$(uuidgen)
set +eu
conda activate bedtools_env
set -eu

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
	input: expand(config["STATSDATADIR"] + "tmp/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.stats.tsv",filtered_product_merge, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODESpluSMERGED)
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
library(ggplot2)
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

plotBase <- \\"p <- ggplot(dat[order(dat\$category), ], aes(x=factor(correctionLevel), y=count, fill=category)) +
geom_bar(stat='identity') +
scale_fill_manual(values=c('HCGM-mono' = '#9ce2bb', 'HCGM-spliced' = '#39c678', 'nonHCGM-mono' = '#fda59b', 'nonHCGM-spliced' = '#fa341e')) +
ylab('# mapped reads') +
xlab('') +
guides(fill = guide_legend(title='Category'))+
#geom_text(position = 'stack', size=geom_textSize, aes(x = factor(correctionLevel), y = count, label = comma(count), hjust = 0.5, vjust = 1))+
scale_y_continuous(labels=scientific)+
{GGPLOT_PUB_QUALITY} + 
theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
\\"

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
 set +eu
conda activate R_env
set -eu
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
	output: "mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.noFilt.tmerge.min{minReadSupport}reads.splicing_status:all.endSupport:all.gff.gz",
	threads:1
	wildcard_constraints:
		barcodes='(?!allTissues).+',
		sizeFrac='[0-9-+\.]+',
		techname='(?!allSeqTechs).+'
	shell:
		'''
uuidTmpOut=$(uuidgen)
uuid=$(uuidgen)
zcat {input} | awk '$1 ~ /SIRV/ || $1=="chrIS"' > {config[TMPDIR]}/$uuid

cat {config[TMPDIR]}/$uuid| tmerge --minReadSupport {wildcards.minReadSupport} --tmPrefix {wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.barcodes}.NAM_ - |sort -T {config[TMPDIR]}  -k1,1 -k4,4n -k5,5n | gzip > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
		'''



rule splitTmsBySplicedStatus:
	input: "mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.tmerge.min{minReadSupport}reads.endSupport:all.gff"
	output: "mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.tmerge.min{minReadSupport}reads.splicing_status:{splicedStatus}.endSupport:all.gff.gz"
	params:
		grepSpliced = lambda wildcards: '| fgrep \'spliced \"1\"\'' if wildcards.splicedStatus == "spliced" else '| fgrep \'spliced \"0\"\'' if wildcards.splicedStatus == "unspliced" else ''
	wildcard_constraints:
		techname='(?!allSeqTechs).+'
	shell:
		'''
uuidTmpOut=$(uuidgen)
cat {input} {params.grepSpliced} | gzip > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
		'''

rule mergeTissuesNonAnchoredMergeReads:
	input: lambda wildcards: expand("mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.tmerge.min{minReadSupport}reads.splicing_status:all.endSupport:all.gff.gz", filtered_product, techname=wildcards.techname, corrLevel=wildcards.corrLevel, capDesign=wildcards.capDesign, sizeFrac=wildcards.sizeFrac, barcodes=BARCODES, minReadSupport=wildcards.minReadSupport)
	output: "mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.tmerge.min{minReadSupport}reads.splicing_status:all.endSupport:all.gff.gz"
	threads:1
	wildcard_constraints:
		barcodes='allTissues',
		sizeFrac='[0-9-+\.]+',
		techname='(?!allSeqTechs).+'
	shell:
		'''
uuidTmpOut=$(uuidgen)
uuid=$(uuidgen)
zcat {input} |sort -T {config[TMPDIR]}  -k1,1 -k4,4n -k5,5n  > {config[TMPDIR]}/$uuid
cat {config[TMPDIR]}/$uuid | tmerge --exonOverhangTolerance {config[exonOverhangTolerance]} --tmPrefix {wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.barcodes}.NAM_ - |sort -T {config[TMPDIR]}  -k1,1 -k4,4n -k5,5n |gzip > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}


		'''


rule mergeAllSeqTechsFracsAndTissuesNonAnchoredMergeReads:
	input:
		lambda wildcards: expand("mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.tmerge.min{minReadSupport}reads.splicing_status:all.endSupport:all.gff.gz", filtered_product_merge, techname=TECHNAMES, corrLevel=wildcards.corrLevel, capDesign=wildcards.capDesign, sizeFrac=wildcards.sizeFrac, barcodes=wildcards.barcodes, minReadSupport=wildcards.minReadSupport),

	output: "mappings/nonAnchoredMergeReads/allSeqTechsCorr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.tmerge.min{minReadSupport}reads.endSupport:all.gff.gz"
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
zcat {input} |sort -T {config[TMPDIR]}  -k1,1 -k4,4n -k5,5n  > {config[TMPDIR]}/$uuid
cat {config[TMPDIR]}/$uuid | tmerge --exonOverhangTolerance {config[exonOverhangTolerance]} --tmPrefix allSeqTechsCorr{wildcards.corrLevel}_{wildcards.capDesign}_{wildcards.sizeFrac}_allTissues.NAM_ - |sort -T {config[TMPDIR]}  -k1,1 -k4,4n -k5,5n | gzip > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}


		'''


rule mergeAllCapDesignsSeqTechsFracsAndTissuesNonAnchoredMergeReads:
	input:
		lambda wildcards: expand("mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.tmerge.min{minReadSupport}reads.splicing_status:all.endSupport:all.gff.gz", filtered_capDesign_product_merge, techname=wildcards.techname, corrLevel=wildcards.corrLevel, capDesign=wildcards.capDesign, sizeFrac=wildcards.sizeFrac, barcodes=wildcards.barcodes, minReadSupport=wildcards.minReadSupport),

	output: "mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.tmerge.min{minReadSupport}reads.splicing_status:all.endSupport:all.gff.gz"
	threads:1
	wildcard_constraints:
		sizeFrac='0+',
		barcodes='allTissues',
		capDesign='|'.join(MERGEDCAPDESIGNS)
	shell:
		'''
uuidTmpOut=$(uuidgen)
uuid=$(uuidgen)
zcat {input} |sort -T {config[TMPDIR]}  -k1,1 -k4,4n -k5,5n  > {config[TMPDIR]}/$uuid
cat {config[TMPDIR]}/$uuid | tmerge --exonOverhangTolerance {config[exonOverhangTolerance]} --tmPrefix allSeqTechsCorr{wildcards.corrLevel}_{wildcards.capDesign}_{wildcards.sizeFrac}_allTissues.NAM_ - |sort -T {config[TMPDIR]}  -k1,1 -k4,4n -k5,5n |gzip > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}


		'''




rule nonAnchoredMergeReadsToBed:
	input: "mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.tmerge.min{minReadSupport}reads.splicing_status:all.endSupport:all.gff.gz"
	output: temp("mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.splicing_status:all.endSupport:all.bed")
	shell:
		'''
uuidTmpOut=$(uuidgen)
zcat {input} | gff2bed_full.pl - > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
		'''



rule getTmStats:
	input: "mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.tmerge.min{minReadSupport}reads.splicing_status:all.endSupport:{endSupport}.gff.gz"
	output: config["STATSDATADIR"] + "tmp/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.splicing_status:all.endSupport:{endSupport}.TmStats.stats.tsv.gz"
	# wildcard_constraints:
	# 	barcodes='(?!allTissues).+',
	# 	sizeFrac='[0-9-+\.]+',
	# 	techname='(?!allSeqTechs).+'
	shell:
		'''
uuidTmpOut=$(uuidgen)
zcat {input} | extractGffAttributeValue.pl transcript_id spliced mature_RNA_length contains_count 3p_dists_to_3p 5p_dists_to_5p meta_3p_dists_to_5p meta_5p_dists_to_5p |sort|uniq | awk -v t={wildcards.techname} -v c={wildcards.corrLevel} -v ca={wildcards.capDesign} -v s={wildcards.sizeFrac} -v b={wildcards.barcodes} '{{print t"Corr"c"\\t"ca"\\t"s"\\t"b"\t"$0}}' | gzip > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
		'''

rule aggTmStats:
	input: lambda wildcards: expand(config["STATSDATADIR"] + "tmp/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.splicing_status:all.endSupport:{endSupport}.TmStats.stats.tsv.gz", filtered_product, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=wildcards.capDesign, sizeFrac=SIZEFRACS, barcodes=BARCODES, minReadSupport=wildcards.minReadSupport, endSupport=wildcards.endSupport)
	output: config["STATSDATADIR"] + "all.{capDesign}.min{minReadSupport}.endSupport:{endSupport}.TmStats.stats.tsv.gz"
	shell:
		'''
uuidTmpOut=$(uuidgen)
echo -e "seqTech\tcorrectionLevel\tcapDesign\tsizeFrac\ttissue\tspliced\tcontains_count\tend\tdistance\tnormDistance" | gzip > {config[TMPDIR]}/$uuidTmpOut

zcat {input} | perl -F"\\t" -slane '@ara=split(",", $F[8]); @arb=split(",", $F[9]); @arc=split(",", $F[10]), @ard=split(",", $F[11]); for ($i=0; $i<=$#ara; $i++){{$threepMinusDist=-$ara[$i];print  "$F[0]\\t$F[1]\\t$F[2]\\t$F[3]\\t$F[5]\\t$F[7]\\t3\\t$threepMinusDist\\t$arc[$i]\\n$F[0]\\t$F[1]\\t$F[2]\\t$F[3]\\t$F[5]\\t$F[7]\\t5\\t$arb[$i]\\t$ard[$i]"}}'| sed 's/Corr0/\\tNo/' | sed 's/Corr{lastK}/\\tYes/' | sort -T {config[TMPDIR]}  | gzip >> {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}

		'''



rule getMergingStats:
	input:
		hcgms = "mappings/highConfidenceReads/HiSS/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.gff.gz",
		pooledMerged = "mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.tmerge.min{minReadSupport}reads.splicing_status:all.endSupport:all.gff.gz"
	output: config["STATSDATADIR"] + "tmp/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.merged.stats.tsv"
	shell:
		'''
uuidTmpOut=$(uuidgen)
hcgms=$(zcat {input.hcgms} | extractGffAttributeValue.pl transcript_id | sort -T {config[TMPDIR]} |uniq|wc -l)
merged=$(zcat {input.pooledMerged} | extractGffAttributeValue.pl transcript_id | sort -T {config[TMPDIR]} |uniq|wc -l)
echo -e "{wildcards.techname}Corr{wildcards.corrLevel}\t{wildcards.capDesign}\t{wildcards.sizeFrac}\t{wildcards.barcodes}\t$hcgms\t$merged" | awk '{{print $0"\t"$6/$5}}' > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}

		'''


rule aggMergingStats:
	input: lambda wildcards: expand(config["STATSDATADIR"] + "tmp/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.merged.stats.tsv",filtered_product_merge, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODESpluSMERGED, minReadSupport=wildcards.minReadSupport)

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
library(ggplot2)

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


plotBase <- \\"p <- ggplot(data=dat, aes(x=factor(correctionLevel), y=count, fill=category)) +
geom_bar(stat='identity', position=position_dodge()) +
geom_text(position = position_dodge(width = 0.9), size=geom_textSize, aes(x = factor(correctionLevel), y = 1, label = comma(count), hjust = 0, vjust = 0.5), angle=90) +
scale_fill_manual(values=c('HCGMreads' = '#d98cb3', 'mergedTMs' = '#cc9966')) +
ylab('# objects') +
xlab('') +
guides(fill = guide_legend(title='Category'))+
scale_y_continuous(labels=comma)+
{params.filterDat[7]}
{GGPLOT_PUB_QUALITY} + 
theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
\\"

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
 set +eu
conda activate R_env
set -eu
cat $(dirname {output[0]})/$(basename {output[0]} .legendOnly.png).r | R --slave


		'''


rule aggTmLengthStats:
	input:
		all=lambda wildcards: expand(config["STATSDATADIR"] + "tmp/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.splicing_status:all.endSupport:all.TmStats.stats.tsv.gz", filtered_product, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES, minReadSupport=wildcards.minReadSupport),
		fl=lambda wildcards: expand(config["STATSDATADIR"] + "tmp/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.splicing_status:all.endSupport:cagePolyASupported.TmStats.stats.tsv.gz", filtered_product, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES, minReadSupport=wildcards.minReadSupport)
	output: config["STATSDATADIR"] + "all.min{minReadSupport}reads.matureRNALength.stats.tsv.gz"
	shell:
		'''
uuidTmpOut=$(uuidgen)
echo -e "seqTech\tcorrectionLevel\tcapDesign\tsizeFrac\ttissue\ttranscript_id\tspliced\tmature_RNA_length\tcategory" |gzip> {config[TMPDIR]}/$uuidTmpOut

zcat {input.all} |cut -f1-7| awk '{{print $0"\\tCLS_TMs"}}' | sed 's/Corr0/\\tNo/' | sed 's/Corr{lastK}/\\tYes/' |sort -T {config[TMPDIR]}  | gzip >> {config[TMPDIR]}/$uuidTmpOut
zcat {input.fl} |cut -f1-7| awk '{{print $0"\\tCLS_FL_TMs"}}' | sed 's/Corr0/\\tNo/' | sed 's/Corr{lastK}/\\tYes/' |sort -T {config[TMPDIR]}  | gzip >> {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
		'''



rule plotHistTmLengthStats:
	input: config["STATSDATADIR"] + "all.min{minReadSupport}reads.matureRNALength.stats.tsv.gz"
	output: 
		hist=returnPlotFilenames(config["PLOTSDIR"] + "matureRNALength.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.splicing_status:{splicedStatus}.matureRNALength.hist.stats")
	params:
		filterDat=lambda wildcards: merge_figures_params(wildcards.capDesign, wildcards.sizeFrac, wildcards.barcodes, wildcards.corrLevel, wildcards.techname, wildcards.splicedStatus)
	shell:
		'''
echo "
library(ggplot2)
library(cowplot)
library(dplyr)
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

wXyPlot = wXyPlot * 1.2

dat %>%
  group_by(seqTech, correctionLevel, sizeFrac, capDesign, tissue, category) %>%
  summarise(med=median(mature_RNA_length)) -> datSumm

summaryStats = transform(datSumm, LabelM = comma(round(med, digits = 0))) 

geom_textSize = geom_textSize + 1


plotBase <- \\"p <- ggplot(dat, aes(x=mature_RNA_length, fill=category)) +
geom_histogram(binwidth=100, alpha=0.35, position='identity') +
geom_vline(data=summaryStats, aes(xintercept=med, color=category), size=lineSize, show.legend=FALSE) +
geom_text(data = summaryStats[summaryStats\$category=='CLS_TMs',], aes(label = LabelM, x = Inf, y = Inf, color=category), hjust=1.2, vjust=1,  size=geom_textSize, show.legend=FALSE) +
geom_text(data = summaryStats[summaryStats\$category=='CLS_FL_TMs',], aes(label = LabelM, x = Inf, y = Inf, color=category), hjust=1.2, vjust=2.1,  size=geom_textSize, show.legend=FALSE) +

coord_cartesian(xlim=c(0, 3000)) +
scale_x_continuous(labels=comma)+
scale_color_manual(values=palette, name='Category', labels = c('GENCODE_protein_coding' = 'GENCODE\nprotein-coding', 'CLS_TMs'='TMs', 'CLS_FL_TMs'='FL TMs')) +
scale_fill_manual(values=palette, name='Category', labels = c('GENCODE_protein_coding' = 'GENCODE\nprotein-coding', 'CLS_TMs'='TMs', 'CLS_FL_TMs'='FL TMs')) +

xlab('Mature RNA length') +
{GGPLOT_PUB_QUALITY} +
theme(axis.text.x = element_text(angle = 45, hjust = 1)) + \\"

{params.filterDat[12]}

wYxPlot = wYxPlot * 1.2
wYxNoLegendPlot<- wYxPlot - wLegendOnly

save_plot('{output.hist[0]}', legendOnly, base_width=wLegendOnly, base_height=hLegendOnly)
save_plot('{output.hist[1]}', legendOnly, base_width=wLegendOnly, base_height=hLegendOnly)

save_plot('{output.hist[2]}', pXy, base_width=wXyPlot, base_height=hXyPlot)
save_plot('{output.hist[3]}', pXy, base_width=wXyPlot, base_height=hXyPlot)

save_plot('{output.hist[4]}', pXyNoLegend, base_width=wXyNoLegendPlot, base_height=hXyNoLegendPlot)
save_plot('{output.hist[5]}', pXyNoLegend, base_width=wXyNoLegendPlot, base_height=hXyNoLegendPlot)

save_plot('{output.hist[6]}', pYx, base_width=wYxPlot, base_height=hYxPlot)
save_plot('{output.hist[7]}', pYx, base_width=wYxPlot, base_height=hYxPlot)

save_plot('{output.hist[8]}', pYxNoLegend, base_width=wYxNoLegendPlot, base_height=hYxNoLegendPlot)
save_plot('{output.hist[9]}', pYxNoLegend, base_width=wYxNoLegendPlot, base_height=hYxNoLegendPlot)

" > $(dirname {output.hist[0]})/$(basename {output.hist[0]} .legendOnly.png).r
 set +eu
conda activate R_env
set -eu
cat $(dirname {output.hist[0]})/$(basename {output.hist[0]} .legendOnly.png).r | R --slave

		'''


rule getGeneReadCoverageStats:
	input: 
		gencode="annotations/simplified/{capDesign}.gencode.simplified_biotypes.gtf",
		bam = "mappings/readMapping/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bam",
		tmerge = "mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.tmerge.min1reads.splicing_status:all.endSupport:all.gff.gz",
		genome=lambda wildcards: config["GENOMESDIR"] + CAPDESIGNTOGENOME[wildcards.capDesign] + ".sorted.genome"
	output: 
		gencode="mappings/geneReadCoverage/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.gencode.coverage.tsv",
		tmerge="mappings/geneReadCoverage/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.coverage.tsv"
	wildcard_constraints:
		barcodes='(?!allTissues).+',
		sizeFrac='[0-9-+\.]+',
		techname='(?!allSeqTechs).+'
	shell:
		'''
uuid=$(uuidgen)
cut -f1 {input.genome} |sort|uniq > {config[TMPDIR]}/$uuid.chr

#gencode
set +eu
conda activate bedtools_env
set -eu

cat {input.gencode} |awk '$3=="exon"' | extract_locus_coords.pl -| fgrep -w -f {config[TMPDIR]}/$uuid.chr |sort -T {config[TMPDIR]}  -k1,1 -k2,2n -k3,3n  > {config[TMPDIR]}/$uuid.1

bedtools coverage -sorted -g {input.genome} -bed -split -nonamecheck -counts -a {config[TMPDIR]}/$uuid.1 -b {input.bam} |sort -k7,7nr | cut -f1-4,7 | awk -v s={wildcards.techname}Corr{wildcards.corrLevel} -v c={wildcards.capDesign} -v si={wildcards.sizeFrac} -v b={wildcards.barcodes} '{{print s"\\t"c"\\t"si"\\t"b"\\t"$0}}'> {config[TMPDIR]}/$uuid.2

mv {config[TMPDIR]}/$uuid.2 {output.gencode}

#tmerge
zcat {input.tmerge} > {config[TMPDIR]}/$uuid.a
bedtools intersect -s -wao -a {config[TMPDIR]}/$uuid.a -b {config[TMPDIR]}/$uuid.a | buildLoci.pl - | extract_locus_coords.pl -| sort -T {config[TMPDIR]}  -k1,1 -k2,2n -k3,3n  > {config[TMPDIR]}/$uuid.3

bedtools coverage -sorted -g {input.genome}  -bed -split -nonamecheck -counts -a {config[TMPDIR]}/$uuid.3 -b {input.bam}  |sort -k7,7nr | cut -f1-4,7 | awk -v s={wildcards.techname}Corr{wildcards.corrLevel} -v c={wildcards.capDesign} -v si={wildcards.sizeFrac} -v b={wildcards.barcodes} '{{print s"\\t"c"\\t"si"\\t"b"\\t"$0}}' > {config[TMPDIR]}/$uuid.4
mv {config[TMPDIR]}/$uuid.4 {output.tmerge}

		'''

rule aggGeneReadCoverageStats:
	input: lambda wildcards: expand("mappings/geneReadCoverage/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.coverage.tsv", filtered_product_merge, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODESpluSMERGED)
	output: 
		config["STATSDATADIR"] + "all.tmerge.GeneReadCoverage.stats.tsv"
	shell:
		'''
uuid=$(uuidgen)
echo -e "seqTech\\tcorrectionLevel\\tcapDesign\\tsizeFrac\\ttissue\\tgene_id\\treadCount" > {config[TMPDIR]}/$uuid
cat {input} | cut -f1-4,8,9| sed 's/Corr0/\tNo/' | sed 's/Corr{lastK}/\tYes/'   >> {config[TMPDIR]}/$uuid
mv {config[TMPDIR]}/$uuid {output}

		'''

rule plotGeneReadCoverageStats:
	input: config["STATSDATADIR"] + "all.tmerge.GeneReadCoverage.stats.tsv"
	output: returnPlotFilenames(config["PLOTSDIR"] + "geneReadCoverage.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.geneReadCoverage.stats")
	params:
		filterDat=lambda wildcards: merge_figures_params(wildcards.capDesign, wildcards.sizeFrac, wildcards.barcodes, wildcards.corrLevel, wildcards.techname)
	shell:
		'''
echo "
library(ggplot2)

library(cowplot)
library(scales)
library(gridExtra)
library(grid)
library(ggplotify)
library(dplyr)
library(data.table)

dat<-fread('{input}', header=T, sep='\\t')

{params.filterDat[10]}
{params.filterDat[0]}
{params.filterDat[1]}
{params.filterDat[2]}
{params.filterDat[3]}
{params.filterDat[4]}
{params.filterDat[5]}
{params.filterDat[8]}

dat <- arrange(dat, desc(readCount))
group_by(dat, seqTech, correctionLevel, capDesign, sizeFrac, tissue) %>% mutate(rank=row_number()) -> dat
mutate(dat, rank=row_number()) -> dat
dat <- mutate(dat, cumSum=cumsum(readCount))
dat <- mutate(dat, cumProp=cumSum/sum(readCount))
filter(dat, rank==10)  -> cumPropTop10Genes


plotBase <- \\"p <- ggplot(dat, aes(x = rank, y = cumProp)) + geom_line(size=lineSize) + scale_x_log10() + ylim(0,1)+ geom_segment(data = cumPropTop10Genes, aes(x=10,xend=10,y=0,yend=cumProp), color='firebrick3',size=lineSize) + geom_segment(data = cumPropTop10Genes, aes(x=0,xend=10,y=cumProp,yend=cumProp,color='Contribution\\nof top 10 genes'),size=lineSize) + labs(y='Proportion of total mapped reads', x='# genes (ranked by expression)', color='') + scale_color_manual(values = c('Contribution\nof top 10 genes' = 'firebrick3')) + facet_grid( seqTech ~ capDesign + tissue) + geom_rect(data = cumPropTop10Genes, aes(xmin=0,xmax=10,ymin=0,ymax=cumProp),fill='firebrick3', alpha=0.3, size=0) +
{GGPLOT_PUB_QUALITY} + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + \\"



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
 set +eu
conda activate R_env
set -eu
cat $(dirname {output[0]})/$(basename {output[0]} .legendOnly.png).r | R --slave

		'''
