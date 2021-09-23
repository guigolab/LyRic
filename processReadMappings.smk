rule integratePolyaAndSjInfo:
	input:
		polyA = "output/mappings/removePolyAERCCs/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.polyAsitesNoErcc.tmp.bed",
		SJs = "output/mappings/getIntronMotif/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.transcripts.tsv.gz"
	output:
		strandInfo="output/mappings/integratePolyaAndSjInfo/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.polyA+SJ.strandInfo.tsv",
		wrongPolyAs="output/mappings/removePolyAERCCs/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.wrongPolyAs.list"
	shell:
		'''
uuid=$(uuidgen)
uuidTmpOutS=$(uuidgen)
uuidTmpOutW=$(uuidgen)
zcat {input.SJs} | skipcomments | cut -f 1,2 | awk '$2!="."' | sort -T {TMPDIR} > {TMPDIR}/$uuid.reads.SJ.strandInfo.tsv
cat {input.polyA} | cut -f4,6 | awk '$2!="."'| sort -T {TMPDIR}  > {TMPDIR}/$uuid.reads.polyA.strandInfo.tsv
join -a1 -a2 -e '.' -o '0,1.2,2.2' {TMPDIR}/$uuid.reads.SJ.strandInfo.tsv  {TMPDIR}/$uuid.reads.polyA.strandInfo.tsv > {TMPDIR}/$uuid.reads.SJ.polyA.strandInfo.tsv

#make list of reads with wrongly called polyA sites (i.e. their strand is different from the one inferred using SJs):
cat {TMPDIR}/$uuid.reads.SJ.polyA.strandInfo.tsv | perl -slane 'if($F[2] ne "." && $F[1] ne "." && $F[1] ne $F[2]){{print join ("\\t", $F[0])}}' |sort -T {TMPDIR} |uniq > {TMPDIR}/$uuidTmpOutW
mv {TMPDIR}/$uuidTmpOutW {output.wrongPolyAs}
#get strand info, prioritizing SJ inference
cat {TMPDIR}/$uuid.reads.SJ.polyA.strandInfo.tsv | perl -slane 'if($F[1] ne "."){{print "$F[0]\\t$F[1]"}} else{{print "$F[0]\\t$F[2]"}}' |sort -T {TMPDIR} |uniq > {TMPDIR}/$uuidTmpOutS
mv {TMPDIR}/$uuidTmpOutS {output.strandInfo}
		'''

rule removeWrongPolyAs:
	input:
		polyA="output/mappings/removePolyAERCCs/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.polyAsitesNoErcc.tmp.bed",
		wrongPolyAs="output/mappings/removePolyAERCCs/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.wrongPolyAs.list"
	output: "output/mappings/removePolyAERCCs/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.polyAsitesNoErcc.bed"
	shell:
		'''
uuidTmpOut=$(uuidgen)
cat {input.polyA} | fgrep -v -w -f {input.wrongPolyAs} > {TMPDIR}/$uuidTmpOut
mv {TMPDIR}/$uuidTmpOut {output}
		'''

rule strandGffs:
	input:
		gff = "output/mappings/readBedToGff/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.gff.gz",
		strandInfo = "output/mappings/integratePolyaAndSjInfo/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.polyA+SJ.strandInfo.tsv"
	output: "output/mappings/strandGffs/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.stranded.gff.gz"
	shell:
		'''
uuid=$(uuidgen)
uuidTmpOut=$(uuidgen)
zcat {input.gff} > {TMPDIR}/$uuid.in.gff
get_right_transcript_strand.pl {TMPDIR}/$uuid.in.gff {input.strandInfo} | fgrep -v ERCC- | sort -T {TMPDIR}  -k1,1 -k4,4n -k5,5n  |gzip> {TMPDIR}/$uuidTmpOut
mv {TMPDIR}/$uuidTmpOut {output}

		'''

rule strandedGffToBed:
	input: strandedGff="output/mappings/strandGffs/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.stranded.gff.gz"
	output: temp("output/mappings/strandGffs/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.stranded.bed")
	conda: "envs/ucsc_env.yml"
	shell:
		'''
uuid=$(uuidgen)
zcat {input.strandedGff} | awk '$3=="exon"'>  {TMPDIR}/$uuid.gff
gtfToGenePred {TMPDIR}/$uuid.gff {TMPDIR}/$uuid.gp
genePredToBed {TMPDIR}/$uuid.gp {TMPDIR}/$uuid.bed
mv {TMPDIR}/$uuid.bed {output}

		'''


rule removeIntraPriming:
	input: 
		strandedBed="output/mappings/strandGffs/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.stranded.bed",
		genome = lambda wildcards: config["GENOMESDIR"] + CAPDESIGNTOGENOME[wildcards.capDesign] + ".sorted.fa"
	output: 
		list="output/mappings/intraPriming/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.list",
		stats="output/statsFiles/" + "tmp/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.intraPriming.stats.tsv"
	conda: "envs/perl_env.yml"
	shell:
		'''
uuid=$(uuidgen)

totalReads=$(cat {input.strandedBed} | wc -l)


cat {input.strandedBed} | findIntraPriming --genomeFa {input.genome} --downSeqLength 20 - > {TMPDIR}/$uuid.bed
cat {TMPDIR}/$uuid.bed | awk '$5>0.6'|cut -f4 | sort|uniq > {TMPDIR}/$uuid
mv {TMPDIR}/$uuid {output.list}
intraPrimed=$(cat {output.list} | wc -l)
let nonIntraPrimed=$totalReads-$intraPrimed || true
echo -e "{wildcards.techname}\t{wildcards.capDesign}\t{wildcards.sizeFrac}\t{wildcards.sampleRep}\t$totalReads\t$intraPrimed" | awk '{{print $0"\t"$6/$5}}' > {TMPDIR}/$uuid.2
mv {TMPDIR}/$uuid.2 {output.stats}
		'''

rule aggIntraPrimingStats:
	input: lambda wildcards: expand("output/statsFiles/" + "tmp/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.intraPriming.stats.tsv",filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, sampleRep=SAMPLEREPS)
	output: "output/statsFiles/" + "all.intraPriming.stats.tsv"
	shell:
		'''
uuid=$(uuidgen)
echo -e "seqTech\tcapDesign\tsizeFrac\tsampleRep\ttotalReads\tintraPrimed\tpercentIntraPrimed" > {TMPDIR}/$uuid
cat {input} | sort -T {TMPDIR}  >> {TMPDIR}/$uuid
mv {TMPDIR}/$uuid {output}
		'''

rule plotIntraPrimingStats:
	input: "output/statsFiles/" + "all.intraPriming.stats.tsv"
	output: returnPlotFilenames("output/plots/" + "intraPriming.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.intraPriming.stats")
	params:
		filterDat=lambda wildcards: multi_figures(wildcards.capDesign, wildcards.sizeFrac, wildcards.sampleRep, wildcards.techname)
	conda: "envs/R_env.yml"
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
{params.filterDat[technameFilterString]}
{params.filterDat[capDesignFilterString]}

{params.filterDat[sizeFracFilterString]}
{params.filterDat[sampleRepFilterString]}
{params.filterDat[substSeqTechString]}
{params.filterDat[substSampleRepString]}
{params.filterDat[graphDimensions]}

plotBase <- \\"p <- ggplot(data=dat, aes(x=1, y=percentIntraPrimed, fill=sizeFrac)) +
geom_bar(width=0.75,stat='identity', position=position_dodge(width=0.9)) +
scale_fill_manual(values={sizeFrac_Rpalette}) +
#geom_hline(aes(yintercept=1), linetype='dashed', alpha=0.7, size=lineSize)+
geom_text(aes(group=sizeFrac, y=0, label = paste(sep='',percent(percentIntraPrimed),'\\n','(',comma(intraPrimed),')')), angle=90, size=geom_textSize, hjust=0, vjust=0.5, position = position_dodge(width=0.9)) +
scale_y_continuous(labels = scales::percent) +
xlab ('') +
ylab ('% intra-primed reads') +

{GGPLOT_PUB_QUALITY} + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) + \\"

{params.filterDat[facetPlotSetup]}


save_plot('{output[0]}', legendOnly, base_width=wLegendOnly, base_height=hLegendOnly)
save_plot('{output[1]}', pXy, base_width=wXyPlot, base_height=hXyPlot)

save_plot('{output[2]}', pXyNoLegend, base_width=wXyNoLegendPlot, base_height=hXyNoLegendPlot)
save_plot('{output[3]}', pYx, base_width=wYxPlot, base_height=hYxPlot)

save_plot('{output[4]}', pYxNoLegend, base_width=wYxNoLegendPlot, base_height=hYxNoLegendPlot)


" > $(dirname {output[0]})/$(basename {output[0]} .legendOnly.png).r

cat $(dirname {output[0]})/$(basename {output[0]} .legendOnly.png).r | R --slave

		'''

rule highConfidenceReads:
	input:
		transcriptStrandInfo = "output/mappings/getIntronMotif/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.transcripts.tsv.gz",
		strandedReads = "output/mappings/strandGffs/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.stranded.gff.gz",
		bam = "output/mappings/longReadMapping/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.bam",
		intraPriming="output/mappings/intraPriming/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.list"
	params: 
		minSJsPhredScore = lambda wildcards: sampleAnnotDict[wildcards.techname + "_" + wildcards.capDesign + "_" + wildcards.sizeFrac + "_" + wildcards.sampleRep]['filter_SJ_Qscore']
	conda: "envs/xtools_env.yml"
	output:
		gff="output/mappings/highConfidenceReads/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.strandedHCGMs.gff.gz",
		stats="output/statsFiles/" + "tmp/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.highConfSplicedReads.stats.tsv"
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
echo -e "{wildcards.techname}\t{wildcards.capDesign}\t{wildcards.sizeFrac}\t{wildcards.sampleRep}\t$totalSplicedReads\t$canonSjReads\t$noFishySjReads\t$noFishyCanonSjReads" | awk '{{print $0"\\t"$6/$5"\\t"$7/$5"\\t"$8/$5}}' > {TMPDIR}/$uuidTmpOutS
mv {TMPDIR}/$uuidTmpOutS {output.stats}


#select read IDs with canonical GT|GC/AG
uuid=$(uuidgen)
zcat {input.transcriptStrandInfo} | tgrep -v -P "^#" | awk '$6==1' | cut -f1 | sort -T {TMPDIR} |uniq > {TMPDIR}/$uuid.reads.hcSJs.list.tmp
wc -l {TMPDIR}/$uuid.reads.hcSJs.list.tmp
zcat {input.strandedReads} > {TMPDIR}/$uuid.str.gff

##high-quality Phred SJs
samtools view -F 256 -F4 -F 2048 {input.bam} | samHQintrons.pl --minQual {params.minSJsPhredScore} - |cut -f1|sort -T {TMPDIR} |uniq > {TMPDIR}/$uuid.HQintrons.reads.list
#tgrep -F -w -f {TMPDIR}/$uuid.HQintrons.reads.list {TMPDIR}/$uuid.reads.hcSJs.list.tmp > {TMPDIR}/$uuid.reads.hcSJs.list
comm -1 -2 {TMPDIR}/$uuid.HQintrons.reads.list {TMPDIR}/$uuid.reads.hcSJs.list.tmp > {TMPDIR}/$uuid.reads.hcSJs.list



tgrep -F -w -f {TMPDIR}/$uuid.reads.hcSJs.list {TMPDIR}/$uuid.str.gff > {TMPDIR}/$uuid.gtag.gff
wc -l {TMPDIR}/$uuid.gtag.gff
cat {TMPDIR}/$uuid.str.gff | extractGffAttributeValue.pl transcript_id | sort -T {TMPDIR} |uniq -u > {TMPDIR}/$uuid.tmp
tgrep -F -w -f {TMPDIR}/$uuid.tmp {TMPDIR}/$uuid.str.gff > {TMPDIR}/$uuid.tmp2
cat {TMPDIR}/$uuid.tmp2  > {TMPDIR}/$uuid.monoPolyA.gff
 echo $?
cat {TMPDIR}/$uuid.gtag.gff {TMPDIR}/$uuid.monoPolyA.gff | fgrep -vw -f {input.intraPriming} -| sort -T {TMPDIR}  -k1,1 -k4,4n -k5,5n  |gzip> {TMPDIR}/$uuidTmpOutG
mv {TMPDIR}/$uuidTmpOutG {output.gff}
		'''

rule aggHighConfSplicedReadsStats:
	input: expand("output/statsFiles/" + "tmp/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.highConfSplicedReads.stats.tsv",filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, sampleRep=SAMPLEREPS)
	output: "output/statsFiles/" + "all.highConfSplicedReads.stats.tsv"
	shell:
		'''
uuidTmpOut=$(uuidgen)
echo -e "seqTech\tcapDesign\tsizeFrac\tsampleRep\ttotalSplicedReads\tcanonSjReads\tnoFishySjReads\tnoFishyCanonSjReads" > {TMPDIR}/$uuidTmpOut
cat {input} | sort -T {TMPDIR}  >> {TMPDIR}/$uuidTmpOut
mv {TMPDIR}/$uuidTmpOut {output}
		'''


rule getHCGMintrons:
	input: "output/mappings/highConfidenceReads/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.strandedHCGMs.gff.gz"
	output: temp("output/mappings/highConfidenceReads/introns/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.strandedHCGMs.introns.tsv")
	shell:
		'''
uuidTmpOut=$(uuidgen)
zcat {input} | makeIntrons.pl -| perl -lane '$F[11]=~s/"|;//g; $start=$F[3]-1; $end=$F[4]+1; print $F[11]."\t".$F[0]."_".$start."_".$end."_".$F[6]' | sort -T {TMPDIR}  -k2,2 > {TMPDIR}/$uuidTmpOut
mv {TMPDIR}/$uuidTmpOut {output}

		'''

rule getHiSeqSupportedHCGMs:
	input:
		hiSeqIntrons=lambda wildcards: "output/mappings/hiSeqIntrons/hiSeq_{capDesign}.canonicalIntrons.list" if sampleAnnotDict[wildcards.techname + "_" + wildcards.capDesign + "_" + wildcards.sizeFrac + "_" + wildcards.sampleRep]['use_matched_HiSeq'] else "/dev/null",
		lrIntrons="output/mappings/highConfidenceReads/introns/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.strandedHCGMs.introns.tsv",
# HCGM: high-confidence genome mappings
		hcgmGTF= "output/mappings/highConfidenceReads/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.strandedHCGMs.gff.gz"
	params:
		useHiSeq = lambda wildcards: 'pleasedo' if sampleAnnotDict[wildcards.techname + "_" + wildcards.capDesign + "_" + wildcards.sizeFrac + "_" + wildcards.sampleRep]['use_matched_HiSeq'] else 'donot'
	output:"output/mappings/highConfidenceReads/HiSS/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.HiSS.gff.gz"
	shell:
		'''

echoerr "#########################"
echoerr "Use HiSeq SJs: {params.useHiSeq}"
echoerr "#########################"
if [ "{params.useHiSeq}" = "pleasedo" ]; then
uuid=$(uuidgen)
uuidTmpOut=$(uuidgen)
join -v1 -1 2 -2 1 {input.lrIntrons} {input.hiSeqIntrons} |awk '{{print $2"\t"$1}}' > {TMPDIR}/$uuid.{wildcards.techname}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.sampleRep}.introns.noHiSeq.tsv
zcat {input.hcgmGTF} > {TMPDIR}/$uuid.$(basename {input.hcgmGTF} .gz)
cut -f1 {TMPDIR}/$uuid.{wildcards.techname}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.sampleRep}.introns.noHiSeq.tsv | sort -T {TMPDIR} |uniq | fgrep -wv -f - {TMPDIR}/$uuid.$(basename {input.hcgmGTF} .gz) |sort -T {TMPDIR}  -k1,1 -k4,4n -k5,5n  |gzip > {TMPDIR}/$uuidTmpOut
mv {TMPDIR}/$uuidTmpOut {output}
else
ln -srf {input.hcgmGTF} {output}
fi

		'''

rule getHiSSStats:
	input:
		reads = "output/mappings/longReadMapping/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.bam",
		HiSSGTF="output/mappings/highConfidenceReads/HiSS/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.HiSS.gff.gz"
	output: "output/statsFiles/" + "tmp/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.HiSS.stats.tsv"
	conda: "envs/xtools_env.yml"
	shell:
		'''
uuid=$(uuidgen)
uuidTmpOut=$(uuidgen)


bedtools bamtobed -i {input.reads} -bed12 > {TMPDIR}/$uuid.{wildcards.techname}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.sampleRep}.merged.bed
zcat {input.HiSSGTF} | gff2bed_full.pl - > {TMPDIR}/$uuid.{wildcards.techname}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.sampleRep}.HiSS.bed

mappedReadsMono=$(cat {TMPDIR}/$uuid.{wildcards.techname}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.sampleRep}.merged.bed | awk '$10<=1'|cut -f4 |sort -T {TMPDIR} |uniq|wc -l)
mappedReadsSpliced=$(cat {TMPDIR}/$uuid.{wildcards.techname}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.sampleRep}.merged.bed | awk '$10>1'|cut -f4 |sort -T {TMPDIR} |uniq|wc -l)

HiSSMono=$(cat {TMPDIR}/$uuid.{wildcards.techname}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.sampleRep}.HiSS.bed| awk '$10<=1'|cut -f4  | sort -T {TMPDIR} |uniq|wc -l)
HiSSSpliced=$(cat {TMPDIR}/$uuid.{wildcards.techname}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.sampleRep}.HiSS.bed| awk '$10>1'|cut -f4  | sort -T {TMPDIR} |uniq|wc -l)

#let totalMapped=$mappedReadsMono+$mappedReadsSpliced || true
let nonHiSSMono=$mappedReadsMono-$HiSSMono || true
let nonHiSSSPliced=$mappedReadsSpliced-$HiSSSpliced || true
echo -e "{wildcards.techname}\t{wildcards.capDesign}\t{wildcards.sizeFrac}\t{wildcards.sampleRep}\t$HiSSMono\t$HiSSSpliced\t$nonHiSSMono\t$nonHiSSSPliced" > {TMPDIR}/$uuidTmpOut
mv {TMPDIR}/$uuidTmpOut {output}

		'''


rule aggHiSSStats:
	input: expand("output/statsFiles/" + "tmp/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.HiSS.stats.tsv",filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, sampleRep=SAMPLEREPS)
	output: "output/statsFiles/" + "all.HiSS.stats.tsv"

	shell:
		'''
uuidTmpOut=$(uuidgen)
echo -e "seqTech\tcapDesign\tsizeFrac\tsampleRep\tcategory\tcount" > {TMPDIR}/$uuidTmpOut
cat {input} | awk '{{print $1"\\t"$2"\\t"$3"\\t"$4"\\tHCGM-mono\\t"$5"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\tHCGM-spliced\\t"$6"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\tnonHCGM-mono\\t"$7"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\tnonHCGM-spliced\\t"$8}}' | sort -T {TMPDIR}  >> {TMPDIR}/$uuidTmpOut
mv {TMPDIR}/$uuidTmpOut {output}

		'''

rule plotAllHiSSStats:
	input: "output/statsFiles/" + "all.HiSS.stats.tsv"
	output:  returnPlotFilenames("output/plots/" + "HiSS.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.HiSS.stats")
	params:
		filterDat=lambda wildcards: multi_figures(wildcards.capDesign, wildcards.sizeFrac, wildcards.sampleRep, wildcards.techname)
	conda: "envs/R_env.yml"
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
{params.filterDat[technameFilterString]}
{params.filterDat[capDesignFilterString]}

{params.filterDat[sizeFracFilterString]}
{params.filterDat[sampleRepFilterString]}
{params.filterDat[substSeqTechString]}
{params.filterDat[substSampleRepString]}
{params.filterDat[graphDimensions]}

plotBase <- \\"p <- ggplot(dat[order(dat\$category), ], aes(x=1, y=count, fill=category)) +
geom_bar(stat='identity') +
scale_fill_manual(values=c('HCGM-mono' = '#9ce2bb', 'HCGM-spliced' = '#39c678', 'nonHCGM-mono' = '#fda59b', 'nonHCGM-spliced' = '#fa341e')) +
ylab('# mapped reads') +
xlab('') +
guides(fill = guide_legend(title='Category'))+
#geom_text(position = 'stack', size=geom_textSize, aes( y = count, label = comma(count), hjust = 0.5, vjust = 1))+
scale_y_continuous(labels=scientific)+
{GGPLOT_PUB_QUALITY} + 
theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
\\"

{params.filterDat[facetPlotSetup]}

save_plot('{output[0]}', legendOnly, base_width=wLegendOnly, base_height=hLegendOnly)
save_plot('{output[1]}', pXy, base_width=wXyPlot, base_height=hXyPlot)

save_plot('{output[2]}', pXyNoLegend, base_width=wXyNoLegendPlot, base_height=hXyNoLegendPlot)
save_plot('{output[3]}', pYx, base_width=wYxPlot, base_height=hYxPlot)

save_plot('{output[4]}', pYxNoLegend, base_width=wYxNoLegendPlot, base_height=hYxNoLegendPlot)


" > $(dirname {output[0]})/$(basename {output[0]} .legendOnly.png).r

cat $(dirname {output[0]})/$(basename {output[0]} .legendOnly.png).r | R --slave


		'''



rule mergedReads:
	input: "output/mappings/highConfidenceReads/HiSS/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.HiSS.gff.gz"
	output: temp("output/mappings/mergedReads/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.HiSS.tmerge.min{minReadSupport}reads.endSupport-all.gff"),
	threads:1
	wildcard_constraints:
		sizeFrac='[0-9-+\.]+',
	shell:
		'''
uuidTmpOut=$(uuidgen)
zcat {input} | tmerge --exonOverhangTolerance {ExonOverhangTolerance} --minReadSupport {wildcards.minReadSupport} --endFuzz {ExonOverhangTolerance} --tmPrefix {wildcards.techname}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.sampleRep}.NAM_ - |sort -T {TMPDIR}  -k1,1 -k4,4n -k5,5n  > {TMPDIR}/$uuidTmpOut
mv {TMPDIR}/$uuidTmpOut {output}
		'''

rule mergedUnfilteredSirvReads:
	input: "output/mappings/strandGffs/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.stranded.gff.gz"
	output: "output/mappings/mergedReads/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.noFilt.tmerge.min{minReadSupport}reads.splicing_status-all.endSupport-all.gff.gz",
	threads:1
	wildcard_constraints:
		sizeFrac='[0-9-+\.]+',
	shell:
		'''
uuidTmpOut=$(uuidgen)
uuid=$(uuidgen)
zcat {input} | awk '$1 ~ /SIRV/ || $1=="chrIS"' > {TMPDIR}/$uuid

cat {TMPDIR}/$uuid| tmerge --exonOverhangTolerance {ExonOverhangTolerance} --minReadSupport {wildcards.minReadSupport} --endFuzz {ExonOverhangTolerance} --tmPrefix {wildcards.techname}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.sampleRep}.NAM_ - |sort -T {TMPDIR}  -k1,1 -k4,4n -k5,5n | gzip > {TMPDIR}/$uuidTmpOut
mv {TMPDIR}/$uuidTmpOut {output}
		'''



rule splitTmsBySplicedStatus:
	input: "output/mappings/mergedReads/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.HiSS.tmerge.min{minReadSupport}reads.endSupport-all.gff"
	output: "output/mappings/mergedReads/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.HiSS.tmerge.min{minReadSupport}reads.splicing_status-{splicedStatus}.endSupport-all.gff.gz"
	params:
		grepSpliced = lambda wildcards: '| fgrep \'spliced \"1\"\'' if wildcards.splicedStatus == "spliced" else '| fgrep \'spliced \"0\"\'' if wildcards.splicedStatus == "unspliced" else ''
	shell:
		'''
uuidTmpOut=$(uuidgen)
cat {input} {params.grepSpliced} | gzip > {TMPDIR}/$uuidTmpOut
mv {TMPDIR}/$uuidTmpOut {output}
		'''


rule catMergedReadsForGroupedSampleReps:
	input: getMergedSampleReps
	output: 
		gff=temp("output/mappings/mergedReads/groupedSampleReps/tmp/{groupedSampleRepBasename}.min{minReadSupport}reads.splicing_status-all.endSupport-all.gff"),
		countFiles=temp("output/mappings/mergedReads/groupedSampleReps/tmp/{groupedSampleRepBasename}.min{minReadSupport}reads.splicing_status-all.endSupport-all.countFiles.txt")
	shell:
		'''
uuid=$(uuidgen)
zcat {input} |sort -T $TMPDIR  -k1,1 -k4,4n -k5,5n > {TMPDIR}/$uuid
mv {TMPDIR}/$uuid {output.gff}

#calculate number of files in input:
echo {input} | awk '{{print NF}}' > {output.countFiles}
		'''

rule mergedReadsGroupedSampleReps:
	input: 
		gff="output/mappings/mergedReads/groupedSampleReps/tmp/{groupedSampleRepBasename}.min{minReadSupport}reads.splicing_status-all.endSupport-all.gff",
		countFiles="output/mappings/mergedReads/groupedSampleReps/tmp/{groupedSampleRepBasename}.min{minReadSupport}reads.splicing_status-all.endSupport-all.countFiles.txt"
	output: 
		gff="output/mappings/mergedReads/groupedSampleReps/{groupedSampleRepBasename}.min{minReadSupport}reads.splicing_status-all.endSupport-all.gff.gz",
		readToTmTsv="output/mappings/mergedReads/groupedSampleReps/{groupedSampleRepBasename}.min{minReadSupport}reads.splicing_status-all.endSupport-all.readsToTm.tsv.gz"
	shell:
		'''
uuid=$(uuidgen)


filesCount=$(cat {input.countFiles})
#min read support is 1 (1 input file) or 2 (more than 1 input files):
if [ $filesCount == 1 ]; then minRS=1; else minRS=2; fi;
echo "Min Read Support: $minRS"
cat {input.gff} | tmerge --minReadSupport $minRS --tmPrefix {wildcards.groupedSampleRepBasename}.NAM_ - |sort -T {TMPDIR}  -k1,1 -k4,4n -k5,5n  > {TMPDIR}/$uuid.gff

echo -e "read_id\ttranscript_id" > {TMPDIR}/$uuid.tsv

tmergeRecoverSampleRepReadIds.pl {input.gff} {TMPDIR}/$uuid.gff | sort -T $TMPDIR | uniq >> {TMPDIR}/$uuid.tsv

gzip {TMPDIR}/$uuid.tsv
gzip {TMPDIR}/$uuid.gff
mv {TMPDIR}/$uuid.gff.gz {output.gff}
mv {TMPDIR}/$uuid.tsv.gz {output.readToTmTsv}

		'''

rule mergedReadsToBed:
	input: "output/mappings/mergedReads/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.HiSS.tmerge.min{minReadSupport}reads.splicing_status-all.endSupport-all.gff.gz"
	output: temp("output/mappings/mergedReads/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.splicing_status-all.endSupport-all.bed")
	shell:
		'''
uuidTmpOut=$(uuidgen)
zcat {input} | gff2bed_full.pl - > {TMPDIR}/$uuidTmpOut
mv {TMPDIR}/$uuidTmpOut {output}
		'''



rule getTmStats:
	input: "output/mappings/mergedReads/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.HiSS.tmerge.min{minReadSupport}reads.splicing_status-all.endSupport-{endSupport}.gff.gz"
	output: "output/statsFiles/" + "tmp/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.min{minReadSupport}reads.splicing_status-all.endSupport-{endSupport}.TmStats.stats.tsv.gz"

	shell:
		'''
uuidTmpOut=$(uuidgen)
zcat {input} | extractGffAttributeValue.pl transcript_id spliced mature_RNA_length contains_count 3p_dists_to_3p 5p_dists_to_5p meta_3p_dists_to_5p meta_5p_dists_to_5p |sort|uniq | awk -v t={wildcards.techname} -v ca={wildcards.capDesign} -v s={wildcards.sizeFrac} -v b={wildcards.sampleRep} '{{print t"\\t"ca"\\t"s"\\t"b"\t"$0}}' | gzip > {TMPDIR}/$uuidTmpOut
mv {TMPDIR}/$uuidTmpOut {output}
		'''

rule aggTmStats:
	input: lambda wildcards: expand("output/statsFiles/" + "tmp/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.min{minReadSupport}reads.splicing_status-all.endSupport-{endSupport}.TmStats.stats.tsv.gz", filtered_product, techname=TECHNAMES, capDesign=wildcards.capDesign, sizeFrac=SIZEFRACS, sampleRep=SAMPLEREPS, minReadSupport=wildcards.minReadSupport, endSupport=wildcards.endSupport)
	output: "output/statsFiles/" + "all.{capDesign}.min{minReadSupport}.endSupport-{endSupport}.TmStats.stats.tsv.gz"
	shell:
		'''
uuidTmpOut=$(uuidgen)
echo -e "seqTech\tcapDesign\tsizeFrac\tsampleRep\tspliced\tcontains_count\tend\tdistance\tnormDistance" | gzip > {TMPDIR}/$uuidTmpOut

zcat {input} | perl -F"\\t" -slane '@ara=split(",", $F[8]); @arb=split(",", $F[9]); @arc=split(",", $F[10]), @ard=split(",", $F[11]); for ($i=0; $i<=$#ara; $i++){{$threepMinusDist=-$ara[$i];print  "$F[0]\\t$F[1]\\t$F[2]\\t$F[3]\\t$F[5]\\t$F[7]\\t3\\t$threepMinusDist\\t$arc[$i]\\n$F[0]\\t$F[1]\\t$F[2]\\t$F[3]\\t$F[5]\\t$F[7]\\t5\\t$arb[$i]\\t$ard[$i]"}}' | sort -T {TMPDIR}  | gzip >> {TMPDIR}/$uuidTmpOut
mv {TMPDIR}/$uuidTmpOut {output}

		'''



rule getMergingStats:
	input:
		hcgms = "output/mappings/highConfidenceReads/HiSS/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.HiSS.gff.gz",
		pooledMerged = "output/mappings/mergedReads/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.HiSS.tmerge.min{minReadSupport}reads.splicing_status-all.endSupport-all.gff.gz"
	output: "output/statsFiles/" + "tmp/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.min{minReadSupport}reads.merged.stats.tsv"
	shell:
		'''
uuidTmpOut=$(uuidgen)
hcgms=$(zcat {input.hcgms} | extractGffAttributeValue.pl transcript_id | sort -T {TMPDIR} |uniq|wc -l)
merged=$(zcat {input.pooledMerged} | extractGffAttributeValue.pl transcript_id | sort -T {TMPDIR} |uniq|wc -l)
echo -e "{wildcards.techname}\t{wildcards.capDesign}\t{wildcards.sizeFrac}\t{wildcards.sampleRep}\t$hcgms\t$merged" | awk '{{print $0"\t"$6/$5}}' > {TMPDIR}/$uuidTmpOut
mv {TMPDIR}/$uuidTmpOut {output}

		'''


rule aggMergingStats:
	input: lambda wildcards: expand("output/statsFiles/" + "tmp/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.min{minReadSupport}reads.merged.stats.tsv",filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, sampleRep=SAMPLEREPS, minReadSupport=wildcards.minReadSupport)

	output: "output/statsFiles/" + "all.min{minReadSupport}reads.merged.stats.tsv"
	shell:
		'''
uuidTmpOut=$(uuidgen)
echo -e "seqTech\tcapDesign\tsizeFrac\tsampleRep\tcategory\tcount" > {TMPDIR}/$uuidTmpOut
cat {input} | awk '{{print $1"\\t"$2"\\t"$3"\\t"$4"\\tHCGMreads\\t"$5"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\tmergedTMs\\t"$6}}' | sort -T {TMPDIR}  >> {TMPDIR}/$uuidTmpOut
mv {TMPDIR}/$uuidTmpOut {output}

		'''
rule plotMergingStats:
	input:  "output/statsFiles/" + "all.min{minReadSupport}reads.merged.stats.tsv"
	output: returnPlotFilenames("output/plots/" + "merged.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.min{minReadSupport}reads.merged.stats")
	params:
		filterDat=lambda wildcards: multi_figures(wildcards.capDesign, wildcards.sizeFrac, wildcards.sampleRep, wildcards.techname)
	conda: "envs/R_env.yml"
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
{params.filterDat[technameFilterString]}
{params.filterDat[capDesignFilterString]}

{params.filterDat[sizeFracFilterString]}
{params.filterDat[sampleRepFilterString]}
{params.filterDat[substSeqTechString]}
{params.filterDat[substSampleRepString]}
{params.filterDat[graphDimensions]}


plotBase <- \\"p <- ggplot(data=dat, aes(x=1, y=count, fill=category)) +
geom_bar(stat='identity', position=position_dodge()) +
geom_text(position = position_dodge(width = 0.9), size=geom_textSize, aes( y = 1, label = comma(count), hjust = 0, vjust = 0.5), angle=90) +
scale_fill_manual(values=c('HCGMreads' = '#d98cb3', 'mergedTMs' = '#cc9966')) +
ylab('# objects') +
xlab('') +
guides(fill = guide_legend(title='Category'))+
scale_y_continuous(labels=comma)+
{params.filterDat[hideXaxisLabels]}
{GGPLOT_PUB_QUALITY} + 
theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
\\"

{params.filterDat[facetPlotSetup]}

save_plot('{output[0]}', legendOnly, base_width=wLegendOnly, base_height=hLegendOnly)
save_plot('{output[1]}', pXy, base_width=wXyPlot, base_height=hXyPlot)

save_plot('{output[2]}', pXyNoLegend, base_width=wXyNoLegendPlot, base_height=hXyNoLegendPlot)
save_plot('{output[3]}', pYx, base_width=wYxPlot, base_height=hYxPlot)

save_plot('{output[4]}', pYxNoLegend, base_width=wYxNoLegendPlot, base_height=hYxNoLegendPlot)


" > $(dirname {output[0]})/$(basename {output[0]} .legendOnly.png).r

cat $(dirname {output[0]})/$(basename {output[0]} .legendOnly.png).r | R --slave


		'''


rule aggTmLengthStats:
	input:
		all=lambda wildcards: expand("output/statsFiles/" + "tmp/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.min{minReadSupport}reads.splicing_status-all.endSupport-all.TmStats.stats.tsv.gz", filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, sampleRep=SAMPLEREPS, minReadSupport=wildcards.minReadSupport),
		fl=lambda wildcards: expand("output/statsFiles/" + "tmp/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.min{minReadSupport}reads.splicing_status-all.endSupport-cagePolyASupported.TmStats.stats.tsv.gz", filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, sampleRep=SAMPLEREPS, minReadSupport=wildcards.minReadSupport)
	output: "output/statsFiles/" + "all.min{minReadSupport}reads.matureRNALength.stats.tsv.gz"
	shell:
		'''
uuidTmpOut=$(uuidgen)
echo -e "seqTech\tcapDesign\tsizeFrac\tsampleRep\ttranscript_id\tspliced\tmature_RNA_length\tcategory" |gzip> {TMPDIR}/$uuidTmpOut

zcat {input.all} |cut -f1-7| awk '{{print $0"\\tCLS_TMs"}}' |sort -T {TMPDIR}  | gzip >> {TMPDIR}/$uuidTmpOut
zcat {input.fl} |cut -f1-7| awk '{{print $0"\\tCLS_FL_TMs"}}' |sort -T {TMPDIR}  | gzip >> {TMPDIR}/$uuidTmpOut
mv {TMPDIR}/$uuidTmpOut {output}
		'''

rule aggHistTmLengthSummaryStats:
	input: "output/statsFiles/" + "all.min{minReadSupport}reads.matureRNALength.stats.tsv.gz"
	output: summary="output/statsFiles/" + "all.min{minReadSupport}reads.matureRNALengthSummary.stats.tsv"
	conda: "envs/R_env.yml"
	shell:
		'''
uuid=$(uuidgen)


echo "
library(tidyverse)
library(data.table)
dat<-fread('{input}', header=T, sep='\\t')
dat %>%
  group_by(seqTech, sizeFrac, capDesign, sampleRep, category) %>%
  summarise(med=median(mature_RNA_length), max=max(mature_RNA_length)) -> datSumm
write_tsv(datSumm, '{TMPDIR}/$uuid')
" | R --slave
mv {TMPDIR}/$uuid {output.summary}

		'''


rule plotHistTmLengthStats:
	input: "output/statsFiles/" + "all.min{minReadSupport}reads.matureRNALength.stats.tsv.gz"
	output: 
		hist=returnPlotFilenames("output/plots/" + "matureRNALength.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.min{minReadSupport}reads.splicing_status-{splicedStatus}.matureRNALength.hist.stats")
	params:
		filterDat=lambda wildcards: multi_figures(wildcards.capDesign, wildcards.sizeFrac, wildcards.sampleRep, wildcards.techname, wildcards.splicedStatus)
	conda: "envs/R_env.yml"
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
{params.filterDat[technameFilterString]}
{params.filterDat[capDesignFilterString]}

{params.filterDat[sizeFracFilterString]}
{params.filterDat[sampleRepFilterString]}
{params.filterDat[substSeqTechString]}
{params.filterDat[substSampleRepString]}
{params.filterDat[graphDimensions]}
{params.filterDat[splicingStatusFilterString]}

wXyPlot = wXyPlot * 1.2

dat %>%
  group_by(seqTech, sizeFrac, capDesign, sampleRep, category) %>%
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

{params.filterDat[facetPlotSetup]}

wYxPlot = wYxPlot * 1.2
wYxNoLegendPlot<- wYxPlot - wLegendOnly

save_plot('{output.hist[0]}', legendOnly, base_width=wLegendOnly, base_height=hLegendOnly)
save_plot('{output.hist[1]}', pXy, base_width=wXyPlot, base_height=hXyPlot)

save_plot('{output.hist[2]}', pXyNoLegend, base_width=wXyNoLegendPlot, base_height=hXyNoLegendPlot)
save_plot('{output.hist[3]}', pYx, base_width=wYxPlot, base_height=hYxPlot)

save_plot('{output.hist[4]}', pYxNoLegend, base_width=wYxNoLegendPlot, base_height=hYxNoLegendPlot)

" > $(dirname {output.hist[0]})/$(basename {output.hist[0]} .legendOnly.png).r

cat $(dirname {output.hist[0]})/$(basename {output.hist[0]} .legendOnly.png).r | R --slave


		'''


rule getGeneReadCoverageStats:
	input: 
		gencode="output/annotations/simplified/{capDesign}.gencode.simplified_biotypes.gtf",
		bam = "output/mappings/longReadMapping/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.bam",
		tmerge = "output/mappings/mergedReads/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.HiSS.tmerge.min2reads.splicing_status-all.endSupport-all.gff.gz",
		genome=lambda wildcards: config["GENOMESDIR"] + CAPDESIGNTOGENOME[wildcards.capDesign] + ".sorted.genome"
	output: 
		gencode="output/mappings/geneReadCoverage/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.gencode.coverage.tsv",
		tmerge="output/mappings/geneReadCoverage/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.coverage.tsv"
	wildcard_constraints:
		sizeFrac='[0-9-+\.]+',
	conda: "envs/xtools_env.yml"
	shell:
		'''
uuid=$(uuidgen)
cut -f1 {input.genome} |sort|uniq > {TMPDIR}/$uuid.chr

#gencode

cat {input.gencode} |awk '$3=="exon"' | extract_locus_coords.pl -| fgrep -w -f {TMPDIR}/$uuid.chr |sort -T {TMPDIR}  -k1,1 -k2,2n -k3,3n  > {TMPDIR}/$uuid.1

bedtools coverage -sorted -g {input.genome} -bed -split -nonamecheck -counts -a {TMPDIR}/$uuid.1 -b {input.bam} |sort -k7,7nr | cut -f1-4,7 | awk -v s={wildcards.techname} -v c={wildcards.capDesign} -v si={wildcards.sizeFrac} -v b={wildcards.sampleRep} '{{print s"\\t"c"\\t"si"\\t"b"\\t"$0}}'> {TMPDIR}/$uuid.2

mv {TMPDIR}/$uuid.2 {output.gencode}

#tmerge
zcat {input.tmerge} > {TMPDIR}/$uuid.a
bedtools intersect -s -wao -a {TMPDIR}/$uuid.a -b {TMPDIR}/$uuid.a | buildLoci.pl - | extract_locus_coords.pl -| sort -T {TMPDIR}  -k1,1 -k2,2n -k3,3n  > {TMPDIR}/$uuid.3

bedtools coverage -sorted -g {input.genome}  -bed -split -nonamecheck -counts -a {TMPDIR}/$uuid.3 -b {input.bam}  |sort -k7,7nr | cut -f1-4,7 | awk -v s={wildcards.techname} -v c={wildcards.capDesign} -v si={wildcards.sizeFrac} -v b={wildcards.sampleRep} '{{print s"\\t"c"\\t"si"\\t"b"\\t"$0}}' > {TMPDIR}/$uuid.4
mv {TMPDIR}/$uuid.4 {output.tmerge}

		'''

rule aggGeneReadCoverageStats:
	input: lambda wildcards: expand("output/mappings/geneReadCoverage/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.coverage.tsv", filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, sampleRep=SAMPLEREPS)
	output: 
		"output/statsFiles/" + "all.tmerge.GeneReadCoverage.stats.tsv"
	shell:
		'''
uuid=$(uuidgen)
echo -e "seqTech\\tcapDesign\\tsizeFrac\\tsampleRep\\tgene_id\\treadCount" > {TMPDIR}/$uuid
cat {input} | cut -f1-4,8,9   >> {TMPDIR}/$uuid
mv {TMPDIR}/$uuid {output}

		'''

rule plotGeneReadCoverageStats:
	input: "output/statsFiles/" + "all.tmerge.GeneReadCoverage.stats.tsv"
	output: returnPlotFilenames("output/plots/" + "geneReadCoverage.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.geneReadCoverage.stats")
	params:
		filterDat=lambda wildcards: multi_figures(wildcards.capDesign, wildcards.sizeFrac, wildcards.sampleRep, wildcards.techname)
	conda: "envs/R_env.yml"
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

{params.filterDat[technameFilterString]}
{params.filterDat[capDesignFilterString]}

{params.filterDat[sizeFracFilterString]}
{params.filterDat[sampleRepFilterString]}
{params.filterDat[substSeqTechString]}
{params.filterDat[substSampleRepString]}
{params.filterDat[graphDimensions]}

dat <- arrange(dat, desc(readCount))
group_by(dat, seqTech, capDesign, sizeFrac, sampleRep) %>% mutate(rank=row_number()) -> dat
mutate(dat, rank=row_number()) -> dat
dat <- mutate(dat, cumSum=cumsum(readCount))
dat <- mutate(dat, cumProp=cumSum/sum(readCount))
filter(dat, rank==10)  -> cumPropTop10Genes


plotBase <- \\"p <- ggplot(dat, aes(x = rank, y = cumProp)) + geom_line(size=lineSize) + scale_x_log10() + ylim(0,1)+ geom_segment(data = cumPropTop10Genes, aes(x=10,xend=10,y=0,yend=cumProp), color='firebrick3',size=lineSize) + geom_segment(data = cumPropTop10Genes, aes(x=0,xend=10,y=cumProp,yend=cumProp,color='Contribution\\nof top 10 genes'),size=lineSize) + labs(y='Proportion of total mapped reads', x='# genes (ranked by expression)', color='') + scale_color_manual(values = c('Contribution\nof top 10 genes' = 'firebrick3')) + facet_grid( seqTech ~ capDesign + sampleRep) + geom_rect(data = cumPropTop10Genes, aes(xmin=0,xmax=10,ymin=0,ymax=cumProp),fill='firebrick3', alpha=0.3, size=0) +
{GGPLOT_PUB_QUALITY} + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + \\"



{params.filterDat[facetPlotSetup]}

save_plot('{output[0]}', legendOnly, base_width=wLegendOnly, base_height=hLegendOnly)
save_plot('{output[1]}', pXy, base_width=wXyPlot, base_height=hXyPlot)

save_plot('{output[2]}', pXyNoLegend, base_width=wXyNoLegendPlot, base_height=hXyNoLegendPlot)
save_plot('{output[3]}', pYx, base_width=wYxPlot, base_height=hYxPlot)

save_plot('{output[4]}', pYxNoLegend, base_width=wYxNoLegendPlot, base_height=hYxNoLegendPlot)

" > $(dirname {output[0]})/$(basename {output[0]} .legendOnly.png).r

cat $(dirname {output[0]})/$(basename {output[0]} .legendOnly.png).r | R --slave

		'''
