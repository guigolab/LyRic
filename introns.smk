rule makeIntrons:
	input: "mappings/readBedToGff/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.gff.gz"
	output: "mappings/makeIntrons/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.introns.gff.gz"
	shell:
		'''
uuidTmpOut=$(uuidgen)
zcat {input} | makeIntrons.pl - | sort -T {config[TMPDIR]}  -k1,1 -k4,4n -k5,5n  |gzip> {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
		'''

rule getIntronMotif:
	input:
		introns = "mappings/makeIntrons/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.introns.gff.gz",
		genome = lambda wildcards: config["GENOMESDIR"] + CAPDESIGNTOGENOME[wildcards.capDesign] + ".sorted.fa"
	output:
		gff = "mappings/getIntronMotif/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.introns.gff.gz",
		tsv = "mappings/getIntronMotif/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.transcripts.tsv.gz"
	shell:
		'''
uuid=$(uuidgen)
mkdir -p {config[TMPDIR]}/$uuid
set +eu
conda activate perl_env
set -eu

zcat {input.introns} | grep -vP "^ERCC"| extract_intron_strand_motif.pl - {input.genome} {config[TMPDIR]}/$uuid/$(basename {output.gff} .introns.gff.gz)

gzip {config[TMPDIR]}/$uuid/*
mv {config[TMPDIR]}/$uuid/* $(dirname {output.gff})
		'''

rule getGencodePcgSpliceSites:
	input: lambda wildcards: CAPDESIGNTOANNOTGTF[wildcards.capDesign]
	output: "annotations/spliceSites/{capDesign}.gencode.PCG.spliceSites.{spliceType}.tsv.gz"
	shell:
		'''
uuidTmpOut=$(uuidgen)
cat {input} | awk '$3=="exon"' | fgrep "transcript_type \\"protein_coding\\";" |sort -T {config[TMPDIR]}  -k1,1 -k4,4n -k5,5n | makeIntrons.pl - | awk '{{print $1"\t"$4-1"\t"$5+1"\t"$7}}' |sort -T {config[TMPDIR]} |uniq | perl -lane 'if($F[3] eq "+"){{$dStart=$F[1]-1; $aStart=$F[2]-2}} elsif($F[3] eq "-"){{$dStart=$F[2]-2; $aStart=$F[1]-1}} else{{die}} $aEnd=$aStart+2; $dEnd=$dStart+2; print "$F[0]"."_$dStart"."_$dEnd"."_$F[3]\t$F[0]"."_"."$F[1]"."_"."$F[2]"."_"."$F[3]:Donor:GENCODE_protein_coding"; print "$F[0]"."_$aStart"."_$aEnd"."_$F[3]\t$F[0]"."_"."$F[1]"."_"."$F[2]"."_"."$F[3]:Acceptor:GENCODE_protein_coding";' | fgrep -w {wildcards.spliceType}| sort -T {config[TMPDIR]}  |uniq| sort -T {config[TMPDIR]}  -k1,1 | gzip > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
		'''

rule getReadsSpliceSites:
	input: "mappings/strandGffs/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.stranded.gff.gz"
	output: "mappings/makeIntrons/readSpliceSites/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.introns.{spliceType}.tsv.gz"
	shell:
		'''
uuidTmpOut=$(uuidgen)
zcat {input} | awk '$3=="exon"' | makeIntrons.pl -| awk '{{print $1"\t"$4-1"\t"$5+1"\t"$7}}' |sort -T {config[TMPDIR]} |uniq | perl -lane 'if($F[3] eq "+"){{$dStart=$F[1]-1; $aStart=$F[2]-2}} elsif($F[3] eq "-"){{$dStart=$F[2]-2; $aStart=$F[1]-1}} else{{die}} $aEnd=$aStart+2; $dEnd=$dStart+2; print "$F[0]"."_$dStart"."_$dEnd"."_$F[3]\t$F[0]"."_"."$F[1]"."_"."$F[2]"."_"."$F[3]:Donor:CLS_reads"; print "$F[0]"."_$aStart"."_$aEnd"."_$F[3]\t$F[0]"."_"."$F[1]"."_"."$F[2]"."_"."$F[3]:Acceptor:CLS_reads";' | fgrep -w {wildcards.spliceType} |sort -T {config[TMPDIR]} |uniq| sort -T {config[TMPDIR]}  -k1,1 |gzip > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
		'''

rule getTmSpliceSites:
	input: "mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.tmerge.min{minReadSupport}reads.splicing_status:spliced.endSupport:all.gff.gz"
	output: "mappings/nonAnchoredMergeReads/spliceSites/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.introns.{spliceType}.tsv.gz"
	shell:
		'''
uuidTmpOut=$(uuidgen)
zcat {input} | awk '$3=="exon"' | makeIntrons.pl - | awk '{{print $1"\t"$4-1"\t"$5+1"\t"$7}}' |sort -T {config[TMPDIR]} |uniq | perl -lane 'if($F[3] eq "+"){{$dStart=$F[1]-1; $aStart=$F[2]-2}} elsif($F[3] eq "-"){{$dStart=$F[2]-2; $aStart=$F[1]-1}} else{{die}} $aEnd=$aStart+2; $dEnd=$dStart+2; print "$F[0]"."_$dStart"."_$dEnd"."_$F[3]\t$F[0]"."_"."$F[1]"."_"."$F[2]"."_"."$F[3]:Donor:CLS_TMs"; print "$F[0]"."_$aStart"."_$aEnd"."_$F[3]\t$F[0]"."_"."$F[1]"."_"."$F[2]"."_"."$F[3]:Acceptor:CLS_TMs";' | fgrep -w {wildcards.spliceType} | sort -T {config[TMPDIR]}  |uniq| sort -T {config[TMPDIR]}  -k1,1 | gzip> {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
		'''

rule getGeneidScores:
	input:
		rawReads="mappings/makeIntrons/readSpliceSites/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.introns.{spliceType}.tsv.gz",
		tmReads="mappings/nonAnchoredMergeReads/spliceSites/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.introns.{spliceType}.tsv.gz",
		geneidScores= lambda wildcards: expand(config["SPLICE_SITE_SCORES_DIR"] + CAPDESIGNTOGENOME[wildcards.capDesign] + ".{spliceType}.geneid.loose.spliceSites.sorted.tsv", spliceType=wildcards.spliceType),
	output:  config["STATSDATADIR"] + "tmp/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.{spliceType}.spliceSites.stats.tsv"
	shell:
		'''
uuidTmpOut=$(uuidgen)
cls=$(uuidgen)
zcat {input.rawReads} {input.tmReads} | sort -T {config[TMPDIR]}   -k1,1 > {config[TMPDIR]}/$cls.SSs.tsv
join -a1 -j1 {config[TMPDIR]}/$cls.SSs.tsv {input.geneidScores} |ssv2tsv | perl -lane '$_=~s/\:/\t/g; @line=split("\\t", $_); unless(defined ($line[5])){{$line[4]=$line[2];$line[5]=-35 + rand(-25 +35)}}; print join("\\t", @line)' | awk -v t={wildcards.techname}Corr{wildcards.corrLevel} -v c={wildcards.capDesign} -v si={wildcards.sizeFrac} -v b={wildcards.barcodes} '{{print t"\t"c"\t"si"\t"b"\t"$1"\t"$3"\t"$4"\t"$6}}'| sed 's/Corr0/\tNo/' | sed 's/Corr{lastK}/\tYes/' > {config[TMPDIR]}/$uuidTmpOut

#verify that all SSs were found:
cut -f1 {config[TMPDIR]}/$cls.SSs.tsv | sort -T {config[TMPDIR]} |uniq > {config[TMPDIR]}/$cls.SSs.list
cut -f6 {config[TMPDIR]}/$uuidTmpOut | sort -T {config[TMPDIR]} |uniq > {config[TMPDIR]}/$cls.SSs.geneid.list
diff=$(diff -q {config[TMPDIR]}/$cls.SSs.list {config[TMPDIR]}/$cls.SSs.geneid.list |wc -l)
if [ ! $diff -eq 0 ]; then echoerr "ERROR: List of SSs differ before/after"; exit 1; fi

mv {config[TMPDIR]}/$uuidTmpOut {output}
		'''

rule getControlGeneidScores:
	input:
		gencode ="annotations/spliceSites/{capDesign}.gencode.PCG.spliceSites.{spliceType}.tsv.gz",
		random =lambda wildcards: "annotations/spliceSites/" + CAPDESIGNTOGENOME[wildcards.capDesign] + ".geneid.loose.spliceSites.OnDetectedRegions.NotDetected.random100K.{spliceType}.tsv",
		geneidScores= lambda wildcards: expand(config["SPLICE_SITE_SCORES_DIR"] + CAPDESIGNTOGENOME[wildcards.capDesign] + ".{spliceType}.geneid.loose.spliceSites.sorted.tsv", spliceType=wildcards.spliceType),

	output: config["STATSDATADIR"] + "tmp/{capDesign}.control.spliceSites.{spliceType}.stats.tsv"
	shell:
		'''
uuid=$(uuidgen)
uuidTmpOut=$(uuidgen)
zcat {input.gencode} | sort -T {config[TMPDIR]}   -k1,1 > {config[TMPDIR]}/$uuid.SSs.tsv
join -a1 -j1 {config[TMPDIR]}/$uuid.SSs.tsv {input.geneidScores} |ssv2tsv | perl -lane '$_=~s/\:/\t/g; @line=split("\\t", $_); unless(defined ($line[5])){{$line[4]=$line[2];$line[5]=-35 + rand(-25 +35)}}; print join("\\t", @line)' | awk '{{print $1"\t"$3"\t"$4"\t"$6}}' > {config[TMPDIR]}/$uuidTmpOut
cat {input.random} | awk '{{print $1"\t"$2"\trandom\t"$3}}' >> {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}


		'''

rule aggControlGeneidScores:
	input: lambda wildcards: expand(config["STATSDATADIR"] + "tmp/{capDesign}.control.spliceSites.{spliceType}.stats.tsv", capDesign=CAPDESIGNS,  spliceType=SPLICE_SITE_TYPES)
	output: config["STATSDATADIR"] + "all.control.spliceSites.stats.tsv"
	shell:
		'''
uuidTmpOut=$(uuidgen)
echo -e "ssCategory\tssScore" > {config[TMPDIR]}/$uuidTmpOut
cat {input} | sort -T {config[TMPDIR]} |uniq| cut -f 3,4 >> {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}

		'''


rule aggGeneidScores:
	input: lambda wildcards: expand(config["STATSDATADIR"] + "tmp/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.{spliceType}.spliceSites.stats.tsv",filtered_product_merge, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODESpluSMERGED, spliceType=SPLICE_SITE_TYPES, minReadSupport=wildcards.minReadSupport)
	output: config["STATSDATADIR"] + "all.min{minReadSupport}reads.spliceSites.stats.tsv"
	threads: 8
	shell:
		'''
uuidTmpOut=$(uuidgen)
echo -e "seqTech\tcorrectionLevel\tcapDesign\tsizeFrac\ttissue\tssCategory\tssScore" > {config[TMPDIR]}/$uuidTmpOut
uuid=$(uuidgen)
#keep non-redundant SSs
cat {input} | cut -f1-6,8,9 > {config[TMPDIR]}/$uuid
sort -T {config[TMPDIR]}  -S 14G --parallel {threads} {config[TMPDIR]}/$uuid | uniq | cut -f1-5,7,8 >> {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}

		'''

rule plotGeneidScores:
	input:
		cls=config["STATSDATADIR"] + "all.min{minReadSupport}reads.spliceSites.stats.tsv",
		control=config["STATSDATADIR"] + "all.control.spliceSites.stats.tsv"
	output: returnPlotFilenames(config["PLOTSDIR"] + "spliceSites.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.spliceSites.stats")
	params:
		filterDat=lambda wildcards: merge_figures_params(wildcards.capDesign, wildcards.sizeFrac, wildcards.barcodes, wildcards.corrLevel, wildcards.techname)
	shell:
		'''
echo "
library(data.table)
library(ggplot2)
library(cowplot)
library(plyr)
library(scales)
library(gridExtra)
library(grid)
library(ggplotify)

palette <- c('random' = '#999999', 'GENCODE_protein_coding' = '#009900', 'CLS_reads' = '#b3d9ff', 'CLS_TMs' = '#cc9966')
dat<-fread('{input.cls}', header=T, sep='\\t')
dat2<-fread('{input.control}', header=T, sep='\\t')
{params.filterDat[10]}
{params.filterDat[0]}
{params.filterDat[1]}
{params.filterDat[2]}
{params.filterDat[3]}
{params.filterDat[4]}
{params.filterDat[5]}
{params.filterDat[8]}

fun_length <- function(x){{
return(data.frame(y=-8.5,label= paste0('N=', comma(length(x)))))
}}

plotBase <- \\"p <- ggplot(dat, aes(x=ssCategory, y=ssScore, color=ssCategory)) +
geom_boxplot(position=position_dodge(0.9), outlier.shape=NA) +
coord_cartesian(ylim=c(-9, 4.5)) +
scale_color_manual(values=palette, name='Category', labels = c(random = 'Random', GENCODE_protein_coding = 'GENCODE\nprotein-coding', CLS_TMs='CLS TMs', CLS_reads='CLS raw reads')) +
geom_boxplot(data=dat2, position=position_dodge(0.9), outlier.shape=NA) +
stat_summary(aes(x=ssCategory, group=ssCategory), position=position_dodge(0.9), fun.data = fun_length, geom = 'text', vjust = +1, hjust=0, angle=90, show.legend=FALSE, color='black', size=geom_textSize) +
stat_summary(data=dat2, aes(x=ssCategory, group=ssCategory), position=position_dodge(0.9), fun.data = fun_length, geom = 'text', vjust = +1, hjust=0, angle=90, show.legend=FALSE, color='black', size=geom_textSize) +
geom_hline(aes(yintercept=0), linetype='dashed', alpha=0.7, size=lineSize)+
ylab('Splice site score') +
xlab('{params.filterDat[6]}') +
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
 set +eu
conda activate R_env
set -eu
cat $(dirname {output[0]})/$(basename {output[0]} .legendOnly.png).r | R --slave




		'''


rule getGencodeSpliceJunctions:
	input: lambda wildcards: CAPDESIGNTOANNOTGTF[wildcards.capDesign]
	output: "annotations/spliceJunctions/{capDesign}.gencode.spliceJunctions.list"
	shell:
		'''
uuidTmpOut=$(uuidgen)
cat {input} | awk '$3=="exon"' |sort -T {config[TMPDIR]}  -k1,1 -k4,4n -k5,5n | makeIntrons.pl - | awk '{{print $1"_"$4"_"$5"_"$7}}' |sort -T {config[TMPDIR]} |uniq > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}

		'''

rule getClsSpliceJunctions:
	input:"mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.tmerge.min{minReadSupport}reads.splicing_status:spliced.endSupport:all.gff.gz"
	output: "mappings/nonAnchoredMergeReads/spliceJunctions/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.spliceJunctions.list"
	shell:
		'''
uuidTmpOut=$(uuidgen)
zcat {input} | awk '$3=="exon"' |sort -T {config[TMPDIR]}  -k1,1 -k4,4n -k5,5n | makeIntrons.pl - | awk '{{print $1"_"$4"_"$5"_"$7}}' |sort -T {config[TMPDIR]} |uniq > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}

		'''

rule getCompareClsGencodeSJsStats:
	input:
		gencodeSJs="annotations/spliceJunctions/{capDesign}.gencode.spliceJunctions.list",
		clsSJs="mappings/nonAnchoredMergeReads/spliceJunctions/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.spliceJunctions.list"
	output: config["STATSDATADIR"] + "tmp/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.vs.Gencode.SJs.stats.tsv"
	shell:
		'''
uuidTmpOut=$(uuidgen)
#annSJs=$(cat {input.gencodeSJs} | wc -l)
clsSJs=$(cat {input.clsSJs} | wc -l)
commonSJs=$(comm -1 -2 {input.gencodeSJs} {input.clsSJs} | wc -l)
novelSJs=$(comm -1 -3 {input.gencodeSJs} {input.clsSJs} | wc -l)
echo -e "{wildcards.techname}Corr{wildcards.corrLevel}\t{wildcards.capDesign}\t{wildcards.sizeFrac}\t{wildcards.barcodes}\t$clsSJs\t$commonSJs\t$novelSJs"  > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}

		'''

rule aggCompareClsGencodeSJsStats:
	input: lambda wildcards: expand(config["STATSDATADIR"] +"tmp/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.vs.Gencode.SJs.stats.tsv",filtered_product_merge, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODESpluSMERGED, minReadSupport=wildcards.minReadSupport)
	output: config["STATSDATADIR"] + "all.tmerge.min{minReadSupport}reads.vs.Gencode.SJs.stats.tsv"
	shell:
		'''
uuidTmpOut=$(uuidgen)
echo -e "seqTech\tcorrectionLevel\tcapDesign\tsizeFrac\ttissue\tcategory\tcount\tpercent" > {config[TMPDIR]}/$uuidTmpOut
cat {input} | awk '{{ print $1"\\t"$2"\\t"$3"\\t"$4"\\tcommon\\t"$6"\t"$6/$5"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\tnovel\\t"$7"\t"$7/$5}}'| sed 's/Corr0/\tNo/' | sed 's/Corr{lastK}/\tYes/' | sort -T {config[TMPDIR]}  >> {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}


		'''

rule plotCompareClsGencodeSJsStats:
	input: config["STATSDATADIR"] + "all.tmerge.min{minReadSupport}reads.vs.Gencode.SJs.stats.tsv"
	output: returnPlotFilenames(config["PLOTSDIR"] + "tmerge.vs.Gencode.SJs.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.vs.Gencode.SJs.stats")
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
library(data.table)

dat <- read.table('{input}', header=T, as.is=T, sep='\\t')
{params.filterDat[10]}
{params.filterDat[0]}
{params.filterDat[1]}
{params.filterDat[2]}
{params.filterDat[3]}
{params.filterDat[4]}
{params.filterDat[5]}
{params.filterDat[8]}

dat\$category<-factor(dat\$category, ordered=TRUE, levels=rev(c('common', 'novel')))

plotBase <- \\"p <- ggplot(dat[order(dat\$category), ], aes(x=factor(correctionLevel), y=count, fill=category)) +
geom_bar(stat='identity') +
scale_fill_manual (values=c(annOnly='#7D96A2',common='#83A458', novel='#B8CF7E'), labels=c(annOnly='Only in GENCODE', common='In sample+GENCODE', novel='Only in sample')) +
ylab('# Splice Junctions')+
xlab('{params.filterDat[6]}') +
geom_text(position = 'stack', size=geom_textSize, aes(x = factor(correctionLevel), y = count, label = paste(sep='',percent(round(percent, digits=2)),' / ','(',comma(count),')'), hjust = 0.5, vjust = 1))+
xlab('Error correction') +
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
