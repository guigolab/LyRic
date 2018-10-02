rule extractPooledTmsFivepEnds:
	input: "mappings/" + "nonAnchoredMergeReads/pooled/{techname}Corr{corrLevel}_{capDesign}_pooled.tmerge.bed"
	output: "mappings/" + "nonAnchoredMergeReads/pooled/5pEnds/{techname}Corr{corrLevel}_{capDesign}_pooled.tmerge.5pEnds.bed"
	shell:
		'''
cat {input} |extractTranscriptEndsFromBed12.pl 5 |sortbed> {output}
		'''

rule cageSupportedfivepEnds:
	input:
		fivePends="mappings/" + "nonAnchoredMergeReads/pooled/5pEnds/{techname}Corr{corrLevel}_{capDesign}_pooled.tmerge.5pEnds.bed",
		tms="mappings/" + "nonAnchoredMergeReads/pooled/{techname}Corr{corrLevel}_{capDesign}_pooled.tmerge.bed",
		cagePeaks=lambda wildcards: CAPDESIGNTOCAGEPEAKS[wildcards.capDesign]
	params: genome = lambda wildcards: config["GENOMESDIR"] + CAPDESIGNTOGENOME[wildcards.capDesign] + ".genome"
	output: "mappings/" + "nonAnchoredMergeReads/pooled/cageSupported/{techname}Corr{corrLevel}_{capDesign}_pooled.tmerge.cageSupported.bed"
	shell:
		'''
cat {input.fivePends} | sortbed | bedtools slop -s -l 50 -r 50 -i stdin -g {params.genome} | bedtools intersect -u -s -a stdin -b {input.cagePeaks} | cut -f4 | fgrep -w -f - {input.tms} > {output}
		'''


rule extractPooledTmsThreepEnds:
	input: "mappings/" + "nonAnchoredMergeReads/pooled/{techname}Corr{corrLevel}_{capDesign}_pooled.tmerge.bed"
	output: "mappings/" + "nonAnchoredMergeReads/pooled/3pEnds/{techname}Corr{corrLevel}_{capDesign}_pooled.tmerge.3pEnds.bed"
	shell:
		'''
cat {input} |extractTranscriptEndsFromBed12.pl 3 |sortbed> {output}
		'''

rule polyASupportedthreepEnds:
	input:
		threePends="mappings/" + "nonAnchoredMergeReads/pooled/3pEnds/{techname}Corr{corrLevel}_{capDesign}_pooled.tmerge.3pEnds.bed",
		tms="mappings/" + "nonAnchoredMergeReads/pooled/{techname}Corr{corrLevel}_{capDesign}_pooled.tmerge.bed",
		polyAsites=lambda wildcards: expand("mappings/" + "removePolyAERCCs/{techname}Corr{corrLevel}_{capDesign}_{barcodes}.polyAsitesNoErcc.bed", filtered_product, techname=wildcards.techname, corrLevel=wildcards.corrLevel, capDesign=wildcards.capDesign, barcodes=BARCODES)
	params: genome = lambda wildcards: config["GENOMESDIR"] + CAPDESIGNTOGENOME[wildcards.capDesign] + ".genome"
	output: "mappings/nonAnchoredMergeReads/pooled/polyASupported/{techname}Corr{corrLevel}_{capDesign}_pooled.tmerge.polyASupported.bed"
	shell:
		'''
cat {input.polyAsites} |sortbed > $TMPDIR/polyAsites.bed
cat {input.threePends} | sortbed | bedtools slop -s -l 5 -r 5 -i stdin -g {params.genome} | bedtools intersect -u -s -a stdin -b $TMPDIR/polyAsites.bed | cut -f4 | fgrep -w -f - {input.tms} > {output}
		'''

rule getCagePolyASupport:
	input:
		polyA="mappings/nonAnchoredMergeReads/pooled/polyASupported/{techname}Corr{corrLevel}_{capDesign}_pooled.tmerge.polyASupported.bed",
		cage="mappings/" + "nonAnchoredMergeReads/pooled/cageSupported/{techname}Corr{corrLevel}_{capDesign}_pooled.tmerge.cageSupported.bed",
		tms="mappings/" + "nonAnchoredMergeReads/pooled/{techname}Corr{corrLevel}_{capDesign}_pooled.tmerge.bed"
	output:
		stats=temp(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}.tmerge.cagePolyASupport.stats.tsv"),
		FLbed="mappings/nonAnchoredMergeReads/pooled/cage+polyASupported/{techname}Corr{corrLevel}_{capDesign}_pooled.tmerge.cage+polyASupported.bed"
	shell:
		'''
cat {input.polyA} | cut -f4 | sort|uniq > $TMPDIR/polyA.list
cat {input.cage} | cut -f4 | sort|uniq > $TMPDIR/cage.list
cat {input.tms} | cut -f4 | sort|uniq > $TMPDIR/all.list
cat $TMPDIR/polyA.list $TMPDIR/cage.list |sort|uniq > $TMPDIR/cageOrPolyA.list
comm -1 -2 $TMPDIR/polyA.list $TMPDIR/cage.list |sort|uniq > $TMPDIR/cage+PolyA.list
noCageNoPolyA=$(comm -2 -3 $TMPDIR/all.list $TMPDIR/cageOrPolyA.list |wc -l)
cageOnly=$(comm -2 -3 $TMPDIR/cage.list $TMPDIR/polyA.list |wc -l)
polyAOnly=$(comm -2 -3 $TMPDIR/polyA.list $TMPDIR/cage.list |wc -l)
cageAndPolyA=$(cat $TMPDIR/cage+PolyA.list | wc -l)
let total=$noCageNoPolyA+$cageOnly+$polyAOnly+$cageAndPolyA
fgrep -w -f $TMPDIR/cage+PolyA.list {input.tms} > {output.FLbed}
echo -e "{wildcards.techname}Corr{wildcards.corrLevel}\t{wildcards.capDesign}\tcageOnly\t$cageOnly
{wildcards.techname}Corr{wildcards.corrLevel}\t{wildcards.capDesign}\tcageAndPolyA\t$cageAndPolyA
{wildcards.techname}Corr{wildcards.corrLevel}\t{wildcards.capDesign}\tpolyAOnly\t$polyAOnly
{wildcards.techname}Corr{wildcards.corrLevel}\t{wildcards.capDesign}\tnoCageNoPolyA\t$noCageNoPolyA" | awk -v t=$total '{{print $0"\t"$4/t}}'> {output.stats}
		'''

rule aggCagePolyAStats:
	input: lambda wildcards: expand(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}.tmerge.cagePolyASupport.stats.tsv", techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS)
	output: config["STATSDATADIR"] + "all.tmerge.cagePolyASupport.stats.tsv"
	shell:
		'''
echo -e "seqTech\tcorrectionLevel\tcapDesign\tcategory\tcount\tpercent" > {output}
cat {input} | sed 's/Corr0/\tNo/' | sed 's/Corr{lastK}/\tYes/' | sort >> {output}
		'''

rule plotCagePolyAStats:
	input: config["STATSDATADIR"] + "all.tmerge.cagePolyASupport.stats.tsv"
	output: config["PLOTSDIR"] + "all.tmerge.cagePolyASupport.stats.{ext}"
	shell:
		'''
echo "library(ggplot2)
library(plyr)
library(scales)
dat <- read.table('{input}', header=T, as.is=T, sep='\\t')
dat\$category<-factor(dat\$category, ordered=TRUE, levels=rev(c('cageOnly', 'cageAndPolyA', 'polyAOnly', 'noCageNoPolyA')))
dat\$seqTech <- gsub(':', '\\n', dat\$seqTech)
horizCats <- length(unique(dat\$correctionLevel)) * length(unique(dat\$capDesign))
vertCats <- length(unique(dat\$seqTech))
plotWidth = horizCats + 3.5
plotHeight = vertCats +2
ggplot(dat[order(dat\$category), ], aes(x=factor(correctionLevel), y=count, fill=category)) +
geom_bar(stat='identity') + ylab('# CLS TMs') +
scale_y_continuous(labels=comma)+ scale_fill_manual (values=c(cageOnly='#66B366', cageAndPolyA='#82865f', polyAOnly = '#D49090', noCageNoPolyA='#a6a6a6'))+ facet_grid( seqTech ~ capDesign)+ xlab('Error correction') + guides(fill = guide_legend(title='Category'))+
geom_text(position = 'stack', aes(x = factor(correctionLevel), y = count, ymax=count, label = paste(sep='',percent(round(percent, digits=2)),' / ','(',comma(count),')'), hjust = 0.5, vjust = 1))+
{GGPLOT_PUB_QUALITY}
ggsave('{output}', width=plotWidth, height=plotHeight)
" > {output}.r
cat {output}.r | R --slave

		'''
