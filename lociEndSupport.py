rule collapseGencode:
	input: "annotations/simplified/{capDesign}.gencode.simplified_biotypes.gtf"
	output: "annotations/simplified/{capDesign}.gencode.collapsed.simplified_biotypes.gtf"
	threads:8
	shell:
		'''
uuid=$(uuidgen)
cat {input} | skipcomments | sortgff | tmerge --cpu {threads} --exonOverhangTolerance 25 - |sortgff > $TMPDIR/$uuid
uuidL=$(uuidgen)
bedtools intersect -s -wao -a $TMPDIR/$uuid -b $TMPDIR/$uuid | buildLoci.pl - |sortgff > $TMPDIR/$uuidL
mergeToRef.pl {input} $TMPDIR/$uuidL | sortgff > {output}

		'''


rule mergeCurrentPreviousPhaseTmsWithGencode:
	input:
		current=lambda wildcards: expand("mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.all.gff", filtered_product_merge, techname=wildcards.techname, corrLevel=wildcards.corrLevel, capDesign=wildcards.capDesign, sizeFrac=wildcards.sizeFrac, barcodes=wildcards.barcodes, minReadSupport=wildcards.minReadSupport),
		previous=lambda wildcards: GENOMETOPREVIOUS[CAPDESIGNTOGENOME[wildcards.capDesign]],
		annot="annotations/simplified/{capDesign}.gencode.collapsed.simplified_biotypes.gtf"
	threads:8
	output: "mappings/nonAnchoredMergeReads/mergeWithPrevious/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.all.gff.gz"
	# wildcard_constraints:
	# 	barcodes='allTissues',
	# 	sizeFrac='allFracs',
	# 	techname='allSeqTechs'
	shell:
		'''
uuid=$(uuidgen)
uuidM=$(uuidgen)

cat {input.previous} {input.annot} {input.current} | skipcomments | sortgff | tmerge --cpu {threads} --exonOverhangTolerance 25 - |sortgff > $TMPDIR/$uuidM

bedtools intersect -s -wao -a $TMPDIR/$uuidM -b $TMPDIR/$uuidM |fgrep -v ERCC| buildLoci.pl - |sortgff | gzip> {output}

		'''

rule mergeCurrentPreviousPhaseTmsWithGencodeBiotypes:
	input:
		clsGencode="mappings/nonAnchoredMergeReads/mergeWithPrevious/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.all.gff.gz",
		gencode="annotations/simplified/{capDesign}.gencode.collapsed.simplified_biotypes.gtf"
	output: "mappings/nonAnchoredMergeReads/mergeWithPrevious+biotypes/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.all.gff.gz"
	shell:
		'''
uuid=$(uuidgen)
zcat  {input.clsGencode} > $TMPDIR/$uuid
mergeToRef.pl {input.gencode} $TMPDIR/$uuid | sortgff |gzip > {output}

		'''


rule getGencodeSupportedEnds:
	input:
		tm="annotations/simplified/{capDesign}.gencode.collapsed.simplified_biotypes.gtf",
		cagePeaks=lambda wildcards: CAPDESIGNTOCAGEPEAKS[wildcards.capDesign],
		genome = lambda wildcards: config["GENOMESDIR"] + CAPDESIGNTOGENOME[wildcards.capDesign] + ".genome",
		PAS=lambda wildcards: GENOMETOPAS[CAPDESIGNTOGENOME[wildcards.capDesign]],

	output:"mappings/nonAnchoredMergeReads/mergeWithPrevious/gencode/{capDesign}.gencode.cage+PASsupported.gff.gz"
	shell:
		'''
uuid5pEnds=$(uuidgen)
cat {input.tm} |gff2bed_full.pl - |extractTranscriptEndsFromBed12.pl 5 |sortbed> $TMPDIR/$uuid5pEnds
uuid3pEnds=$(uuidgen)
cat {input.tm} |gff2bed_full.pl - |extractTranscriptEndsFromBed12.pl 3 |sortbed> $TMPDIR/$uuid3pEnds

uuidCageSupported=$(uuidgen)
cat $TMPDIR/$uuid5pEnds | sortbed | bedtools slop -s -l 50 -r 50 -i stdin -g {input.genome} | bedtools intersect -u -s -a stdin -b {input.cagePeaks} | cut -f4  |sort|uniq > $TMPDIR/$uuidCageSupported

uuidPASsupported=$(uuidgen)
cat $TMPDIR/$uuid3pEnds | sortbed | bedtools slop -s -l 50 -r -10 -i stdin -g {input.genome} | bedtools intersect -u -s -a stdin -b {input.PAS} |cut -f4 |sort|uniq > $TMPDIR/$uuidPASsupported

uuidcagePASsupported=$(uuidgen)
comm -1 -2 $TMPDIR/$uuidCageSupported $TMPDIR/$uuidPASsupported |sort|uniq > $TMPDIR/$uuidcagePASsupported

fgrep -w -f $TMPDIR/$uuidcagePASsupported {input.tm} |sortgff | gzip> {output}



		'''

rule getCurrentPreviousPhaseTmsWithGencodeSupportedEnds:
	input:
		tm="mappings/nonAnchoredMergeReads/mergeWithPrevious+biotypes/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.all.gff.gz",
		cagePeaks=lambda wildcards: CAPDESIGNTOCAGEPEAKS[wildcards.capDesign],
		genome = lambda wildcards: config["GENOMESDIR"] + CAPDESIGNTOGENOME[wildcards.capDesign] + ".genome",
		PAS=lambda wildcards: GENOMETOPAS[CAPDESIGNTOGENOME[wildcards.capDesign]],

	output:"mappings/nonAnchoredMergeReads/mergeWithPrevious/cls/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.cage+PASsupported.gff.gz"
	shell:
		'''
uuid5pEnds=$(uuidgen)
zcat {input.tm} |gff2bed_full.pl - |extractTranscriptEndsFromBed12.pl 5 |sortbed> $TMPDIR/$uuid5pEnds
uuid3pEnds=$(uuidgen)
zcat {input.tm} |gff2bed_full.pl - |extractTranscriptEndsFromBed12.pl 3 |sortbed> $TMPDIR/$uuid3pEnds

uuidCageSupported=$(uuidgen)
cat $TMPDIR/$uuid5pEnds | sortbed | bedtools slop -s -l 50 -r 50 -i stdin -g {input.genome} | bedtools intersect -u -s -a stdin -b {input.cagePeaks} | cut -f4  |sort|uniq > $TMPDIR/$uuidCageSupported

uuidPASsupported=$(uuidgen)
cat $TMPDIR/$uuid3pEnds | sortbed | bedtools slop -s -l 50 -r -10 -i stdin -g {input.genome} | bedtools intersect -u -s -a stdin -b {input.PAS} |cut -f4 |sort|uniq > $TMPDIR/$uuidPASsupported

uuidcagePASsupported=$(uuidgen)
comm -1 -2 $TMPDIR/$uuidCageSupported $TMPDIR/$uuidPASsupported |sort|uniq > $TMPDIR/$uuidcagePASsupported

zcat {input.tm} | fgrep -w -f $TMPDIR/$uuidcagePASsupported - |sortgff | gzip> {output}

		'''


rule getFlLocusStats:
	input:
		gencode="annotations/simplified/{capDesign}.gencode.collapsed.simplified_biotypes.gtf",
		clsGencode="mappings/nonAnchoredMergeReads/mergeWithPrevious+biotypes/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.all.gff.gz",
		gencodeFL="mappings/nonAnchoredMergeReads/mergeWithPrevious/gencode/{capDesign}.gencode.cage+PASsupported.gff.gz",
		clsGencodeFL="mappings/nonAnchoredMergeReads/mergeWithPrevious/cls/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.cage+PASsupported.gff.gz"
	output: temp(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.FLloci.stats.tsv")
	shell:
		'''
#PCG stats
gencodePcg=$(cat {input.gencode} | fgrep "gene_type \\"protein_coding\\";" | extractGffAttributeValue.pl gene_id|sort|uniq|wc -l)
gencodePcgFL=$(zcat {input.gencodeFL} | fgrep "gene_type \\"protein_coding\\";" | extractGffAttributeValue.pl gene_id|sort|uniq|wc -l)
clsGencodePcg=$(zcat {input.clsGencode} | fgrep "gene_type \\"protein_coding\\";" | extractGffAttributeValue.pl gene_id|sort|uniq|wc -l)
clsGencodePcgFL=$(zcat {input.clsGencodeFL} | fgrep "gene_type \\"protein_coding\\";" | extractGffAttributeValue.pl gene_id|sort|uniq|wc -l)

#lncRNA stats
uuidgencodeLnc=$(uuidgen)
uuidgencodeLncFL=$(uuidgen)
uuidclsGencodeLnc=$(uuidgen)
uuidclsGencodeLncFL=$(uuidgen)

cat {input.gencode} | fgrep "gene_type \\"lncRNA\\";" |awk '$3=="exon"' > $TMPDIR/$uuidgencodeLnc
zcat {input.gencodeFL} | fgrep "gene_type \\"lncRNA\\";" |awk '$3=="exon"'> $TMPDIR/$uuidgencodeLncFL
zcat {input.clsGencode} | fgrep "gene_type \\"lncRNA\\";" |awk '$3=="exon"'> $TMPDIR/$uuidclsGencodeLnc
zcat {input.clsGencodeFL} | fgrep "gene_type \\"lncRNA\\";" |awk '$3=="exon"'> $TMPDIR/$uuidclsGencodeLncFL

gencodeLncSpliced=$(cat $TMPDIR/$uuidgencodeLnc | extractGffAttributeValue.pl transcript_id |sort|uniq -d | fgrep -w -f - $TMPDIR/$uuidgencodeLnc | extractGffAttributeValue.pl gene_id|sort|uniq|wc -l)
gencodeLncFLSpliced=$(cat $TMPDIR/$uuidgencodeLncFL| extractGffAttributeValue.pl transcript_id |sort|uniq -d | fgrep -w -f - $TMPDIR/$uuidgencodeLncFL | extractGffAttributeValue.pl gene_id|sort|uniq|wc -l)
clsGencodeLncSpliced=$(cat $TMPDIR/$uuidclsGencodeLnc | extractGffAttributeValue.pl transcript_id |sort|uniq -d | fgrep -w -f - $TMPDIR/$uuidclsGencodeLnc | extractGffAttributeValue.pl gene_id|sort|uniq|wc -l)
clsGencodeLncFLSpliced=$(cat $TMPDIR/$uuidclsGencodeLncFL | extractGffAttributeValue.pl transcript_id |sort|uniq -d | fgrep -w -f - $TMPDIR/$uuidclsGencodeLncFL | extractGffAttributeValue.pl gene_id|sort|uniq|wc -l)

clsGencodeLncKnownSpliced=$(cat $TMPDIR/$uuidclsGencodeLnc | fgrep "gene_ref_status \\"known\\";" |extractGffAttributeValue.pl transcript_id |sort|uniq -d | fgrep -w -f - $TMPDIR/$uuidclsGencodeLnc | extractGffAttributeValue.pl gene_id|sort|uniq|wc -l)
clsGencodeLncFLKnownSpliced=$(cat $TMPDIR/$uuidclsGencodeLncFL | fgrep "gene_ref_status \\"known\\";" | extractGffAttributeValue.pl transcript_id |sort|uniq -d | fgrep -w -f - $TMPDIR/$uuidclsGencodeLncFL | extractGffAttributeValue.pl gene_id|sort|uniq|wc -l)


echo -e "{wildcards.techname}Corr{wildcards.corrLevel}\t{wildcards.capDesign}\t{wildcards.sizeFrac}\t{wildcards.barcodes}\t$gencodePcg\t$gencodePcgFL\t$clsGencodePcg\t$clsGencodePcgFL\t$gencodeLncSpliced\t$gencodeLncFLSpliced\t$clsGencodeLncSpliced\t$clsGencodeLncFLSpliced\t$clsGencodeLncKnownSpliced\t$clsGencodeLncFLKnownSpliced"  > {output}


		'''


rule aggFlLocusStats:
	input: lambda wildcards: expand(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.FLloci.stats.tsv", techname='allSeqTechs', corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNSplusMERGED, sizeFrac='allFracs', barcodes='allTissues', minReadSupport=wildcards.minReadSupport)
	output: config["STATSDATADIR"] + "all.min{minReadSupport}reads.FLloci.stats.tsv"
	shell:
		'''
echo -e "seqTech\tcorrectionLevel\tcapDesign\tsizeFrac\ttissue\tannotation_set\tbiotype\tcategory\tcount\tpercent" > {output}

cat {input} | awk -v p='{config[PROJECT_NAME]}' '{{print $1"\\t"$2"\\t"$3"\\t"$4"\\tGENCODE\\tprotein-coding\\tnon-FL\\t"$5-$6"\\t"($5-$6)/$5"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\tGENCODE\\tprotein-coding\\tFL\\t"$6"\\t"$6/$5"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\tGENCODE++"p"\\tprotein-coding\\tnon-FL\\t"$7-$8"\\t"($7-$8)/$7"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\tGENCODE++"p"\\tprotein-coding\\tFL\\t"$8"\\t"$8/$7"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\tGENCODE\\tlncRNA\\tnon-FL\\t"$9-$10"\\t"($9-$10)/$9"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\tGENCODE\\tlncRNA\\tFL\\t"$10"\\t"$10/$9"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\tGENCODE++"p"\\tlncRNA\\tnon-FL\\t"$11-$12"\\t"($11-$12)/$11"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\tGENCODE++"p"\\tlncRNA\\tFL\\t"$12"\\t"$12/$11}}' | sed 's/Corr0/\tNo/' | sed 's/Corr{lastK}/\tYes/' | sort >> {output}
		'''

rule plotFlLocusGencodeOnlyStats:
	input: config["STATSDATADIR"] + "all.min{minReadSupport}reads.FLloci.stats.tsv"
	output: config["PLOTSDIR"] + "FLloci.gencodeOnly.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.FLloci.gencodeOnly.stats.{ext}"
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


dat <- subset(dat, annotation_set=='GENCODE')

dat\$category<-factor(dat\$category, ordered=TRUE, levels=rev(c('FL', 'non-FL')))

dat\$annotbiotype <- paste(sep='', dat\$biotype, '\n(', dat\$annotation_set, ')')
plotHeight = plotHeight +1
plotWidth = plotWidth +0.5

ggplot(dat[order(dat\$category), ], aes(x=factor(annotbiotype), y=count, fill=category)) +
geom_bar(stat='identity') +
ylab('# gene loci') +
scale_y_continuous(labels=comma)+
scale_fill_manual (values=c('FL'='#C453C4', 'non-FL'='#a6a6a6'))+
theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
facet_grid( seqTech + sizeFrac ~ capDesign + tissue)+
xlab('Biotype (Annotation set)') +
guides(fill = guide_legend(title='Category'))+
geom_text(position = 'stack', size=geom_textSize, aes(x = factor(annotbiotype), y = count, label = paste(sep='',percent(round(percent, digits=2)),' / ','(',comma(count),')'), hjust = 0.5, vjust = 1))+
{GGPLOT_PUB_QUALITY}
ggsave('{output}', width=plotWidth, height=plotHeight)
" > {output}.r
cat {output}.r | R --slave

		'''


rule plotFlLocusStats:
	input: config["STATSDATADIR"] + "all.min{minReadSupport}reads.FLloci.stats.tsv"
	output: config["PLOTSDIR"] + "FLloci.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.FLloci.stats.{ext}"
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



dat\$category<-factor(dat\$category, ordered=TRUE, levels=rev(c('FL', 'non-FL')))

dat\$annotbiotype <- paste(sep='', dat\$biotype, '\n(', dat\$annotation_set, ')')

plotHeight = plotHeight +1
plotWidth = plotWidth +0.5

ggplot(dat[order(dat\$category), ], aes(x=factor(annotbiotype), y=count, fill=category)) +
geom_bar(stat='identity') +
ylab('# gene loci') +
scale_y_continuous(labels=comma)+
scale_fill_manual (values=c('FL'='#C453C4', 'non-FL'='#a6a6a6'))+
theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
facet_grid( seqTech + sizeFrac ~ capDesign + tissue)+
xlab('Biotype (Annotation set)') +
guides(fill = guide_legend(title='Category'))+
geom_text(position = 'stack', size=geom_textSize, aes(x = factor(annotbiotype), y = count, label = paste(sep='',percent(round(percent, digits=2)),' / ','(',comma(count),')'), hjust = 0.5, vjust = 1))+
{GGPLOT_PUB_QUALITY}
ggsave('{output}', width=plotWidth, height=plotHeight)
" > {output}.r
cat {output}.r | R --slave

		'''








