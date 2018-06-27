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
		transcriptStrandInfo = "mappings/" + "getIntronMotif/{techname}_{capDesign}_{barcodes}.transcripts.tsv",
		strandedReads = "mappings/" + "strandGffs/{techname}_{capDesign}_{barcodes}.stranded.gff"
	output:
		"mappings/" + "highConfidenceReads/{techname}_{capDesign}_{barcodes}.strandedHCGMs.gff.gz"
	shell:
		'''
#select read IDs with canonical GT|GC/AG and high-confidence SJs
cat {input.transcriptStrandInfo} | skipcomments | awk '$6==1 && $7==1' | cut -f1 | sort|uniq > $TMPDIR/reads.hcSJs.list
wc -l $TMPDIR/reads.hcSJs.list
set +e
fgrep -w -f $TMPDIR/reads.hcSJs.list {input.strandedReads} > $TMPDIR/gtag.gff
set -e
wc -l $TMPDIR/gtag.gff
cat {input.strandedReads} | extractGffAttributeValue.pl transcript_id | sort|uniq -u > $TMPDIR/tmp
set +e
fgrep -w -f $TMPDIR/tmp {input.strandedReads} > $TMPDIR/tmp2
set -e
cat $TMPDIR/tmp2 | fgrep -v ERCC- > $TMPDIR/monoPolyA.gff
 echo $?
cat $TMPDIR/gtag.gff $TMPDIR/monoPolyA.gff | sortgff |gzip> {output}
 echo $?
		'''


rule getHCGMintrons:
	input: "mappings/" + "highConfidenceReads/{techname}_{capDesign}_{barcodes}.strandedHCGMs.gff.gz"
	output: "mappings/" + "highConfidenceReads/introns/{techname}_{capDesign}_{barcodes}.strandedHCGMs.introns.tsv"
	shell:
		'''
zcat {input} | sort -k12,12 -k4,4n -k5,5n | awk -v fldgn=10 -v fldtr=12 -f ~/julien_utils/make_introns.awk | perl -lane '$F[11]=~s/"|;//g; $start=$F[3]-1; $end=$F[4]+1; print $F[11]."\t".$F[0]."_".$start."_".$end."_".$F[6]' | sort -k2,2 > {output}

		'''

rule getHiSeqSupportedHCGMs:
	input:
		hiSeqIntrons="mappings/hiSeqIntrons/" + "hiSeq_{capDesign}.canonicalIntrons.list",
		lrIntrons="mappings/" + "highConfidenceReads/introns/{techname}_{capDesign}_{barcodes}.strandedHCGMs.introns.tsv",
		#mergedGTF="mappings/" + "nonAnchoredMergeReads/pooled/{techname}_{capDesign}.tmerge.gff"
		hcgmGTF= "mappings/" + "highConfidenceReads/{techname}_{capDesign}_{barcodes}.strandedHCGMs.gff.gz"
	output:"mappings/" + "highConfidenceReads/HiSS/{techname}_{capDesign}_{barcodes}.HiSS.gff.gz"
	shell:
		'''
join -v1 -1 2 -2 1 {input.lrIntrons} {input.hiSeqIntrons} |awk '{{print $2"\t"$1}}' > $TMPDIR/{wildcards.techname}_{wildcards.capDesign}_{wildcards.barcodes}.introns.noHiSeq.tsv
zcat {input.hcgmGTF} > $TMPDIR/$(basename {input.hcgmGTF} .gz)
cut -f1 $TMPDIR/{wildcards.techname}_{wildcards.capDesign}_{wildcards.barcodes}.introns.noHiSeq.tsv | sort|uniq | fgrep -wv -f - $TMPDIR/$(basename {input.hcgmGTF} .gz) |sortgff |gzip > {output}

		'''

rule getHiSSStats:
	input:
		reads = "mappings/" + "mergeSizeFracBams/{techname}_{capDesign}_{barcodes}.merged.bam",
		HiSSGTF="mappings/" + "highConfidenceReads/HiSS/{techname}_{capDesign}_{barcodes}.HiSS.gff.gz"
	output: config["STATSDATADIR"] + "{techname}_{capDesign}_{barcodes}.HiSS.stats.tsv"
	shell:
		'''
bedtools bamtobed -i {input.reads} -bed12 > $TMPDIR/{wildcards.techname}_{wildcards.capDesign}_{wildcards.barcodes}.merged.bed
zcat {input.HiSSGTF} | gff2bed_full.pl - > $TMPDIR/{wildcards.techname}_{wildcards.capDesign}_{wildcards.barcodes}.HiSS.bed

mappedReadsMono=$(cat $TMPDIR/{wildcards.techname}_{wildcards.capDesign}_{wildcards.barcodes}.merged.bed | awk '$10<=1'|cut -f4 |sort|uniq|wc -l)
mappedReadsSpliced=$(cat $TMPDIR/{wildcards.techname}_{wildcards.capDesign}_{wildcards.barcodes}.merged.bed | awk '$10>1'|cut -f4 |sort|uniq|wc -l)

HiSSMono=$(cat $TMPDIR/{wildcards.techname}_{wildcards.capDesign}_{wildcards.barcodes}.HiSS.bed| awk '$10<=1'|cut -f4  | sort|uniq|wc -l)
HiSSSpliced=$(cat $TMPDIR/{wildcards.techname}_{wildcards.capDesign}_{wildcards.barcodes}.HiSS.bed| awk '$10>1'|cut -f4  | sort|uniq|wc -l)

#let totalMapped=$mappedReadsMono+$mappedReadsSpliced || true
let nonHiSSMono=$mappedReadsMono-$HiSSMono || true
let nonHiSSSPliced=$mappedReadsSpliced-$HiSSSpliced || true
echo -e "{wildcards.techname}\t{wildcards.capDesign}\t{wildcards.barcodes}\t$HiSSMono\t$HiSSSpliced\t$nonHiSSMono\t$nonHiSSSPliced" > {output}

		'''

rule aggHiSSCapDesignStats:
	input: lambda wildcards: expand(config["STATSDATADIR"] + "{techname}_{capDesign}_{barcodes}.HiSS.stats.tsv", filtered_product, techname=wildcards.techname, capDesign=wildcards.capDesign, barcodes=BARCODES)
	output: config["STATSDATADIR"] + "{techname}_{capDesign}.HiSS.stats.tsv"
	shell:
		'''
totalHiSSMono=$(cat {input} | cut -f4 | sum.sh)
totalHiSSSpliced=$(cat {input} | cut -f5 | sum.sh)
totalNonHiSSMono=$(cat {input} | cut -f6 | sum.sh)
totalNonHiSSSpliced=$(cat {input} | cut -f7 | sum.sh)

echo -e "\
{wildcards.techname}\t{wildcards.capDesign}\tHCGM-mono\t$totalHiSSMono
{wildcards.techname}\t{wildcards.capDesign}\tHCGM-spliced\t$totalHiSSSpliced
{wildcards.techname}\t{wildcards.capDesign}\tnonHCGM-mono\t$totalNonHiSSMono
{wildcards.techname}\t{wildcards.capDesign}\tnonHCGM-spliced\t$totalNonHiSSSpliced" > {output}
		'''


rule aggHiSSStats:
	input: lambda wildcards: expand(config["STATSDATADIR"] + "{techname}_{capDesign}.HiSS.stats.tsv",techname=TECHNAMES, capDesign=CAPDESIGNS)
	output: config["STATSDATADIR"] + "all.pooled.merged.HiSS.stats.tsv"
	shell:
		'''
echo -e "seqTech\tcorrectionLevel\tcapDesign\tcategory\tcount" > {output}
cat {input} | sed 's/Corr0/\tNo/' | sed 's/Corr90/\tYes/' | sort >> {output}

		'''

rule plotHiSSStats:
	input: config["STATSDATADIR"] + "all.pooled.merged.HiSS.stats.tsv"
	output:  config["PLOTSDIR"] + "all.pooled.merged.HiSS.stats.{ext}"
	shell:
		'''
echo "library(ggplot2)
library(plyr)
library(scales)
dat <- read.table('{input}', header=T, as.is=T, sep='\\t')
dat\$category<-factor(dat\$category, ordered=TRUE, levels=rev(c('HCGM-mono', 'HCGM-spliced', 'nonHCGM-mono', 'nonHCGM-spliced')))

ggplot(dat[order(dat\$category), ], aes(x=factor(correctionLevel), y=count, fill=category)) +
geom_bar(stat='identity') + scale_fill_manual(values=c('HCGM-mono' = '#9ce2bb', 'HCGM-spliced' = '#39c678', 'nonHCGM-mono' = '#fda59b', 'nonHCGM-spliced' = '#fa341e')) + facet_grid( seqTech ~ capDesign)+ ylab('# mapped reads') + xlab('Error correction') + guides(fill = guide_legend(title='Category'))+
geom_text(position = 'stack', aes(x = factor(correctionLevel), y = count, ymax=count, label = comma(count), hjust = 0.5, vjust = 1), size=2)+ scale_y_continuous(labels=scientific)+
{GGPLOT_PUB_QUALITY}
ggsave('{output}', width=7, height=3)
" > {output}.r
cat {output}.r | R --slave

		'''



rule nonAnchoredMergeReads:
	input: "mappings/" + "highConfidenceReads/HiSS/{techname}_{capDesign}_{barcodes}.HiSS.gff.gz"
	output: "mappings/" + "nonAnchoredMergeReads/{techname}_{capDesign}_{barcodes}.tmerge.gff"
	threads:8
	shell:
		'''
zcat {input} |sortgff | tmerge --cpu {threads} --tmPrefix {wildcards.techname}_{wildcards.capDesign}_{wildcards.barcodes}.NAM_ - |sortgff > {output}
		'''


rule checkNonAnchoredMerging:
	input:
		before="mappings/" + "highConfidenceReads/HiSS/{techname}_{capDesign}_{barcodes}.HiSS.gff.gz",
		after="mappings/" + "nonAnchoredMergeReads/{techname}_{capDesign}_{barcodes}.tmerge.gff"
	output: "mappings/" + "nonAnchoredMergeReads/qc/{techname}_{capDesign}_{barcodes}.tmerge.qc.txt"
	shell:
		'''
#nt coverage should be equal before/after
zcat {input.before} |awk '$3=="exon"' |sortgff| bedtools merge >$TMPDIR/before.merged.bed
cat {input.after} |awk '$3=="exon"' |sortgff| bedtools merge > $TMPDIR/after.merged.bed
beforeNts=$(cat $TMPDIR/before.merged.bed| awk '{{print $3-$2}}' |sum.sh)
afterNts=$(cat $TMPDIR/after.merged.bed| awk '{{print $3-$2}}' |sum.sh)
let diff=$beforeNts-$afterNts || true # "true" serves to avoid "let" exiting with status > 0 when its output is = 0
echo "nts before: $beforeNts" > {output}
echo "nts after: $afterNts" >> {output}
echo "diff nts: $diff" >> {output}
if [ ! $diff -eq 0 ]; then echo "ERROR: Nucleotide coverage differ before/after merging" >> {output}; \mv {input.after} {input.after}.bkp; exit 1; fi

#all transcript_ids should be present before/after
zcat {input.before} | extractGffAttributeValue.pl transcript_id | sort|uniq > $TMPDIR/before.list
cat {input.after} |extractGffAttributeValue.pl contains |sed 's/,/\\n/g' | sort|uniq > $TMPDIR/after.list
diffLines=$(diff -q $TMPDIR/before.list $TMPDIR/after.list |wc -l)

echo "difflines: $diffLines"  >> {output}
if [ ! $diffLines -eq 0 ]; then echo "ERROR: List of transcripts differ before/after" >> {output}; \mv {input.after} {input.after}.bkp; exit 1; fi

#intron list should be the same before/after:
zcat {input.before} |sort -k12,12 -k4,4n -k5,5n | awk -v fldgn=10 -v fldtr=12 -f ~/julien_utils/make_introns.awk |awk '{{print $1"_"$4"_"$5"_"$7}}' |sort|uniq > $TMPDIR/before.introns.list
cat {input.after} |sort -k12,12 -k4,4n -k5,5n | awk -v fldgn=10 -v fldtr=12 -f ~/julien_utils/make_introns.awk |awk '{{print $1"_"$4"_"$5"_"$7}}' |sort|uniq > $TMPDIR/after.introns.list
diffIntronLines=$(diff -q $TMPDIR/before.introns.list $TMPDIR/after.introns.list |wc -l)

echo "diffintronlines: $diffIntronLines" >> {output}
if [ ! $diffIntronLines -eq 0 ]; then echo "ERROR: List of introns differ before/after" >> {output}; \mv {input.after} {input.after}.bkp; exit 1; fi

echo XXdoneXX  >> {output}

		'''

rule poolNonAnchoredMergeReads:
	input: lambda wildcards: expand("mappings/" + "nonAnchoredMergeReads/{techname}_{capDesign}_{barcodes}.tmerge.gff", filtered_product, techname=wildcards.techname, capDesign=wildcards.capDesign, barcodes=BARCODES)
	output: "mappings/" + "nonAnchoredMergeReads/pooled/{techname}_{capDesign}_pooled.tmerge.gff"
	threads: 8
	shell:
		'''
cat {input} |sortgff | tmerge --cpu {threads} --tmPrefix {wildcards.techname}_{wildcards.capDesign}_pooled.NAM_ - |sortgff > {output}
		'''

# rule poolNonAnchoredMergeReadsMergeByChr:
# 	input: lambda wildcards: expand("mappings/" + "nonAnchoredMergeReads/chr/pooled/{techname}_{capDesign}.tmerge.{chrom}.gff", filtered_product, techname=wildcards.techname, capDesign=wildcards.capDesign, chrom=GENOMECHROMS)
# 	output: "mappings/" + "nonAnchoredMergeReads/pooled/{techname}_{capDesign}.tmerge.gff"
# 	shell:
# 		'''
# cat {input} |sortgff > {output}
# 		'''

rule checkPooledNonAnchoredMerging:
	input:
		before=lambda wildcards: expand("mappings/" + "nonAnchoredMergeReads/{techname}_{capDesign}_{barcodes}.tmerge.gff", filtered_product, techname=wildcards.techname, capDesign=wildcards.capDesign, barcodes=BARCODES),
		after="mappings/" + "nonAnchoredMergeReads/pooled/{techname}_{capDesign}_pooled.tmerge.gff"
	output: "mappings/" + "nonAnchoredMergeReads/pooled/qc/{techname}_{capDesign}.tmerge.qc.txt"
	shell:
		'''
#nt coverage should be equal before/after
cat {input.before} |awk '$3=="exon"' |sortgff| bedtools merge >$TMPDIR/before.merged.bed
cat {input.after} |awk '$3=="exon"' |sortgff| bedtools merge > $TMPDIR/after.merged.bed
beforeNts=$(cat $TMPDIR/before.merged.bed| awk '{{print $3-$2}}' |sum.sh)
afterNts=$(cat $TMPDIR/after.merged.bed| awk '{{print $3-$2}}' |sum.sh)
let diff=$beforeNts-$afterNts || true # "true" serves to avoid "let" exiting with status > 0 when its output is = 0
echo "nts before: $beforeNts" > {output}
echo "nts after: $afterNts" >> {output}
echo "diff nts: $diff" >> {output}

if [ ! $diff -eq 0 ]; then echo "ERROR: Nucleotide coverage differ before/after merging" >> {output};\mv {input.after} {input.after}.bkp;  exit 1; fi

#all transcript_ids should be present before/after
cat {input.before} | extractGffAttributeValue.pl transcript_id | sort|uniq > $TMPDIR/before.list
cat {input.after} |extractGffAttributeValue.pl contains |sed 's/,/\\n/g' | sort|uniq > $TMPDIR/after.list
diffLines=$(diff -q $TMPDIR/before.list $TMPDIR/after.list |wc -l)

echo "difflines: $diffLines"  >> {output}
if [ ! $diffLines -eq 0 ]; then echo "ERROR: List of transcripts differ before/after" >> {output}; \mv {input.after} {input.after}.bkp; exit 1; fi

#intron list should be the same before/after:
cat {input.before} |sort -k12,12 -k4,4n -k5,5n | awk -v fldgn=10 -v fldtr=12 -f ~/julien_utils/make_introns.awk |awk '{{print $1"_"$4"_"$5"_"$7}}' |sort|uniq > $TMPDIR/before.introns.list
cat {input.after} |sort -k12,12 -k4,4n -k5,5n | awk -v fldgn=10 -v fldtr=12 -f ~/julien_utils/make_introns.awk |awk '{{print $1"_"$4"_"$5"_"$7}}' |sort|uniq > $TMPDIR/after.introns.list
diffIntronLines=$(diff -q $TMPDIR/before.introns.list $TMPDIR/after.introns.list |wc -l)

echo "diffintronlines: $diffIntronLines" >> {output}
if [ ! $diffIntronLines -eq 0 ]; then echo "ERROR: List of introns differ before/after" >> {output}; \mv {input.after} {input.after}.bkp; exit 1; fi

echo "XXdoneXX" >> {output}

		'''

rule getPooledMergingStats:
	input:
		hcgms = lambda wildcards: expand("mappings/" + "highConfidenceReads/HiSS/{techname}_{capDesign}_{barcodes}.HiSS.gff.gz", filtered_product, techname=wildcards.techname, capDesign=wildcards.capDesign, barcodes=BARCODES),
		pooledMerged = "mappings/" + "nonAnchoredMergeReads/pooled/{techname}_{capDesign}_pooled.tmerge.gff"
	output: config["STATSDATADIR"] + "{techname}_{capDesign}.pooled.merged.stats.tsv"
	shell:
		'''
hcgms=$(zcat {input.hcgms} | extractGffAttributeValue.pl transcript_id | sort|uniq|wc -l)
merged=$(cat {input.pooledMerged} | extractGffAttributeValue.pl transcript_id | sort|uniq|wc -l)
echo -e "{wildcards.techname}\t{wildcards.capDesign}\t$hcgms\t$merged" | awk '{{print $0"\t"$4/$3}}' > {output}

		'''

rule aggPooledMergingStats:
	input: lambda wildcards: expand(config["STATSDATADIR"] + "{techname}_{capDesign}.pooled.merged.stats.tsv", techname=TECHNAMES, capDesign=CAPDESIGNS)
	output: config["STATSDATADIR"] + "all.pooled.merged.stats.tsv"
	shell:
		'''
echo -e "seqTech\tcorrectionLevel\tcapDesign\tcategory\tcount" > {output}
cat {input} | awk '{{print $1"\\t"$2"\\tHCGMreads\\t"$3"\\n"$1"\\t"$2"\\tmergedTMs\\t"$4}}' | sed 's/Corr0/\tNo/' | sed 's/Corr90/\tYes/' | sort >> {output}

		'''
rule plotPooledMergingStats:
	input:  config["STATSDATADIR"] + "all.pooled.merged.stats.tsv"
	output: config["PLOTSDIR"] + "all.pooled.merged.stats.{ext}"
	shell:
		'''
echo "library(ggplot2)
library(plyr)
library(scales)
dat <- read.table('{input}', header=T, as.is=T, sep='\\t')
ggplot(data=dat, aes(x=factor(correctionLevel), y=count, fill=category)) +
geom_bar(stat='identity', position=position_dodge()) + geom_text(position = position_dodge(width = 0.9), aes(x = factor(correctionLevel), y = 1, ymax=count, label = comma(count), hjust = 0, vjust = 0.5), angle=90, size=4) + scale_fill_manual(values=c('HCGMreads' = '#d98cb3', 'mergedTMs' = '#cc9966')) + facet_grid( seqTech ~ capDesign)+ ylab('# objects') + xlab('Error correction') + guides(fill = guide_legend(title='Category'))+ scale_y_continuous(labels=comma)+
{GGPLOT_PUB_QUALITY}
ggsave('{output}', width=7, height=3)
" > {output}.r
cat {output}.r | R --slave

		'''




# rule getHiSSSplicedLengthStats:
# 	input: "mappings/" + "nonAnchoredMergeReads/pooled/HiSS/{techname}_{capDesign}.tmerge.HiSS.gff"
# 	output: config["STATSDATADIR"] + "{techname}_{capDesign}.pooled.merged.HiSS.splicedLength.stats.tsv"
# 	shell:
# 		'''
# cat {input} | gff2bed_full.pl -| perl -slane '@blocksizes=split(",",$F[10]); $splicedLength=0; foreach $exonLength (@blocksizes){{$splicedLength+=$exonLength}}; print "$F[3]\t$splicedLength"' | awk -v t={wildcards.techname} -v c={wildcards.capDesign} '{{print t\"\\t\"c"\\t\"$0}}'> {output}
# 		'''

# rule aggHiSSSplicedLengthStats:
# 	input: lambda wildcards: expand(config["STATSDATADIR"] + "{techname}_{capDesign}.pooled.merged.HiSS.splicedLength.stats.tsv",techname=TECHNAMES, capDesign=CAPDESIGNS)
# 	output: config["STATSDATADIR"] + "all.pooled.merged.HiSS.splicedLength.stats.tsv"
# 	shell:
# 		'''
# cat {input}  | sed 's/Corr/\t/' > {output}
# 		'''

# rule plotHiSSSplicedLengthStats:
# 	input: config["STATSDATADIR"] + "all.pooled.merged.HiSS.splicedLength.stats.tsv"
# 	output: config["PLOTSDIR"] + "all.pooled.merged.HiSS.splicedLength.stats.{ext}"
# 	shell:
# 		'''
# echo "library(ggplot2)
# library(plyr)
# library(scales)
# dat <- read.table('{input}', header=F, as.is=T, sep='\\t')
# colnames(dat)<-c('seqTech','correctionLevel','capDesign', 'TMid', 'length')
# ggplot(dat, aes(x=factor(correctionLevel), y=length)) +
# geom_boxplot() +
# facet_grid( seqTech ~ capDesign)+
# coord_cartesian(ylim=c(500, 3000)) +
# xlab('Error correction') +
# ylab('Spliced length (nts)') +
# scale_y_continuous(labels=comma)+
# {GGPLOT_PUB_QUALITY}
# ggsave('{output}', width=10, height=3)
# " > {output}.r
# cat {output}.r | R --slave

# 		'''