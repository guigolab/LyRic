rule integratePolyaAndSjInfo:
	input:
		polyA = "mappings/removePolyAERCCs/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.polyAsitesNoErcc.tmp.bed",
		SJs = "mappings/getIntronMotif/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.transcripts.tsv"
	output:
		strandInfo="mappings/integratePolyaAndSjInfo/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.polyA+SJ.strandInfo.tsv",
		wrongPolyAs="mappings/" + "removePolyAERCCs/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.wrongPolyAs.list"
	shell:
		'''
cat {input.SJs} | skipcomments | cut -f 1,2 | awk '$2!="."' | sort> $TMPDIR/reads.SJ.strandInfo.tsv
cat {input.polyA} | cut -f4,6 | awk '$2!="."'| sort > $TMPDIR/reads.polyA.strandInfo.tsv
join -a1 -a2 -e '.' -o '0,1.2,2.2' $TMPDIR/reads.SJ.strandInfo.tsv  $TMPDIR/reads.polyA.strandInfo.tsv > $TMPDIR/reads.SJ.polyA.strandInfo.tsv

#make list of reads with wrongly called polyA sites (i.e. their strand is different from the one inferred using SJs):
cat $TMPDIR/reads.SJ.polyA.strandInfo.tsv | perl -slane 'if($F[2] ne "." && $F[1] ne "." && $F[1] ne $F[2]){{print join ("\\t", $F[0])}}' |sort|uniq > {output.wrongPolyAs}
#get strand info, prioritizing SJ inference
cat $TMPDIR/reads.SJ.polyA.strandInfo.tsv | perl -slane 'if($F[1] ne "."){{print "$F[0]\\t$F[1]"}} else{{print "$F[0]\t$F[2]"}}' |sort|uniq > {output.strandInfo}
		'''

rule removeWrongPolyAs:
	input:
		polyA="mappings/removePolyAERCCs/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.polyAsitesNoErcc.tmp.bed",
		wrongPolyAs="mappings/" + "removePolyAERCCs/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.wrongPolyAs.list"
	output: "mappings/" + "removePolyAERCCs/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.polyAsitesNoErcc.bed"
	shell:
		'''
cat {input.polyA} | fgrep -v -w -f {input.wrongPolyAs} > {output}
		'''

rule strandGffs:
	input:
		gff = "mappings/" + "readBedToGff/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.gff.gz",
		strandInfo = "mappings/" + "integratePolyaAndSjInfo/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.polyA+SJ.strandInfo.tsv"
	output: "mappings/" + "strandGffs/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.stranded.gff.gz"
	shell:
		'''
zcat {input.gff} > $TMPDIR/in.gff
get_right_transcript_strand.pl $TMPDIR/in.gff {input.strandInfo} | fgrep -v ERCC- | sortgff |gzip> {output}

		'''

rule highConfidenceReads:
	input:
		transcriptStrandInfo = "mappings/" + "getIntronMotif/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.transcripts.tsv",
		strandedReads = "mappings/" + "strandGffs/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.stranded.gff.gz"
	output:
		"mappings/" + "highConfidenceReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.strandedHCGMs.gff.gz"
	shell:
		'''
#select read IDs with canonical GT|GC/AG and high-confidence SJs
cat {input.transcriptStrandInfo} | skipcomments | awk '$6==1 && $7==1' | cut -f1 | sort|uniq > $TMPDIR/reads.hcSJs.list
wc -l $TMPDIR/reads.hcSJs.list
zcat {input.strandedReads} > $TMPDIR/str.gff
set +e
fgrep -w -f $TMPDIR/reads.hcSJs.list $TMPDIR/str.gff > $TMPDIR/gtag.gff
set -e
wc -l $TMPDIR/gtag.gff
cat $TMPDIR/str.gff | extractGffAttributeValue.pl transcript_id | sort|uniq -u > $TMPDIR/tmp
set +e
fgrep -w -f $TMPDIR/tmp $TMPDIR/str.gff > $TMPDIR/tmp2
set -e
cat $TMPDIR/tmp2 | fgrep -v ERCC- > $TMPDIR/monoPolyA.gff
 echo $?
cat $TMPDIR/gtag.gff $TMPDIR/monoPolyA.gff | sortgff |gzip> {output}
 echo $?
		'''


rule getHCGMintrons:
	input: "mappings/" + "highConfidenceReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.strandedHCGMs.gff.gz"
	output: "mappings/" + "highConfidenceReads/introns/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.strandedHCGMs.introns.tsv"
	shell:
		'''
zcat {input} | sort -k12,12 -k4,4n -k5,5n | awk -v fldgn=10 -v fldtr=12 -f ~/julien_utils/make_introns.awk | perl -lane '$F[11]=~s/"|;//g; $start=$F[3]-1; $end=$F[4]+1; print $F[11]."\t".$F[0]."_".$start."_".$end."_".$F[6]' | sort -k2,2 > {output}

		'''

rule getHiSeqSupportedHCGMs:
	input:
		hiSeqIntrons="mappings/hiSeqIntrons/" + "hiSeq_{capDesign}.canonicalIntrons.list" if config["USE_MATCHED_ILLUMINA"] and config["DEMULTIPLEX"] else "mappings/hiSeqIntrons/hiSeq_{techname}_{capDesign}.{barcodes}.canonicalIntrons.list" if config["USE_MATCHED_ILLUMINA"] and not config["DEMULTIPLEX"] else config["SUPPORT_INTRONS_DB"],
		lrIntrons="mappings/" + "highConfidenceReads/introns/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.strandedHCGMs.introns.tsv",
		#mergedGTF="mappings/" + "nonAnchoredMergeReads/pooled/{techname}Corr{corrLevel}_{capDesign}.tmerge.gff"
		hcgmGTF= "mappings/" + "highConfidenceReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.strandedHCGMs.gff.gz"
	output:"mappings/" + "highConfidenceReads/HiSS/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.gff.gz"
	shell:
		'''
join -v1 -1 2 -2 1 {input.lrIntrons} {input.hiSeqIntrons} |awk '{{print $2"\t"$1}}' > $TMPDIR/{wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.barcodes}.introns.noHiSeq.tsv
zcat {input.hcgmGTF} > $TMPDIR/$(basename {input.hcgmGTF} .gz)
cut -f1 $TMPDIR/{wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.barcodes}.introns.noHiSeq.tsv | sort|uniq | fgrep -wv -f - $TMPDIR/$(basename {input.hcgmGTF} .gz) |sortgff |gzip > {output}

		'''

rule getHiSSStats:
	input:
		reads = "mappings/" + "readMapping/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bam",
		HiSSGTF="mappings/" + "highConfidenceReads/HiSS/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.gff.gz"
	output: temp(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.stats.tsv")
	shell:
		'''
bedtools bamtobed -i {input.reads} -bed12 > $TMPDIR/{wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.barcodes}.merged.bed
zcat {input.HiSSGTF} | gff2bed_full.pl - > $TMPDIR/{wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.barcodes}.HiSS.bed

mappedReadsMono=$(cat $TMPDIR/{wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.barcodes}.merged.bed | awk '$10<=1'|cut -f4 |sort|uniq|wc -l)
mappedReadsSpliced=$(cat $TMPDIR/{wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.barcodes}.merged.bed | awk '$10>1'|cut -f4 |sort|uniq|wc -l)

HiSSMono=$(cat $TMPDIR/{wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.barcodes}.HiSS.bed| awk '$10<=1'|cut -f4  | sort|uniq|wc -l)
HiSSSpliced=$(cat $TMPDIR/{wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.barcodes}.HiSS.bed| awk '$10>1'|cut -f4  | sort|uniq|wc -l)

#let totalMapped=$mappedReadsMono+$mappedReadsSpliced || true
let nonHiSSMono=$mappedReadsMono-$HiSSMono || true
let nonHiSSSPliced=$mappedReadsSpliced-$HiSSSpliced || true
echo -e "{wildcards.techname}Corr{wildcards.corrLevel}\t{wildcards.capDesign}\t{wildcards.sizeFrac}\t{wildcards.barcodes}\t$HiSSMono\t$HiSSSpliced\t$nonHiSSMono\t$nonHiSSSPliced" |sed 's/{wildcards.capDesign}_//' > {output}

		'''


rule aggHiSSStats:
	input: expand(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.stats.tsv",filtered_product_merge, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACSpluSMERGED, barcodes=BARCODESpluSMERGED)
	output: config["STATSDATADIR"] + "all.HiSS.stats.tsv"

	shell:
		'''
echo -e "seqTech\tcorrectionLevel\tcapDesign\tsizeFrac\ttissue\tcategory\tcount" > {output}
cat {input} | awk '{{print $1"\\t"$2"\\t"$3"\\t"$4"\\tHCGM-mono\\t"$5"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\tHCGM-spliced\\t"$6"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\tnonHCGM-mono\\t"$7"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\tnonHCGM-spliced\\t"$8"\\n"}}'| sed 's/Corr0/\tNo/' | sed 's/Corr{lastK}/\tYes/' | sort >> {output}

		'''

rule plotAllHiSSStats:
	input: config["STATSDATADIR"] + "all.HiSS.stats.tsv"
	output:  config["PLOTSDIR"] + "{capDesign}_{byFrac}_{byTissue}.HiSS.stats.{ext}"
	params:
		filterDat=lambda wildcards: merge_figures_params(wildcards.capDesign, wildcards.byFrac, wildcards.byTissue)
	shell:
		'''
echo "library(ggplot2)
library(plyr)
library(scales)
dat <- read.table('{input}', header=T, as.is=T, sep='\\t')
dat\$category<-factor(dat\$category, ordered=TRUE, levels=rev(c('HCGM-mono', 'HCGM-spliced', 'nonHCGM-mono', 'nonHCGM-spliced')))
{params.filterDat}

ggplot(dat[order(dat\$category), ], aes(x=factor(correctionLevel), y=count, fill=category)) +
geom_bar(stat='identity') + scale_fill_manual(values=c('HCGM-mono' = '#9ce2bb', 'HCGM-spliced' = '#39c678', 'nonHCGM-mono' = '#fda59b', 'nonHCGM-spliced' = '#fa341e')) + facet_grid( seqTech + sizeFrac ~ capDesign + tissue)+ ylab('# mapped reads') + xlab('Error correction') + guides(fill = guide_legend(title='Category'))+
geom_text(position = 'stack', size=geom_textSize, aes(x = factor(correctionLevel), y = count, ymax=count, label = comma(count), hjust = 0.5, vjust = 1))+ scale_y_continuous(labels=scientific)+
{GGPLOT_PUB_QUALITY}
ggsave('{output}', width=plotWidth, height=plotHeight)
" > {output}.r
cat {output}.r | R --slave

		'''



rule nonAnchoredMergeReads:
	input: "mappings/" + "highConfidenceReads/HiSS/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.gff.gz"
	output: "mappings/" + "nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.gff"
	threads:8
	wildcard_constraints:
		barcodes='(?!allTissues).+',
		sizeFrac='[0-9-+\.]+'
	shell:
		'''
zcat {input} |sortgff | tmerge --cpu {threads} --tmPrefix {wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.barcodes}.NAM_ - |sortgff > {output}
		'''

rule mergeTissuesNonAnchoredMergeReads:
	input: lambda wildcards: expand("mappings/" + "nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.gff", filtered_product, techname=wildcards.techname, corrLevel=wildcards.corrLevel, capDesign=wildcards.capDesign, sizeFrac=wildcards.sizeFrac, barcodes=BARCODES)
	output: "mappings/" + "nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.gff"
	threads:8
	wildcard_constraints:
		barcodes='allTissues',
		sizeFrac='[0-9-+\.]+'
	shell:
		'''
cat {input} |sortgff | tmerge --cpu {threads} --tmPrefix {wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.barcodes}.NAM_ - |sortgff > {output}
		'''

rule mergeFracsNonAnchoredMergeReads:
	input: lambda wildcards: expand("mappings/" + "nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.gff", filtered_product, techname=wildcards.techname, corrLevel=wildcards.corrLevel, capDesign=wildcards.capDesign, sizeFrac=SIZEFRACS, barcodes=wildcards.barcodes)
	output: "mappings/" + "nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.gff"
	threads:8
	wildcard_constraints:
		sizeFrac='allFracs',
		barcodes='(?!allTissues).+'
	shell:
		'''
cat {input} |sortgff | tmerge --cpu {threads} --tmPrefix {wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.barcodes}.NAM_ - |sortgff > {output}
		'''


rule mergeFracsAndTissuesNonAnchoredMergeReads:
	input:
		lambda wildcards: expand("mappings/" + "nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.gff", filtered_product, techname=wildcards.techname, corrLevel=wildcards.corrLevel, capDesign=wildcards.capDesign, sizeFrac=SIZEFRACS, barcodes=BARCODES),

	output: "mappings/" + "nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_allFracs_allTissues.tmerge.gff"
	threads:8

	shell:
		'''
cat {input} |sortgff | tmerge --cpu {threads} --tmPrefix {wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}_allFracs_allTissues.NAM_ - |sortgff > {output}
		'''



rule checkNonAnchoredMerging:
	input:
		before="mappings/" + "highConfidenceReads/HiSS/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.gff.gz",
		after="mappings/" + "nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.gff"
	output: "mappings/" + "nonAnchoredMergeReads/qc/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.qc.txt"
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


rule nonAnchoredMergeReadsToBed:
	input: "mappings/" + "nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.gff"
	output: "mappings/" + "nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.bed"
	shell:
		'''
cat {input} | gff2bed_full.pl - > {output}
		'''


rule getMergingStats:
	input:
		hcgms = "mappings/" + "highConfidenceReads/HiSS/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.gff.gz",
		pooledMerged = "mappings/" + "nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.gff"
	output: temp(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.merged.stats.tsv")
	shell:
		'''
hcgms=$(zcat {input.hcgms} | extractGffAttributeValue.pl transcript_id | sort|uniq|wc -l)
merged=$(cat {input.pooledMerged} | extractGffAttributeValue.pl transcript_id | sort|uniq|wc -l)
echo -e "{wildcards.techname}Corr{wildcards.corrLevel}\t{wildcards.capDesign}\t{wildcards.sizeFrac}\t{wildcards.barcodes}\t$hcgms\t$merged" | awk '{{print $0"\t"$6/$5}}' |sed 's/{wildcards.capDesign}_//' > {output}

		'''


rule aggMergingStats:
	input: expand(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.merged.stats.tsv",filtered_product_merge, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACSpluSMERGED, barcodes=BARCODESpluSMERGED)

	output: config["STATSDATADIR"] + "all.merged.stats.tsv"
	shell:
		'''
echo -e "seqTech\tcorrectionLevel\tcapDesign\tsizeFrac\ttissue\tcategory\tcount" > {output}
cat {input} | awk '{{print $1"\\t"$2"\\t"$3"\\t"$4"\\tHCGMreads\\t"$5"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\tmergedTMs\\t"$6}}' | sed 's/Corr0/\tNo/' | sed 's/Corr{lastK}/\tYes/' | sort >> {output}

		'''
rule plotMergingStats:
	input:  config["STATSDATADIR"] + "all.merged.stats.tsv"
	output: config["PLOTSDIR"] + "{capDesign}_{byFrac}_{byTissue}.merged.stats.{ext}"
	params:
		filterDat=lambda wildcards: merge_figures_params(wildcards.capDesign, wildcards.byFrac, wildcards.byTissue)
	shell:
		'''
echo "library(ggplot2)
library(plyr)
library(scales)
dat <- read.table('{input}', header=T, as.is=T, sep='\\t')
{params.filterDat}
ggplot(data=dat, aes(x=factor(correctionLevel), y=count, fill=category)) +
geom_bar(stat='identity', position=position_dodge()) + geom_text(position = position_dodge(width = 0.9), size=geom_textSize, aes(x = factor(correctionLevel), y = 1, ymax=count, label = comma(count), hjust = 0, vjust = 0.5), angle=90) + scale_fill_manual(values=c('HCGMreads' = '#d98cb3', 'mergedTMs' = '#cc9966')) + facet_grid( seqTech + sizeFrac ~ capDesign + tissue)+ ylab('# objects') + xlab('Error correction') + guides(fill = guide_legend(title='Category'))+ scale_y_continuous(labels=comma)+
{GGPLOT_PUB_QUALITY}
ggsave('{output}', width=plotWidth, height=plotHeight)
" > {output}.r
cat {output}.r | R --slave

		'''

rule getTmLengthStats:
	input:
		gencode=lambda wildcards: CAPDESIGNTOANNOTGTF[wildcards.capDesign],
		tms="mappings/" + "nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.bed",
		flTms="mappings/nonAnchoredMergeReads/cage+polyASupported/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.cage+polyASupported.bed"
	output:temp(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.splicedLength.stats.tsv")
	shell:
		'''
cat {input.gencode} | awk '$3=="exon"' | fgrep "transcript_type \\"protein_coding\\";" |gff2bed_full.pl - | bed12ToTranscriptLength.pl - | awk -v t={wildcards.techname}Corr{wildcards.corrLevel} -v c={wildcards.capDesign} -v si={wildcards.sizeFrac} -v b={wildcards.barcodes} '{{print t"\t"c"\t"si"\t"b"\tGENCODE_protein_coding\t"$2}}'| sed 's/Corr0/\tNo/' | sed 's/Corr{lastK}/\tYes/' |sed 's/{wildcards.capDesign}_//' > $TMPDIR/gencode.pcg.tsv

cat {input.tms} | bed12ToTranscriptLength.pl - | awk -v t={wildcards.techname}Corr{wildcards.corrLevel} -v c={wildcards.capDesign} -v si={wildcards.sizeFrac} -v b={wildcards.barcodes} '{{print t"\t"c"\t"si"\t"b"\tCLS_TMs\t"$2}}'| sed 's/Corr0/\tNo/' | sed 's/Corr{lastK}/\tYes/' |sed 's/{wildcards.capDesign}_//' > $TMPDIR/tms.tsv

cat {input.flTms} | bed12ToTranscriptLength.pl - | awk -v t={wildcards.techname}Corr{wildcards.corrLevel} -v c={wildcards.capDesign} -v si={wildcards.sizeFrac} -v b={wildcards.barcodes} '{{print t"\t"c"\t"si"\t"b"\tCLS_FL_TMs\t"$2}}'| sed 's/Corr0/\tNo/' | sed 's/Corr{lastK}/\tYes/' |sed 's/{wildcards.capDesign}_//' > $TMPDIR/flTms.tsv

cat $TMPDIR/gencode.pcg.tsv $TMPDIR/tms.tsv $TMPDIR/flTms.tsv | sed 's/Corr0/\tNo/' | sed 's/Corr{lastK}/\tYes/' |sed 's/{wildcards.capDesign}_//' > {output}
		'''


rule aggTmLengthStats:
	input: expand(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.splicedLength.stats.tsv",filtered_product_merge, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACSpluSMERGED, barcodes=BARCODESpluSMERGED)
	#input: lambda wildcards: expand(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}.{spliceType}.splice.sites.stats.tsv", techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, spliceType=SPLICE_SITE_TYPES)
	output: config["STATSDATADIR"] + "all.splicedLength.stats.tsv"
	shell:
		'''
echo -e "seqTech\tcorrectionLevel\tcapDesign\tsizeFrac\ttissue\tcategory\tsplicedLength" > {output}
cat {input} |sort >> {output}
		'''


rule plotTmLengthStats:
	input: config["STATSDATADIR"] + "all.splicedLength.stats.tsv"
	output: config["PLOTSDIR"] + "{capDesign}_{byFrac}_{byTissue}.splicedLength.stats.{ext}"
	params:
		filterDat=lambda wildcards: merge_figures_params(wildcards.capDesign, wildcards.byFrac, wildcards.byTissue)
	shell:
		'''
echo "
library(data.table)
library(ggplot2)
library(scales)
library(plyr)
palette <- c('GENCODE_protein_coding' = '#009900', 'CLS_TMs' = '#cc9966', 'CLS_FL_TMs' = '#cc00cc')
dat<-fread('{input}', header=T, sep='\\t')
{params.filterDat}

fun_length <- function(x){{
return(data.frame(y=-8.5,label= paste0('N=', comma(length(x)))))
}}
fun_median <- function(x){{
return(data.frame(y=-8.5,label= paste0('Median=', comma(median(x)))))
}}
ggplot(dat, aes(x=factor(correctionLevel), y=splicedLength, color=category)) +
geom_boxplot(position=position_dodge(0.9), outlier.shape=NA) +
coord_cartesian(ylim=c(100, 3500)) +
scale_y_continuous(labels=comma)+
scale_color_manual(values=palette, name='Category', labels = c(GENCODE_protein_coding = 'GENCODE\nprotein-coding', CLS_TMs='CLS TMs', CLS_FL_TMs='CLS FL TMs')) +
facet_grid( seqTech + sizeFrac ~ capDesign + tissue)+
stat_summary(aes(x=factor(correctionLevel)), position=position_dodge(0.9), fun.data = fun_length, geom = 'text', vjust = +1, hjust=0, show.legend=FALSE, size=geom_textSize) +
#stat_summary(aes(x=factor(correctionLevel)), position=position_dodge(0.9), fun.data = fun_median, geom = 'text', vjust = 0, hjust=0, show.legend=FALSE, color='black', size=geom_textSize) +
stat_summary(aes(x=factor(correctionLevel)), position=position_dodge(0.9), fun.data = fun_median, geom = 'text', vjust = 0, hjust=0, show.legend=FALSE, size=geom_textSize) +
#geom_hline(aes(yintercept=0), linetype='dashed', alpha=0.7)+
ylab('Spliced length') + xlab('Error correction') + coord_flip(ylim=c(100, 3500)) +
{GGPLOT_PUB_QUALITY}
ggsave('{output}', width=plotWidth, height=plotHeight)

" > {output}.r
cat {output}.r | R --slave

		'''
