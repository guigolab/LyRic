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
cat {input.SJs} | skipcomments | cut -f 1,2 | awk '$2!="."' | sort> $TMPDIR/$uuid.reads.SJ.strandInfo.tsv
cat {input.polyA} | cut -f4,6 | awk '$2!="."'| sort > $TMPDIR/$uuid.reads.polyA.strandInfo.tsv
join -a1 -a2 -e '.' -o '0,1.2,2.2' $TMPDIR/$uuid.reads.SJ.strandInfo.tsv  $TMPDIR/$uuid.reads.polyA.strandInfo.tsv > $TMPDIR/$uuid.reads.SJ.polyA.strandInfo.tsv

#make list of reads with wrongly called polyA sites (i.e. their strand is different from the one inferred using SJs):
cat $TMPDIR/$uuid.reads.SJ.polyA.strandInfo.tsv | perl -slane 'if($F[2] ne "." && $F[1] ne "." && $F[1] ne $F[2]){{print join ("\\t", $F[0])}}' |sort|uniq > {output.wrongPolyAs}
#get strand info, prioritizing SJ inference
cat $TMPDIR/$uuid.reads.SJ.polyA.strandInfo.tsv | perl -slane 'if($F[1] ne "."){{print "$F[0]\\t$F[1]"}} else{{print "$F[0]\\t$F[2]"}}' |sort|uniq > {output.strandInfo}
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
zcat {input.gff} > $TMPDIR/$uuid.in.gff
get_right_transcript_strand.pl $TMPDIR/$uuid.in.gff {input.strandInfo} | fgrep -v ERCC- | sortgff |gzip> {output}

		'''

rule highConfidenceReads:
	input:
		transcriptStrandInfo = "mappings/getIntronMotif/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.transcripts.tsv",
		strandedReads = "mappings/strandGffs/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.stranded.gff.gz",
		bam = "mappings/readMapping/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bam"
	output:
		"mappings/highConfidenceReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.strandedHCGMs.gff.gz"
	shell:
		'''
#select read IDs with canonical GT|GC/AG and high-confidence SJs
uuid=$(uuidgen)
set +e
cat {input.transcriptStrandInfo} | skipcomments | awk '$6==1 && $7==1' | cut -f1 | sort|uniq > $TMPDIR/$uuid.reads.hcSJs.list.tmp
set -e
wc -l $TMPDIR/$uuid.reads.hcSJs.list.tmp
zcat {input.strandedReads} > $TMPDIR/$uuid.str.gff

##high-quality Phred SJs
samtools view -F 256 -F4 -F 2048 {input.bam} | samHQintrons.pl --minQual {config[minPhredQualAroundSJ]} - |cut -f1|sort|uniq > $TMPDIR/$uuid.HQintrons.reads.list
set +e
fgrep -w -f $TMPDIR/$uuid.HQintrons.reads.list $TMPDIR/$uuid.reads.hcSJs.list.tmp > $TMPDIR/$uuid.reads.hcSJs.list
set -e

set +e
fgrep -w -f $TMPDIR/$uuid.reads.hcSJs.list $TMPDIR/$uuid.str.gff > $TMPDIR/$uuid.gtag.gff
set -e
wc -l $TMPDIR/$uuid.gtag.gff
cat $TMPDIR/$uuid.str.gff | extractGffAttributeValue.pl transcript_id | sort|uniq -u > $TMPDIR/$uuid.tmp
set +e
fgrep -w -f $TMPDIR/$uuid.tmp $TMPDIR/$uuid.str.gff > $TMPDIR/$uuid.tmp2
set -e
cat $TMPDIR/$uuid.tmp2  > $TMPDIR/$uuid.monoPolyA.gff
 echo $?
cat $TMPDIR/$uuid.gtag.gff $TMPDIR/$uuid.monoPolyA.gff | sortgff |gzip> {output}
 echo $?
		'''


rule getHCGMintrons:
	input: "mappings/highConfidenceReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.strandedHCGMs.gff.gz"
	output: temp("mappings/highConfidenceReads/introns/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.strandedHCGMs.introns.tsv")
	shell:
		'''
zcat {input} | makeIntrons.pl -| perl -lane '$F[11]=~s/"|;//g; $start=$F[3]-1; $end=$F[4]+1; print $F[11]."\t".$F[0]."_".$start."_".$end."_".$F[6]' | sort -k2,2 > {output}

		'''

rule getHiSeqSupportedHCGMs:
	input:
		hiSeqIntrons="mappings/hiSeqIntrons/hiSeq_{capDesign}.canonicalIntrons.list" if config["USE_MATCHED_ILLUMINA"] and config["DEMULTIPLEX"] else "mappings/hiSeqIntrons/hiSeq_{techname}_{capDesign}.{barcodes}.canonicalIntrons.list" if config["USE_MATCHED_ILLUMINA"] and not config["DEMULTIPLEX"] else config["SUPPORT_INTRONS_DB"],
		lrIntrons="mappings/highConfidenceReads/introns/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.strandedHCGMs.introns.tsv",
		hcgmGTF= "mappings/highConfidenceReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.strandedHCGMs.gff.gz"
	output:"mappings/highConfidenceReads/HiSS/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.gff.gz"
	shell:
		'''
uuid=$(uuidgen)
join -v1 -1 2 -2 1 {input.lrIntrons} {input.hiSeqIntrons} |awk '{{print $2"\t"$1}}' > $TMPDIR/$uuid.{wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.barcodes}.introns.noHiSeq.tsv
zcat {input.hcgmGTF} > $TMPDIR/$uuid.$(basename {input.hcgmGTF} .gz)
cut -f1 $TMPDIR/$uuid.{wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.barcodes}.introns.noHiSeq.tsv | sort|uniq | fgrep -wv -f - $TMPDIR/$uuid.$(basename {input.hcgmGTF} .gz) |sortgff |gzip > {output}

		'''

rule getHiSSStats:
	input:
		reads = "mappings/readMapping/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bam",
		HiSSGTF="mappings/highConfidenceReads/HiSS/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.gff.gz"
	output: temp(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.stats.tsv")
	shell:
		'''
uuid=$(uuidgen)
bedtools bamtobed -i {input.reads} -bed12 > $TMPDIR/$uuid.{wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.barcodes}.merged.bed
zcat {input.HiSSGTF} | gff2bed_full.pl - > $TMPDIR/$uuid.{wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.barcodes}.HiSS.bed

mappedReadsMono=$(cat $TMPDIR/$uuid.{wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.barcodes}.merged.bed | awk '$10<=1'|cut -f4 |sort|uniq|wc -l)
mappedReadsSpliced=$(cat $TMPDIR/$uuid.{wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.barcodes}.merged.bed | awk '$10>1'|cut -f4 |sort|uniq|wc -l)

HiSSMono=$(cat $TMPDIR/$uuid.{wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.barcodes}.HiSS.bed| awk '$10<=1'|cut -f4  | sort|uniq|wc -l)
HiSSSpliced=$(cat $TMPDIR/$uuid.{wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.barcodes}.HiSS.bed| awk '$10>1'|cut -f4  | sort|uniq|wc -l)

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
cat {input} | awk '{{print $1"\\t"$2"\\t"$3"\\t"$4"\\tHCGM-mono\\t"$5"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\tHCGM-spliced\\t"$6"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\tnonHCGM-mono\\t"$7"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\tnonHCGM-spliced\\t"$8}}'| sed 's/Corr0/\tNo/' | sed 's/Corr{lastK}/\tYes/' | sort >> {output}

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
zcat {input} | tmerge --exonOverhangTolerance {config[exonOverhangTolerance]} --minReadSupport {wildcards.minReadSupport} --tmPrefix {wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.barcodes}.NAM_ - |sortgff > {output}
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
cat {input} |sortgff > $TMPDIR/$uuid
cat $TMPDIR/$uuid | tmerge --exonOverhangTolerance {config[exonOverhangTolerance]} --tmPrefix {wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.barcodes}.NAM_ - |sortgff > {output}
#cp {output} {output}.bkp
#checkTmergeOutput.sh $TMPDIR/$uuid {output} &> {output}.qc.txt
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
cat {input} |sortgff > $TMPDIR/$uuid
cat $TMPDIR/$uuid | tmerge --exonOverhangTolerance {config[exonOverhangTolerance]} --tmPrefix {wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}_{wildcards.sizeFrac}_{wildcards.barcodes}.NAM_ - |sortgff > {output}
#cp {output} {output}.bkp
#checkTmergeOutput.sh $TMPDIR/$uuid {output} &> {output}.qc.txt
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
cat {input} |sortgff > $TMPDIR/$uuid
cat $TMPDIR/$uuid | tmerge --exonOverhangTolerance {config[exonOverhangTolerance]} --tmPrefix {wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}_allFracs_allTissues.NAM_ - |sortgff > {output}
#cp {output} {output}.bkp
#checkTmergeOutput.sh $TMPDIR/$uuid {output} &> {output}.qc.txt
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
cat {input} |sortgff > $TMPDIR/$uuid
cat $TMPDIR/$uuid | tmerge --exonOverhangTolerance {config[exonOverhangTolerance]} --tmPrefix allSeqTechsCorr{wildcards.corrLevel}_{wildcards.capDesign}_allFracs_allTissues.NAM_ - |sortgff > {output}
#cp {output} {output}.bkp
#checkTmergeOutput.sh $TMPDIR/$uuid {output} &> {output}.qc.txt
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
cat {input} |sortgff > $TMPDIR/$uuid
cat $TMPDIR/$uuid | tmerge --exonOverhangTolerance {config[exonOverhangTolerance]} --tmPrefix allSeqTechsCorr{wildcards.corrLevel}_{wildcards.capDesign}_allFracs_allTissues.NAM_ - |sortgff > {output}
#cp {output} {output}.bkp
#checkTmergeOutput.sh $TMPDIR/$uuid {output} &> {output}.qc.txt
#rm {output}.bkp

		'''




rule nonAnchoredMergeReadsToBed:
	input: "mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.all.gff"
	output: temp("mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.all.bed")
	shell:
		'''
cat {input} | gff2bed_full.pl - > {output}
		'''


rule getMergingStats:
	input:
		hcgms = "mappings/highConfidenceReads/HiSS/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.gff.gz",
		pooledMerged = "mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.all.gff"
	output: temp(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.merged.stats.tsv")
	shell:
		'''
hcgms=$(zcat {input.hcgms} | extractGffAttributeValue.pl transcript_id | sort|uniq|wc -l)
merged=$(cat {input.pooledMerged} | extractGffAttributeValue.pl transcript_id | sort|uniq|wc -l)
echo -e "{wildcards.techname}Corr{wildcards.corrLevel}\t{wildcards.capDesign}\t{wildcards.sizeFrac}\t{wildcards.barcodes}\t$hcgms\t$merged" | awk '{{print $0"\t"$6/$5}}' > {output}

		'''


rule aggMergingStats:
	input: lambda wildcards: expand(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.merged.stats.tsv",filtered_product_merge, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACSpluSMERGED, barcodes=BARCODESpluSMERGED, minReadSupport=wildcards.minReadSupport)

	output: config["STATSDATADIR"] + "all.min{minReadSupport}reads.merged.stats.tsv"
	shell:
		'''
echo -e "seqTech\tcorrectionLevel\tcapDesign\tsizeFrac\ttissue\tcategory\tcount" > {output}
cat {input} | awk '{{print $1"\\t"$2"\\t"$3"\\t"$4"\\tHCGMreads\\t"$5"\\n"$1"\\t"$2"\\t"$3"\\t"$4"\\tmergedTMs\\t"$6}}' | sed 's/Corr0/\tNo/' | sed 's/Corr{lastK}/\tYes/' | sort >> {output}

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

rule getTmLengthStats:
	input:
		gencode=lambda wildcards: CAPDESIGNTOANNOTGTF[wildcards.capDesign],
		tms="mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.all.bed",
		flTms="mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.cage+polyASupported.bed"
	output:temp(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.splicedLength.stats.tsv")
	shell:
		'''
uuid=$(uuidgen)
cat {input.gencode} | awk '$3=="exon"' | fgrep "transcript_type \\"protein_coding\\";" |gff2bed_full.pl - | bed12ToTranscriptLength.pl - | awk -v t={wildcards.techname}Corr{wildcards.corrLevel} -v c={wildcards.capDesign} -v si={wildcards.sizeFrac} -v b={wildcards.barcodes} '{{print t"\t"c"\t"si"\t"b"\tGENCODE_protein_coding\t"$2}}'| sed 's/Corr0/\tNo/' | sed 's/Corr{lastK}/\tYes/' > $TMPDIR/$uuid.gencode.pcg.tsv

cat {input.tms} | bed12ToTranscriptLength.pl - | awk -v t={wildcards.techname}Corr{wildcards.corrLevel} -v c={wildcards.capDesign} -v si={wildcards.sizeFrac} -v b={wildcards.barcodes} '{{print t"\t"c"\t"si"\t"b"\tCLS_TMs\t"$2}}'| sed 's/Corr0/\tNo/' | sed 's/Corr{lastK}/\tYes/'  > $TMPDIR/$uuid.tms.tsv

cat {input.flTms} | bed12ToTranscriptLength.pl - | awk -v t={wildcards.techname}Corr{wildcards.corrLevel} -v c={wildcards.capDesign} -v si={wildcards.sizeFrac} -v b={wildcards.barcodes} '{{print t"\t"c"\t"si"\t"b"\tCLS_FL_TMs\t"$2}}'| sed 's/Corr0/\tNo/' | sed 's/Corr{lastK}/\tYes/'  > $TMPDIR/$uuid.flTms.tsv

cat $TMPDIR/$uuid.gencode.pcg.tsv $TMPDIR/$uuid.tms.tsv $TMPDIR/$uuid.flTms.tsv | sed 's/Corr0/\tNo/' | sed 's/Corr{lastK}/\tYes/'  > {output}
		'''


rule aggTmLengthStats:
	input: lambda wildcards: expand(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.splicedLength.stats.tsv",filtered_product_merge, techname=TECHNAMESplusMERGED, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACSpluSMERGED, barcodes=BARCODESpluSMERGED, minReadSupport=wildcards.minReadSupport)
	#input: lambda wildcards: expand(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}.{spliceType}.splice.sites.stats.tsv", techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, spliceType=SPLICE_SITE_TYPES)
	output: config["STATSDATADIR"] + "all.min{minReadSupport}reads.splicedLength.stats.tsv"
	shell:
		'''
echo -e "seqTech\tcorrectionLevel\tcapDesign\tsizeFrac\ttissue\tcategory\tsplicedLength" > {output}
cat {input} |sort >> {output}
		'''


rule plotTmLengthStats:
	input: config["STATSDATADIR"] + "all.min{minReadSupport}reads.splicedLength.stats.tsv"
	output: config["PLOTSDIR"] + "splicedLength.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.splicedLength.stats.{ext}"
	params:
		filterDat=lambda wildcards: merge_figures_params(wildcards.capDesign, wildcards.sizeFrac, wildcards.barcodes, wildcards.corrLevel, wildcards.techname)
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

stat_summary(aes(x=factor(correctionLevel), group=category), position=position_dodge(0.9), fun.data = fun_length, geom = 'text', vjust = +1, hjust=0, show.legend=FALSE, size=geom_textSize, color='#666666') +
stat_summary(aes(x=factor(correctionLevel), group=category), position=position_dodge(0.9), fun.data = fun_median, geom = 'text', vjust = 0, hjust=0, show.legend=FALSE, size=geom_textSize, color='#666666') +
ylab('Spliced length') +
xlab('{params.filterDat[6]}') +
coord_flip(ylim=c(100, 3500)) +
{params.filterDat[9]}
{GGPLOT_PUB_QUALITY}
ggsave('{output}', width=plotWidth, height=plotHeight)

" > {output}.r
cat {output}.r | R --slave

		'''