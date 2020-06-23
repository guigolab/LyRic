
rule makeReadsBlastDB:
	input: FQ_CORR_PATH + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.fastq.gz"
	params:
		dbprefix=DEMULTIPLEX_DIR + "blastdb/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.fastq"
	output: DEMULTIPLEX_DIR + "blastdb/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.fastq.nal"
	shell:
		'''

zcat {input} | fastq2fasta.pl | makeblastdb -in - -dbtype nucl -out {params.dbprefix} -title {params.dbprefix}
if [ ! -f {output} ]; then
echoerr "File {output} not created. Making one"
echo -e "#
# Alias file created: $(date)
#
TITLE {params.dbprefix}\nDBLIST {wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}_{wildcards.sizeFrac}.fastq" > {output}
fi

		'''

rule blastDemultiplex:
	input:
		reads =DEMULTIPLEX_DIR + "blastdb/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.fastq.nal" ,
		adapters =DEMULTIPLEX_DIR + "all_adapters.seq_ids.fa"
	params:
		dbprefix=DEMULTIPLEX_DIR + "blastdb/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.fastq"
	threads: 12
	output:
		DEMULTIPLEX_DIR + "blastDemultiplex/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.fastq.all_adapters.blast.tsv.gz"
	shell:
		'''
#-dust no
blastn -query {input.adapters} -task blastn -soft_masking false -word_size 4 -qcov_hsp_perc 70 -evalue 1 -max_target_seqs 100000000 -num_threads {threads} -db {params.dbprefix} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen score gaps" | gzip > {output}

		'''

rule processBlastDemultiplex:
	input:
		blastOut=DEMULTIPLEX_DIR + "blastDemultiplex/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.fastq.all_adapters.blast.tsv",
		adaptersTsv=adaptersTSV
	output: DEMULTIPLEX_DIR + "processBlastDemultiplex/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.fastq.demultiplex.tsv"
	shell:
		'''
zcat {input.blastOut}| demultiplexFromBlastTab.pl - | barcodeSeqTobarcodeId.pl {input.adaptersTsv} - {wildcards.capDesign} > {output}

		'''

rule discardErroneousandBarcodesAndUP:
#remove records with barcodes from wrong capture design
#this will also exclude "UP" records
	input: DEMULTIPLEX_DIR + "processBlastDemultiplex/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.fastq.demultiplex.tsv"
	output: DEMULTIPLEX_DIR + "discardErroneousandBarcodesAndUP/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.demultiplexSameCapDesignBarcode.tsv"
	shell:
		'''
cat {input} | fgrep "{wildcards.capDesign}_" |sort -T {config[TMPDIR]} |uniq > {output}

		'''

rule getBasicDemultiplexingStats:
	input:
		totals = config["STATSDATADIR"] + "{techname}Corr{corrLevel}.fastq.readCounts.tsv",
		demul = DEMULTIPLEX_DIR + "processBlastDemultiplex/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.fastq.demultiplex.tsv",
		demulFiltered = DEMULTIPLEX_DIR + "discardErroneousandBarcodesAndUP/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.demultiplexSameCapDesignBarcode.tsv"
	output: config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.fastq.basicDemultiplexing.stats.tsv"
	shell:
		'''
total=$(cat {input.totals} | awk -v f={wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}_{wildcards.sizeFrac}.fastq.gz '$1==f' | cut -f2)
withUP=$(cat {input.demul} | awk '$5=="UP"'| cut -f1|sort -T {config[TMPDIR]} |uniq|wc -l)
withBarcode=$(cat {input.demul} | fgrep -v -w "UP" | cut -f1|sort -T {config[TMPDIR]} |uniq|wc -l)
withSameCapDesignBarcode=$(cat {input.demulFiltered}| cut -f1|sort -T {config[TMPDIR]} |uniq|wc -l)
let withDifferentCapDesignBarcode=$withBarcode-$withSameCapDesignBarcode
echo -e "{wildcards.capDesign}\t{wildcards.sizeFrac}\t$total\t$withUP\t$withBarcode\t$withSameCapDesignBarcode\t$withDifferentCapDesignBarcode" | awk '{{print $0"\t"$4/$3"\t"$5/$3"\t"$6/$3"\t"$7/$5}}' > {output}

		'''


rule aggBasicDemultiplexingStats:
	input:
		lambda wildcards: expand(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.fastq.basicDemultiplexing.stats.tsv", filtered_product, techname=wildcards.techname, corrLevel=wildcards.corrLevel, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS)
	output:config["STATSDATADIR"] + "{techname}Corr{corrLevel}.fastq.basicDemultiplexing.stats.tsv"
	shell:
		'''
echo -e "capDesign\tsizeFraction\treads\treadsWithUP\treadsWithAnyBarcode\treadsWithSameCapDesignBarcode\treadsWithDifferentCapDesignBarcode\tpcReadsWithUP\tpcReadsWithAnyBarcode\tpcReadsWithSameCapDesignBarcode\tpcReadsWithWrongBarcode"> {output};
cat {input} | sort -T {config[TMPDIR]}  >> {output}
		'''

rule plotUpStats:
	input: config["STATSDATADIR"] + "{techname}Corr{corrLevel}.fastq.basicDemultiplexing.stats.tsv"
	output: config["PLOTSDIR"] + "{techname}Corr{corrLevel}.fastq.UP.stats.{ext}" if config["DEMULTIPLEX"] else "dummy1.txt"
	shell:
		'''
echo "library(ggplot2)
library(plyr)
library(scales)
dat <- read.table('{input}', header=T, as.is=T, sep='\\t')
ggplot(data=dat, aes(x=capDesign, y=pcReadsWithUP, fill=sizeFraction)) +
geom_bar(width=0.75,stat='identity', position=position_dodge(width=0.9)) +
scale_fill_manual(values={sizeFrac_Rpalette}) +
geom_text(aes(group=sizeFraction, y=pcReadsWithUP-0.0005, label = paste(sep='',percent(pcReadsWithUP),'\n','(',comma(readsWithUP),')')), size=3, position = position_dodge(width=0.9)) +
scale_y_continuous(labels = scales::percent) +
xlab('Capture design') +
ylab('% reads with UP') +
{GGPLOT_PUB_QUALITY}
ggsave('{output}', width=7, height=3)
" > {output}.r
 set +eu
conda activate R_env
set -eu
cat {output}.r | R --slave
		'''

rule plotBcStats:
	input: config["STATSDATADIR"] + "{techname}Corr{corrLevel}.fastq.basicDemultiplexing.stats.tsv"
	output: config["PLOTSDIR"] + "{techname}Corr{corrLevel}.fastq.BC.stats.{ext}"
	shell:
		'''
echo "library(ggplot2)
library(plyr)
library(scales)
dat <- read.table('{input}', header=T, as.is=T, sep='\\t')
ggplot(data=dat, aes(x=capDesign, y=pcReadsWithAnyBarcode, fill=sizeFraction)) +
geom_bar(width=0.75,stat='identity', position=position_dodge(width=0.9)) +
scale_fill_manual(values={sizeFrac_Rpalette}) +
geom_text(aes(group=sizeFraction, y=pcReadsWithAnyBarcode-0.0005, label = paste(sep='',percent(pcReadsWithAnyBarcode),'\n','(',comma(readsWithAnyBarcode),')')), size=3, position = position_dodge(width=0.9)) +
scale_y_continuous(labels = scales::percent) +
xlab('Capture design') +
ylab('% reads with any barcode') +
{GGPLOT_PUB_QUALITY}
ggsave('{output}', width=7, height=3)
" > {output}.r
 set +eu
conda activate R_env
set -eu
cat {output}.r | R --slave

		'''


rule plotForeignBcStats:
	input: config["STATSDATADIR"] +  "{techname}Corr{corrLevel}.fastq.basicDemultiplexing.stats.tsv"
	output: config["PLOTSDIR"] + "{techname}Corr{corrLevel}.fastq.foreignBC.stats.{ext}"
	shell:
		'''
echo "library(ggplot2)
library(plyr)
library(scales)
dat <- read.table('{input}', header=T, as.is=T, sep='\\t')
ggplot(data=dat, aes(x=capDesign, y=pcReadsWithWrongBarcode, fill=sizeFraction)) +
geom_bar(width=0.75,stat='identity', position=position_dodge(width=0.9)) +
scale_fill_manual(values={sizeFrac_Rpalette}) +
geom_text(aes(group=sizeFraction, y=pcReadsWithWrongBarcode-0.0005, label = paste(sep='',percent(pcReadsWithWrongBarcode),'\n','(',comma(readsWithDifferentCapDesignBarcode),')')), size=3, position = position_dodge(width=0.9)) +
scale_y_continuous(labels = scales::percent) +
xlab('Capture design') +
ylab('% reads with foreign barcode') +
{GGPLOT_PUB_QUALITY}
ggsave('{output}', width=7, height=3)
" > {output}.r
  set +eu
conda activate R_env
set -eu
cat {output}.r | R --slave

		'''



rule discardAmbiguousBarcodes:
#make list of reads with multiple distinct barcodes
	input: DEMULTIPLEX_DIR + "processBlastDemultiplex/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.fastq.demultiplex.tsv"
	output: DEMULTIPLEX_DIR + "discardAmbiguousBarcodes/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.ambiguousDemultiplex.list"
	shell:
		'''
cat {input} | awk '$5!="UP"'|cut -f1,5 |sort -T {config[TMPDIR]} |uniq | cut -f1|sort -T {config[TMPDIR]} |uniq -d > {output}

		'''

rule getAmbiguousBarcodeStats:
	input:
		demul=DEMULTIPLEX_DIR + "processBlastDemultiplex/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.fastq.demultiplex.tsv",
		ambig=DEMULTIPLEX_DIR + "discardAmbiguousBarcodes/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.ambiguousDemultiplex.list"
	output: config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.ambiguousBarcodes.reads.stats.tsv"
	shell:
		'''
totalDemul=$(cat {input.demul} | cut -f1|sort -T {config[TMPDIR]} |uniq | wc -l)
totalAmbigBarcode=$(cat {input.ambig} | wc -l)
echo -e "{wildcards.capDesign}\t{wildcards.sizeFrac}\t$totalDemul\t$totalAmbigBarcode" | awk '{{print $0"\t"$4/$3}}' > {output}

		'''

rule aggAmbiguousBarcodeStats:
	input: #expand(config["STATSDATADIR"] + "{{techname}Corr{corrLevel}}_{capDesign}_{sizeFrac}.ambiguousBarcodes.reads.stats.tsv", capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS)
		lambda wildcards: expand(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.ambiguousBarcodes.reads.stats.tsv", filtered_product, techname=wildcards.techname , corrLevel=wildcards.corrLevel, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS)
	output:config["STATSDATADIR"] + "{techname}Corr{corrLevel}.ambiguousBarcodes.reads.stats.tsv"
	shell:
		'''
cat {input} | sort -T {config[TMPDIR]}  > {output}
		'''

rule plotAmbiguousBarcodeStats:
	input: config["STATSDATADIR"] + "{techname}Corr{corrLevel}.ambiguousBarcodes.reads.stats.tsv"
	output: config["PLOTSDIR"] + "{techname}Corr{corrLevel}.ambiguousBarcodes.reads.stats.{ext}"
	shell:
		'''
echo "library(ggplot2)
library(plyr)
library(scales)
dat <- read.table('{input}', header=F, as.is=T, sep='\\t')
colnames(dat)<-c('capDesign', 'fraction','total','ambig', 'percentAmbig')
ggplot(data=dat, aes(x=capDesign, y=percentAmbig, fill=fraction)) +
geom_bar(stat='identity', position=position_dodge()) +
scale_fill_manual(values={sizeFrac_Rpalette}) +
geom_text(aes(group=fraction, y=percentAmbig+0.0005, label = comma(ambig)), size=3.5, position = position_dodge(width=0.9)) +
scale_y_continuous(labels = scales::percent) +
xlab('Capture design') +
ylab('% reads with multiple\ndistinct barcodes') +
{GGPLOT_PUB_QUALITY}
ggsave('{output}', width=7, height=3)
" > {output}.r
  set +eu
conda activate R_env
set -eu

 cat {output}.r | R --slave
		'''

rule getAdaptersLocationOverReads:
	input: DEMULTIPLEX_DIR + "processBlastDemultiplex/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.fastq.demultiplex.tsv"
	output: config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.adapters.location.stats.tsv"
	shell:
		'''
cat {input} | perl -slane '@line=split "\t"; if($line[4] ne "UP"){{$line[4]="BC"}}; print "$frac\t$line[1]\t$line[4]"' -- -frac={wildcards.sizeFrac} > {output}

		'''

rule aggAdaptersLocationOverReads:
	input:
		lambda wildcards: expand(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.adapters.location.stats.tsv", filtered_product, techname=wildcards.techname, corrLevel=wildcards.corrLevel, capDesign=wildcards.capDesign, sizeFrac=SIZEFRACS)
	output: config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}.all.adapters.location.stats.tsv"
	shell:
		'''
cat {input} | sort -T {config[TMPDIR]}  > {output}
		'''

rule plotAdaptersLocationOverReads:
	input: config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}.all.adapters.location.stats.tsv"
	output: config["PLOTSDIR"] + "{techname}Corr{corrLevel}_{capDesign}.adapters.location.stats.{ext}"
	shell:
		'''
echo "library(ggplot2)
library(plyr)
dat <- read.table('{input}', header=F, as.is=T, sep='\\t')
palette <- c('BC' ='#9f3ca4','UP' ='#80cc28')
colnames(dat)<-c('fraction','location','adapter')
ggplot(dat, aes(x=location, fill=adapter)) +
geom_histogram(binwidth=0.01, position='dodge') +
facet_grid( fraction ~ ., scales='free') +
xlab ('Mapped location of adapter on ROI') +
scale_fill_manual(values=palette) +
scale_x_continuous(breaks=seq(0, 1, 0.1),labels = scales::percent) +
theme_bw(base_size=17) +
{GGPLOT_PUB_QUALITY}
ggsave('{output}', width=9, height=3)
" > {output}.r
 set +eu
conda activate R_env
set -eu
cat {output}.r | R --slave

		'''

rule detectReadChimeras:
	input: DEMULTIPLEX_DIR + "processBlastDemultiplex/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.fastq.demultiplex.tsv"
	output: DEMULTIPLEX_DIR + "detectReadChimeras/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.fastq.chimeras.list"
	shell:
		'''
cat {input} | awk '$2>0.1 && $2<0.9' | cut -f1 |sort -T {config[TMPDIR]} |uniq > {output}

		'''

rule getRateOfReadChimeras:
	input:
		totals=config["STATSDATADIR"] + "{techname}Corr{corrLevel}.fastq.readCounts.tsv",
		chim=DEMULTIPLEX_DIR + "detectReadChimeras/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.fastq.chimeras.list"
	output:config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.chimeric.reads.stats.tsv"
	shell:
		'''
total=$(cat {input.totals} | awk -v f={wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}_{wildcards.sizeFrac}.fastq.gz '$1==f' | cut -f2)
chimeric=$(cat {input.chim} | wc -l)
echo -e "{wildcards.capDesign}\t{wildcards.sizeFrac}\t$total\t$chimeric" | awk '{{print $0"\t"$4/$3}}' > {output}
		'''

rule aggRateOfReadChimeras:
	input:
		lambda wildcards: expand(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.chimeric.reads.stats.tsv", filtered_product, techname=wildcards.techname, corrLevel=wildcards.corrLevel, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS)
	output: config["STATSDATADIR"] + "{techname}Corr{corrLevel}.chimeric.reads.stats.tsv"
	shell:
		'''
cat {input} |sort -T {config[TMPDIR]}  > {output}
		'''

rule plotRateOfReadChimeras:
	input: config["STATSDATADIR"] + "{techname}Corr{corrLevel}.chimeric.reads.stats.tsv"
	output: config["PLOTSDIR"] + "{techname}Corr{corrLevel}.chimeric.reads.stats.{ext}"
	shell:
		'''
echo "library(ggplot2)
library(plyr)
library(scales)
dat <- read.table('{input}', header=F, as.is=T, sep='\\t')
colnames(dat)<-c('capDesign', 'fraction','total','chimeric', 'percentChimeric')
ggplot(data=dat, aes(x=capDesign, y=percentChimeric, fill=fraction)) +
geom_bar(stat='identity', position=position_dodge()) +
scale_fill_manual(values={sizeFrac_Rpalette}) +
geom_text(aes(group=fraction, y=percentChimeric+0.005, label = comma(chimeric)), size=3.5, position = position_dodge(width=0.9)) +
scale_y_continuous(labels = scales::percent) +
xlab('Capture design') +
ylab('% chimeric reads') +
{GGPLOT_PUB_QUALITY}
ggsave('{output}', width=7, height=3)
" > {output}.r
 set +eu
conda activate R_env
set -eu

cat {output}.r | R --slave

		'''

rule selectNonAmbiguousNonChimericReads:
	input:
		demulFiltered = DEMULTIPLEX_DIR + "discardErroneousandBarcodesAndUP/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.demultiplexSameCapDesignBarcode.tsv",
		chim = DEMULTIPLEX_DIR + "detectReadChimeras/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.fastq.chimeras.list",
		ambig = DEMULTIPLEX_DIR + "discardAmbiguousBarcodes/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.ambiguousDemultiplex.list"
	output: DEMULTIPLEX_DIR + "selectNonAmbiguousNonChimericReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.demultiplex.noAmbig.noChim.tsv"
	shell:
		'''
cat {input.chim} {input.ambig} | sort -T {config[TMPDIR]} |uniq | fgrep -w -v -f - {input.demulFiltered} |sort -T {config[TMPDIR]} |uniq> {output}

		'''


rule getNonAmbiguousNonChimericReadsStats:
	input:
		totals=config["STATSDATADIR"] + "{techname}Corr{corrLevel}.fastq.readCounts.tsv",
		noAmbigNoChim= DEMULTIPLEX_DIR + "selectNonAmbiguousNonChimericReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.demultiplex.noAmbig.noChim.tsv"
	output: config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.finalDemul.reads.stats.tsv"
	shell:
		'''
total=$(cat {input.totals} | awk -v f={wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}_{wildcards.sizeFrac}.fastq.gz '$1==f' | cut -f2)
finalDemul=$(cat {input.noAmbigNoChim} | cut -f1 |sort -T {config[TMPDIR]} |uniq|wc -l)
echo -e "{wildcards.capDesign}\t{wildcards.sizeFrac}\t$total\t$finalDemul" | awk '{{print $0"\t"$4/$3}}' > {output}

		'''

rule aggNonAmbiguousNonChimericReadsStats:
	input:
		lambda wildcards: expand(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.finalDemul.reads.stats.tsv", filtered_product, techname=wildcards.techname, corrLevel=wildcards.corrLevel, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS)
	output: config["STATSDATADIR"] + "{techname}Corr{corrLevel}.finalDemul.reads.stats.tsv"
	shell:
		'''
cat {input} | sort -T {config[TMPDIR]}  > {output}
		'''

rule plotNonAmbiguousNonChimericReadsStats:
	input: config["STATSDATADIR"] + "{techname}Corr{corrLevel}.finalDemul.reads.stats.tsv"
	output: config["PLOTSDIR"] + "{techname}Corr{corrLevel}.finalDemul.reads.stats.{ext}"
	shell:
		'''
echo "library(ggplot2)
library(plyr)
library(scales)
dat <- read.table('{input}', header=F, as.is=T, sep='\\t')
colnames(dat)<-c('capDesign', 'fraction','total','finalDemul', 'percentFinalDemul')
ggplot(data=dat, aes(x=capDesign, y=percentFinalDemul, fill=fraction)) +
geom_bar(width=0.75,stat='identity', position=position_dodge(width=0.9)) +
scale_fill_manual(values={sizeFrac_Rpalette}) +
geom_text(aes(group=fraction, y=percentFinalDemul-0.01, label = paste(sep='',percent(percentFinalDemul),'\n','(',comma(finalDemul),')')), size=3, position = position_dodge(width=0.9)) +
scale_y_continuous(labels = scales::percent) +
xlab('Capture design') +
ylab('% unambiguously demultiplexed,\nnon-chimeric reads') +
{GGPLOT_PUB_QUALITY}
ggsave('{output}', width=7, height=3)
" > {output}.r
 set +eu
conda activate R_env
set -eu
cat {output}.r | R --slave

		'''


rule checkForOnlyOneBarcodePerRead:
	input:
		DEMULTIPLEX_DIR + "selectNonAmbiguousNonChimericReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.demultiplex.noAmbig.noChim.tsv"
	output: DEMULTIPLEX_DIR + "qc/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.demul.QC1.txt"
	shell:
		'''
cut -f1,5 {input} | sort -T {config[TMPDIR]} |uniq | cut -f1 |sort -T {config[TMPDIR]} |uniq -d |wc -l > {output}
cat {output} | while read count; do if [ $count -gt 0 ]; then echo "$count duplicates found";mv {output} {output}.tmp;  exit 1; fi; done
		'''

rule convertFastqToTsv:
	input: FQ_CORR_PATH + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.fastq.gz"
	output: temp(DEMULTIPLEX_DIR + "convertFastqToTsv/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.tsv")
	shell:
		'''
zcat  {input} | fastq2tsv.pl > {output}
		'''

rule demultiplexFastqs:
	input:
		tsvFastq = DEMULTIPLEX_DIR + "convertFastqToTsv/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.tsv",
		noAmbigNoChim= DEMULTIPLEX_DIR + "selectNonAmbiguousNonChimericReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.demultiplex.noAmbig.noChim.tsv"
	output: DEMULTIPLEXED_FASTQS + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.{barcodes}.fastq.gz"
	shell:
		'''
tmpUuid=$(uuidgen)
cat {input.noAmbigNoChim} | awk -v b={wildcards.barcodes} '$5==b'| cut -f1|sort -T {config[TMPDIR]} |uniq > {config[TMPDIR]}/$tmpUuid
tgrep -F -w -f {config[TMPDIR]}/$tmpUuid {input.tsvFastq} > {config[TMPDIR]}/$tmpUuid.fastq
cat {config[TMPDIR]}/$tmpUuid.fastq | tsv2fastq.pl | gzip> {output}

		'''

rule checkTotalsDemultiplexed:
	input:
		demul= lambda wildcards: expand(DEMULTIPLEXED_FASTQS + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.{barcodes}.fastq.gz", filtered_product, techname=wildcards.techname, corrLevel=wildcards.corrLevel, capDesign=wildcards.capDesign, sizeFrac=wildcards.sizeFrac, barcodes=BARCODES),
		totals=config["STATSDATADIR"] + "{techname}Corr{corrLevel}.fastq.readCounts.tsv"
	output: DEMULTIPLEX_DIR + "qc/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.demul.QC2.txt"
	shell:
		'''
total=$(cat {input.totals} | awk -v f={wildcards.techname}Corr{wildcards.corrLevel}_{wildcards.capDesign}_{wildcards.sizeFrac}.fastq.gz '$1==f' | cut -f2)
uuid=$(uuidgen)
for file in `echo {input.demul}`; do zcat $file; done > {config[TMPDIR]}/$uuid.demul.fastq
demul=$(cat {config[TMPDIR]}/$uuid.demul.fastq | fastq2tsv.pl | wc -l)
echo -e "{wildcards.capDesign}.{wildcards.sizeFrac}\t$total\t$demul" | awk '{{print $2-$3}}' > {output}
cat {output}| while read diff; do if [ $diff -lt 1 ]; then echo "Number of demultiplexed reads is greater than number of input reads (diff: $diff)";mv {output} {output}.tmp;  exit 1; fi; done

		'''



rule getDemultiplexingStatsPerSample:
	input: DEMULTIPLEXED_FASTQS + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.{barcodes}.fastq.gz"
	output: config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.{barcodes}.demultiplexing.perSample.stats.tsv"
	shell:
		'''
demul=$(zcat {input} | fastq2tsv.pl | wc -l)
echo -e "{wildcards.capDesign}\t{wildcards.sizeFrac}\t{wildcards.barcodes}\t$demul" > {output}

		'''

rule aggDemultiplexingStatsPerSample:
	input: #expand(config["STATSDATADIR"] + "{{techname}Corr{corrLevel}}_{capDesign}_{sizeFrac}.{barcodes}.demultiplexing.perSample.stats.tsv", capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES)
		lambda wildcards: expand(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.{barcodes}.demultiplexing.perSample.stats.tsv", filtered_product, techname=wildcards.techname, corrLevel=wildcards.corrLevel, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES)
	output: config["STATSDATADIR"] + "{techname}Corr{corrLevel}.demultiplexing.perSample.stats.tsv"
	shell:
		'''
cat {input} | sort -T {config[TMPDIR]}  | perl -slane 'if($F[2]=~/$F[0]/ || $F[2] eq "Undeter"){{print}}'> {output}

		'''
rule plotDemultiplexingStatsPerSample:
	input: config["STATSDATADIR"] + "{techname}Corr{corrLevel}.demultiplexing.perSample.stats.tsv"
	output: config["PLOTSDIR"] + "{techname}Corr{corrLevel}.demultiplexing.perSample.stats.{ext}"
	shell:
		'''
echo "library(ggplot2)
library(plyr)
library(scales)
dat <- read.table('{input}', header=F, as.is=T, sep='\\t')
colnames(dat)<-c('capDesign', 'sizeFraction','barcode','count')
ggplot(dat, aes(x=barcode, y=count, fill=sizeFraction)) +
geom_bar(width=0.75,stat='identity') +
#geom_text(aes(y=count-0.01, label = comma(count)), size=3, angle=90) +
scale_fill_manual(values={sizeFrac_Rpalette}) +
facet_grid(sizeFraction ~ capDesign, scales='free') +
xlab ('Sample (barcode)') +
theme_bw(base_size=17) +
{GGPLOT_PUB_QUALITY} + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('{output}', width=13, height=9)
" > {output}.r
 set +eu
conda activate R_env
set -eu
cat {output}.r | R --slave

		'''
