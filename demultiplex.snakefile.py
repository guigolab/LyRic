
rule makeReadsBlastDB:
	input: config["FQPATH"] + "{techname}_{capDesign}_{sizeFrac}.fastq.gz"
	params:
		dbprefix=config["DEMULTIPLEX_DIR"] + "blastdb/" + "{techname}_{capDesign}_{sizeFrac}.fastq"
	output: config["DEMULTIPLEX_DIR"] + "blastdb/" + "{techname}_{capDesign}_{sizeFrac}.fastq.nal"
	shell:
		'''

zcat {input} | fastq2fasta.pl | makeblastdb -in - -dbtype nucl -out {params.dbprefix} -title {params.dbprefix}
if [ ! -f {output} ]; then
echoerr "File {output} not created. Making one"
echo -e "#
# Alias file created: $(date)
#
TITLE {params.dbprefix}\nDBLIST {wildcards.techname}_{wildcards.capDesign}_{wildcards.sizeFrac}.fastq" > {output}
fi

		'''

rule blastDemultiplex:
	input:
		reads =config["DEMULTIPLEX_DIR"] + "blastdb/" + "{techname}_{capDesign}_{sizeFrac}.fastq.nal" ,
		adapters =config["DEMULTIPLEX_DIR"] + "all_adapters.seq_ids.fa"
	params:
		dbprefix=config["DEMULTIPLEX_DIR"] + "blastdb/" + "{techname}_{capDesign}_{sizeFrac}.fastq"
	threads: 12
	output:
		temp(config["DEMULTIPLEX_DIR"] + "blastDemultiplex/{techname}_{capDesign}_{sizeFrac}.fastq.all_adapters.blast.tsv")
	shell:
		'''
#-dust no
blastn -query {input.adapters} -task blastn -soft_masking false -word_size 4 -qcov_hsp_perc 70 -evalue 1 -max_target_seqs 100000000 -num_threads {threads} -db {params.dbprefix} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen score gaps" > {output}

		'''

rule processBlastDemultiplex:
	input:
		blastOut=config["DEMULTIPLEX_DIR"] + "blastDemultiplex/{techname}_{capDesign}_{sizeFrac}.fastq.all_adapters.blast.tsv",
		adaptersTsv=config["DEMULTIPLEX_DIR"] + "all_adapters.tsv"
	output: temp(config["DEMULTIPLEX_DIR"] + "processBlastDemultiplex/{techname}_{capDesign}_{sizeFrac}.fastq.demultiplex.tsv")
	shell:
		'''
cat {input.blastOut}| demultiplexFromBlastTab.pl - | barcodeSeqTobarcodeId.pl {input.adaptersTsv} - {wildcards.capDesign} > {output}

		'''

rule discardErroneousandBarcodesAndUP:
#remove records with barcodes from wrong capture design
#this will also exclude "UP" records
	input: config["DEMULTIPLEX_DIR"] + "processBlastDemultiplex/{techname}_{capDesign}_{sizeFrac}.fastq.demultiplex.tsv"
	output: temp(config["DEMULTIPLEX_DIR"] + "discardErroneousandBarcodesAndUP/{techname}_{capDesign}_{sizeFrac}.demultiplexSameCapDesignBarcode.tsv")
	shell:
		'''
cat {input} | fgrep "{wildcards.capDesign}_" |sort|uniq > {output}

		'''

rule getBasicDemultiplexingStats:
	input:
		totals = config["STATSDATADIR"] + "{techname}.fastq.readCounts.tsv",
		demul = config["DEMULTIPLEX_DIR"] + "processBlastDemultiplex/{techname}_{capDesign}_{sizeFrac}.fastq.demultiplex.tsv",
		demulFiltered = config["DEMULTIPLEX_DIR"] + "discardErroneousandBarcodesAndUP/{techname}_{capDesign}_{sizeFrac}.demultiplexSameCapDesignBarcode.tsv"
	output: temp(config["STATSDATADIR"] + "{techname}_{capDesign}_{sizeFrac}.fastq.basicDemultiplexing.stats.tsv")
	shell:
		'''
total=$(cat {input.totals} | awk -v f={wildcards.techname}_{wildcards.capDesign}_{wildcards.sizeFrac}.fastq.gz '$1==f' | cut -f2)
withUP=$(cat {input.demul} | awk '$5=="UP"'| cut -f1|sort|uniq|wc -l)
withBarcode=$(cat {input.demul} | fgrep -v -w "UP" | cut -f1|sort|uniq|wc -l)
withSameCapDesignBarcode=$(cat {input.demulFiltered}| cut -f1|sort|uniq|wc -l)
let withDifferentCapDesignBarcode=$withBarcode-$withSameCapDesignBarcode
echo -e "{wildcards.capDesign}\t{wildcards.sizeFrac}\t$total\t$withUP\t$withBarcode\t$withSameCapDesignBarcode\t$withDifferentCapDesignBarcode" | awk '{{print $0"\t"$4/$3"\t"$5/$3"\t"$6/$3"\t"$7/$5}}' > {output}

		'''


rule aggBasicDemultiplexingStats:
	input: #expand(config["STATSDATADIR"] + "{{techname}}_{capDesign}_{sizeFrac}.fastq.basicDemultiplexing.stats.tsv", capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS)
		lambda wildcards: expand(config["STATSDATADIR"] + "{techname}_{capDesign}_{sizeFrac}.fastq.basicDemultiplexing.stats.tsv", filtered_product, techname=wildcards.techname, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS)
	output:config["STATSDATADIR"] + "{techname}.fastq.basicDemultiplexing.stats.tsv"
	shell:
		'''
echo -e "capDesign\tsizeFraction\treads\treadsWithUP\treadsWithAnyBarcode\treadsWithSameCapDesignBarcode\treadsWithDifferentCapDesignBarcode\tpcReadsWithUP\tpcReadsWithAnyBarcode\tpcReadsWithSameCapDesignBarcode\tpcReadsWithWrongBarcode"> {output};
cat {input} | sort >> {output}
		'''

rule plotUpStats:
	input: config["STATSDATADIR"] + "{techname}.fastq.basicDemultiplexing.stats.tsv"
	output: config["PLOTSDIR"] + "{techname}.fastq.UP.stats.{ext}"
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
cat {output}.r | R --slave
#dropbox_uploader.sh upload {output} {DROPBOX_PLOTS};
		'''

rule plotBcStats:
	input: config["STATSDATADIR"] + "{techname}.fastq.basicDemultiplexing.stats.tsv"
	output: config["PLOTSDIR"] + "{techname}.fastq.BC.stats.{ext}"
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
cat {output}.r | R --slave
#dropbox_uploader.sh upload {output} {DROPBOX_PLOTS};

		'''


rule plotForeignBcStats:
	input: config["STATSDATADIR"] +  "{techname}.fastq.basicDemultiplexing.stats.tsv"
	output: config["PLOTSDIR"] + "{techname}.fastq.foreignBC.stats.{ext}"
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
cat {output}.r | R --slave
#dropbox_uploader.sh upload {output} {DROPBOX_PLOTS};

		'''



rule discardAmbiguousBarcodes:
#make list of reads with multiple distinct barcodes
	input: config["DEMULTIPLEX_DIR"] + "processBlastDemultiplex/{techname}_{capDesign}_{sizeFrac}.fastq.demultiplex.tsv"
	output: config["DEMULTIPLEX_DIR"] + "discardAmbiguousBarcodes/{techname}_{capDesign}_{sizeFrac}.ambiguousDemultiplex.list"
	shell:
		'''
cat {input} | awk '$5!="UP"'|cut -f1,5 |sort|uniq | cut -f1|sort|uniq -d > {output}

		'''

rule getAmbiguousBarcodeStats:
	input:
		demul=config["DEMULTIPLEX_DIR"] + "processBlastDemultiplex/{techname}_{capDesign}_{sizeFrac}.fastq.demultiplex.tsv",
		ambig=config["DEMULTIPLEX_DIR"] + "discardAmbiguousBarcodes/{techname}_{capDesign}_{sizeFrac}.ambiguousDemultiplex.list"
	output: temp(config["STATSDATADIR"] + "{techname}_{capDesign}_{sizeFrac}.ambiguousBarcodes.reads.stats.tsv")
	shell:
		'''
totalDemul=$(cat {input.demul} | cut -f1|sort|uniq | wc -l)
totalAmbigBarcode=$(cat {input.ambig} | wc -l)
echo -e "{wildcards.capDesign}\t{wildcards.sizeFrac}\t$totalDemul\t$totalAmbigBarcode" | awk '{{print $0"\t"$4/$3}}' > {output}

		'''

rule aggAmbiguousBarcodeStats:
	input: #expand(config["STATSDATADIR"] + "{{techname}}_{capDesign}_{sizeFrac}.ambiguousBarcodes.reads.stats.tsv", capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS)
		lambda wildcards: expand(config["STATSDATADIR"] + "{techname}_{capDesign}_{sizeFrac}.ambiguousBarcodes.reads.stats.tsv", filtered_product, techname=wildcards.techname , capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS)
	output:config["STATSDATADIR"] + "{techname}.ambiguousBarcodes.reads.stats.tsv"
	shell:
		'''
cat {input} | sort > {output}
		'''

rule plotAmbiguousBarcodeStats:
	input: config["STATSDATADIR"] + "{techname}.ambiguousBarcodes.reads.stats.tsv"
	output: config["PLOTSDIR"] + "{techname}.ambiguousBarcodes.reads.stats.{ext}"
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
cat {output}.r | R --slave
#dropbox_uploader.sh upload {output} {DROPBOX_PLOTS};
		'''

rule getAdaptersLocationOverReads:
#	params:
#		sizeFrac=expand(SIZEFRACS)
	input: config["DEMULTIPLEX_DIR"] + "processBlastDemultiplex/{techname}_{capDesign}_{sizeFrac}.fastq.demultiplex.tsv"
	output: temp(config["STATSDATADIR"] + "{techname}_{capDesign}_{sizeFrac}.adapters.location.stats.tsv")
	shell:
		'''
cat {input} | perl -slane '@line=split "\t"; if($line[4] ne "UP"){{$line[4]="BC"}}; print "$frac\t$line[1]\t$line[4]"' -- -frac={wildcards.sizeFrac} > {output}

		'''

rule aggAdaptersLocationOverReads:
	input: #expand(config["STATSDATADIR"] + "{{techname}}_{{capDesign}}_{sizeFrac}.adapters.location.stats.tsv", sizeFrac=SIZEFRACS)
		lambda wildcards: expand(config["STATSDATADIR"] + "{techname}_{capDesign}_{sizeFrac}.adapters.location.stats.tsv", filtered_product, techname=wildcards.techname, capDesign=wildcards.capDesign, sizeFrac=SIZEFRACS)
	output: config["STATSDATADIR"] + "{techname}_{capDesign}.all.adapters.location.stats.tsv"
	shell:
		'''
cat {input} | sort > {output}
		'''

rule plotAdaptersLocationOverReads:
	input: config["STATSDATADIR"] + "{techname}_{capDesign}.all.adapters.location.stats.tsv"
	output: config["PLOTSDIR"] + "{techname}_{capDesign}.adapters.location.stats.{ext}"
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
cat {output}.r | R --slave
#dropbox_uploader.sh upload {output} {DROPBOX_PLOTS};

		'''

rule detectReadChimeras:
	input: config["DEMULTIPLEX_DIR"] + "processBlastDemultiplex/{techname}_{capDesign}_{sizeFrac}.fastq.demultiplex.tsv"
	output: temp(config["DEMULTIPLEX_DIR"] + "detectReadChimeras/{techname}_{capDesign}_{sizeFrac}.fastq.chimeras.list")
	shell:
		'''
cat {input} | awk '$2>0.1 && $2<0.9' | cut -f1 |sort|uniq > {output}

		'''

rule getRateOfReadChimeras:
	input:
		totals=config["STATSDATADIR"] + "{techname}.fastq.readCounts.tsv",
		chim=config["DEMULTIPLEX_DIR"] + "detectReadChimeras/{techname}_{capDesign}_{sizeFrac}.fastq.chimeras.list"
	output:temp(config["STATSDATADIR"] + "{techname}_{capDesign}_{sizeFrac}.chimeric.reads.stats.tsv")
	shell:
		'''
total=$(cat {input.totals} | awk -v f={wildcards.techname}_{wildcards.capDesign}_{wildcards.sizeFrac}.fastq.gz '$1==f' | cut -f2)
chimeric=$(cat {input.chim} | wc -l)
echo -e "{wildcards.capDesign}\t{wildcards.sizeFrac}\t$total\t$chimeric" | awk '{{print $0"\t"$4/$3}}' > {output}
		'''

rule aggRateOfReadChimeras:
	input: #expand(config["STATSDATADIR"] + "{{techname}}_{capDesign}_{sizeFrac}.chimeric.reads.stats.tsv", capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS)
		lambda wildcards: expand(config["STATSDATADIR"] + "{techname}_{capDesign}_{sizeFrac}.chimeric.reads.stats.tsv", filtered_product, techname=wildcards.techname, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS)
	output: config["STATSDATADIR"] + "{techname}.chimeric.reads.stats.tsv"
	shell:
		'''
cat {input} |sort > {output}
		'''

rule plotRateOfReadChimeras:
	input: config["STATSDATADIR"] + "{techname}.chimeric.reads.stats.tsv"
	output: config["PLOTSDIR"] + "{techname}.chimeric.reads.stats.{ext}"
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
cat {output}.r | R --slave
#dropbox_uploader.sh upload {output} {DROPBOX_PLOTS};

		'''

rule selectNonAmbiguousNonChimericReads:
	input:
		demulFiltered = config["DEMULTIPLEX_DIR"] + "discardErroneousandBarcodesAndUP/{techname}_{capDesign}_{sizeFrac}.demultiplexSameCapDesignBarcode.tsv",
		chim = config["DEMULTIPLEX_DIR"] + "detectReadChimeras/{techname}_{capDesign}_{sizeFrac}.fastq.chimeras.list",
		ambig = config["DEMULTIPLEX_DIR"] + "discardAmbiguousBarcodes/{techname}_{capDesign}_{sizeFrac}.ambiguousDemultiplex.list"
	output: temp(config["DEMULTIPLEX_DIR"] + "selectNonAmbiguousNonChimericReads/{techname}_{capDesign}_{sizeFrac}.demultiplex.noAmbig.noChim.tsv")
	shell:
		'''
cat {input.chim} {input.ambig} | sort|uniq | fgrep -w -v -f - {input.demulFiltered} |sort|uniq> {output}

		'''


rule getNonAmbiguousNonChimericReadsStats:
	input:
		totals=config["STATSDATADIR"] + "{techname}.fastq.readCounts.tsv",
		noAmbigNoChim= config["DEMULTIPLEX_DIR"] + "selectNonAmbiguousNonChimericReads/{techname}_{capDesign}_{sizeFrac}.demultiplex.noAmbig.noChim.tsv"
	output: temp(config["STATSDATADIR"] + "{techname}_{capDesign}_{sizeFrac}.finalDemul.reads.stats.tsv")
	shell:
		'''
total=$(cat {input.totals} | awk -v f={wildcards.techname}_{wildcards.capDesign}_{wildcards.sizeFrac}.fastq.gz '$1==f' | cut -f2)
finalDemul=$(cat {input.noAmbigNoChim} | cut -f1 |sort|uniq|wc -l)
echo -e "{wildcards.capDesign}\t{wildcards.sizeFrac}\t$total\t$finalDemul" | awk '{{print $0"\t"$4/$3}}' > {output}

		'''

rule aggNonAmbiguousNonChimericReadsStats:
	input: #expand(config["STATSDATADIR"] + "{{techname}}_{capDesign}_{sizeFrac}.finalDemul.reads.stats.tsv", capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS)
		lambda wildcards: expand(config["STATSDATADIR"] + "{techname}_{capDesign}_{sizeFrac}.finalDemul.reads.stats.tsv", filtered_product, techname=wildcards.techname, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS)
	output: config["STATSDATADIR"] + "{techname}.finalDemul.reads.stats.tsv"
	shell:
		'''
cat {input} | sort > {output}
		'''

rule plotNonAmbiguousNonChimericReadsStats:
	input: config["STATSDATADIR"] + "{techname}.finalDemul.reads.stats.tsv"
	output: config["PLOTSDIR"] + "{techname}.finalDemul.reads.stats.{ext}"
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
cat {output}.r | R --slave
#dropbox_uploader.sh upload {output} {DROPBOX_PLOTS};

		'''


rule checkForOnlyOneBarcodePerRead:
	input:
		config["DEMULTIPLEX_DIR"] + "selectNonAmbiguousNonChimericReads/{techname}_{capDesign}_{sizeFrac}.demultiplex.noAmbig.noChim.tsv"
	output: config["DEMULTIPLEX_DIR"] + "qc/{techname}_{capDesign}_{sizeFrac}.demul.QC1.txt"
	shell:
		'''
cut -f1,5 {input} | sort|uniq | cut -f1 |sort|uniq -d |wc -l > {output}
cat {output} | while read count; do if [ $count -gt 0 ]; then echo "$count duplicates found";mv {output} {output}.tmp;  exit 1; fi; done
		'''

rule convertFastqToTsv:
	input: config["FQPATH"] + "{techname}_{capDesign}_{sizeFrac}.fastq.gz"
	output: temp(config["DEMULTIPLEX_DIR"] + "convertFastqToTsv/{techname}_{capDesign}_{sizeFrac}.tsv")
	shell:
		'''
zcat  {input} | fastq2tsv.pl > {output}
		'''

rule demultiplexFastqs:
	input:
		tsvFastq = config["DEMULTIPLEX_DIR"] + "convertFastqToTsv/{techname}_{capDesign}_{sizeFrac}.tsv",
		noAmbigNoChim= config["DEMULTIPLEX_DIR"] + "selectNonAmbiguousNonChimericReads/{techname}_{capDesign}_{sizeFrac}.demultiplex.noAmbig.noChim.tsv"
	output: config["DEMULTIPLEX_DIR"] + "demultiplexFastqs/{techname}_{capDesign}_{sizeFrac}.{barcodes}.fastq.gz"
	shell:
		'''
cat {input.noAmbigNoChim} | awk -v b={wildcards.barcodes} '$5==b'| cut -f1|sort|uniq > $TMPDIR/tmp
set +e
fgrep -w -f $TMPDIR/tmp {input.tsvFastq} > $TMPDIR/tmp.fastq
set -e
cat $TMPDIR/tmp.fastq | tsv2fastq.pl | gzip> {output}

		'''

rule getUndeterminedReads:
	input:
		tsvFastq = config["DEMULTIPLEX_DIR"] + "convertFastqToTsv/{techname}_{capDesign}_{sizeFrac}.tsv",
		noAmbigNoChim= config["DEMULTIPLEX_DIR"] + "selectNonAmbiguousNonChimericReads/{techname}_{capDesign}_{sizeFrac}.demultiplex.noAmbig.noChim.tsv"
	output: config["DEMULTIPLEX_DIR"] + "demultiplexFastqs/{techname}_{capDesign}_{sizeFrac}.Undeter.fastq.gz"
	shell:
		'''
cat {input.noAmbigNoChim} | cut -f1 | sort | uniq | fgrep -v -w -f - {input.tsvFastq} | tsv2fastq.pl | gzip > {output}
		'''

rule checkTotalsDemultiplexed:
	input:
		demul= lambda wildcards: expand(config["DEMULTIPLEX_DIR"] + "demultiplexFastqs/{techname}_{capDesign}_{sizeFrac}.{barcodes}.fastq.gz", filtered_product, techname=wildcards.techname, capDesign=wildcards.capDesign, sizeFrac=wildcards.sizeFrac, barcodes=BARCODES),
			#expand(config["DEMULTIPLEX_DIR"] + "{{techname}}_{{capDesign}}_{{sizeFrac}}.{barcodesU}.fastq", barcodesU=BARCODESUNDETER),
		totals=config["STATSDATADIR"] + "{techname}.fastq.readCounts.tsv"
	output: config["DEMULTIPLEX_DIR"] + "qc/{techname}_{capDesign}_{sizeFrac}.demul.QC2.txt"
	shell:
		'''
total=$(cat {input.totals} | awk -v f={wildcards.techname}_{wildcards.capDesign}_{wildcards.sizeFrac}.fastq.gz '$1==f' | cut -f2)
for file in `echo {input.demul}`; do zcat $file; done > $TMPDIR/{wildcards.techname}_{wildcards.capDesign}_{wildcards.sizeFrac}.demul.fastq
demul=$(cat $TMPDIR/{wildcards.techname}_{wildcards.capDesign}_{wildcards.sizeFrac}.demul.fastq | fastq2tsv.pl | wc -l)
echo -e "{wildcards.capDesign}.{wildcards.sizeFrac}\t$total\t$demul" | awk '{{print $2-$3}}' > {output}
cat {output}| while read diff; do if [ $diff -lt 1 ]; then echo "Number of demultiplexed reads is greater than number of input reads (diff: $diff)";mv {output} {output}.tmp;  exit 1; fi; done

		'''



rule getDemultiplexingStatsPerSample:
	input: config["DEMULTIPLEX_DIR"] + "demultiplexFastqs/{techname}_{capDesign}_{sizeFrac}.{barcodes}.fastq.gz"
	output: temp(config["STATSDATADIR"] + "{techname}_{capDesign}_{sizeFrac}.{barcodes}.demultiplexing.perSample.stats.tsv")
	shell:
		'''
demul=$(zcat {input} | fastq2tsv.pl | wc -l)
echo -e "{wildcards.capDesign}\t{wildcards.sizeFrac}\t{wildcards.barcodes}\t$demul" > {output}

		'''

rule aggDemultiplexingStatsPerSample:
	input: #expand(config["STATSDATADIR"] + "{{techname}}_{capDesign}_{sizeFrac}.{barcodes}.demultiplexing.perSample.stats.tsv", capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES)
		lambda wildcards: expand(config["STATSDATADIR"] + "{techname}_{capDesign}_{sizeFrac}.{barcodes}.demultiplexing.perSample.stats.tsv", filtered_product, techname=wildcards.techname, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES)
	output: config["STATSDATADIR"] + "{techname}.demultiplexing.perSample.stats.tsv"
	shell:
		'''
cat {input} | sort | perl -slane 'if($F[2]=~/$F[0]/ || $F[2] eq "Undeter"){{print}}'> {output}

		'''
rule plotDemultiplexingStatsPerSample:
	input: config["STATSDATADIR"] + "{techname}.demultiplexing.perSample.stats.tsv"
	output: config["PLOTSDIR"] + "{techname}.demultiplexing.perSample.stats.{ext}"
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
cat {output}.r | R --slave
#dropbox_uploader.sh upload {output} {DROPBOX_PLOTS};

		'''




