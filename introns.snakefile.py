rule makeIntrons:
	input: "mappings/" + "readBedToGff/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.gff.gz"
	output: "mappings/" + "makeIntrons/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.introns.gff.gz"
	shell:
		'''
zcat {input} | sort -k12,12 -k4,4n -k5,5n | awk -v fldgn=10 -v fldtr=12 -f ~/julien_utils/make_introns.awk | sortgff |gzip> {output}
		'''

rule getIntronMotif:
	input:
		introns = "mappings/" + "makeIntrons/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.introns.gff.gz",
		genome = lambda wildcards: config["GENOMESDIR"] + CAPDESIGNTOGENOME[wildcards.capDesign] + ".fa"
	output:
		gff = temp("mappings/" + "getIntronMotif/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.introns.gff"),
		tsv = temp("mappings/" + "getIntronMotif/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.transcripts.tsv")
	shell:
		'''
zcat {input.introns} | grep -vP "^ERCC"| extract_intron_strand_motif.pl - {input.genome} $(dirname {output.gff})/$(basename {output.gff} .introns.gff)

		'''

rule getGencodePcgSpliceSites:
	input: lambda wildcards: CAPDESIGNTOANNOTGTF[wildcards.capDesign]
	output: "annotations/" + "spliceSites/{capDesign}.gencode.PCG.spliceSites.{spliceType}.tsv.gz"
	shell:
		'''
cat {input} | awk '$3=="exon"' | fgrep "transcript_type \\"protein_coding\\";" | sort -k12,12 -k4,4n -k5,5n | awk -v fldgn=10 -v fldtr=12 -f ~/julien_utils/make_introns.awk | awk '{{print $1"\t"$4-1"\t"$5+1"\t"$7}}' |sort|uniq | perl -lane 'if($F[3] eq "+"){{$dStart=$F[1]-1; $aStart=$F[2]-2}} elsif($F[3] eq "-"){{$dStart=$F[2]-2; $aStart=$F[1]-1}} else{{die}} $aEnd=$aStart+2; $dEnd=$dStart+2; print "$F[0]"."_$dStart"."_$dEnd"."_$F[3]\t$F[0]"."_"."$F[1]"."_"."$F[2]"."_"."$F[3]:Donor:GENCODE_protein_coding"; print "$F[0]"."_$aStart"."_$aEnd"."_$F[3]\t$F[0]"."_"."$F[1]"."_"."$F[2]"."_"."$F[3]:Acceptor:GENCODE_protein_coding";' | fgrep -w {wildcards.spliceType}| sort |uniq| sort -k1,1 | gzip > {output}
		'''

rule getReadsSpliceSites:
#	input: lambda wildcards: expand("mappings/" + "strandGffs/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.stranded.gff", filtered_product, techname=wildcards.techname,  corrLevel=wildcards.corrLevel,capDesign=wildcards.capDesign, barcodes=BARCODES)
	input: "mappings/" + "strandGffs/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.stranded.gff.gz"
	output: "mappings/" + "makeIntrons/readSpliceSites/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.introns.{spliceType}.tsv.gz"
	shell:
		'''
zcat {input} | awk '$3=="exon"' | sort -k12,12 -k4,4n -k5,5n | awk -v fldgn=10 -v fldtr=12 -f ~/julien_utils/make_introns.awk| awk '{{print $1"\t"$4-1"\t"$5+1"\t"$7}}' |sort|uniq | perl -lane 'if($F[3] eq "+"){{$dStart=$F[1]-1; $aStart=$F[2]-2}} elsif($F[3] eq "-"){{$dStart=$F[2]-2; $aStart=$F[1]-1}} else{{die}} $aEnd=$aStart+2; $dEnd=$dStart+2; print "$F[0]"."_$dStart"."_$dEnd"."_$F[3]\t$F[0]"."_"."$F[1]"."_"."$F[2]"."_"."$F[3]:Donor:CLS_reads"; print "$F[0]"."_$aStart"."_$aEnd"."_$F[3]\t$F[0]"."_"."$F[1]"."_"."$F[2]"."_"."$F[3]:Acceptor:CLS_reads";' | fgrep -w {wildcards.spliceType} |sort|uniq| sort -k1,1 |gzip > {output}
		'''

rule getTmSpliceSites:
	input: "mappings/" + "nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.gff"
	output: "mappings/" + "nonAnchoredMergeReads/spliceSites/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.introns.{spliceType}.tsv.gz"
	shell:
		'''
cat {input} | awk '$3=="exon"' | sort -k12,12 -k4,4n -k5,5n | awk -v fldgn=10 -v fldtr=12 -f ~/julien_utils/make_introns.awk | awk '{{print $1"\t"$4-1"\t"$5+1"\t"$7}}' |sort|uniq | perl -lane 'if($F[3] eq "+"){{$dStart=$F[1]-1; $aStart=$F[2]-2}} elsif($F[3] eq "-"){{$dStart=$F[2]-2; $aStart=$F[1]-1}} else{{die}} $aEnd=$aStart+2; $dEnd=$dStart+2; print "$F[0]"."_$dStart"."_$dEnd"."_$F[3]\t$F[0]"."_"."$F[1]"."_"."$F[2]"."_"."$F[3]:Donor:CLS_TMs"; print "$F[0]"."_$aStart"."_$aEnd"."_$F[3]\t$F[0]"."_"."$F[1]"."_"."$F[2]"."_"."$F[3]:Acceptor:CLS_TMs";' | fgrep -w {wildcards.spliceType} | sort |uniq| sort -k1,1 | gzip> {output}
		'''

rule getGeneidScores:
	input:
		gencode="annotations/" + "spliceSites/{capDesign}.gencode.PCG.spliceSites.{spliceType}.tsv.gz",
		rawReads="mappings/" + "makeIntrons/readSpliceSites/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.introns.{spliceType}.tsv.gz",
		tmReads="mappings/" + "nonAnchoredMergeReads/spliceSites/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.introns.{spliceType}.tsv.gz",
		geneidScores= lambda wildcards: expand(config["SPLICE_SITE_SCORES_DIR"] + CAPDESIGNTOGENOME[wildcards.capDesign] + ".{spliceType}.geneid.loose.spliceSites.sorted.tsv", spliceType=wildcards.spliceType),
		random=lambda wildcards: "annotations/" + "spliceSites/" + CAPDESIGNTOGENOME[wildcards.capDesign] + ".geneid.loose.spliceSites.OnDetectedRegions.NotDetected.random100K.{spliceType}.tsv"
	output:  temp(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.{spliceType}.spliceSites.stats.tsv")
	shell:
		'''
zcat {input.gencode} {input.rawReads} {input.tmReads} | sort  -k1,1 > $TMPDIR/cls.SSs.tsv
join -a1 -j1 $TMPDIR/cls.SSs.tsv {input.geneidScores} |ssv2tsv | perl -lane '$_=~s/\:/\t/g; @line=split("\\t", $_); unless(defined ($line[5])){{$line[4]=$line[2];$line[5]=-35 + rand(-25 +35)}}; print join("\\t", @line)' | awk -v t={wildcards.techname}Corr{wildcards.corrLevel} -v c={wildcards.capDesign} -v si={wildcards.sizeFrac} -v b={wildcards.barcodes} '{{print t"\t"c"\t"si"\t"b"\t"$1"\t"$3"\t"$4"\t"$6}}'| sed 's/Corr0/\tNo/' | sed 's/Corr{lastK}/\tYes/' |sed 's/{wildcards.capDesign}_//' > {output}

#verify that all SSs were found:
cut -f1 $TMPDIR/cls.SSs.tsv | sort|uniq > $TMPDIR/cls.SSs.list
cut -f6 {output} | sort|uniq > $TMPDIR/cls.SSs.geneid.list
diff=$(diff -q $TMPDIR/cls.SSs.list $TMPDIR/cls.SSs.geneid.list |wc -l)
if [ ! $diff -eq 0 ]; then echoerr "ERROR: List of SSs differ before/after"; exit 1; fi

cat {input.random} | awk -v t={wildcards.techname}Corr{wildcards.corrLevel} -v c={wildcards.capDesign} -v si={wildcards.sizeFrac} -v b={wildcards.barcodes} '{{print t"\t"c"\t"si"\t"b"\t"$1"\t"$2"\trandom\t"$3}}' | sed 's/Corr0/\tNo/' | sed 's/Corr{lastK}/\tYes/' >> {output}
#cat {input.random} | fgrep -w {wildcards.spliceType} |shuf | head -n100000 | awk -v t={wildcards.techname}Corr{wildcards.corrLevel} -v c={wildcards.capDesign} -v si={wildcards.sizeFrac} -v b={wildcards.barcodes} '{{print t"\t"c"\t"si"\t"b"\t"$1"\t"$2"\trandom\t"$3}}' | sed 's/Corr0/\tNo/' | sed 's/Corr{lastK}/\tYes/' >> {output}


		'''

rule aggGeneidScores:
	input: expand(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.{spliceType}.spliceSites.stats.tsv",filtered_product_merge, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACSpluSMERGED, barcodes=BARCODESpluSMERGED, spliceType=SPLICE_SITE_TYPES)
	#input: lambda wildcards: expand(config["STATSDATADIR"] + "{techname}Corr{corrLevel}_{capDesign}.{spliceType}.splice.sites.stats.tsv", techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, spliceType=SPLICE_SITE_TYPES)
	output: config["STATSDATADIR"] + "all.spliceSites.stats.tsv"
	threads: 8
	shell:
		'''
echo -e "seqTech\tcorrectionLevel\tcapDesign\tsizeFrac\ttissue\tssId\tssCategory\tssScore" > {output}
uuid=$(uuidgen)
#keep non-redundant SSs
cat {input} | cut -f1-6,8,9 > $TMPDIR/$uuid
sort -S 28G --parallel {threads} $TMPDIR/$uuid | uniq >> {output}
		'''

rule plotGeneidScores:
	input: config["STATSDATADIR"] + "all.spliceSites.stats.tsv"
	output: config["PLOTSDIR"] + "{capDesign}_{byFrac}_{byTissue}.spliceSites.stats.{ext}"
	params:
		filterDat=lambda wildcards: merge_figures_params(wildcards.capDesign, wildcards.byFrac, wildcards.byTissue)
	shell:
		'''
echo "
library(data.table)
library(ggplot2)
library(scales)
library(plyr)
palette <- c('random' = '#666666', 'GENCODE_protein_coding' = '#009900', 'CLS_reads' = '#b3d9ff', 'CLS_TMs' = '#cc9966')
dat<-fread('{input}', header=T, sep='\\t')
{params.filterDat}

fun_length <- function(x){{
return(data.frame(y=-8.5,label= paste0('N=', comma(length(x)))))
}}

ggplot(dat, aes(x=factor(correctionLevel), y=ssScore, color=ssCategory)) +
geom_boxplot(position=position_dodge(0.9), outlier.shape=NA) +
coord_cartesian(ylim=c(-9, 4.5)) +
scale_color_manual(values=palette, name='Category', labels = c(random = 'Random', GENCODE_protein_coding = 'GENCODE\nprotein-coding', CLS_TMs='CLS TMs', CLS_reads='CLS raw reads')) +
facet_grid( seqTech + sizeFrac ~ capDesign + tissue)+
stat_summary(aes(x=factor(correctionLevel)), position=position_dodge(0.9), fun.data = fun_length, geom = 'text', vjust = +1, hjust=0, angle=90, show.legend=FALSE) +
geom_hline(aes(yintercept=0), linetype='dashed', alpha=0.7)+
ylab('Splice site score') + xlab('Error correction') +
{GGPLOT_PUB_QUALITY}
ggsave('{output}', width=plotWidth, height=plotHeight)

" > {output}.r
cat {output}.r | R --slave

		'''
