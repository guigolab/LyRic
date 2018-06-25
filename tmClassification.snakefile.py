rule compareTargetsToTms:
	input:
		tms= "mappings/" + "nonAnchoredMergeReads/pooled/{techname}_{capDesign}_pooled.tmerge.gff",
		targetedSegments=config["TARGETSDIR"] + "{capDesign}_primary_targets.exons.reduced.gene_type.segments.gtf"
	output: "mappings/" + "nonAnchoredMergeReads/vsTargets/{techname}_{capDesign}_pooled.tmerge.gfftsv.gz"
	shell:
		'''
bedtools intersect -wao -a {input.targetedSegments} -b {input.tms} |gzip > {output}

		'''

rule getTargetCoverageStats:
	input: "mappings/" + "nonAnchoredMergeReads/vsTargets/{techname}_{capDesign}_pooled.tmerge.gfftsv.gz"
	output: config["STATSDATADIR"] + "{techname}_{capDesign}_pooled.targetCoverage.stats.tsv"
	shell:
		'''
for type in `zcat {input} | extractGffAttributeValue.pl gene_type | sort|uniq`; do
all=$(zcat {input} | fgrep "gene_type \\"$type\\";" | extractGffAttributeValue.pl transcript_id | sort|uniq|wc -l)
detected=$(zcat {input} | fgrep "gene_type \\"$type\\";" | awk '$NF>0' | extractGffAttributeValue.pl transcript_id | sort|uniq|wc -l)
let undetected=$all-$detected || true
echo -e "{wildcards.techname}\t{wildcards.capDesign}\t$type\t$all\t$detected" | awk '{{print $0"\t"$5/$4}}'
done > {output}
		'''

rule aggTargetCoverageStats:
	input: lambda wildcards: expand(config["STATSDATADIR"] + "{techname}_{capDesign}_pooled.targetCoverage.stats.tsv", techname=TECHNAMES, capDesign=CAPDESIGNS)
	output: config["STATSDATADIR"] + "all_pooled.targetCoverage.stats.tsv"
	shell:
		'''
echo -e "seqTech\tcorrectionLevel\tcapDesign\ttargetType\ttotalTargets\tdetectedTargets\tpercentDetectedTargets" > {output}
cat {input} | grep -v erccSpikein | sed 's/Corr0/\tNo/' | sed 's/Corr90/\tYes/'| sort >> {output}
		'''

rule plotTargetCoverageStats:
	input: config["STATSDATADIR"] + "all_pooled.targetCoverage.stats.tsv"
	output: config["PLOTSDIR"] + "all_pooled.targetCoverage.stats.{ext}"
	shell:
		'''
echo "library(ggplot2)
library(plyr)
library(scales)
dat <- read.table('{input}', header=T, as.is=T, sep='\\t')
#nrTypes <- length(unique(dat\$targetType))
#darkcols <- brewer.pal(nrTypes, "Dark2")
ggplot(dat, aes(x=factor(correctionLevel), y=percentDetectedTargets, fill=targetType)) +
geom_bar(width=0.75,stat='identity', position=position_dodge(width=0.9)) +
scale_fill_brewer(palette='Set3') +
 facet_grid( seqTech ~ capDesign) +
 geom_hline(aes(yintercept=1), linetype='dashed', alpha=0.7) +
 geom_text(aes(group=targetType, y=0.01, label = paste(sep='',percent(percentDetectedTargets),' / ','(',comma(detectedTargets),')')), angle=90, size=2.5, hjust=0, vjust=0.5, position = position_dodge(width=0.9)) +
 ylab('% targeted regions detected') + xlab('Error correction') + scale_y_continuous(limits = c(0, 1), labels = scales::percent)+
{GGPLOT_PUB_QUALITY}
ggsave('{output}', width=8, height=4)
" > {output}.r
cat {output}.r | R --slave


		'''


rule gffcompareToAnnotation:
	input:
		annot=lambda wildcards: CAPDESIGNTOANNOTGTF[wildcards.capDesign],
		tm="mappings/" + "nonAnchoredMergeReads/pooled/{techname}_{capDesign}_pooled.tmerge.gff"
	output: "mappings/" + "nonAnchoredMergeReads/pooled/gffcompare/{techname}_{capDesign}_pooled.tmerge.vs.gencode.simple.tsv"
	shell:
		'''
pref=$(basename {output} .simple.tsv)
annotFullPath=$(fullpath {input.annot})
tmFullPath=$(fullpath {input.tm})
cd $(dirname {output})
gffcompare -T -o $pref -r $annotFullPath $tmFullPath
cat $pref.tracking | simplifyGffCompareClasses.pl - > $(basename {output})

		'''

rule getGffCompareStats:
	input: "mappings/" + "nonAnchoredMergeReads/pooled/gffcompare/{techname}_{capDesign}_pooled.tmerge.vs.gencode.simple.tsv"
	output: config["STATSDATADIR"] + "{techname}_{capDesign}_pooled.tmerge.vs.gencode.stats.tsv"
	shell:
		'''
cat {input} |cut -f4 | sort|uniq -c | awk -v s={wildcards.techname} -v c={wildcards.capDesign} '{{print s"\t"c"\t"$2"\t"$1}}' | sed 's/Corr0/\tNo/' | sed 's/Corr90/\tYes/'>> {output}
		'''


#echo -e "seqTech\tcorrectionLevel\tcapDesign\tcategory\tcount" > {output}

rule aggGffCompareStats:
	input: lambda wildcards: expand(config["STATSDATADIR"] + "{techname}_{capDesign}_pooled.tmerge.vs.gencode.stats.tsv", techname=TECHNAMES, capDesign=CAPDESIGNS)
	output: config["STATSDATADIR"] + "all.pooled.tmerge.vs.gencode.stats.tsv"
	shell:
		'''
echo -e "seqTech\tcorrectionLevel\tcapDesign\tcategory\tcount" > {output}
cat {input} >> {output}
		'''

rule plotGffCompareStats:
	input: config["STATSDATADIR"] + "all.pooled.tmerge.vs.gencode.stats.tsv"
	output: config["PLOTSDIR"] + "all.pooled.tmerge.vs.gencode.stats.{ext}"
	shell:
		'''
echo "library(ggplot2)
library(plyr)
library(scales)

dat <- read.table('{input}', header=T, as.is=T, sep='\\t')
dat\$category<-factor(dat\$category, ordered=TRUE, levels=rev(c('Intergenic', 'Extends', 'Intronic', 'Overlaps', 'Antisense', 'Equal', 'Included')))
palette <- c('Intergenic' = '#0099cc', 'Extends' ='#00bfff', 'Intronic' = '#4dd2ff', 'Overlaps' = '#80dfff', 'Antisense' = '#ccf2ff', 'Equal' = '#c65353', 'Included' ='#d98c8c')
ggplot(dat[order(dat\$category), ], aes(x=factor(correctionLevel), y=count, fill=category)) +
geom_bar(stat='identity') +
scale_fill_manual(values=palette) +
facet_grid( seqTech ~ capDesign)+ ylab('# TMs') + xlab('Error correction') + guides(fill = guide_legend(title='Category'))+
scale_y_continuous(labels=scientific)+
{GGPLOT_PUB_QUALITY}
ggsave('{output}', width=7, height=3)
" > {output}.r
cat {output}.r | R --slave

		'''

#calculate amount of novel nucleotides
