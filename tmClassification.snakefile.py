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
cat {input} | grep -v erccSpikein | sed 's/Corr/\t/'| sort >> {output}
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
 geom_text(aes(group=targetType, y=0.01, label = paste(sep='',percent(percentDetectedTargets),' / ','(',comma(detectedTargets),')')), angle=90, size=1, hjust=0, vjust=0.5, position = position_dodge(width=0.9)) +
 ylab('% targeted regions detected') + xlab('Correction level (k-mer size)') + scale_y_continuous(limits = c(0, 1), labels = scales::percent)+
{GGPLOT_PUB_QUALITY}
ggsave('{output}', width=8, height=3)
" > {output}.r
cat {output}.r | R --slave


		'''


#calculate amount of novel nucleotides