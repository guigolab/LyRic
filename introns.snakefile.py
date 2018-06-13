rule makeIntrons:
	input: "mappings/" + "readBedToGff/{techname}_{capDesign}_{barcodes}.merged.gff"
	output: "mappings/" + "makeIntrons/{techname}_{capDesign}_{barcodes}.introns.gff.gz"
	shell:
		'''
cat {input} | sort -k12,12 -k4,4n -k5,5n | awk -v fldgn=10 -v fldtr=12 -f ~/julien_utils/make_introns.awk | sortgff |gzip> {output}
		'''

rule getIntronMotif:
	input:
		introns = "mappings/" + "makeIntrons/{techname}_{capDesign}_{barcodes}.introns.gff.gz",
		genome = lambda wildcards: config["GENOMESDIR"] + CAPDESIGNTOGENOME[wildcards.capDesign] + ".fa"
	output:
		gff = "mappings/" + "getIntronMotif/{techname}_{capDesign}_{barcodes}.introns.gff",
		tsv = "mappings/" + "getIntronMotif/{techname}_{capDesign}_{barcodes}.transcripts.tsv"
	shell:
		'''
zcat {input.introns} | grep -vP "^ERCC"| extract_intron_strand_motif.pl - {input.genome} $(dirname {output.gff})/$(basename {output.gff} .introns.gff)

		'''

rule getGencodePcgSpliceSites:
	input: lambda wildcards: CAPDESIGNTOANNOTGTF[wildcards.capDesign]
	output: "annotations/" + "spliceSites/{capDesign}.gencode.PCG.spliceSites.tsv.gz"
	shell:
		'''
cat {input} | awk '$3=="exon"' | fgrep "transcript_type \\"protein_coding\\";" | sort -k12,12 -k4,4n -k5,5n | awk -v fldgn=10 -v fldtr=12 -f ~/julien_utils/make_introns.awk | awk '{{print $1"\t"$4-1"\t"$5+1"\t"$7}}' |sort|uniq | perl -lane 'if($F[3] eq "+"){{$dStart=$F[1]-1; $aStart=$F[2]-2}} elsif($F[3] eq "-"){{$dStart=$F[2]-2; $aStart=$F[1]-1}} else{{die}} $aEnd=$aStart+2; $dEnd=$dStart+2; print "$F[0]"."_$dStart"."_$dEnd"."_$F[3]\t$F[0]"."_"."$F[1]"."_"."$F[2]"."_"."$F[3]:Donor:GENCODE_protein_coding"; print "$F[0]"."_$aStart"."_$aEnd"."_$F[3]\t$F[0]"."_"."$F[1]"."_"."$F[2]"."_"."$F[3]:Acceptor:GENCODE_protein_coding";' | sort |uniq| sort -k1,1 | gzip > {output}
		'''

rule getReadsSpliceSites:
	input: lambda wildcards: expand("mappings/" + "strandGffs/{techname}_{capDesign}_{barcodes}.stranded.gff", filtered_product, techname=wildcards.techname, capDesign=wildcards.capDesign, barcodes=BARCODES)
	output: "mappings/" + "makeIntrons/readSpliceSites/{techname}_{capDesign}.introns.tsv.gz"
	shell:
		'''
cat {input} | awk '$3=="exon"' | sort -k12,12 -k4,4n -k5,5n | awk -v fldgn=10 -v fldtr=12 -f ~/julien_utils/make_introns.awk| awk '{{print $1"\t"$4-1"\t"$5+1"\t"$7}}' |sort|uniq | perl -lane 'if($F[3] eq "+"){{$dStart=$F[1]-1; $aStart=$F[2]-2}} elsif($F[3] eq "-"){{$dStart=$F[2]-2; $aStart=$F[1]-1}} else{{die}} $aEnd=$aStart+2; $dEnd=$dStart+2; print "$F[0]"."_$dStart"."_$dEnd"."_$F[3]\t$F[0]"."_"."$F[1]"."_"."$F[2]"."_"."$F[3]:Donor:CLS_reads"; print "$F[0]"."_$aStart"."_$aEnd"."_$F[3]\t$F[0]"."_"."$F[1]"."_"."$F[2]"."_"."$F[3]:Acceptor:CLS_reads";' |sort|uniq| sort -k1,1 |gzip > {output}
		'''

rule getTmSpliceSites:
	input: "mappings/" + "nonAnchoredMergeReads/pooled/{techname}_{capDesign}_pooled.tmerge.gff"
	output: "mappings/" + "nonAnchoredMergeReads/pooled/spliceSites/{techname}_{capDesign}.introns.tsv.gz"
	shell:
		'''
cat {input} | awk '$3=="exon"' | sort -k12,12 -k4,4n -k5,5n | awk -v fldgn=10 -v fldtr=12 -f ~/julien_utils/make_introns.awk | awk '{{print $1"\t"$4-1"\t"$5+1"\t"$7}}' |sort|uniq | perl -lane 'if($F[3] eq "+"){{$dStart=$F[1]-1; $aStart=$F[2]-2}} elsif($F[3] eq "-"){{$dStart=$F[2]-2; $aStart=$F[1]-1}} else{{die}} $aEnd=$aStart+2; $dEnd=$dStart+2; print "$F[0]"."_$dStart"."_$dEnd"."_$F[3]\t$F[0]"."_"."$F[1]"."_"."$F[2]"."_"."$F[3]:Donor:CLS_TMs"; print "$F[0]"."_$aStart"."_$aEnd"."_$F[3]\t$F[0]"."_"."$F[1]"."_"."$F[2]"."_"."$F[3]:Acceptor:CLS_TMs";' | sort |uniq| sort -k1,1 | gzip> {output}
		'''

rule getGeneidScores:
	input:
		gencode="annotations/" + "spliceSites/{capDesign}.gencode.PCG.spliceSites.tsv.gz",
		rawReads="mappings/" + "makeIntrons/readSpliceSites/{techname}_{capDesign}.introns.tsv.gz",
		tmReads="mappings/" + "nonAnchoredMergeReads/pooled/spliceSites/{techname}_{capDesign}.introns.tsv.gz",
		geneidScores= lambda wildcards: config["SPLICE_SITE_SCORES_DIR"] + CAPDESIGNTOGENOME[wildcards.capDesign] + ".geneid.loose.spliceSites.sorted.tsv.gz",
		random=lambda wildcards: "annotations/" + "spliceSites/" + CAPDESIGNTOGENOME[wildcards.capDesign] + ".geneid.loose.spliceSites.OnDetectedRegions.NotDetected.random1M.tsv"
	output:  config["STATSDATADIR"] + "{techname}_{capDesign}.splice.sites.stats.tsv"
	shell:
		'''
zcat {input.gencode} {input.rawReads} {input.tmReads} | sort  -k1,1 > $TMPDIR/cls.SSs.tsv
zcat {input.geneidScores} > $TMPDIR/all.SSs.tsv
join -j1 $TMPDIR/cls.SSs.tsv $TMPDIR/all.SSs.tsv |ssv2tsv | grep -P "(Acceptor.+Acceptor)|(Donor.+Donor)" | awk -v t={wildcards.techname} -v c={wildcards.capDesign} '{{print t"\t"c"\t"$1"\t"$2"\t"$4}}'| sed 's/Corr/\t/' | sed 's/:/\t/g'|cut -f1-3,5-7> {output}
cat {input.random} | awk -v t={wildcards.techname} -v c={wildcards.capDesign} '{{print t"\t"c"\t"$1"\t"$2"\trandom\t"$3}}' | sed 's/Corr/\t/' >> {output}

#verify that all SSs were found:
cut -f1 $TMPDIR/cls.SSs.tsv | sort|uniq > $TMPDIR/cls.SSs.list
cut -f1 {output} | sort|uniq > $TMPDIR/cls.SSs.geneid.list
diff=$(diff -q $TMPDIR/cls.SSs.list $TMPDIR/cls.SSs.geneid.list |wc -l)
if [ ! $diff -eq 0 ]; then echo "ERROR: List of SSs differ before/after"; exit 1; fi

		'''

rule aggGeneidScores:
	input: lambda wildcards: expand(config["STATSDATADIR"] + "{techname}_{capDesign}.splice.sites.stats.tsv", techname=TECHNAMES, capDesign=CAPDESIGNS)
	output: config["STATSDATADIR"] + "all.splice.sites.stats.tsv"
	shell:
		'''
echo -e "seqTech\tcorrectionLevel\tcapDesign\tssId\tsjId\tssType\tssCategory\tssScore" > {output}
cat {input} >> {output}
		'''
