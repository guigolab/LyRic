import glob
from collections import defaultdict
import os.path
from itertools import product
import sys


print ("TODO:\n ## for UCSC, filter out non-chr lines and gene_id attributes (some are ultra long and lead to buffer overflow): | fgrep chr | perl -slane ' $_=~s/gene_id \S+//g; print' \n ## include contents of README_correction_error_rate.sh into the snakemake workflow (see plotPolyAreadsStats for inspiration)")


GGPLOT_PUB_QUALITY=config["GGPLOT_PUB_QUALITY"]
CAPDESIGNTOGENOME=config["capDesignToGenome"]
CAPDESIGNTOANNOTGTF=config["capDesignToAnnotGtf"]
sizeFrac_Rpalette=config["SIZEFRACTION_RPALETTE"]
# no underscores allowed in wildcards, to avoid greedy matching since we use them as separators
wildcard_constraints:
 	capDesign = "[^_/]+",
 	sizeFrac = "[^_/]+",
 	techname = "[^_/]+",
# 	barcodes = "[^\.]",
# 	barcodesU ="^(merged)"

# get CAPDESIGNS (capture designs, i.e. Hv1, Hv2, Mv1) and SIZEFRACS (size fractions) variables from FASTQ file names (warning: this will generate duplicate entries so "set" them afterwards):
(TECHNAMES, CAPDESIGNS, SIZEFRACS) = glob_wildcards(config["FQPATH"] + "{techname}_{capDesign}_{sizeFrac}.fastq.gz")
# remove duplicate entries:
CAPDESIGNS=set(CAPDESIGNS)
SIZEFRACS=set(SIZEFRACS)
TECHNAMES=set(TECHNAMES)
adaptersTSV = config["adaptersTSV"]
f = open(adaptersTSV, 'r')
BARCODES = []
BARCODESUNDETER = []
for line in f:
	columns = line.split("\t")
	barcodeId = columns[0].split("_")
	if columns[0] != 'UP':
		BARCODES.append(columns[0])
		BARCODESUNDETER.append(columns[0])

BARCODESUNDETER.append("Undeter")
BARCODES=set(BARCODES)
BARCODESUNDETER=set(BARCODESUNDETER)

SPLICE_SITE_TYPES=["Donor", "Acceptor"]


# ### list of chromosomes for each genome
# GENOMES=[]

# for k, v in config["capDesignToGenome"].items():
# 	GENOMES.append(v)
# GENOMES=set(GENOMES)

# GENOMECHROMS=[]
# for genome in GENOMES:
# 	print(genome)
# 	f = open(config["GENOMESDIR"] + genome + ".genome", 'r')
# 	for line in f:
# 		columns = line.split("\t")
# 		#print("\t", columns[0])
# 		GENOMECHROMS.append(columns[0])

# GENOMECHROMS=set(GENOMECHROMS)
#GENCOMECHROMS contains full list of chr for all genomes. this is not optimal.
#print(GENOMECHROMS)

########################
### make list of authorized wildcard combinations
### inspired by https://stackoverflow.com/questions/41185567/how-to-use-expand-in-snakemake-when-some-particular-combinations-of-wildcards-ar
###

AUTHORIZEDCOMBINATIONS = []

for comb in product(TECHNAMES,CAPDESIGNS,SIZEFRACS,BARCODES):
	if(os.path.isfile(config["FQPATH"] + comb[0] + "_" + comb[1] + "_" + comb[2] + ".fastq.gz")): #allow only combinations corresponding to existing FASTQs
		tup=(("techname", comb[0]),("capDesign", comb[1]),("sizeFrac", comb[2]))
		AUTHORIZEDCOMBINATIONS.append(tup)
		for ext in config["PLOTFORMATS"]:
			tup2=(("techname", comb[0]),("capDesign", comb[1]),("sizeFrac", comb[2]),("ext",ext))
			AUTHORIZEDCOMBINATIONS.append(tup2)
#		print(comb[3], comb[1], sep=",")
		if (comb[3]).find(comb[1]) == 0: # check that barcode ID's prefix matches capDesign (i.e. ignore demultiplexed FASTQs with unmatched barcode/capDesign)
#			print ("MATCH")
			tup3=(("techname", comb[0]),("capDesign", comb[1]),("sizeFrac", comb[2]),("barcodes",comb[3]))
			tup4=(("techname", comb[0]),("capDesign", comb[1]),("barcodes",comb[3]))
			AUTHORIZEDCOMBINATIONS.append(tup3)
			AUTHORIZEDCOMBINATIONS.append(tup4)
			for strand in config["STRANDS"]:
				tup5=(("techname", comb[0]),("capDesign", comb[1]),("barcodes",comb[3]),("strand", strand))
				AUTHORIZEDCOMBINATIONS.append(tup5)
			# for chrom in GENOMECHROMS:
			# 	tup6=(("techname", comb[0]),("capDesign", comb[1]),("barcodes",comb[3]),("chrom", chrom))
			# 	AUTHORIZEDCOMBINATIONS.append(tup6)
			# 	tup7=(("techname", comb[0]),("capDesign", comb[1]),("chrom", chrom))
			# 	AUTHORIZEDCOMBINATIONS.append(tup7)


AUTHORIZEDCOMBINATIONS=set(AUTHORIZEDCOMBINATIONS)

#print ("AUTHORIZEDCOMBINATIONS:", AUTHORIZEDCOMBINATIONS)

def filtered_product(*args):
	for wc_comb in product(*args):
		found=False
#		print (wc_comb)
		if wc_comb in AUTHORIZEDCOMBINATIONS:
#			print ("AUTH")
			found=True
			yield(wc_comb)
#		if not found:
#			print ("NONAUTH")

include: "demultiplex.snakefile.py"
include: "fastqStats.snakefile.py"
include: "lrMapping.snakefile.py"
include: "srMapping.snakefile.py"
include: "polyAmapping.snakefile.py"
include: "introns.snakefile.py"
include: "processReadMappings.snakefile.py"
include: "tmClassification.snakefile.py"

#pseudo-rule specifying the target files we ultimately want.
rule all:
	input:
		expand(config["PLOTSDIR"] + "{techname}_{capDesign}_all.readlength.{ext}", techname=TECHNAMES, capDesign=CAPDESIGNS, ext=config["PLOTFORMATS"]), # facetted histograms of read length
		expand(config["FQPATH"] + "qc/{techname}_{capDesign}_{sizeFrac}.dupl.txt", filtered_product,techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS),
 		expand(config["PLOTSDIR"] + "{techname}.fastq.UP.stats.{ext}", techname=TECHNAMES, ext=config["PLOTFORMATS"]), # UP reads plots
 		expand(config["PLOTSDIR"] + "{techname}.fastq.BC.stats.{ext}", techname=TECHNAMES, ext=config["PLOTFORMATS"]), # barcode reads plots
 		expand(config["PLOTSDIR"] + "{techname}.fastq.foreignBC.stats.{ext}", techname=TECHNAMES, ext=config["PLOTFORMATS"]), #foreign barcode reads plots
 		expand ("mappings/" + "hiSeq_{capDesign}.bam", capDesign=CAPDESIGNS),  # mapped short reads
 #		expand ("mappings/" + "mergeSizeFracBams/{techname}_{capDesign}_{barcodes}.merged.bam", filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, barcodes=BARCODES),  # mapped reads
 		expand(config["DEMULTIPLEX_DIR"] + "demultiplexFastqs/{techname}_{capDesign}_{sizeFrac}.{barcodes}.fastq.gz", filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES),
		expand ("mappings/" + "clusterPolyAsites/{techname}_{capDesign}_{barcodes}.polyAsites.clusters.bed", filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, barcodes=BARCODES),
#		expand("mappings/" + "getPolyAreadsList/{techname}_{capDesign}_{barcodes}.polyAreads.list", filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, barcodes=BARCODES),
#		expand("mappings/" + "getIntronMotif/{techname}_{capDesign}_{barcodes}.introns.gff", filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, barcodes=BARCODES),
#		expand("mappings/" + "getIntronMotif/{techname}_{capDesign}_{barcodes}.transcripts.tsv", filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, barcodes=BARCODES),
#		expand("mappings/" + "highConfidenceReads/{techname}_{capDesign}_{barcodes}.strandedHCGMs.gff", filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, barcodes=BARCODES),
		#expand("mappings/" + "nonAnchoredMergeReads/{techname}_{capDesign}_{barcodes}.compmerge.gff", filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, barcodes=BARCODES),
		#expand("mappings/" + "nonAnchoredMergeReads/{techname}_{capDesign}_{barcodes}.tmerge.gff", filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, barcodes=BARCODES),
		#expand("mappings/" + "highConfidenceReads/tmp/{techname}_{capDesign}_{barcodes}.strandedHCGMs.{chrom}.gff", filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, barcodes=BARCODES, chrom=GENOMECHROMS),
		#expand("mappings/" + "nonAnchoredMergeReads/chr/{techname}_{capDesign}_{barcodes}.tmerge.{chrom}.gff", filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, barcodes=BARCODES, chrom=GENOMECHROMS),
		#expand("mappings/" + "nonAnchoredMergeReads/{techname}_{capDesign}_{barcodes}.tmerge.gff", filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, barcodes=BARCODES),
		expand(config["PLOTSDIR"] + "all.pooled.merged.HiSS.stats.{ext}", ext=config["PLOTFORMATS"]),
		expand("mappings/" + "mergeCapDesignBams/{techname}_{capDesign}.merged2.bam", techname=TECHNAMES, capDesign=CAPDESIGNS),
		expand("mappings/" + "nonAnchoredMergeReads/pooled/{techname}_{capDesign}_pooled.tmerge.gff", techname=TECHNAMES, capDesign=CAPDESIGNS),
		expand("mappings/" + "nonAnchoredMergeReads/qc/{techname}_{capDesign}_{barcodes}.tmerge.qc.txt",filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, barcodes=BARCODES),
		expand( "mappings/" + "nonAnchoredMergeReads/pooled/qc/{techname}_{capDesign}.tmerge.qc.txt", techname=TECHNAMES, capDesign=CAPDESIGNS),
		expand(config["PLOTSDIR"] + "all.pooled.merged.stats.{ext}", ext=config["PLOTFORMATS"]),
#		expand(config["PLOTSDIR"] + "all.pooled.merged.HiSS.stats.{ext}", ext=config["PLOTFORMATS"]),
		expand(config["PLOTSDIR"] + "all.polyAreads.stats.{ext}", ext=config["PLOTFORMATS"]),
		#expand(config["PLOTSDIR"] + "all.HCGMs.stats.{ext}", ext=config["PLOTFORMATS"]),
#		expand(config["PLOTSDIR"] + "all.pooled.merged.HiSS.splicedLength.stats.{ext}", ext=config["PLOTFORMATS"]),
		 expand("mappings/" + "makePolyABigWigs/{techname}_{capDesign}_{barcodes}.polyAsitesNoErcc.{strand}.bw", filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, barcodes=BARCODES, strand=config["STRANDS"]),

 		expand("mappings/" + "qc/{techname}_{capDesign}_{barcodes}.merged.bam.dupl.txt",filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, barcodes=BARCODES),
# 		expand ("mappings/" + "readBedToGff/{techname}_{capDesign}_{barcodes}.merged.gff",filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, barcodes=BARCODES),
 		expand("mappings/qualimap_reports/" + "{techname}_{capDesign}.merged2/genome_results.txt", techname=TECHNAMES, capDesign=CAPDESIGNS),
 		expand(config["PLOTSDIR"] + "{techname}.ambiguousBarcodes.reads.stats.{ext}", techname=TECHNAMES, ext=config["PLOTFORMATS"]), # ambiguous barcodes plots
 		expand(config["PLOTSDIR"] + "{techname}_{capDesign}.adapters.location.stats.{ext}",techname=TECHNAMES, capDesign=CAPDESIGNS, ext=config["PLOTFORMATS"]), #location of adapters over reads
 		expand(config["PLOTSDIR"] + "{techname}.chimeric.reads.stats.{ext}", techname=TECHNAMES, ext=config["PLOTFORMATS"]), # stats on chimeric reads
 		expand(config["PLOTSDIR"] + "{techname}.finalDemul.reads.stats.{ext}", techname=TECHNAMES, ext=config["PLOTFORMATS"]), #final demultiplexing stats
 		expand(config["DEMULTIPLEX_DIR"] + "qc/{techname}_{capDesign}_{sizeFrac}.demul.QC1.txt",filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS), # QC on demultiplexing (checks that there is only one barcode assigned per read
 		expand(config["DEMULTIPLEX_DIR"] + "qc/{techname}_{capDesign}_{sizeFrac}.demul.QC2.txt", filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS), # QC on demultiplexing (checks that there is only one barcode assigned per read
 		expand(config["PLOTSDIR"] + "{techname}.demultiplexing.perSample.stats.{ext}", techname=TECHNAMES, ext=config["PLOTFORMATS"]),
 		expand(config["PLOTSDIR"] + "{techname}.mapping.perSample.perFraction.stats.{ext}", techname=TECHNAMES, ext=config["PLOTFORMATS"]),
 		expand("mappings/" + "nonAnchoredMergeReads/vsTargets/{techname}_{capDesign}_pooled.tmerge.gfftsv.gz", techname=TECHNAMES, capDesign=CAPDESIGNS),
 		expand(config["PLOTSDIR"] + "all_pooled.targetCoverage.stats.{ext}", ext=config["PLOTFORMATS"]),
 		expand(config["PLOTSDIR"] + "all.splice.sites.stats.{ext}", ext=config["PLOTFORMATS"]),
 		expand(config["PLOTSDIR"] + "all.pooled.tmerge.vs.gencode.stats.{ext}",  ext=config["PLOTFORMATS"]),
