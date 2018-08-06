import glob
from collections import defaultdict
import os.path
from itertools import product
import sys


print ("TODO:\n ## include contents of README_correction_error_rate.sh into the snakemake workflow (see plotPolyAreadsStats for inspiration)")


GGPLOT_PUB_QUALITY=config["GGPLOT_PUB_QUALITY"]
CAPDESIGNTOGENOME=config["capDesignToGenome"]
CAPDESIGNTOANNOTGTF=config["capDesignToAnnotGtf"]
sizeFrac_Rpalette=config["SIZEFRACTION_RPALETTE"]
long_Rpalette=config["LONG_RPALETTE"]
CAPDESIGNTOCAGEPEAKS=config["capDesignToCAGEpeaks"]
CAPDESIGNTOPAS=config["capDesignToGenomePAS"]
DUMMY_DIR="dummy/"
# no underscores allowed in wildcards, to avoid greedy matching since we use them as separators
wildcard_constraints:
 	capDesign = "[^_/]+",
 	sizeFrac = "[^_/]+",
 	techname = "[^_/]+",
# 	barcodes = "[^\.]",
# 	barcodesU ="^(merged)"


if config["DEMULTIPLEX"]:
	# get CAPDESIGNS (capture designs, i.e. Hv1, Hv2, Mv1) and SIZEFRACS (size fractions) variables from FASTQ file names (warning: this will generate duplicate entries so "set" them afterwards):
	DEMULTIPLEX_DIR=config["DEMULTIPLEX_DIR"]
	DEMULTIPLEXED_FASTQS= DEMULTIPLEX_DIR + "demultiplexFastqs/"
	(TECHNAMES, CAPDESIGNS, SIZEFRACS) = glob_wildcards(config["FQPATH"] + "{techname}_{capDesign}_{sizeFrac}.fastq.gz")
	# remove duplicate entries:
	adaptersTSV = DEMULTIPLEX_DIR + "all_adapters.tsv"
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

else:
#	DEMULTIPLEX_DIR=config["FQPATH"]
	DEMULTIPLEXED_FASTQS=config["FQPATH"]
	(TECHNAMES, CAPDESIGNS, SIZEFRACS, BARCODES) = glob_wildcards(config["FQPATH"] + "{techname}_{capDesign}_{sizeFrac}.{barcodes}.fastq.gz")
	BARCODES=set(BARCODES)
	BARCODESUNDETER=BARCODES


CAPDESIGNS=set(CAPDESIGNS)
SIZEFRACS=set(SIZEFRACS)
TECHNAMES=set(TECHNAMES)

#print (DEMULTIPLEXED_FASTQS, TECHNAMES, CAPDESIGNS, SIZEFRACS,BARCODES)

SPLICE_SITE_TYPES=["Donor", "Acceptor"]

# for polyA calling calibration:
minPolyAlength = []
for x in range(5, 30, 1):
	minPolyAlength.append(x)


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
	if(os.path.isfile(config["FQPATH"] + comb[0] + "_" + comb[1] + "_" + comb[2] + ".fastq.gz") or os.path.isfile(config["FQPATH"] + comb[0] + "_" + comb[1] + "_" + comb[2]  + "." + comb[3] + ".fastq.gz")): #allow only combinations corresponding to existing FASTQs
		tup=(("techname", comb[0]),("capDesign", comb[1]),("sizeFrac", comb[2]))
		AUTHORIZEDCOMBINATIONS.append(tup)
		for ext in config["PLOTFORMATS"]:
			tup2=(("techname", comb[0]),("capDesign", comb[1]),("sizeFrac", comb[2]),("ext",ext))
			AUTHORIZEDCOMBINATIONS.append(tup2)
#		print(comb[3], comb[1], sep=",")
		if (comb[3]).find(comb[1]) == 0 or config["DEMULTIPLEX"] is False: # check that barcode ID's prefix matches capDesign (i.e. ignore demultiplexed FASTQs with unmatched barcode/capDesign)
#			print ("MATCH")
			tup3=(("techname", comb[0]),("capDesign", comb[1]),("sizeFrac", comb[2]),("barcodes",comb[3]))
			tup4=(("techname", comb[0]),("capDesign", comb[1]),("barcodes",comb[3]))
			AUTHORIZEDCOMBINATIONS.append(tup3)
			AUTHORIZEDCOMBINATIONS.append(tup4)
			for strand in config["STRANDS"]:
				tup5=(("techname", comb[0]),("capDesign", comb[1]),("barcodes",comb[3]),("strand", strand))
				AUTHORIZEDCOMBINATIONS.append(tup5)
			for ext in config["PLOTFORMATS"]:
				tup6=(("techname", comb[0]),("capDesign", comb[1]),("barcodes",comb[3]),("ext", ext))
				AUTHORIZEDCOMBINATIONS.append(tup6)
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
		#print (wc_comb)
		if wc_comb in AUTHORIZEDCOMBINATIONS:
			#print ("AUTH")
			found=True
			yield(wc_comb)
		#if not found:
			#print ("NONAUTH")

#def fastqStatsDemul:


if config["DEMULTIPLEX"]:
	include: "demultiplex.snakefile.py"
include: "fastqStats.snakefile.py"
include: "lrMapping.snakefile.py"
include: "srMapping.snakefile.py"
include: "polyAmapping.snakefile.py"
include: "introns.snakefile.py"
include: "processReadMappings.snakefile.py"
include: "tmClassification.snakefile.py"
include: "tmEndSupport.py"

#pseudo-rule specifying the target files we ultimately want.
rule all:
	input:
		expand(config["PLOTSDIR"] + "{techname}_{capDesign}_all.readlength.{ext}", techname=TECHNAMES, capDesign=CAPDESIGNS, ext=config["PLOTFORMATS"]) if config["DEMULTIPLEX"] else expand(config["PLOTSDIR"] + "{techname}_{capDesign}.{barcodes}_all.readlength.{ext}",filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, barcodes=BARCODES, ext=config["PLOTFORMATS"]), # facetted histograms of read length
 		expand(config["FQPATH"] + "qc/{techname}_{capDesign}_{sizeFrac}.dupl.txt", filtered_product,techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS) if config["DEMULTIPLEX"] else expand(config["FQPATH"] + "qc/{techname}_{capDesign}_{sizeFrac}.{barcodes}.dupl.txt", filtered_product,techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES),
  		expand(config["PLOTSDIR"] + "{techname}.fastq.UP.stats.{ext}", techname=TECHNAMES, ext=config["PLOTFORMATS"]) if config["DEMULTIPLEX"] else expand( DUMMY_DIR + "dummy{number}.txt", number='1'), # UP reads plots
  		expand(config["PLOTSDIR"] + "{techname}.fastq.BC.stats.{ext}", techname=TECHNAMES, ext=config["PLOTFORMATS"]) if config["DEMULTIPLEX"]  else expand( DUMMY_DIR + "dummy{number}.txt", number='2') , # barcode reads plots
  		expand(config["PLOTSDIR"] + "{techname}.fastq.foreignBC.stats.{ext}", techname=TECHNAMES, ext=config["PLOTFORMATS"]) if config["DEMULTIPLEX"]  else expand( DUMMY_DIR + "dummy{number}.txt", number='3') , #foreign barcode reads plots
  		expand ("mappings/" + "hiSeq_{capDesign}.bam", capDesign=CAPDESIGNS) if config["USE_MATCHED_ILLUMINA"] else expand( DUMMY_DIR + "dummy{number}.txt", number='5') ,  # mapped short reads
  		expand(DEMULTIPLEX_DIR + "demultiplexFastqs/{techname}_{capDesign}_{sizeFrac}.{barcodes}.fastq.gz", filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES) if config["DEMULTIPLEX"]  else expand( DUMMY_DIR + "dummy{number}.txt", number='4'),
 		expand ("mappings/" + "clusterPolyAsites/{techname}_{capDesign}_{barcodes}.polyAsites.clusters.bed", filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, barcodes=BARCODES),
 		expand(config["PLOTSDIR"] + "all.pooled.merged.HiSS.stats.{ext}", ext=config["PLOTFORMATS"]),
 		expand("mappings/" + "mergeCapDesignBams/{techname}_{capDesign}.merged2.bam", techname=TECHNAMES, capDesign=CAPDESIGNS),
 		expand("mappings/" + "nonAnchoredMergeReads/pooled/{techname}_{capDesign}_pooled.tmerge.gff", techname=TECHNAMES, capDesign=CAPDESIGNS),
 		expand("mappings/" + "nonAnchoredMergeReads/qc/{techname}_{capDesign}_{barcodes}.tmerge.qc.txt",filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, barcodes=BARCODES),
 		expand( "mappings/" + "nonAnchoredMergeReads/pooled/qc/{techname}_{capDesign}.tmerge.qc.txt", techname=TECHNAMES, capDesign=CAPDESIGNS),
 		expand(config["PLOTSDIR"] + "all.pooled.merged.stats.{ext}", ext=config["PLOTFORMATS"]),
 		expand(config["PLOTSDIR"] + "all.polyAreads.stats.{ext}", ext=config["PLOTFORMATS"]),
 		 expand("mappings/" + "makePolyABigWigs/{techname}_{capDesign}_{barcodes}.polyAsitesNoErcc.{strand}.bw", filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, barcodes=BARCODES, strand=config["STRANDS"]),

  		expand("mappings/" + "qc/{techname}_{capDesign}_{barcodes}.merged.bam.dupl.txt",filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, barcodes=BARCODES),
  		expand("mappings/qualimap_reports/" + "{techname}_{capDesign}.merged2/genome_results.txt", techname=TECHNAMES, capDesign=CAPDESIGNS),
  		expand(config["PLOTSDIR"] + "{techname}.ambiguousBarcodes.reads.stats.{ext}", techname=TECHNAMES, ext=config["PLOTFORMATS"]) if config["DEMULTIPLEX"]  else expand( DUMMY_DIR + "dummy{number}.txt", number='6'), # ambiguous barcodes plots
  		expand(config["PLOTSDIR"] + "{techname}_{capDesign}.adapters.location.stats.{ext}",techname=TECHNAMES, capDesign=CAPDESIGNS, ext=config["PLOTFORMATS"]) if config["DEMULTIPLEX"]  else expand( DUMMY_DIR + "dummy{number}.txt", number='7'), #location of adapters over reads
  		expand(config["PLOTSDIR"] + "{techname}.chimeric.reads.stats.{ext}", techname=TECHNAMES, ext=config["PLOTFORMATS"]) if config["DEMULTIPLEX"]  else expand( DUMMY_DIR + "dummy{number}.txt", number='8'), # stats on chimeric reads
  		expand(config["PLOTSDIR"] + "{techname}.finalDemul.reads.stats.{ext}", techname=TECHNAMES, ext=config["PLOTFORMATS"]) if config["DEMULTIPLEX"]  else expand( DUMMY_DIR + "dummy{number}.txt", number='8'), #final demultiplexing stats
  		expand(DEMULTIPLEX_DIR + "qc/{techname}_{capDesign}_{sizeFrac}.demul.QC1.txt",filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS) if config["DEMULTIPLEX"]  else expand( DUMMY_DIR + "dummy{number}.txt", number='9'), # QC on demultiplexing (checks that there is only one barcode assigned per read
  		expand(DEMULTIPLEX_DIR + "qc/{techname}_{capDesign}_{sizeFrac}.demul.QC2.txt", filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS) if config["DEMULTIPLEX"]  else expand( DUMMY_DIR + "dummy{number}.txt", number='10'), # QC on demultiplexing (checks that there is only one barcode assigned per read
  		expand(config["PLOTSDIR"] + "{techname}.demultiplexing.perSample.stats.{ext}", techname=TECHNAMES, ext=config["PLOTFORMATS"])  if config["DEMULTIPLEX"]  else expand( DUMMY_DIR + "dummy{number}.txt", number='11'),
  		expand(config["PLOTSDIR"] + "{techname}.mapping.perSample.perFraction.stats.{ext}", techname=TECHNAMES, ext=config["PLOTFORMATS"]),
  		expand("mappings/" + "nonAnchoredMergeReads/vsTargets/{techname}_{capDesign}_pooled.tmerge.gfftsv.gz", techname=TECHNAMES, capDesign=CAPDESIGNS),
  		expand(config["PLOTSDIR"] + "all_pooled.targetCoverage.stats.{ext}", ext=config["PLOTFORMATS"]),
  		expand(config["PLOTSDIR"] + "all.splice.sites.stats.{ext}", ext=config["PLOTFORMATS"]),
  		expand(config["PLOTSDIR"] + "all.pooled.tmerge.vs.gencode.stats.{ext}",  ext=config["PLOTFORMATS"]),
  		expand("mappings/nonAnchoredMergeReads/pooled/cageSupported/{techname}_{capDesign}_pooled.tmerge.cageSupported.bed", techname=TECHNAMES, capDesign=CAPDESIGNS),
  		expand("mappings/nonAnchoredMergeReads/pooled/polyASupported/{techname}_{capDesign}_pooled.tmerge.polyASupported.bed", techname=TECHNAMES, capDesign=CAPDESIGNS),
  		expand("mappings/nonAnchoredMergeReads/pooled/cage+polyASupported/{techname}_{capDesign}_pooled.tmerge.cage+polyASupported.bed", techname=TECHNAMES, capDesign=CAPDESIGNS),
   		expand(config["PLOTSDIR"] + "all.tmerge.cagePolyASupport.stats.{ext}", ext=config["PLOTFORMATS"]),
  		expand(config["PLOTSDIR"] + "all.polyA.vs.PAS.precisionRecall.stats.{ext}", ext=config["PLOTFORMATS"]),




# "Dummy" rule to skip undesired targets depending on config

rule dummy:
	input: "/dev/null"
	output: DUMMY_DIR + "dummy{number}.txt"
	shell:
		'''
touch {output}
		'''
