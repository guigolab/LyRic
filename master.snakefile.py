import glob
from collections import defaultdict
import os
from itertools import product
import sys


print ("TODO:\n ## include contents of README_correction_error_rate.sh into the snakemake workflow (see plotPolyAreadsStats for inspiration)\n\ninclude basic HiSeq mapping stats\n\nWhat is the intersection between PacBio TMs and ONT TMs???")


GGPLOT_PUB_QUALITY=config["GGPLOT_PUB_QUALITY"]
CAPDESIGNTOGENOME=config["capDesignToGenome"]
CAPDESIGNTOANNOTGTF=config["capDesignToAnnotGtf"]
sizeFrac_Rpalette=config["SIZEFRACTION_RPALETTE"]
long_Rpalette=config["LONG_RPALETTE"]
CAPDESIGNTOCAGEPEAKS=config["capDesignToCAGEpeaks"]
CAPDESIGNTOPAS=config["capDesignToGenomePAS"]
DUMMY_DIR="dummy/"
FQPATH=config["FQPATH"]
FQ_CORR_PATH=FQPATH + "corr/"

print ("FASTQs are in: " + FQPATH)

# no underscores allowed in wildcards, to avoid greedy matching since we use them as separators
wildcard_constraints:
 	capDesign = "[^_/]+",
 	sizeFrac = "[^_/]+",
 	techname = "[^_/]+",
 	corrLevel = "[^_/]+",
# 	barcodes = "[^\.]",
# 	barcodesU ="^(merged)"

if config["DEMULTIPLEX"]:
	# get CAPDESIGNS (capture designs, i.e. Hv1, Hv2, Mv1) and SIZEFRACS (size fractions) variables from FASTQ file names (warning: this will generate duplicate entries so "set" them afterwards):
	DEMULTIPLEX_DIR=config["DEMULTIPLEX_DIR"]
	DEMULTIPLEXED_FASTQS= DEMULTIPLEX_DIR + "demultiplexFastqs/"
	(TECHNAMES, CAPDESIGNS, SIZEFRACS) = glob_wildcards(FQPATH + "{techname, [^_/]+}_{capDesign}_{sizeFrac}.fastq.gz")
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
#	DEMULTIPLEX_DIR=FQPATH
	DEMULTIPLEXED_FASTQS=FQ_CORR_PATH
	(TECHNAMES, CAPDESIGNS, SIZEFRACS, BARCODES) = glob_wildcards(FQPATH + "{techname, [^_/]+}_{capDesign}_{sizeFrac}.{barcodes}.fastq.gz")
	BARCODES=set(BARCODES)
	BARCODESUNDETER=list(BARCODES)


CAPDESIGNS=set(CAPDESIGNS)
SIZEFRACS=set(SIZEFRACS)
TECHNAMES=set(TECHNAMES)


BARCODESpluSMERGED=set(BARCODES)
BARCODESpluSMERGED.add("allTissues")
SIZEFRACSpluSMERGED=set(SIZEFRACS)
SIZEFRACSpluSMERGED.add("allFracs")

print ("NONMERGED: ", TECHNAMES, CAPDESIGNS, SIZEFRACS,BARCODES)

print ("MERGED: ", TECHNAMES, CAPDESIGNS, SIZEFRACSpluSMERGED,BARCODESpluSMERGED)


SPLICE_SITE_TYPES=["Donor", "Acceptor"]

# for polyA calling calibration:
minPolyAlength = []
for x in range(5, 30, 1):
	minPolyAlength.append(x)



if config["LORDEC_CORRECT"]:
	graph_kmers=["17", "18", "19", "20", "21", "30", "40", "50", "60", "70", "80", "90"]
	lastK=graph_kmers[-1]
	FINALCORRECTIONLEVELS=["0", lastK]
	#make string to interpolate into bash script in "lordecCorrectLr" rule
	graph_kmers_string= "(" + " ".join(graph_kmers) + ")"
	solid_kmer_abundance_threshold="3"
	splitFastqsInto=99
	splitFasta =[]
	for x in range(0, splitFastqsInto, 1):
		splitFasta.append(str(x).zfill(4))
else:
	FINALCORRECTIONLEVELS=["0"]
	lastK='NULL'

########################
### make list of authorized wildcard combinations
### inspired by https://stackoverflow.com/questions/41185567/how-to-use-expand-in-snakemake-when-some-particular-combinations-of-wildcards-ar
###

AUTHORIZEDCOMBINATIONS = []
AUTHORIZEDCOMBINATIONSMERGE = []

for comb in product(TECHNAMES,CAPDESIGNS,SIZEFRACS,BARCODES):
	if(os.path.isfile(FQPATH + comb[0] + "_" + comb[1] + "_" + comb[2] + ".fastq.gz") or os.path.isfile(FQPATH + comb[0] + "_" + comb[1] + "_" + comb[2]  + "." + comb[3] + ".fastq.gz")): #allow only combinations corresponding to existing FASTQs
		AUTHORIZEDCOMBINATIONS.append((("techname", comb[0]),("capDesign", comb[1]),("sizeFrac", comb[2])))
		AUTHORIZEDCOMBINATIONSMERGE.append((("techname", comb[0]),("capDesign", comb[1]),("sizeFrac", "allFracs")))
		for corrLevel in FINALCORRECTIONLEVELS:
			AUTHORIZEDCOMBINATIONS.append((("techname", comb[0]),("corrLevel", corrLevel),("capDesign", comb[1]),("sizeFrac", comb[2])))
			AUTHORIZEDCOMBINATIONSMERGE.append((("techname", comb[0]),("corrLevel", corrLevel),("capDesign", comb[1]),("sizeFrac", "allFracs")))
		for ext in config["PLOTFORMATS"]:
			for corrLevel in FINALCORRECTIONLEVELS:
				AUTHORIZEDCOMBINATIONS.append((("techname", comb[0]),("corrLevel", corrLevel),("capDesign", comb[1]),("sizeFrac", comb[2]),("ext",ext)))
				AUTHORIZEDCOMBINATIONSMERGE.append((("techname", comb[0]),("corrLevel", corrLevel),("capDesign", comb[1]),("sizeFrac", "allFracs"),("ext",ext)))
#		print(comb[3], comb[1], sep=",")
		if (comb[3]).find(comb[1]) == 0 or config["DEMULTIPLEX"] is False: # check that barcode ID's prefix matches capDesign (i.e. ignore demultiplexed FASTQs with unmatched barcode/capDesign)
#			print ("MATCH")
			for corrLevel in FINALCORRECTIONLEVELS:
				AUTHORIZEDCOMBINATIONS.append((("techname", comb[0]),("corrLevel", corrLevel),("capDesign", comb[1]),("sizeFrac", comb[2]),("barcodes",comb[3])))
				AUTHORIZEDCOMBINATIONSMERGE.append((("techname", comb[0]),("corrLevel", corrLevel),("capDesign", comb[1]),("sizeFrac", "allFracs"),("barcodes",comb[3])))
				AUTHORIZEDCOMBINATIONSMERGE.append((("techname", comb[0]),("corrLevel", corrLevel),("capDesign", comb[1]),("sizeFrac", "allFracs"),("barcodes","allTissues")))
				AUTHORIZEDCOMBINATIONSMERGE.append((("techname", comb[0]),("corrLevel", corrLevel),("capDesign", comb[1]),("sizeFrac", comb[2]),("barcodes","allTissues")))
				AUTHORIZEDCOMBINATIONS.append((("techname", comb[0]),("corrLevel", corrLevel),("capDesign", comb[1]),("barcodes",comb[3])))
				AUTHORIZEDCOMBINATIONSMERGE.append((("techname", comb[0]),("corrLevel", corrLevel),("capDesign", comb[1]),("barcodes","allTissues")))

			for strand in config["STRANDS"]:
				for corrLevel in FINALCORRECTIONLEVELS:
					AUTHORIZEDCOMBINATIONS.append((("techname", comb[0]),("corrLevel", corrLevel),("capDesign", comb[1]),("barcodes",comb[3]),("strand", strand)))
					AUTHORIZEDCOMBINATIONSMERGE.append((("techname", comb[0]),("corrLevel", corrLevel),("capDesign", comb[1]),("barcodes","allTissues"),("strand", strand)))
			for ext in config["PLOTFORMATS"]:
				for corrLevel in FINALCORRECTIONLEVELS:
					AUTHORIZEDCOMBINATIONS.append((("techname", comb[0]),("corrLevel", corrLevel),("capDesign", comb[1]),("barcodes",comb[3]),("ext", ext)))
					AUTHORIZEDCOMBINATIONSMERGE.append((("techname", comb[0]),("corrLevel", corrLevel),("capDesign", comb[1]),("barcodes","allTissues"),("ext", ext)))



AUTHORIZEDCOMBINATIONS=set(AUTHORIZEDCOMBINATIONS)
AUTHORIZEDCOMBINATIONSMERGE=set(AUTHORIZEDCOMBINATIONSMERGE)

#print ("AUTHORIZEDCOMBINATIONS:", AUTHORIZEDCOMBINATIONS)
#print ("AUTHORIZEDCOMBINATIONSMERGE:", AUTHORIZEDCOMBINATIONSMERGE)

def filtered_product(*args):
	for wc_comb in product(*args):
		found=False
		print (wc_comb)
		if wc_comb in AUTHORIZEDCOMBINATIONS:
			#print ("AUTH")
			found=True
			yield(wc_comb)
		#if not found:
			#print ("NONAUTH")

def filtered_product_merge(*args):
	for wc_comb in product(*args):
		found=False
		#print (wc_comb)
		if wc_comb in AUTHORIZEDCOMBINATIONS or wc_comb in AUTHORIZEDCOMBINATIONSMERGE:
			#print ("AUTH")
			found=True
			yield(wc_comb)



include: "lrCorrection.snakefile.py"
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
#		expand(FQPATH + "corr/{techname}Corr0_{capDesign}_{sizeFrac}.fastq.gz", filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS),
		expand(FQ_CORR_PATH + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.fastq.gz", filtered_product, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS) if config["DEMULTIPLEX"] else expand(FQ_CORR_PATH + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.{barcodes}.fastq.gz", filtered_product, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES),

		expand(config["PLOTSDIR"] + "{techname}Corr{corrLevel}_{capDesign}_all.readlength.{ext}", techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, ext=config["PLOTFORMATS"]) if config["DEMULTIPLEX"] else expand(config["PLOTSDIR"] + "{techname}Corr{corrLevel}_{capDesign}.{barcodes}_all.readlength.{ext}",filtered_product, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, barcodes=BARCODES, ext=config["PLOTFORMATS"]), # facetted histograms of read length
 		expand(FQ_CORR_PATH + "qc/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.dupl.txt", filtered_product,techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS) if config["DEMULTIPLEX"] else expand(FQ_CORR_PATH + "qc/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.{barcodes}.dupl.txt", filtered_product,techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES),
  		expand(config["PLOTSDIR"] + "{techname}Corr{corrLevel}.fastq.UP.stats.{ext}", techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, ext=config["PLOTFORMATS"]) if config["DEMULTIPLEX"] else expand( DUMMY_DIR + "dummy{number}.txt", number='1'), # UP reads plots
  		expand(config["PLOTSDIR"] + "{techname}Corr{corrLevel}.fastq.BC.stats.{ext}", techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, ext=config["PLOTFORMATS"]) if config["DEMULTIPLEX"]  else expand( DUMMY_DIR + "dummy{number}.txt", number='2') , # barcode reads plots
  		expand(config["PLOTSDIR"] + "{techname}Corr{corrLevel}.fastq.foreignBC.stats.{ext}", techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, ext=config["PLOTFORMATS"]) if config["DEMULTIPLEX"]  else expand( DUMMY_DIR + "dummy{number}.txt", number='3') , #foreign barcode reads plots
  		expand ("mappings/" + "hiSeq_{capDesign}.bam", capDesign=CAPDESIGNS) if config["USE_MATCHED_ILLUMINA"] else expand( DUMMY_DIR + "dummy{number}.txt", number='5') ,  # mapped short reads
  		expand(DEMULTIPLEX_DIR + "demultiplexFastqs/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.{barcodes}.fastq.gz", filtered_product, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES) if config["DEMULTIPLEX"]  else expand( DUMMY_DIR + "dummy{number}.txt", number='4'),


 		expand ("mappings/" + "clusterPolyAsites/{techname}Corr{corrLevel}_{capDesign}_{barcodes}.polyAsites.clusters.bed", filtered_product, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, barcodes=BARCODES),
 		expand(config["PLOTSDIR"] + "all.pooled.merged.HiSS.stats.{ext}", ext=config["PLOTFORMATS"]),
 		expand("mappings/" + "mergeCapDesignBams/{techname}Corr{corrLevel}_{capDesign}.merged2.bam", techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS),
 		expand("mappings/" + "nonAnchoredMergeReads/pooled/{techname}Corr{corrLevel}_{capDesign}_pooled.tmerge.gff", techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS),
 		expand("mappings/" + "nonAnchoredMergeReads/qc/{techname}Corr{corrLevel}_{capDesign}_{barcodes}.tmerge.qc.txt",filtered_product, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, barcodes=BARCODES),
 		expand( "mappings/" + "nonAnchoredMergeReads/pooled/qc/{techname}Corr{corrLevel}_{capDesign}.tmerge.qc.txt", techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS),
 		expand(config["PLOTSDIR"] + "all.pooled.merged.stats.{ext}", ext=config["PLOTFORMATS"]),
 		expand(config["PLOTSDIR"] + "all.polyAreads.stats.{ext}", ext=config["PLOTFORMATS"]),
 		expand("mappings/" + "makePolyABigWigs/{techname}Corr{corrLevel}_{capDesign}_{barcodes}.polyAsitesNoErcc.{strand}.bw", filtered_product, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, barcodes=BARCODES, strand=config["STRANDS"]),

  		expand("mappings/" + "qc/{techname}Corr{corrLevel}_{capDesign}_{barcodes}.merged.bam.dupl.txt",filtered_product, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, barcodes=BARCODES),
  		expand("mappings/qualimap_reports/" + "{techname}Corr{corrLevel}_{capDesign}.merged2/genome_results.txt", techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS),
  		expand(config["PLOTSDIR"] + "{techname}Corr{corrLevel}.ambiguousBarcodes.reads.stats.{ext}", techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, ext=config["PLOTFORMATS"]) if config["DEMULTIPLEX"]  else expand( DUMMY_DIR + "dummy{number}.txt", number='6'), # ambiguous barcodes plots
  		expand(config["PLOTSDIR"] + "{techname}Corr{corrLevel}_{capDesign}.adapters.location.stats.{ext}",techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, ext=config["PLOTFORMATS"]) if config["DEMULTIPLEX"]  else expand( DUMMY_DIR + "dummy{number}.txt", number='7'), #location of adapters over reads
  		expand(config["PLOTSDIR"] + "{techname}Corr{corrLevel}.chimeric.reads.stats.{ext}", techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, ext=config["PLOTFORMATS"]) if config["DEMULTIPLEX"]  else expand( DUMMY_DIR + "dummy{number}.txt", number='8'), # stats on chimeric reads
  		expand(config["PLOTSDIR"] + "{techname}Corr{corrLevel}.finalDemul.reads.stats.{ext}", techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, ext=config["PLOTFORMATS"]) if config["DEMULTIPLEX"]  else expand( DUMMY_DIR + "dummy{number}.txt", number='9'), #final demultiplexing stats
  		expand(DEMULTIPLEX_DIR + "qc/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.demul.QC1.txt",filtered_product, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS) if config["DEMULTIPLEX"]  else expand( DUMMY_DIR + "dummy{number}.txt", number='10'), # QC on demultiplexing (checks that there is only one barcode assigned per read
  		expand(DEMULTIPLEX_DIR + "qc/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.demul.QC2.txt", filtered_product, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS) if config["DEMULTIPLEX"]  else expand( DUMMY_DIR + "dummy{number}.txt", number='11'), # QC on demultiplexing (checks that there is only one barcode assigned per read
  		expand(config["PLOTSDIR"] + "{techname}Corr{corrLevel}.demultiplexing.perSample.stats.{ext}", techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, ext=config["PLOTFORMATS"])  if config["DEMULTIPLEX"]  else expand( DUMMY_DIR + "dummy{number}.txt", number='12'),
  		expand(config["PLOTSDIR"] + "{techname}Corr{corrLevel}.mapping.perSample.perFraction.stats.{ext}", techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, ext=config["PLOTFORMATS"]),
  		expand("mappings/" + "nonAnchoredMergeReads/vsTargets/{techname}Corr{corrLevel}_{capDesign}_pooled.tmerge.gfftsv.gz", techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS),
  		expand(config["PLOTSDIR"] + "all_pooled.targetCoverage.stats.{ext}", ext=config["PLOTFORMATS"]),
  		expand(config["PLOTSDIR"] + "all.splice.sites.stats.{ext}", ext=config["PLOTFORMATS"]),
  		expand(config["PLOTSDIR"] + "all.pooled.tmerge.vs.gencode.stats.{ext}",  ext=config["PLOTFORMATS"]),
  		expand("mappings/nonAnchoredMergeReads/pooled/cageSupported/{techname}Corr{corrLevel}_{capDesign}_pooled.tmerge.cageSupported.bed", techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS),
  		expand("mappings/nonAnchoredMergeReads/pooled/polyASupported/{techname}Corr{corrLevel}_{capDesign}_pooled.tmerge.polyASupported.bed", techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS),
  		expand("mappings/nonAnchoredMergeReads/pooled/cage+polyASupported/{techname}Corr{corrLevel}_{capDesign}_pooled.tmerge.cage+polyASupported.bed", techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS),
   		expand(config["PLOTSDIR"] + "all.tmerge.cagePolyASupport.stats.{ext}", ext=config["PLOTFORMATS"]),
  		expand(config["PLOTSDIR"] + "all.polyA.vs.PAS.precisionRecall.stats.{ext}", ext=config["PLOTFORMATS"]),

  		#expand("mappings/" + "readBamToBed/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bed", filtered_product_merge, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACSpluSMERGED, barcodes=BARCODESpluSMERGED)
		expand("mappings/readMapping/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bam", filtered_product, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES),
		expand("mappings/readMapping/{techname}Corr{corrLevel}_{capDesign}_allFracs_{barcodes}.bam", filtered_product, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, barcodes=BARCODES),
		#expand("mappings/readMapping/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_allTissues.bam", filtered_product, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS)




# "Dummy" rule to skip undesired targets depending on config

rule dummy:
	input: "/dev/null"
	output: DUMMY_DIR + "dummy{number}.txt"
	shell:
		'''
touch {output}
		'''
