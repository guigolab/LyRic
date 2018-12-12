import glob
from collections import defaultdict
import os
import itertools
import sys
#import pprint #pretty printing

print ("TODO:\n ## include basic HiSeq mapping stats\n\n## Calculate amount of novel captured nts\n\n ## What is the intersection between PacBio TMs and ONT TMs???\n\n")


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
TRACK_HUB_DATA_URL=config["TRACK_HUB_BASE_URL"] + "dataFiles/"
print ("Input LR FASTQs are in: " + FQPATH)

GENOMES=[]
GENOMETOCAPDESIGNS=defaultdict(list)
for capD in CAPDESIGNTOGENOME:
	genome=CAPDESIGNTOGENOME[capD]
	GENOMES.append(genome)
	GENOMETOCAPDESIGNS[genome].append(capD)
GENOMES=set(GENOMES)


# no underscores allowed in wildcards, to avoid greedy matching since we use them as separators
wildcard_constraints:
 	capDesign = "[^_/]+",
 	sizeFrac = "[^_/]+",
 	techname = "[^_/]+",
 	corrLevel = "[^_/]+",

SIRVpresent=None
try:
	config["SIRVgff"] and config["SIRVinfo"]
except KeyError:
	SIRVpresent=False
	print ("Will NOT compare TMs to SIRVs because config['SIRVgff'] and/or config['SIRVinfo'] not found.")
else:
	SIRVpresent=True
	print ("Will compare TMs to SIRVs")


if config["DEMULTIPLEX"]:
	print ("Demultiplexing is ON")
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
	print ("Demultiplexing is OFF")
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
CAPDESIGNSplusMERGED=set(CAPDESIGNS)
CAPDESIGNSplusMERGED.add("allCapDesigns")
TECHNAMESplusMERGED=set(TECHNAMES)
TECHNAMESplusMERGED.add("allTechs")


SPLICE_SITE_TYPES=["Donor", "Acceptor"]

# for polyA calling calibration:
minPolyAlength = []
for x in range(5, 30, 1):
	minPolyAlength.append(x)



if config["LORDEC_CORRECT"]:
	graph_kmers=["17", "18", "19", "20", "21", "30", "40", "50", "60", "70", "80", "90"]
	lastK=graph_kmers[-1]
	FINALCORRECTIONLEVELS=["0", lastK]
	print ("LR error correction is ON. Correction levels: ")
	print (FINALCORRECTIONLEVELS)

	#make string to interpolate into bash script in "lordecCorrectLr" rule
	graph_kmers_string= "(" + " ".join(graph_kmers) + ")"
	solid_kmer_abundance_threshold="3"
	splitFastqsInto=19
	splitFasta =[]
	for x in range(0, splitFastqsInto, 1):
		splitFasta.append(str(x).zfill(4))
else:
	print ("LR error correction is OFF")
	FINALCORRECTIONLEVELS=["0"]
	lastK='NULL'



FINALCORRECTIONLEVELSplusMERGED=set(FINALCORRECTIONLEVELS)
FINALCORRECTIONLEVELSplusMERGED.add("byCorr")

########################
### make list of authorized wildcard combinations
### inspired by https://stackoverflow.com/questions/41185567/how-to-use-expand-in-snakemake-when-some-particular-combinations-of-wildcards-ar
###

AUTHORIZEDCOMBINATIONS = []
AUTHORIZEDCOMBINATIONSMERGE = []

for comb in itertools.product(TECHNAMES,CAPDESIGNS,SIZEFRACS,BARCODES):
	if(os.path.isfile(FQPATH + comb[0] + "_" + comb[1] + "_" + comb[2] + ".fastq.gz") or
	os.path.isfile(FQPATH + comb[0] + "_" + comb[1] + "_" + comb[2]  + "." + comb[3] + ".fastq.gz")): #allow only combinations corresponding to existing FASTQs
		AUTHORIZEDCOMBINATIONS.append((("techname", comb[0]),("capDesign", comb[1]),("sizeFrac", comb[2])))
		# for mergedTerm in (MERGEDTERMS):
		for mergedComb in itertools.product([("techname", comb[0]), ("techname", "allTechs")], [("capDesign", comb[1]), ("capDesign", "allCapDesigns")], [("sizeFrac", comb[2]), ("sizeFrac", "allFracs")]):
			AUTHORIZEDCOMBINATIONSMERGE.append((mergedComb))

		for corrLevel in FINALCORRECTIONLEVELS:
			AUTHORIZEDCOMBINATIONS.append((("techname", comb[0]),("corrLevel", corrLevel),("capDesign", comb[1]),("sizeFrac", comb[2])))
			for mergedComb in itertools.product([("techname", comb[0]), ("techname", "allTechs")], [("corrLevel", corrLevel), ("corrLevel", "byCorr")] , [("capDesign", comb[1]), ("capDesign", "allCapDesigns")], [("sizeFrac", comb[2]), ("sizeFrac", "allFracs")]):
				AUTHORIZEDCOMBINATIONSMERGE.append((mergedComb))
		for ext in config["PLOTFORMATS"]:
			for corrLevel in FINALCORRECTIONLEVELS:
				AUTHORIZEDCOMBINATIONS.append((("techname", comb[0]),("corrLevel", corrLevel),("capDesign", comb[1]),("sizeFrac", comb[2]),("ext",ext)))
				for mergedComb in itertools.product([("techname", comb[0]), ("techname", "allTechs")], [("corrLevel", corrLevel), ("corrLevel", "byCorr")] , [("capDesign", comb[1]), ("capDesign", "allCapDesigns")], [("sizeFrac", comb[2]), ("sizeFrac", "allFracs")], [("ext",ext)]):
					AUTHORIZEDCOMBINATIONSMERGE.append((mergedComb))

		if (comb[3]).find(comb[1]) == 0 or config["DEMULTIPLEX"] is False: # check that barcode ID's prefix matches capDesign (i.e. ignore demultiplexed FASTQs with unmatched barcode/capDesign)
			AUTHORIZEDCOMBINATIONS.append((("techname", comb[0]),("capDesign", comb[1]),("barcodes", comb[3])))

			for split in splitFasta:
				AUTHORIZEDCOMBINATIONS.append((("techname", comb[0]),("capDesign", comb[1]),("sizeFrac", comb[2]),("barcodes",comb[3]), ("split", split)))
			for corrLevel in FINALCORRECTIONLEVELS:
				AUTHORIZEDCOMBINATIONS.append((("techname", comb[0]),("corrLevel", corrLevel),("capDesign", comb[1]),("sizeFrac", comb[2]),("barcodes",comb[3])))
				for mergedComb in itertools.product([("techname", comb[0]), ("techname", "allTechs")], [("corrLevel", corrLevel), ("corrLevel", "byCorr")] , [("capDesign", comb[1]), ("capDesign", "allCapDesigns")], [("sizeFrac", comb[2]), ("sizeFrac", "allFracs")], [("barcodes",comb[3]), ("barcodes", "allTissues")]):
					AUTHORIZEDCOMBINATIONSMERGE.append((mergedComb))
				AUTHORIZEDCOMBINATIONS.append((("techname", comb[0]),("corrLevel", corrLevel),("capDesign", comb[1]),("barcodes",comb[3])))
				for mergedComb in itertools.product([("techname", comb[0]), ("techname", "allTechs")], [("corrLevel", corrLevel), ("corrLevel", "byCorr")] , [("capDesign", comb[1]), ("capDesign", "allCapDesigns")], [("barcodes",comb[3]), ("barcodes", "allTissues")]):
					AUTHORIZEDCOMBINATIONSMERGE.append((mergedComb))

			for strand in config["STRANDS"]:
				for corrLevel in FINALCORRECTIONLEVELS:
					AUTHORIZEDCOMBINATIONS.append((("techname", comb[0]),("corrLevel", corrLevel),("capDesign", comb[1]),("barcodes",comb[3]),("strand", strand)))
					AUTHORIZEDCOMBINATIONS.append((("techname", comb[0]),("corrLevel", corrLevel),("capDesign", comb[1]),("sizeFrac", comb[2]),("barcodes",comb[3]),("strand", strand)))
					for mergedComb in itertools.product([("techname", comb[0]), ("techname", "allTechs")], [("corrLevel", corrLevel), ("corrLevel", "byCorr")] , [("capDesign", comb[1]), ("capDesign", "allCapDesigns")], [("barcodes",comb[3]), ("barcodes", "allTissues")], [("strand", strand)]):
						AUTHORIZEDCOMBINATIONSMERGE.append((mergedComb))
					for mergedComb in itertools.product([("techname", comb[0]), ("techname", "allTechs")], [("corrLevel", corrLevel), ("corrLevel", "byCorr")] , [("capDesign", comb[1]), ("capDesign", "allCapDesigns")], [("sizeFrac", comb[2]), ("sizeFrac", "allFracs")], [("barcodes",comb[3]), ("barcodes", "allTissues")], [("strand", strand)]):
						AUTHORIZEDCOMBINATIONSMERGE.append((mergedComb))
			for ext in config["PLOTFORMATS"]:
				for corrLevel in FINALCORRECTIONLEVELS:
					AUTHORIZEDCOMBINATIONS.append((("techname", comb[0]),("corrLevel", corrLevel),("capDesign", comb[1]),("barcodes",comb[3]),("ext", ext)))
					for mergedComb in itertools.product([("techname", comb[0]), ("techname", "allTechs")], [("corrLevel", corrLevel), ("corrLevel", "byCorr")] , [("capDesign", comb[1]), ("capDesign", "allCapDesigns")], [("barcodes",comb[3]), ("barcodes", "allTissues")], [("ext", ext)]):
						AUTHORIZEDCOMBINATIONSMERGE.append((mergedComb))

				#for merged figures:
				for mergedComb in itertools.product([("capDesign", comb[1]), ("capDesign", "allCapDesigns")], [("sizeFrac", comb[2]), ("sizeFrac", "allFracs")], [("barcodes",comb[3]), ("barcodes", "allTissues")], [("ext", ext)]):
					AUTHORIZEDCOMBINATIONSMERGE.append((mergedComb))

			for splice in SPLICE_SITE_TYPES:
				for corrLevel in FINALCORRECTIONLEVELS:
					AUTHORIZEDCOMBINATIONS.append((("techname", comb[0]),("corrLevel", corrLevel),("capDesign", comb[1]),("sizeFrac", comb[2]),("barcodes",comb[3]),("spliceType",splice)))
					for mergedComb in itertools.product([("techname", comb[0]), ("techname", "allTechs")], [("corrLevel", corrLevel), ("corrLevel", "byCorr")] , [("capDesign", comb[1]), ("capDesign", "allCapDesigns")], [("sizeFrac", comb[2]), ("sizeFrac", "allFracs")], [("barcodes",comb[3]), ("barcodes", "allTissues")], [("spliceType",splice)]):
						AUTHORIZEDCOMBINATIONSMERGE.append((mergedComb))




AUTHORIZEDCOMBINATIONS=set(AUTHORIZEDCOMBINATIONS)
AUTHORIZEDCOMBINATIONSMERGE=set(AUTHORIZEDCOMBINATIONSMERGE)

PLOTSbyTISSUE=set(BARCODESpluSMERGED)
PLOTSbyTISSUE.add("byTissue")
PLOTSbySIZEFRAC=set(SIZEFRACSpluSMERGED)
PLOTSbySIZEFRAC.add("byFrac")
#PLOTSbyTISSUE=['allTissues', 'byTissue']
#PLOTSbySIZEFRAC=['allFracs','byFrac']



def filtered_product(*args):
	for wc_comb in itertools.product(*args):
		found=False
		if wc_comb in AUTHORIZEDCOMBINATIONS:
			found=True
			yield(wc_comb)

def filtered_product_merge(*args):
	for wc_comb in itertools.product(*args):
		found=False

		if wc_comb in AUTHORIZEDCOMBINATIONSMERGE:
			#print("FPM ")
			#print (wc_comb)
		#	print ("AUTH")
			found=True
			yield(wc_comb)

def filtered_merge_figures(*args):
	for wc_comb in itertools.product(*args):
		#print (wc_comb)
		if wc_comb[0] in (('capDesign', 'allCapDesigns'),):
			if wc_comb[3] in (('barcodes', 'allTissues'),):
#			if wc_comb[1] in (('sizeFrac', 'allFracs'),) and wc_comb[2] in (('barcodes', 'allTissues'),):
				#print ("AUTH 2")
				yield(wc_comb)
		else:
			if wc_comb in AUTHORIZEDCOMBINATIONSMERGE:
				#print ("AUTH 1")
				yield(wc_comb)
			else:
				if wc_comb[2] not in (('sizeFrac', 'allFracs'),) and wc_comb[3] in (('barcodes', 'allTissues'),):
					#print ("AUTH 3")
					yield(wc_comb)
				if wc_comb[2] in (('sizeFrac', 'allFracs'),) and (wc_comb[3] not in (('barcodes', 'allTissues'),) and wc_comb[3][1].find(wc_comb[0][1]) == 0) or wc_comb[3] in (('barcodes', 'byTissue'),):
						#print ("AUTH 4")
						yield(wc_comb)

def merge_figures_params(c,bf,bt,cl):
	clBool=None
	if cl == FINALCORRECTIONLEVELS[1]:
		clBool="Yes"
	elif cl == FINALCORRECTIONLEVELS[0]:
		clBool="No"
	else:
		clBool=cl

	substSeqTechString="""
dat\$seqTech <- gsub(':', '\\n', dat\$seqTech)
"""
	substTissueString="dat\$tissue <- gsub('" + c + "_', '', dat\$tissue)\n"
	capDesignFilterString=''
	corrLevelFilterString=''
	sizeFracFilterString=''
	tissueFilterString=''
	xlabString='Error correction'
	hideXaxisLabels=''
	hideYaxisLabels=''
	graphDimensions="""
horizCats <- length(unique(dat\$correctionLevel)) * length(unique(dat\$capDesign)) * length(unique(dat\$tissue))
vertCats <- length(unique(dat\$seqTech)) * length(unique(dat\$sizeFrac))
plotWidth = horizCats + 3.5
plotHeight = vertCats + 1
geom_textSize=2.4
"""
	if c not in ('allCapDesigns'):
		capDesignFilterString+="dat <- subset(dat, capDesign=='" + c + "')\n"
	if clBool not in ('byCorr'):
		corrLevelFilterString+="dat <- subset(dat, correctionLevel=='" + clBool + "')\n"
		xlabString=''
		hideXaxisLabels="theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +"
		hideYaxisLabels="theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) +"
	if bf in ('allFracs'):
		sizeFracFilterString+="dat <- subset(dat, sizeFrac=='" + bf + "')\n"
	else:
		if bf in ('byFrac'):
			sizeFracFilterString+="dat <- subset(dat, sizeFrac!='allFracs')\n"
		else:
			sizeFracFilterString+="dat <- subset(dat, sizeFrac=='" + bf + "')\n"
	if bt in ('allTissues'):
		tissueFilterString+="dat <- subset(dat, tissue=='" + bt + "')\n"
	else:
		if bt in ('byTissue'):
			tissueFilterString+="dat <- subset(dat, tissue!='allTissues')\n"
		else:
			tissueFilterString+="dat <- subset(dat, tissue=='" + bt + "')\n"
	return(capDesignFilterString, corrLevelFilterString, sizeFracFilterString, tissueFilterString, substSeqTechString, substTissueString,  xlabString, hideXaxisLabels , graphDimensions, hideYaxisLabels)


def trackHubSubGroupString(tn, cd, sf, bc, cl):
	techname=(("techname", tn),)
	capDesign=(("capDesign", cd),)
	sizeFracs=[]
	for sizeF in sf:
		sizeFracs.append(("sizeFrac", sizeF))
	barCodes=[]
	for barC in bc:
		barCodes.append(("barcodes", barC))
	corrLevels=[]
	for corrL in cl:
		corrLevels.append(("corrLevel", corrL))

	returnSubGroup1String="subGroup1 sample Sample"
	returnSubGroup2String="subGroup2 sizeFraction Size_fraction"
	returnSubGroup3String="subGroup3 corrLevel LoRDEC_sequence_correction"
	returnSubGroup1StringData=[]
	returnSubGroup2StringData=[]
	returnSubGroup3StringData=[]
	for wc_comb in itertools.product(techname, corrLevels, capDesign, sizeFracs, barCodes):
		#print("THSS ")
		#print(wc_comb)
		if wc_comb in AUTHORIZEDCOMBINATIONSMERGE:
		#	print ("AUTH")
			returnSubGroup1StringData.append(wc_comb[4][1] + "=" + wc_comb[4][1])
			returnSubGroup2StringData.append(wc_comb[3][1] + "=" + wc_comb[3][1])
			returnSubGroup3StringData.append(wc_comb[1][1] + "=" + wc_comb[1][1])
	#print ("returnSubGroup2StringData list")
	#print (returnSubGroup2StringData)
	returnSubGroup1StringData=set(returnSubGroup1StringData)
	returnSubGroup2StringData=set(returnSubGroup2StringData)
	returnSubGroup3StringData=set(returnSubGroup3StringData)
	#print ("returnSubGroup2StringData set")
	#print (returnSubGroup2StringData)
	return(returnSubGroup1String + " " + " ".join(returnSubGroup1StringData) + "\n" + returnSubGroup2String + " " + " ".join(returnSubGroup2StringData) + "\n" + returnSubGroup3String + " " + " ".join(returnSubGroup3StringData))

def trackChecked(t,c,cd,s,b):
	checked=''
	if t.find('pacBio') == 0 and c == FINALCORRECTIONLEVELS[1]:
		checked='off'
	elif t.find('ont') == 0 and c == FINALCORRECTIONLEVELS[1]:
		checked='off'
	else:
		if s == 'allFracs' and b == 'allTissues':
			checked='on'
		else:
			checked='off'
	return(checked)

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
include: "trackHub.py"

#pseudo-rule specifying the target files we ultimately want.
rule all:
	input:
		expand(FQ_CORR_PATH + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.fastq.gz", filtered_product, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS) if config["DEMULTIPLEX"] else expand(FQ_CORR_PATH + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.{barcodes}.fastq.gz", filtered_product, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES),

		expand(config["PLOTSDIR"] + "{techname}Corr{corrLevel}_{capDesign}_all.readlength.{ext}", techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, ext=config["PLOTFORMATS"]) if config["DEMULTIPLEX"] else expand(config["PLOTSDIR"] + "{techname}Corr{corrLevel}_{capDesign}.{barcodes}_all.readlength.{ext}",filtered_product, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, barcodes=BARCODES, ext=config["PLOTFORMATS"]), # facetted histograms of read length
 		expand(FQ_CORR_PATH + "qc/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.dupl.txt", filtered_product,techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS) if config["DEMULTIPLEX"] else expand(FQ_CORR_PATH + "qc/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.{barcodes}.dupl.txt", filtered_product,techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES),
  		expand(config["PLOTSDIR"] + "{techname}Corr{corrLevel}.fastq.UP.stats.{ext}", techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, ext=config["PLOTFORMATS"]) if config["DEMULTIPLEX"] else expand( DUMMY_DIR + "dummy{number}.txt", number='1'), # UP reads plots
  		expand(config["PLOTSDIR"] + "{techname}Corr{corrLevel}.fastq.BC.stats.{ext}", techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, ext=config["PLOTFORMATS"]) if config["DEMULTIPLEX"]  else expand( DUMMY_DIR + "dummy{number}.txt", number='2') , # barcode reads plots
  		expand(config["PLOTSDIR"] + "{techname}Corr{corrLevel}.fastq.foreignBC.stats.{ext}", techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, ext=config["PLOTFORMATS"]) if config["DEMULTIPLEX"]  else expand( DUMMY_DIR + "dummy{number}.txt", number='3') , #foreign barcode reads plots

   		expand(DEMULTIPLEX_DIR + "demultiplexFastqs/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.{barcodes}.fastq.gz", filtered_product, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES) if config["DEMULTIPLEX"]  else expand( DUMMY_DIR + "dummy{number}.txt", number='4'),
  		expand ("mappings/clusterPolyAsites/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.polyAsites.clusters.bed", filtered_product_merge, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACSpluSMERGED, barcodes=BARCODESpluSMERGED),
 		expand(config["PLOTSDIR"] + "HiSS.stats/{capDesign}_Corr{corrLevel}_{sizeFrac}_{barcodes}.HiSS.stats.{ext}", filtered_merge_figures, capDesign=CAPDESIGNSplusMERGED, corrLevel=FINALCORRECTIONLEVELSplusMERGED, sizeFrac=PLOTSbySIZEFRAC, barcodes=PLOTSbyTISSUE, ext=config["PLOTFORMATS"]),
  		expand("mappings/nonAnchoredMergeReads/colored/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.bed", filtered_product_merge, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACSpluSMERGED, barcodes=BARCODESpluSMERGED),

 		expand(config["PLOTSDIR"] + "merged.stats/{capDesign}_Corr{corrLevel}_{sizeFrac}_{barcodes}.merged.stats.{ext}", filtered_merge_figures, capDesign=CAPDESIGNSplusMERGED, corrLevel=FINALCORRECTIONLEVELSplusMERGED, sizeFrac=PLOTSbySIZEFRAC, barcodes=PLOTSbyTISSUE, ext=config["PLOTFORMATS"]),

  		expand("mappings/makePolyABigWigs/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.polyAsitesNoErcc.{strand}.bw", filtered_product_merge, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACSpluSMERGED, barcodes=BARCODESpluSMERGED, strand=config["STRANDS"]),

   		expand("mappings/readMapping/qc/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bam.dupl.txt",filtered_product_merge, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACSpluSMERGED, barcodes=BARCODESpluSMERGED),
   		expand(config["PLOTSDIR"] + "{techname}Corr{corrLevel}.ambiguousBarcodes.reads.stats.{ext}", techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, ext=config["PLOTFORMATS"]) if config["DEMULTIPLEX"]  else expand( DUMMY_DIR + "dummy{number}.txt", number='6'), # ambiguous barcodes plots
   		expand(config["PLOTSDIR"] + "{techname}Corr{corrLevel}_{capDesign}.adapters.location.stats.{ext}",techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, ext=config["PLOTFORMATS"]) if config["DEMULTIPLEX"]  else expand( DUMMY_DIR + "dummy{number}.txt", number='7'), #location of adapters over reads
   		expand(config["PLOTSDIR"] + "{techname}Corr{corrLevel}.chimeric.reads.stats.{ext}", techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, ext=config["PLOTFORMATS"]) if config["DEMULTIPLEX"]  else expand( DUMMY_DIR + "dummy{number}.txt", number='8'), # stats on chimeric reads
   		expand(config["PLOTSDIR"] + "{techname}Corr{corrLevel}.finalDemul.reads.stats.{ext}", techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, ext=config["PLOTFORMATS"]) if config["DEMULTIPLEX"]  else expand( DUMMY_DIR + "dummy{number}.txt", number='9'), #final demultiplexing stats
   		expand(DEMULTIPLEX_DIR + "qc/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.demul.QC1.txt",filtered_product, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS) if config["DEMULTIPLEX"]  else expand( DUMMY_DIR + "dummy{number}.txt", number='10'), # QC on demultiplexing (checks that there is only one barcode assigned per read
   		expand(DEMULTIPLEX_DIR + "qc/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.demul.QC2.txt", filtered_product, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS) if config["DEMULTIPLEX"]  else expand( DUMMY_DIR + "dummy{number}.txt", number='11'), # QC on demultiplexing (checks that there is only one barcode assigned per read
   		expand(config["PLOTSDIR"] + "{techname}Corr{corrLevel}.demultiplexing.perSample.stats.{ext}", techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, ext=config["PLOTFORMATS"])  if config["DEMULTIPLEX"]  else expand( DUMMY_DIR + "dummy{number}.txt", number='12'),
  		expand(config["PLOTSDIR"] + "{techname}Corr{corrLevel}.mapping.perSample.perFraction.stats.{ext}", techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, ext=config["PLOTFORMATS"]),
  		expand("mappings/nonAnchoredMergeReads/cage+polyASupported/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.cage+polyASupported.bed", filtered_product_merge, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACSpluSMERGED, barcodes=BARCODESpluSMERGED),
		expand(config["PLOTSDIR"] + "cagePolyASupport.stats/{capDesign}_Corr{corrLevel}_{sizeFrac}_{barcodes}.cagePolyASupport.stats.{ext}", filtered_merge_figures, capDesign=CAPDESIGNSplusMERGED, corrLevel=FINALCORRECTIONLEVELSplusMERGED, sizeFrac=PLOTSbySIZEFRAC, barcodes=PLOTSbyTISSUE, ext=config["PLOTFORMATS"]),
		expand(config["PLOTSDIR"] + "tmerge.vs.gencode.stats/{capDesign}_Corr{corrLevel}_{sizeFrac}_{barcodes}.tmerge.vs.gencode.stats.{ext}", filtered_merge_figures, capDesign=CAPDESIGNSplusMERGED, corrLevel=FINALCORRECTIONLEVELSplusMERGED, sizeFrac=PLOTSbySIZEFRAC, barcodes=PLOTSbyTISSUE, ext=config["PLOTFORMATS"]),
		expand(config["PLOTSDIR"] + "targetCoverage.stats/{capDesign}_Corr{corrLevel}_{sizeFrac}_{barcodes}.targetCoverage.stats.{ext}", filtered_merge_figures, capDesign=CAPDESIGNSplusMERGED, corrLevel=FINALCORRECTIONLEVELSplusMERGED, sizeFrac=PLOTSbySIZEFRAC, barcodes=PLOTSbyTISSUE, ext=config["PLOTFORMATS"]),
		expand(config["PLOTSDIR"] + "spliceSites.stats/{capDesign}_Corr{corrLevel}_{sizeFrac}_{barcodes}.spliceSites.stats.{ext}", filtered_merge_figures, capDesign=CAPDESIGNSplusMERGED, corrLevel=FINALCORRECTIONLEVELSplusMERGED, sizeFrac=PLOTSbySIZEFRAC, barcodes=PLOTSbyTISSUE, ext=config["PLOTFORMATS"]),
		expand(config["PLOTSDIR"] + "splicedLength.stats/{capDesign}_Corr{corrLevel}_{sizeFrac}_{barcodes}.splicedLength.stats.{ext}", filtered_merge_figures, capDesign=CAPDESIGNSplusMERGED, corrLevel=FINALCORRECTIONLEVELSplusMERGED, sizeFrac=PLOTSbySIZEFRAC, barcodes=PLOTSbyTISSUE, ext=config["PLOTFORMATS"]),
		expand(config["PLOTSDIR"] + "polyAreads.stats/{capDesign}_Corr{corrLevel}_{sizeFrac}_{barcodes}.polyAreads.stats.{ext}", filtered_merge_figures, capDesign=CAPDESIGNSplusMERGED, corrLevel=FINALCORRECTIONLEVELSplusMERGED, sizeFrac=PLOTSbySIZEFRAC, barcodes=PLOTSbyTISSUE, ext=config["PLOTFORMATS"]),
		expand(config["PLOTSDIR"] + "tmerge.novelLoci.stats/{capDesign}_Corr{corrLevel}_{sizeFrac}_{barcodes}.tmerge.novelLoci.stats.{ext}", filtered_merge_figures, capDesign=CAPDESIGNSplusMERGED, corrLevel=FINALCORRECTIONLEVELSplusMERGED, sizeFrac=PLOTSbySIZEFRAC, barcodes=PLOTSbyTISSUE, ext=config["PLOTFORMATS"]),
		expand(config["PLOTSDIR"] + "tmerge.vs.Gencode.SJs.stats/{capDesign}_Corr{corrLevel}_{sizeFrac}_{barcodes}.tmerge.vs.Gencode.SJs.stats.{ext}", filtered_merge_figures, capDesign=CAPDESIGNSplusMERGED, corrLevel=FINALCORRECTIONLEVELSplusMERGED, sizeFrac=PLOTSbySIZEFRAC, barcodes=PLOTSbyTISSUE, ext=config["PLOTFORMATS"]),
		expand(config["PLOTSDIR"] + "tmerge.vs.SIRVs.stats/{capDesign}_Corr{corrLevel}_{sizeFrac}_{barcodes}.tmerge.vs.SIRVs.stats.{ext}", filtered_merge_figures, capDesign=CAPDESIGNSplusMERGED, corrLevel=FINALCORRECTIONLEVELSplusMERGED, sizeFrac=PLOTSbySIZEFRAC, barcodes=PLOTSbyTISSUE, ext=config["PLOTFORMATS"]) if SIRVpresent  else expand( DUMMY_DIR + "dummy{number}.txt", number='13'),
		expand(config["PLOTSDIR"] + "tmerge.vs.SIRVs.detection.stats/{capDesign}_Corr{corrLevel}_{sizeFrac}_{barcodes}.tmerge.vs.SIRVs.detection.stats.{ext}", filtered_merge_figures, capDesign=CAPDESIGNSplusMERGED, corrLevel=FINALCORRECTIONLEVELSplusMERGED, sizeFrac=PLOTSbySIZEFRAC, barcodes=PLOTSbyTISSUE, ext=config["PLOTFORMATS"]) if SIRVpresent  else expand( DUMMY_DIR + "dummy{number}.txt", number='14'),

		config["TRACK_HUB_DIR"] + "hub.txt",
		config["TRACK_HUB_DIR"] + "genomes.txt",
		expand(config["TRACK_HUB_DIR"] + "{genome}/trackDb.txt", genome=GENOMES),


# "Dummy" rule to skip undesired targets depending on config

rule dummy:
	input: "/dev/null"
	output: DUMMY_DIR + "dummy{number}.txt"
	shell:
		'''
touch {output}
		'''
