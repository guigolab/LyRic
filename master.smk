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
GENOMETOPAS=config["genomeToPAS"]
DUMMY_DIR="dummy/"
FQPATH=config["FQPATH"]
FQ_CORR_PATH=FQPATH + "corr/"
TRACK_HUB_DATA_URL=config["TRACK_HUB_BASE_URL"] + "dataFiles/"
GENOMETOPREVIOUS=config["genomeToPreviousPhaseTMs"]
print ("Input LR FASTQs are in: " + FQPATH)

ENDSUPPORTcategories=["all", "cage+polyASupported"]

GENOMES=[]
GENOMETOCAPDESIGNS=defaultdict(list)
for capD in CAPDESIGNTOGENOME:
	genome=CAPDESIGNTOGENOME[capD]
	GENOMES.append(genome)
	GENOMETOCAPDESIGNS[genome].append(capD)
GENOMES=set(GENOMES)

MERGEDCAPDESIGNS=defaultdict(list)
for genome in GENOMETOCAPDESIGNS:
	mergedCapDString=''
	capDperGenome=0
	for capD in sorted(GENOMETOCAPDESIGNS[genome]):
		if capD.find('preCap') == -1:  #avoid merging pre-capture with post-capture
			mergedCapDString+=capD
			capDperGenome+=1
	if capDperGenome > 1:
		for capD in GENOMETOCAPDESIGNS[genome]:
			MERGEDCAPDESIGNS[mergedCapDString].append(capD)
			CAPDESIGNTOCAGEPEAKS[mergedCapDString]=CAPDESIGNTOCAGEPEAKS[capD]
			CAPDESIGNTOANNOTGTF[mergedCapDString]=CAPDESIGNTOANNOTGTF[capD]
		CAPDESIGNTOGENOME[mergedCapDString]=genome


GENOMETOCAPDESIGNSplusMERGED=defaultdict(list)
for capD in CAPDESIGNTOGENOME:
	genome=CAPDESIGNTOGENOME[capD]
	GENOMETOCAPDESIGNSplusMERGED[genome].append(capD)

print ("Genomes:")
print (GENOMES)
# no underscores allowed in wildcards, to avoid greedy matching since we use them as separators

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

wildcard_constraints:
 	capDesign = "[^_/]+",
 	sizeFrac = "[^_/]+",
 	techname = "[^_/]+",
 	corrLevel = "[^_/]+",


BARCODESpluSMERGED=set(BARCODES)
BARCODESpluSMERGED.add("allTissues")
SIZEFRACSpluSMERGED=set(SIZEFRACS)
SIZEFRACSpluSMERGED.add("allFracs")
CAPDESIGNSplusMERGED=set(CAPDESIGNS)
for capD in MERGEDCAPDESIGNS:
	CAPDESIGNSplusMERGED.add(capD)
TECHNAMESplusMERGED=set(TECHNAMES)
TECHNAMESplusMERGED.add("allSeqTechs")


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
	print ("LR error correction is OFF. Correction levels: ")
	FINALCORRECTIONLEVELS=["0",]
	lastK='NULL'
	print (FINALCORRECTIONLEVELS)





########################
### make list of authorized wildcard combinations
### inspired by https://stackoverflow.com/questions/41185567/how-to-use-expand-in-snakemake-when-some-particular-combinations-of-wildcards-ar
###

AUTHORIZEDCOMBINATIONS = []
AUTHORIZEDCOMBINATIONSMERGE = []

for comb in itertools.product(TECHNAMES,CAPDESIGNS,SIZEFRACS,BARCODES):
	if(os.path.isfile(FQPATH + comb[0] + "_" + comb[1] + "_" + comb[2] + ".fastq.gz") or
	os.path.isfile(FQPATH + comb[0] + "_" + comb[1] + "_" + comb[2]  + "." + comb[3] + ".fastq.gz")): #allow only combinations corresponding to existing FASTQs
		AUTHORIZEDCOMBINATIONS.append((("techname", comb[0]),("capDesign", comb[1])))
		AUTHORIZEDCOMBINATIONS.append((("techname", comb[0]),("capDesign", comb[1]),("sizeFrac", comb[2])))
		# for mergedTerm in (MERGEDTERMS):
		for mergedComb in itertools.product([("techname", comb[0]), ("techname", "byTech"), ("techname", "allSeqTechs")], [("capDesign", comb[1]), ("capDesign", "byCapDesign")], [("sizeFrac", comb[2]), ("sizeFrac", "allFracs"), ("sizeFrac", "byFrac")]):
			AUTHORIZEDCOMBINATIONSMERGE.append((mergedComb))

		for corrLevel in FINALCORRECTIONLEVELS:
			AUTHORIZEDCOMBINATIONS.append((("techname", comb[0]),("corrLevel", corrLevel),("capDesign", comb[1]),("sizeFrac", comb[2])))
			for mergedComb in itertools.product([("techname", comb[0]), ("techname", "byTech"), ("techname", "allSeqTechs")], [("corrLevel", corrLevel), ("corrLevel", "byCorr")] , [("capDesign", comb[1]), ("capDesign", "byCapDesign")], [("sizeFrac", comb[2]), ("sizeFrac", "allFracs"), ("sizeFrac", "byFrac")]):
				AUTHORIZEDCOMBINATIONSMERGE.append((mergedComb))
		for ext in config["PLOTFORMATS"]:
			for corrLevel in FINALCORRECTIONLEVELS:
				AUTHORIZEDCOMBINATIONS.append((("techname", comb[0]),("corrLevel", corrLevel),("capDesign", comb[1]),("sizeFrac", comb[2]),("ext",ext)))
				for mergedComb in itertools.product([("techname", comb[0]), ("techname", "byTech"), ("techname", "allSeqTechs")], [("corrLevel", corrLevel), ("corrLevel", "byCorr")] , [("capDesign", comb[1]), ("capDesign", "byCapDesign")], [("sizeFrac", comb[2]), ("sizeFrac", "allFracs"), ("sizeFrac", "byFrac")], [("ext",ext)]):
					AUTHORIZEDCOMBINATIONSMERGE.append((mergedComb))

		if (comb[3]).find(comb[1]) == 0 or config["DEMULTIPLEX"] is False: # check that barcode ID's prefix matches capDesign (i.e. ignore demultiplexed FASTQs with unmatched barcode/capDesign)
			AUTHORIZEDCOMBINATIONS.append((("techname", comb[0]),("capDesign", comb[1]),("barcodes", comb[3])))

			if config["LORDEC_CORRECT"]:
				for split in splitFasta:
					AUTHORIZEDCOMBINATIONS.append((("techname", comb[0]),("capDesign", comb[1]),("sizeFrac", comb[2]),("barcodes",comb[3]), ("split", split)))
			for corrLevel in FINALCORRECTIONLEVELS:
				AUTHORIZEDCOMBINATIONS.append((("techname", comb[0]),("corrLevel", corrLevel),("capDesign", comb[1]),("sizeFrac", comb[2]),("barcodes",comb[3])))
				AUTHORIZEDCOMBINATIONS.append((("techname", comb[0]),("corrLevel", corrLevel),("capDesign", comb[1]),("barcodes",comb[3])))
				for mergedComb in itertools.product([("techname", comb[0]), ("techname", "byTech"), ("techname", "allSeqTechs")], [("corrLevel", corrLevel), ("corrLevel", "byCorr")] , [("capDesign", comb[1]), ("capDesign", "byCapDesign")], [("sizeFrac", comb[2]), ("sizeFrac", "allFracs"), ("sizeFrac", "byFrac")], [("barcodes",comb[3]), ("barcodes", "allTissues"), ("barcodes", "byTissue")]):
					AUTHORIZEDCOMBINATIONSMERGE.append((mergedComb))

				for mergedComb in itertools.product([("techname", comb[0]), ("techname", "byTech"), ("techname", "allSeqTechs")], [("corrLevel", corrLevel), ("corrLevel", "byCorr")] , [("capDesign", comb[1]), ("capDesign", "byCapDesign")], [("barcodes",comb[3]), ("barcodes", "allTissues"), ("barcodes", "byTissue")]):
					AUTHORIZEDCOMBINATIONSMERGE.append((mergedComb))

				for minReadSupport in config["MINIMUM_TMERGE_READ_SUPPORT"]:
					AUTHORIZEDCOMBINATIONS.append((("techname", comb[0]),("corrLevel", corrLevel),("capDesign", comb[1]),("sizeFrac", comb[2]),("barcodes",comb[3]), ("minReadSupport", minReadSupport)))
					for mergedComb in itertools.product([("techname", comb[0]), ("techname", "byTech"), ("techname", "allSeqTechs")], [("corrLevel", corrLevel), ("corrLevel", "byCorr")] , [("capDesign", comb[1]), ("capDesign", "byCapDesign")], [("barcodes",comb[3]), ("barcodes", "allTissues"), ("barcodes", "byTissue")], [("minReadSupport", minReadSupport),]):
						AUTHORIZEDCOMBINATIONSMERGE.append((mergedComb))
					for mergedComb in itertools.product([("techname", comb[0]), ("techname", "byTech"), ("techname", "allSeqTechs")], [("corrLevel", corrLevel), ("corrLevel", "byCorr")] , [("capDesign", comb[1]), ("capDesign", "byCapDesign")], [("sizeFrac", comb[2]), ("sizeFrac", "allFracs"), ("sizeFrac", "byFrac")], [("barcodes",comb[3]), ("barcodes", "allTissues"), ("barcodes", "byTissue")], [("minReadSupport", minReadSupport),]):
						AUTHORIZEDCOMBINATIONSMERGE.append((mergedComb))
				for endSupport in ENDSUPPORTcategories:
					for mergedComb in itertools.product([("techname", comb[0]), ("techname", "byTech"), ("techname", "allSeqTechs")], [("corrLevel", corrLevel), ("corrLevel", "byCorr")] , [("capDesign", comb[1]), ("capDesign", "byCapDesign")], [("barcodes",comb[3]), ("barcodes", "allTissues"), ("barcodes", "byTissue")], [("endSupport", endSupport),]):
						AUTHORIZEDCOMBINATIONSMERGE.append((mergedComb))
					for mergedComb in itertools.product([("techname", comb[0]), ("techname", "byTech"), ("techname", "allSeqTechs")], [("corrLevel", corrLevel), ("corrLevel", "byCorr")] , [("capDesign", comb[1]), ("capDesign", "byCapDesign")], [("sizeFrac", comb[2]), ("sizeFrac", "allFracs"), ("sizeFrac", "byFrac")], [("barcodes",comb[3]), ("barcodes", "allTissues"), ("barcodes", "byTissue")], [("endSupport", endSupport),]):
						AUTHORIZEDCOMBINATIONSMERGE.append((mergedComb))
					for minReadSupport in config["MINIMUM_TMERGE_READ_SUPPORT"]:
						for mergedComb in itertools.product([("techname", comb[0]), ("techname", "byTech"), ("techname", "allSeqTechs")], [("corrLevel", corrLevel), ("corrLevel", "byCorr")] , [("capDesign", comb[1]), ("capDesign", "byCapDesign")], [("barcodes",comb[3]), ("barcodes", "allTissues"), ("barcodes", "byTissue")], [("endSupport", endSupport),],  [("minReadSupport", minReadSupport),]):
							AUTHORIZEDCOMBINATIONSMERGE.append((mergedComb))
						for mergedComb in itertools.product([("techname", comb[0]), ("techname", "byTech"), ("techname", "allSeqTechs")], [("corrLevel", corrLevel), ("corrLevel", "byCorr")] , [("capDesign", comb[1]), ("capDesign", "byCapDesign")], [("sizeFrac", comb[2]), ("sizeFrac", "allFracs"), ("sizeFrac", "byFrac")], [("barcodes",comb[3]), ("barcodes", "allTissues"), ("barcodes", "byTissue")], [("endSupport", endSupport),],  [("minReadSupport", minReadSupport),]):
							AUTHORIZEDCOMBINATIONSMERGE.append((mergedComb))

			# for strand in config["STRANDS"]:
			# 	for corrLevel in FINALCORRECTIONLEVELS:
			# 		AUTHORIZEDCOMBINATIONS.append((("techname", comb[0]),("corrLevel", corrLevel),("capDesign", comb[1]),("barcodes",comb[3]),("strand", strand)))
			# 		AUTHORIZEDCOMBINATIONS.append((("techname", comb[0]),("corrLevel", corrLevel),("capDesign", comb[1]),("sizeFrac", comb[2]),("barcodes",comb[3]),("strand", strand)))
			# 		for mergedComb in itertools.product([("techname", comb[0]), ("techname", "byTech"), ("techname", "allSeqTechs")], [("corrLevel", corrLevel), ("corrLevel", "byCorr")] , [("capDesign", comb[1]), ("capDesign", "byCapDesign")], [("barcodes",comb[3]), ("barcodes", "allTissues"), ("barcodes", "byTissue")], [("strand", strand)]):
			# 			AUTHORIZEDCOMBINATIONSMERGE.append((mergedComb))
			# 		for mergedComb in itertools.product([("techname", comb[0]), ("techname", "byTech"), ("techname", "allSeqTechs")], [("corrLevel", corrLevel), ("corrLevel", "byCorr")] , [("capDesign", comb[1]), ("capDesign", "byCapDesign")], [("sizeFrac", comb[2]), ("sizeFrac", "allFracs"), ("sizeFrac", "byFrac")], [("barcodes",comb[3]), ("barcodes", "allTissues"), ("barcodes", "byTissue")], [("strand", strand)]):
			# 			AUTHORIZEDCOMBINATIONSMERGE.append((mergedComb))
			for ext in config["PLOTFORMATS"]:
				for corrLevel in FINALCORRECTIONLEVELS:
					AUTHORIZEDCOMBINATIONS.append((("techname", comb[0]),("corrLevel", corrLevel),("capDesign", comb[1]),("barcodes",comb[3]),("ext", ext)))
					for mergedComb in itertools.product([("techname", comb[0]), ("techname", "byTech"), ("techname", "allSeqTechs")], [("corrLevel", corrLevel), ("corrLevel", "byCorr")] , [("capDesign", comb[1]), ("capDesign", "byCapDesign")], [("barcodes",comb[3]), ("barcodes", "allTissues"), ("barcodes", "byTissue")], [("ext", ext)]):
						AUTHORIZEDCOMBINATIONSMERGE.append((mergedComb))
					for mergedComb in itertools.product([("techname", comb[0]), ("techname", "byTech"), ("techname", "allSeqTechs")], [("corrLevel", corrLevel), ("corrLevel", "byCorr")] , [("capDesign", comb[1]), ("capDesign", "byCapDesign")], [("sizeFrac", comb[2]), ("sizeFrac", "allFracs"), ("sizeFrac", "byFrac")], [("barcodes",comb[3]), ("barcodes", "allTissues"), ("barcodes", "byTissue")], [("ext", ext)]):
						AUTHORIZEDCOMBINATIONSMERGE.append((mergedComb))
					for minReadSupport in config["MINIMUM_TMERGE_READ_SUPPORT"]:
						for mergedComb in itertools.product([("techname", comb[0]), ("techname", "byTech"), ("techname", "allSeqTechs")], [("corrLevel", corrLevel), ("corrLevel", "byCorr")] , [("capDesign", comb[1]), ("capDesign", "byCapDesign")], [("barcodes",comb[3]), ("barcodes", "allTissues"), ("barcodes", "byTissue")], [("ext", ext)], [("minReadSupport", minReadSupport),]):
							AUTHORIZEDCOMBINATIONSMERGE.append((mergedComb))
						for mergedComb in itertools.product([("techname", comb[0]), ("techname", "byTech"), ("techname", "allSeqTechs")], [("corrLevel", corrLevel), ("corrLevel", "byCorr")] , [("capDesign", comb[1]), ("capDesign", "byCapDesign")], [("sizeFrac", comb[2]), ("sizeFrac", "allFracs"), ("sizeFrac", "byFrac")], [("barcodes",comb[3]), ("barcodes", "allTissues"), ("barcodes", "byTissue")], [("ext", ext)], [("minReadSupport", minReadSupport),]):
							AUTHORIZEDCOMBINATIONSMERGE.append((mergedComb))
					for endSupport in ENDSUPPORTcategories:
						for mergedComb in itertools.product([("techname", comb[0]), ("techname", "byTech"), ("techname", "allSeqTechs")], [("corrLevel", corrLevel), ("corrLevel", "byCorr")] , [("capDesign", comb[1]), ("capDesign", "byCapDesign")], [("barcodes",comb[3]), ("barcodes", "allTissues"), ("barcodes", "byTissue")], [("endSupport", endSupport),], [("ext", ext)]):
							AUTHORIZEDCOMBINATIONSMERGE.append((mergedComb))
						for mergedComb in itertools.product([("techname", comb[0]), ("techname", "byTech"), ("techname", "allSeqTechs")], [("corrLevel", corrLevel), ("corrLevel", "byCorr")] , [("capDesign", comb[1]), ("capDesign", "byCapDesign")], [("sizeFrac", comb[2]), ("sizeFrac", "allFracs"), ("sizeFrac", "byFrac")], [("barcodes",comb[3]), ("barcodes", "allTissues"), ("barcodes", "byTissue")], [("endSupport", endSupport),], [("ext", ext)]):
							AUTHORIZEDCOMBINATIONSMERGE.append((mergedComb))
						for minReadSupport in config["MINIMUM_TMERGE_READ_SUPPORT"]:
							for mergedComb in itertools.product([("techname", comb[0]), ("techname", "byTech"), ("techname", "allSeqTechs")], [("corrLevel", corrLevel), ("corrLevel", "byCorr")] , [("capDesign", comb[1]), ("capDesign", "byCapDesign")], [("sizeFrac", comb[2]), ("sizeFrac", "allFracs"), ("sizeFrac", "byFrac")], [("barcodes",comb[3]), ("barcodes", "allTissues"), ("barcodes", "byTissue")], [("endSupport", endSupport),], [("ext", ext)], [("minReadSupport", minReadSupport),]):
								AUTHORIZEDCOMBINATIONSMERGE.append((mergedComb))

				#for merged figures:
#				for mergedComb in itertools.product([("capDesign", comb[1]), ("capDesign", "byCapDesign")], [("sizeFrac", comb[2]), ("sizeFrac", "allFracs"), ("sizeFrac", "byFrac")], [("barcodes",comb[3]), ("barcodes", "allTissues"), ("barcodes", "byTissue")], [("ext", ext)]):
#					AUTHORIZEDCOMBINATIONSMERGE.append((mergedComb))

			for splice in SPLICE_SITE_TYPES:
				for corrLevel in FINALCORRECTIONLEVELS:
					AUTHORIZEDCOMBINATIONS.append((("techname", comb[0]),("corrLevel", corrLevel),("capDesign", comb[1]),("sizeFrac", comb[2]),("barcodes",comb[3]),("spliceType",splice)))
					for mergedComb in itertools.product([("techname", comb[0]), ("techname", "byTech"), ("techname", "allSeqTechs")], [("corrLevel", corrLevel), ("corrLevel", "byCorr")] , [("capDesign", comb[1]), ("capDesign", "byCapDesign")], [("sizeFrac", comb[2]), ("sizeFrac", "allFracs"), ("sizeFrac", "byFrac")], [("barcodes",comb[3]), ("barcodes", "allTissues"), ("barcodes", "byTissue")], [("spliceType",splice)]):
						AUTHORIZEDCOMBINATIONSMERGE.append((mergedComb))
					for minReadSupport in config["MINIMUM_TMERGE_READ_SUPPORT"]:
						AUTHORIZEDCOMBINATIONS.append((("techname", comb[0]),("corrLevel", corrLevel),("capDesign", comb[1]),("sizeFrac", comb[2]),("barcodes",comb[3]),("spliceType",splice),("minReadSupport", minReadSupport)))
						for mergedComb in itertools.product([("techname", comb[0]), ("techname", "byTech"), ("techname", "allSeqTechs")], [("corrLevel", corrLevel), ("corrLevel", "byCorr")] , [("capDesign", comb[1]), ("capDesign", "byCapDesign")], [("sizeFrac", comb[2]), ("sizeFrac", "allFracs"), ("sizeFrac", "byFrac")], [("barcodes",comb[3]), ("barcodes", "allTissues"), ("barcodes", "byTissue")], [("spliceType",splice)], [("minReadSupport",minReadSupport)]):
							AUTHORIZEDCOMBINATIONSMERGE.append((mergedComb))





AUTHORIZEDCOMBINATIONS=set(AUTHORIZEDCOMBINATIONS)
AUTHORIZEDCOMBINATIONSMERGE=set(AUTHORIZEDCOMBINATIONSMERGE)


PLOTSbyTISSUE=set(BARCODESpluSMERGED)
PLOTSbyTISSUE.add("byTissue")
PLOTSbyTISSUEnoALL=set(BARCODES)
PLOTSbyTISSUEnoALL.add("byTissue")
PLOTSbySIZEFRAC=set(SIZEFRACSpluSMERGED)
PLOTSbySIZEFRAC.add("byFrac")
PLOTSbySIZEFRACnoALL=set(SIZEFRACS)
PLOTSbySIZEFRACnoALL.add("byFrac")

PLOTSbySEQTECH=set(TECHNAMESplusMERGED)
PLOTSbySEQTECH.add("byTech")
PLOTSbySEQTECHnoALLTECHS=set(TECHNAMES)
PLOTSbySEQTECHnoALLTECHS.add("byTech")

PLOTSbyCORRLEVEL=set(FINALCORRECTIONLEVELS)
PLOTSbyCORRLEVEL.add("byCorr")
PLOTSbyCAPDESIGN=set(CAPDESIGNS)
PLOTSbyCAPDESIGN.add("byCapDesign")
PLOTSbyCAPDESIGNplusMERGED=set(PLOTSbyCAPDESIGN)
for capD in MERGEDCAPDESIGNS:
	PLOTSbyCAPDESIGNplusMERGED.add(capD)

def filtered_product(*args):
	found=False
	for wc_comb in itertools.product(*args):
		#print(wc_comb)
		if wc_comb in AUTHORIZEDCOMBINATIONS:
			found=True
			yield(wc_comb)
		else:
			if len(wc_comb) > 4:
				if wc_comb[2][0] in ('capDesign') and wc_comb[2][1] in MERGEDCAPDESIGNS and wc_comb[3] in (('sizeFrac', 'allFracs'),) and wc_comb[4] in (('barcodes', 'allTissues'),):
					found=True
					yield(wc_comb)
	if not found:
		print(" Error in function filtered_product. Args were:")
		print((args))
		quit(" Error. Could not yield any input file.")

def filtered_product_merge(*args):
	found=False
	#print("NEW COMB")
	#print(args)

	for wc_comb in itertools.product(*args):
		#print(wc_comb)
		if wc_comb in AUTHORIZEDCOMBINATIONSMERGE:
			#print("FPM ")
			#print (wc_comb)
		#	print ("AUTH")
			if wc_comb[0] in (('techname', 'allSeqTechs'),):
				if wc_comb[3] in (('sizeFrac', 'allFracs'),) and wc_comb[4] in (('barcodes', 'allTissues'),):
					found=True
					#print("YIELD 1")
					yield(wc_comb)

			else:
				found=True
				#print("YIELD 2")
				yield(wc_comb)
		else:
			if wc_comb[2][0] in ('capDesign') and wc_comb[2][1] in MERGEDCAPDESIGNS and wc_comb[3] in (('sizeFrac', 'allFracs'),) and wc_comb[4] in (('barcodes', 'allTissues'),):
				found=True
				#print("YIELD 3")
				yield(wc_comb)
	if not found:
		print(" Error in function filtered_product_merge. Args were:")
		print((args))
		quit(" Error. Could not yield any input file.")


def filtered_capDesign_product_merge(*args):
	found=False
	for wc_comb in itertools.product(*args):
		#print(wc_comb)
		if wc_comb in AUTHORIZEDCOMBINATIONSMERGE:
			#print("FPM ")
			#print (wc_comb)
		#	print ("AUTH")
			if wc_comb[0] in (('techname', 'allSeqTechs'),):
				if wc_comb[3] in (('sizeFrac', 'allFracs'),) and wc_comb[4] in (('barcodes', 'allTissues'),):
					found=True
					yield(wc_comb)

			else:
				found=True
				yield(wc_comb)
		else:
			if wc_comb[2][0] in ('capDesign') and wc_comb[2][1] in MERGEDCAPDESIGNS and wc_comb[3] in (('sizeFrac', 'allFracs'),) and wc_comb[4] in (('barcodes', 'allTissues'),):
				#print ("OKcapDesign_product_merge")
				for capD in MERGEDCAPDESIGNS[wc_comb[2][1]]:
					wc_comb2=list(wc_comb)
					wc_comb2[2]=('capDesign', capD)
					wc_comb2=tuple(wc_comb2)
					#print (wc_comb2)
					found=True
					yield(wc_comb2)
	if not found:
		print(" Error in function filtered_capDesign_product_merge. Args were:")
		print((args))
		quit(" Error. Could not yield any input file.")



def filtered_merge_figures(*args):
	found=False
	for wc_comb in itertools.product(*args):
		#print (wc_comb)
		yieldWC=False
		if wc_comb in AUTHORIZEDCOMBINATIONSMERGE:
			if wc_comb[1] not in (('corrLevel', '90'),):
		#	print ("PRE AUTH ")
				if wc_comb[2] in (('capDesign', 'byCapDesign'),):
					if wc_comb[4] not in (('barcodes', 'byTissue'),): #DO NOT plot all tissues and capDesigns on the same axis (too big)
		#			print ("AUTH 2")
						yieldWC=True
					else:
					#print ("NO AUTH 1")
						continue
				if wc_comb[0] in (('techname', 'allSeqTechs'),):
					if wc_comb[3] in (('sizeFrac', 'allFracs'),) and wc_comb[4] in (('barcodes', 'allTissues'),):
						yieldWC=True
					else:
					#print ("NO AUTH 2")
						continue
				else:
		#		print ("AUTH 1")
					yieldWC=True
		else:
			if wc_comb[2][0] in ('capDesign') and wc_comb[2][1] in MERGEDCAPDESIGNS and wc_comb[3] in (('sizeFrac', 'allFracs'),) and wc_comb[4] in (('barcodes', 'allTissues'),):
				yieldWC=True
		if yieldWC == True:
			#print("AUTH")
			yield(wc_comb)
			found=True
	if not found:
		print(" Error in function filtered_merge_figures. Args were:")
		print((args))
		quit(" Error. Could not yield any input file.")


def merge_figures_params(c,bf,bt,cl, tn):
	clBool=None
	if len(FINALCORRECTIONLEVELS)>1 and cl == FINALCORRECTIONLEVELS[1]:
		clBool="Yes"
	elif cl == FINALCORRECTIONLEVELS[0]:
		clBool="No"
	else:
		clBool=cl
	seqTechFilterString=''
	if tn not in ('byTech'):
		seqTechFilterString="dat <- subset(dat, seqTech=='" + tn + "')\n"
	else:
		seqTechFilterString="dat <- subset(dat, seqTech!='allSeqTechs')\n"
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
	if c not in ('byCapDesign'):
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
	return(capDesignFilterString, corrLevelFilterString, sizeFracFilterString, tissueFilterString, substSeqTechString, substTissueString,  xlabString, hideXaxisLabels , graphDimensions, hideYaxisLabels, seqTechFilterString)


def trackHubSubGroupString(tn, cd, sf, bc, cl, m):
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
	minRS=[]
	for minReadSupport in m:
		minRS.append(("minReadSupport", minReadSupport))
	returnSubGroup1String="subGroup1 sample Sample"
	returnSubGroup2String="subGroup2 sizeFraction Size_fraction"
	returnSubGroup3String="subGroup3 corrLevel LoRDEC_sequence_correction"
	returnSubGroup4String="subGroup4 minReadSupport Min_read_support_per_TM"

	returnSubGroup1StringData=[]
	returnSubGroup2StringData=[]
	returnSubGroup3StringData=[]
	returnSubGroup4StringData=[]
	for wc_comb in itertools.product(techname, corrLevels, capDesign, sizeFracs, barCodes, minRS):
		#print("THSS ")
		#print(wc_comb)
		if wc_comb in AUTHORIZEDCOMBINATIONSMERGE:
		#	print ("AUTH")
			returnSubGroup1StringData.append(wc_comb[4][1] + "=" + wc_comb[4][1])
			returnSubGroup2StringData.append(wc_comb[3][1] + "=" + wc_comb[3][1])
			returnSubGroup3StringData.append(wc_comb[1][1] + "=" + wc_comb[1][1])
			returnSubGroup4StringData.append(wc_comb[5][1] + "=" + wc_comb[5][1])

	#print ("returnSubGroup2StringData list")
	#print (returnSubGroup2StringData)
	returnSubGroup1StringData=set(returnSubGroup1StringData)
	returnSubGroup2StringData=set(returnSubGroup2StringData)
	returnSubGroup3StringData=set(returnSubGroup3StringData)
	returnSubGroup4StringData=set(returnSubGroup4StringData)
	#print ("returnSubGroup2StringData set")
	#print (returnSubGroup2StringData)
	return(returnSubGroup1String + " " + " ".join(returnSubGroup1StringData) + "\n" + returnSubGroup2String + " " + " ".join(returnSubGroup2StringData) + "\n" + returnSubGroup3String + " " + " ".join(returnSubGroup3StringData)+ "\n" + returnSubGroup4String + " " + " ".join(returnSubGroup4StringData))

def trackChecked(t,c,cd,s,b,m):
	checked=''
	if t.find('pacBio') == 0 and (len(FINALCORRECTIONLEVELS)>1 and c == FINALCORRECTIONLEVELS[1]):
		checked='off'
	elif t.find('ont') == 0 and (len(FINALCORRECTIONLEVELS)>1 and c == FINALCORRECTIONLEVELS[1]):
		checked='off'
	else:
		if s == 'allFracs' and b == 'allTissues':
			checked='on'
		else:
			checked='off'
	if m == '1':
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
include: "lociEndSupport.py"

#pseudo-rule specifying the target files we ultimately want.
rule all:
	input:
		#expand("mappings/nonAnchoredMergeReads/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.all.gff", filtered_product_merge, techname=TECHNAMESplusMERGED, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNSplusMERGED, sizeFrac=SIZEFRACSpluSMERGED, barcodes=BARCODESpluSMERGED), #needed only for testing
		expand(FQ_CORR_PATH + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.fastq.gz", filtered_product, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS) if config["DEMULTIPLEX"] else expand(FQ_CORR_PATH + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.{barcodes}.fastq.gz", filtered_product, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES),

		expand(config["PLOTSDIR"] + "readLength.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.readLength.stats.{ext}", filtered_merge_figures, techname=PLOTSbySEQTECHnoALLTECHS, corrLevel=PLOTSbyCORRLEVEL, capDesign=PLOTSbyCAPDESIGN, sizeFrac=PLOTSbySIZEFRACnoALL, ext=config["PLOTFORMATS"]) if config["DEMULTIPLEX"] else expand(config["PLOTSDIR"] + "readLength.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.readLength.stats.{ext}", filtered_merge_figures, techname=PLOTSbySEQTECHnoALLTECHS, corrLevel=PLOTSbyCORRLEVEL, capDesign=PLOTSbyCAPDESIGN, sizeFrac=PLOTSbySIZEFRACnoALL, barcodes=PLOTSbyTISSUEnoALL, ext=config["PLOTFORMATS"]), # facetted histograms of read length
		expand(FQ_CORR_PATH + "qc/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.dupl.txt", filtered_product,techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS) if config["DEMULTIPLEX"] else expand(FQ_CORR_PATH + "qc/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.{barcodes}.dupl.txt", filtered_product,techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES),
		expand(config["PLOTSDIR"] + "{techname}Corr{corrLevel}.fastq.UP.stats.{ext}", techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, ext=config["PLOTFORMATS"]) if config["DEMULTIPLEX"] else expand( DUMMY_DIR + "dummy{number}.txt", number='1'), # UP reads plots
		expand(config["PLOTSDIR"] + "{techname}Corr{corrLevel}.fastq.BC.stats.{ext}", techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, ext=config["PLOTFORMATS"]) if config["DEMULTIPLEX"]  else expand( DUMMY_DIR + "dummy{number}.txt", number='2') , # barcode reads plots
		expand(config["PLOTSDIR"] + "{techname}Corr{corrLevel}.fastq.foreignBC.stats.{ext}", techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, ext=config["PLOTFORMATS"]) if config["DEMULTIPLEX"]  else expand( DUMMY_DIR + "dummy{number}.txt", number='3') , #foreign barcode reads plots
		expand(DEMULTIPLEX_DIR + "demultiplexFastqs/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.{barcodes}.fastq.gz", filtered_product, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES) if config["DEMULTIPLEX"]  else expand( DUMMY_DIR + "dummy{number}.txt", number='4'),
		expand("mappings/readMapping/qc/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bam.dupl.txt",filtered_product_merge, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACSpluSMERGED, barcodes=BARCODESpluSMERGED),
		expand(config["PLOTSDIR"] + "{techname}Corr{corrLevel}.ambiguousBarcodes.reads.stats.{ext}", techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, ext=config["PLOTFORMATS"]) if config["DEMULTIPLEX"]  else expand( DUMMY_DIR + "dummy{number}.txt", number='6'), # ambiguous barcodes plots
		expand(config["PLOTSDIR"] + "{techname}Corr{corrLevel}_{capDesign}.adapters.location.stats.{ext}",techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, ext=config["PLOTFORMATS"]) if config["DEMULTIPLEX"]  else expand( DUMMY_DIR + "dummy{number}.txt", number='7'), #location of adapters over reads
		expand(config["PLOTSDIR"] + "{techname}Corr{corrLevel}.chimeric.reads.stats.{ext}", techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, ext=config["PLOTFORMATS"]) if config["DEMULTIPLEX"]  else expand( DUMMY_DIR + "dummy{number}.txt", number='8'), # stats on chimeric reads
		expand(config["PLOTSDIR"] + "{techname}Corr{corrLevel}.finalDemul.reads.stats.{ext}", techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, ext=config["PLOTFORMATS"]) if config["DEMULTIPLEX"]  else expand( DUMMY_DIR + "dummy{number}.txt", number='9'), #final demultiplexing stats
		expand(DEMULTIPLEX_DIR + "qc/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.demul.QC1.txt",filtered_product, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS) if config["DEMULTIPLEX"]  else expand( DUMMY_DIR + "dummy{number}.txt", number='10'), # QC on demultiplexing (checks that there is only one barcode assigned per read
		expand(DEMULTIPLEX_DIR + "qc/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.demul.QC2.txt", filtered_product, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS) if config["DEMULTIPLEX"]  else expand( DUMMY_DIR + "dummy{number}.txt", number='11'), # QC on demultiplexing (checks that there is only one barcode assigned per read
		expand(config["PLOTSDIR"] + "{techname}Corr{corrLevel}.demultiplexing.perSample.stats.{ext}", techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, ext=config["PLOTFORMATS"])  if config["DEMULTIPLEX"]  else expand( DUMMY_DIR + "dummy{number}.txt", number='12'),
		expand(config["PLOTSDIR"] + "{techname}Corr{corrLevel}.mapping.perSample.perFraction.stats.{ext}", techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, ext=config["PLOTFORMATS"]),
		expand ("mappings/clusterPolyAsites/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.polyAsites.clusters.bed", filtered_product_merge, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACSpluSMERGED, barcodes=BARCODESpluSMERGED),
		expand(config["PLOTSDIR"] + "HiSS.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.stats.{ext}", filtered_merge_figures, techname=PLOTSbySEQTECHnoALLTECHS, corrLevel=PLOTSbyCORRLEVEL,  capDesign=PLOTSbyCAPDESIGN, sizeFrac=PLOTSbySIZEFRAC, barcodes=PLOTSbyTISSUE, ext=config["PLOTFORMATS"]),
		expand("mappings/nonAnchoredMergeReads/colored/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.{endSupport}.bed", filtered_product_merge, techname=TECHNAMESplusMERGED, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNSplusMERGED, sizeFrac=SIZEFRACSpluSMERGED, barcodes=BARCODESpluSMERGED, endSupport=ENDSUPPORTcategories, minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"]),

		expand(config["PLOTSDIR"] + "merged.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.merged.stats.{ext}", filtered_merge_figures, techname=PLOTSbySEQTECHnoALLTECHS, corrLevel=PLOTSbyCORRLEVEL,  capDesign=PLOTSbyCAPDESIGN, sizeFrac=PLOTSbySIZEFRAC, barcodes=PLOTSbyTISSUE, ext=config["PLOTFORMATS"], minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"]),
		expand(config["PLOTSDIR"] + "cagePolyASupport.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.cagePolyASupport.stats.{ext}", filtered_merge_figures, techname=PLOTSbySEQTECH, corrLevel=PLOTSbyCORRLEVEL,  capDesign=PLOTSbyCAPDESIGNplusMERGED, sizeFrac=PLOTSbySIZEFRAC, barcodes=PLOTSbyTISSUE, ext=config["PLOTFORMATS"], minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"]),
		expand(config["PLOTSDIR"] + "tmerge.vs.gencode.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.{endSupport}.vs.gencode.stats.{ext}", filtered_merge_figures, techname=PLOTSbySEQTECH, corrLevel=PLOTSbyCORRLEVEL,  capDesign=PLOTSbyCAPDESIGNplusMERGED, sizeFrac=PLOTSbySIZEFRAC, barcodes=PLOTSbyTISSUE, endSupport=ENDSUPPORTcategories, ext=config["PLOTFORMATS"], minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"]),
		expand(config["PLOTSDIR"] + "targetCoverage.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.targetCoverage.stats.{ext}", filtered_merge_figures, techname=PLOTSbySEQTECH, corrLevel=PLOTSbyCORRLEVEL,  capDesign=PLOTSbyCAPDESIGN, sizeFrac=PLOTSbySIZEFRAC, barcodes=PLOTSbyTISSUE, ext=config["PLOTFORMATS"], minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"]),
		expand(config["PLOTSDIR"] + "spliceSites.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.spliceSites.stats.{ext}", filtered_merge_figures, techname=PLOTSbySEQTECHnoALLTECHS, corrLevel=PLOTSbyCORRLEVEL,  capDesign=PLOTSbyCAPDESIGN, sizeFrac=PLOTSbySIZEFRAC, barcodes=PLOTSbyTISSUE, ext=config["PLOTFORMATS"], minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"]),
		expand(config["PLOTSDIR"] + "splicedLength.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.splicedLength.stats.{ext}", filtered_merge_figures, techname=PLOTSbySEQTECH, corrLevel=PLOTSbyCORRLEVEL,  capDesign=PLOTSbyCAPDESIGN, sizeFrac=PLOTSbySIZEFRAC, barcodes=PLOTSbyTISSUE, ext=config["PLOTFORMATS"], minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"]),
		expand(config["PLOTSDIR"] + "polyAreads.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.polyAreads.stats.{ext}", filtered_merge_figures, techname=PLOTSbySEQTECHnoALLTECHS, corrLevel=PLOTSbyCORRLEVEL,  capDesign=PLOTSbyCAPDESIGN, sizeFrac=PLOTSbySIZEFRAC, barcodes=PLOTSbyTISSUE, ext=config["PLOTFORMATS"]),
		expand(config["PLOTSDIR"] + "tmerge.novelLoci.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.{endSupport}.novelLoci.stats.{ext}", filtered_merge_figures, techname=PLOTSbySEQTECH, corrLevel=PLOTSbyCORRLEVEL,  capDesign=PLOTSbyCAPDESIGNplusMERGED, sizeFrac=PLOTSbySIZEFRAC, barcodes=PLOTSbyTISSUE, endSupport=ENDSUPPORTcategories, ext=config["PLOTFORMATS"], minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"]),
		expand(config["PLOTSDIR"] + "tmerge.vs.Gencode.SJs.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.vs.Gencode.SJs.stats.{ext}", filtered_merge_figures, techname=PLOTSbySEQTECH, corrLevel=PLOTSbyCORRLEVEL,  capDesign=PLOTSbyCAPDESIGN, sizeFrac=PLOTSbySIZEFRAC, barcodes=PLOTSbyTISSUE, ext=config["PLOTFORMATS"], minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"]),
		expand(config["PLOTSDIR"] + "tmerge.vs.SIRVs.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.vs.SIRVs.stats.{ext}", filtered_merge_figures, techname=PLOTSbySEQTECH, corrLevel=PLOTSbyCORRLEVEL,  capDesign=PLOTSbyCAPDESIGN, sizeFrac=PLOTSbySIZEFRAC, barcodes=PLOTSbyTISSUE, ext=config["PLOTFORMATS"], minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"]) if SIRVpresent  else expand( DUMMY_DIR + "dummy{number}.txt", number='13'),
		expand(config["PLOTSDIR"] + "tmerge.vs.SIRVs.detection.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.vs.SIRVs.detection.stats.{ext}", filtered_merge_figures, techname=PLOTSbySEQTECH, corrLevel=PLOTSbyCORRLEVEL,  capDesign=PLOTSbyCAPDESIGN, sizeFrac=PLOTSbySIZEFRAC, barcodes=PLOTSbyTISSUE, ext=config["PLOTFORMATS"], minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"]) if SIRVpresent  else expand( DUMMY_DIR + "dummy{number}.txt", number='14'),
		expand(config["PLOTSDIR"] + "FLloci.gencodeOnly.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.FLloci.gencodeOnly.stats.{ext}", filtered_merge_figures, techname='allSeqTechs', corrLevel=FINALCORRECTIONLEVELS,  capDesign=PLOTSbyCAPDESIGNplusMERGED, sizeFrac='allFracs', barcodes='allTissues', ext=config["PLOTFORMATS"], minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"]),
		expand(config["PLOTSDIR"] + "FLloci.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.FLloci.stats.{ext}", filtered_merge_figures, techname='allSeqTechs', corrLevel=FINALCORRECTIONLEVELS,  capDesign=PLOTSbyCAPDESIGNplusMERGED, sizeFrac='allFracs', barcodes='allTissues', ext=config["PLOTFORMATS"], minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"]),
		#################
		### Track hub ###
		#################
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
