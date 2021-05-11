import glob
from collections import defaultdict
import os
import itertools
import sys
import pandas as pd
from pprint import pprint


#prefix all shell scripts with the following command
shell.prefix("source ~/.bashrc; set +eu; conda deactivate;  set -euo pipefail; ")



GGPLOT_PUB_QUALITY="theme(axis.text= element_text(size=themeSize*1.8), axis.ticks = element_line(size=lineSize), axis.line = element_line(colour = '#595959', size=lineSize), axis.title=element_text(size = themeSize*2), panel.grid.major = element_line(colour='#d9d9d9', size=lineSize),panel.grid.minor = element_line(colour='#e6e6e6', size=minorLineSize),panel.border = element_blank(),panel.background = element_blank(), strip.background = element_rect(colour='#737373',fill='white'), legend.key.size=unit(0.5,'line'), legend.title=element_text(size=themeSize*1.2), legend.text=element_text(size=themeSize), strip.text = element_text(size = themeSize))"

sizeFrac_Rpalette="c('0+'='#b3b3b3', '0-1' ='#f765ac','1+' ='#b370f9')"
long_Rpalette="c('#8dd3c7','#ffffb3','#80ff80','#ff4dff','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd','#ccebc5','#ffed6f','#c4ff4d','#ff66a3')"
sampleAnnot_Rpalette="list(seqPlatform = c(ONT = '#8dd3c7', pacBioSI = '#ffffb3', pacBioSII = '#5d513a'), libraryPrep = c(CapTrap = '#ff4dff', SMARTer = '#bebada', Teloprime = '#fb8072', directRNA = '#80b1d3', 'Rt' = '#8c8c8c', PcrOnt='#3333ff'), tissue = c(Brain = '#b3de69', Heart = '#fccde5', K562 ='#00bfff', A549 ='#009933', EmbryoBrain ='#77a725',  EmbryoHeart ='#f76eb5', EmbryoLiver ='#4d94ff', EmbryoSC='#ff9900', WBlood='#e6e600', Liver='#005ce6', iPSC='#995c00', HEK293T='#cc4400', BLaER1='#ffd9b3'))"
simpleBiotypes_Rpalette="c('lncRNA' = '#66ccff',  'nonExonic' = '#6666ff',  'protein_coding' = '#ff8c66', miRNA = '#808000', 'misc_RNA' = '#99ff99','pseudogene' = '#d98cb3', 'rRNA' = '#d9d9d9', 'ERCC' = '#8dd3c7', 'SIRV' = '#399384')"


CAPDESIGNTOGENOME=config["capDesignToGenome"]
CAPDESIGNTOANNOTGTF=config["capDesignToAnnotGtf"]
GENOMETOCAGEPEAKS=config["genomeToCAGEpeaks"]
GENOMETODHSPEAKS=config["genomeToDHSpeaks"]
GENOMETOPAS=config["genomeToPAS"]
DUMMY_DIR="dummy/"
FQPATH=config["FQPATH"]
FQ_CORR_PATH=FQPATH + "corr/"
TRACK_HUB_DATA_URL=config["TRACK_HUB_BASE_URL"] + "dataFiles/"
if 'genomeToPreviousPhaseTMs' in config:
	GENOMETOPREVIOUS=config["genomeToPreviousPhaseTMs"]
else:
	GENOMETOPREVIOUS=None

print ("Input LR FASTQs are in: " + FQPATH, file=sys.stderr)





ENDSUPPORTcategories=["all", "cagePolyASupported"]
TMSPLICEDSTATUScategories=["all", "spliced", "unspliced"]
READFILTERS=['HiSS', 'noFilt']

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
			#CAPDESIGNTOCAGEPEAKS[mergedCapDString]=CAPDESIGNTOCAGEPEAKS[capD]
			CAPDESIGNTOANNOTGTF[mergedCapDString]=CAPDESIGNTOANNOTGTF[capD]
		CAPDESIGNTOGENOME[mergedCapDString]=genome


GENOMETOCAPDESIGNSplusMERGED=defaultdict(list)
for capD in CAPDESIGNTOGENOME:
	genome=CAPDESIGNTOGENOME[capD]
	GENOMETOCAPDESIGNSplusMERGED[genome].append(capD)

print ("Genomes:", file=sys.stderr)
print (GENOMES, file=sys.stderr)
# no underscores allowed in wildcards, to avoid greedy matching since we use them as separators

SIRVpresent=None
try:
	config["SIRVinfo"]
except KeyError:
	SIRVpresent=False
	print ("Will NOT compare TMs to SIRVs because config['SIRVinfo'] not found.", file=sys.stderr)
else:
	SIRVpresent=True
	print ("Will compare TMs to SIRVs", file=sys.stderr)


if config["DEMULTIPLEX"]:
	print ("Demultiplexing is ON", file=sys.stderr)
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
	print ("Demultiplexing is OFF", file=sys.stderr)
#	DEMULTIPLEX_DIR=FQPATH
	DEMULTIPLEXED_FASTQS=FQ_CORR_PATH
	(TECHNAMES, CAPDESIGNS, SIZEFRACS, BARCODES) = glob_wildcards(FQPATH + "{techname, [^_/]+}_{capDesign}_{sizeFrac}_{barcodes}.fastq.gz")
	BARCODES=set(BARCODES)
	BARCODESUNDETER=list(BARCODES)


CAPDESIGNS=set(CAPDESIGNS)
SIZEFRACS=set(SIZEFRACS)
TECHNAMES=set(TECHNAMES)

SIZEFRACSnoSIZESELECTONLY=["0+"]

wildcard_constraints:
 	capDesign = "[^_/]+",
 	sizeFrac = "[^_/]+",
 	techname = "[^_/]+",
 	corrLevel = "[^_/]+",
 	barcodes = "[^_/]+",
	minReadSupport = '|'.join(config["MINIMUM_TMERGE_READ_SUPPORT"]),
	endSupport  = '|'.join(ENDSUPPORTcategories),
	splicedStatus = '|'.join(TMSPLICEDSTATUScategories),
	filt= '|'.join(READFILTERS),
	sample_name= "[^/]+",

## read sample metadata annotations into a pandas dataframe:
sampleAnnot = pd.read_table(config['SAMPLE_ANNOT'], header=0, sep='\t')

## check that samples are properly annotated:
SAMPLES, = glob_wildcards(FQPATH + "{sample_name, [^/]+}.fastq.gz")
for sample_name in SAMPLES:
	sampleRow=sampleAnnot[sampleAnnot.sample_name == sample_name]
	assert len(sampleRow) < 2, "Duplicate found for sample " + sample_name + " in " + config['SAMPLE_ANNOT']
	assert len(sampleRow) > 0, "Sample " + sample_name + " not found in " + config['SAMPLE_ANNOT']

sampleAnnot.set_index('sample_name', inplace=True)
sampleAnnotDict = sampleAnnot.to_dict('index')





BARCODESpluSMERGED=set(BARCODES)

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
	print ("LR error correction is ON. Correction levels: ", file=sys.stderr)
	print (FINALCORRECTIONLEVELS, file=sys.stderr)

	#make string to interpolate into bash script in "lordecCorrectLr" rule
	graph_kmers_string= "(" + " ".join(graph_kmers) + ")"
	solid_kmer_abundance_threshold="3"
	splitFastqsInto=19
	splitFasta =[]
	for x in range(0, splitFastqsInto, 1):
		splitFasta.append(str(x).zfill(4))
else:
	print ("LR error correction is OFF. Correction levels: ", file=sys.stderr)
	FINALCORRECTIONLEVELS=["0",]
	lastK='NULL'
	print (FINALCORRECTIONLEVELS, file=sys.stderr)





########################
### make list of authorized wildcard combinations
### inspired by https://stackoverflow.com/questions/41185567/how-to-use-expand-in-snakemake-when-some-particular-combinations-of-wildcards-ar
###

AUTHORIZEDCOMBINATIONS = []
AUTHORIZEDCOMBINATIONSMERGE = []

for comb in itertools.product(TECHNAMES,CAPDESIGNS,SIZEFRACS,BARCODES):
	if(os.path.isfile(FQPATH + comb[0] + "_" + comb[1] + "_" + comb[2] + ".fastq.gz") or
	os.path.isfile(FQPATH + comb[0] + "_" + comb[1] + "_" + comb[2]  + "_" + comb[3] + ".fastq.gz")): #allow only combinations corresponding to existing FASTQs
		AUTHORIZEDCOMBINATIONS.append((("techname", comb[0]),("capDesign", comb[1])))
		AUTHORIZEDCOMBINATIONS.append((("techname", comb[0]),("capDesign", comb[1]),("sizeFrac", comb[2])))
		AUTHORIZEDCOMBINATIONS.append((("techname", comb[0]),("capDesign", comb[1]),("sizeFrac", comb[2]),("barcodes", comb[3])))

		for mergedComb in itertools.product([("techname", comb[0]), ("techname", "byTech"), ("techname", "allSeqTechs")], [("capDesign", comb[1]), ("capDesign", "byCapDesign")], [("sizeFrac", comb[2]),]):
			AUTHORIZEDCOMBINATIONSMERGE.append((mergedComb))
		for corrLevel in FINALCORRECTIONLEVELS:
			AUTHORIZEDCOMBINATIONS.append((("techname", comb[0]),("corrLevel", corrLevel),("capDesign", comb[1]),("sizeFrac", comb[2])))
			for mergedComb in itertools.product([("techname", comb[0]), ("techname", "byTech"), ("techname", "allSeqTechs")], [("corrLevel", corrLevel), ("corrLevel", "byCorr")] , [("capDesign", comb[1]), ("capDesign", "byCapDesign")], [("sizeFrac", comb[2]),]):
				AUTHORIZEDCOMBINATIONSMERGE.append((mergedComb))

		if (comb[3]).find(comb[1]) == 0 or config["DEMULTIPLEX"] is False: # check that barcode ID's prefix matches capDesign (i.e. ignore demultiplexed FASTQs with unmatched barcode/capDesign)
			AUTHORIZEDCOMBINATIONS.append((("techname", comb[0]),("capDesign", comb[1]),("barcodes", comb[3])))

			if config["LORDEC_CORRECT"]:
				for split in splitFasta:
					AUTHORIZEDCOMBINATIONS.append((("techname", comb[0]),("capDesign", comb[1]),("sizeFrac", comb[2]),("barcodes",comb[3]), ("split", split)))
			for corrLevel in FINALCORRECTIONLEVELS:
				AUTHORIZEDCOMBINATIONS.append((("techname", comb[0]),("corrLevel", corrLevel),("capDesign", comb[1]),("sizeFrac", comb[2]),("barcodes",comb[3])))
				AUTHORIZEDCOMBINATIONS.append((("techname", comb[0]),("corrLevel", corrLevel),("capDesign", comb[1]),("barcodes",comb[3])))
				AUTHORIZEDCOMBINATIONS.append((("corrLevel", corrLevel),("capDesign", comb[1]),("sizeFrac", comb[2]),("barcodes",comb[3])))
				for mergedComb in itertools.product([("techname", comb[0]), ("techname", "byTech"), ("techname", "allSeqTechs")], [("corrLevel", corrLevel), ("corrLevel", "byCorr")] , [("capDesign", comb[1]), ("capDesign", "byCapDesign")], [("sizeFrac", comb[2]), ], [("barcodes",comb[3]), ("barcodes", "allTissues"), ("barcodes", "byTissue")]):
					AUTHORIZEDCOMBINATIONSMERGE.append((mergedComb))

				for mergedComb in itertools.product([("techname", comb[0]), ("techname", "byTech"), ("techname", "allSeqTechs")], [("corrLevel", corrLevel), ("corrLevel", "byCorr")] , [("capDesign", comb[1]), ("capDesign", "byCapDesign")], [("barcodes",comb[3]), ("barcodes", "allTissues"), ("barcodes", "byTissue")]):
					AUTHORIZEDCOMBINATIONSMERGE.append((mergedComb))

				for minReadSupport in config["MINIMUM_TMERGE_READ_SUPPORT"]:
					AUTHORIZEDCOMBINATIONS.append((("techname", comb[0]),("corrLevel", corrLevel),("capDesign", comb[1]),("sizeFrac", comb[2]),("barcodes",comb[3]), ("minReadSupport", minReadSupport)))
					for mergedComb in itertools.product([("techname", comb[0]), ("techname", "byTech"), ("techname", "allSeqTechs")], [("corrLevel", corrLevel), ("corrLevel", "byCorr")] , [("capDesign", comb[1]), ("capDesign", "byCapDesign")], [("barcodes",comb[3]), ("barcodes", "allTissues"), ("barcodes", "byTissue")], [("minReadSupport", minReadSupport),]):
						AUTHORIZEDCOMBINATIONSMERGE.append((mergedComb))
					for mergedComb in itertools.product([("techname", comb[0]), ("techname", "byTech"), ("techname", "allSeqTechs")], [("corrLevel", corrLevel), ("corrLevel", "byCorr")] , [("capDesign", comb[1]), ("capDesign", "byCapDesign")], [("sizeFrac", comb[2]), ], [("barcodes",comb[3]), ("barcodes", "allTissues"), ("barcodes", "byTissue")], [("minReadSupport", minReadSupport),]):
						AUTHORIZEDCOMBINATIONSMERGE.append((mergedComb))
					for filt in READFILTERS:
						for mergedComb in itertools.product([("techname", comb[0]), ("techname", "byTech"), ("techname", "allSeqTechs")], [("corrLevel", corrLevel), ("corrLevel", "byCorr")] , [("capDesign", comb[1]), ("capDesign", "byCapDesign")], [("barcodes",comb[3]), ("barcodes", "allTissues"), ("barcodes", "byTissue")], [("minReadSupport", minReadSupport),], [("filt", filt),]):
							AUTHORIZEDCOMBINATIONSMERGE.append((mergedComb))
						for mergedComb in itertools.product([("techname", comb[0]), ("techname", "byTech"), ("techname", "allSeqTechs")], [("corrLevel", corrLevel), ("corrLevel", "byCorr")] , [("capDesign", comb[1]), ("capDesign", "byCapDesign")], [("sizeFrac", comb[2]), ], [("barcodes",comb[3]), ("barcodes", "allTissues"), ("barcodes", "byTissue")], [("minReadSupport", minReadSupport),], [("filt", filt),]):
							AUTHORIZEDCOMBINATIONSMERGE.append((mergedComb))
				for endSupport in ENDSUPPORTcategories:
					for mergedComb in itertools.product([("techname", comb[0]), ("techname", "byTech"), ("techname", "allSeqTechs")], [("corrLevel", corrLevel), ("corrLevel", "byCorr")] , [("capDesign", comb[1]), ("capDesign", "byCapDesign")], [("barcodes",comb[3]), ("barcodes", "allTissues"), ("barcodes", "byTissue")], [("endSupport", endSupport),]):
						AUTHORIZEDCOMBINATIONSMERGE.append((mergedComb))
					for mergedComb in itertools.product([("techname", comb[0]), ("techname", "byTech"), ("techname", "allSeqTechs")], [("corrLevel", corrLevel), ("corrLevel", "byCorr")] , [("capDesign", comb[1]), ("capDesign", "byCapDesign")], [("sizeFrac", comb[2]), ], [("barcodes",comb[3]), ("barcodes", "allTissues"), ("barcodes", "byTissue")], [("endSupport", endSupport),]):
						AUTHORIZEDCOMBINATIONSMERGE.append((mergedComb))
					for minReadSupport in config["MINIMUM_TMERGE_READ_SUPPORT"]:
						AUTHORIZEDCOMBINATIONS.append((("techname", comb[0]),("corrLevel", corrLevel),("capDesign", comb[1]),("sizeFrac", comb[2]),("barcodes",comb[3]), ("minReadSupport", minReadSupport), ("endSupport", endSupport)))
						for mergedComb in itertools.product([("techname", comb[0]), ("techname", "byTech"), ("techname", "allSeqTechs")], [("corrLevel", corrLevel), ("corrLevel", "byCorr")] , [("capDesign", comb[1]), ("capDesign", "byCapDesign")], [("barcodes",comb[3]), ("barcodes", "allTissues"), ("barcodes", "byTissue")], [("endSupport", endSupport),],  [("minReadSupport", minReadSupport),]):
							AUTHORIZEDCOMBINATIONSMERGE.append((mergedComb))
						for mergedComb in itertools.product([("techname", comb[0]), ("techname", "byTech"), ("techname", "allSeqTechs")], [("corrLevel", corrLevel), ("corrLevel", "byCorr")] , [("capDesign", comb[1]), ("capDesign", "byCapDesign")], [("sizeFrac", comb[2]), ], [("barcodes",comb[3]), ("barcodes", "allTissues"), ("barcodes", "byTissue")], [("endSupport", endSupport),],  [("minReadSupport", minReadSupport),]):
							AUTHORIZEDCOMBINATIONSMERGE.append((mergedComb))
						for splicedStatus in TMSPLICEDSTATUScategories:
							for mergedComb in itertools.product([("techname", comb[0]), ("techname", "byTech"), ("techname", "allSeqTechs")], [("corrLevel", corrLevel), ("corrLevel", "byCorr")] , [("capDesign", comb[1]), ("capDesign", "byCapDesign")], [("sizeFrac", comb[2]), ], [("barcodes",comb[3]), ("barcodes", "allTissues"), ("barcodes", "byTissue")], [("minReadSupport", minReadSupport),], [("splicedStatus", splicedStatus),]):
								AUTHORIZEDCOMBINATIONSMERGE.append((mergedComb))
							for mergedComb in itertools.product([("techname", comb[0]), ("techname", "byTech"), ("techname", "allSeqTechs")], [("corrLevel", corrLevel), ("corrLevel", "byCorr")] , [("capDesign", comb[1]), ("capDesign", "byCapDesign")], [("sizeFrac", comb[2]), ], [("barcodes",comb[3]), ("barcodes", "allTissues"), ("barcodes", "byTissue")], [("endSupport", endSupport),], [("minReadSupport", minReadSupport),], [("splicedStatus", splicedStatus),]):
								AUTHORIZEDCOMBINATIONSMERGE.append((mergedComb))


			for splice in SPLICE_SITE_TYPES:
				for corrLevel in FINALCORRECTIONLEVELS:
					AUTHORIZEDCOMBINATIONS.append((("techname", comb[0]),("corrLevel", corrLevel),("capDesign", comb[1]),("sizeFrac", comb[2]),("barcodes",comb[3]),("spliceType",splice)))
					for mergedComb in itertools.product([("techname", comb[0]), ("techname", "byTech"), ("techname", "allSeqTechs")], [("corrLevel", corrLevel), ("corrLevel", "byCorr")] , [("capDesign", comb[1]), ("capDesign", "byCapDesign")], [("sizeFrac", comb[2]), ], [("barcodes",comb[3]), ("barcodes", "allTissues"), ("barcodes", "byTissue")], [("spliceType",splice)]):
						AUTHORIZEDCOMBINATIONSMERGE.append((mergedComb))
					for minReadSupport in config["MINIMUM_TMERGE_READ_SUPPORT"]:
						AUTHORIZEDCOMBINATIONS.append((("techname", comb[0]),("corrLevel", corrLevel),("capDesign", comb[1]),("sizeFrac", comb[2]),("barcodes",comb[3]),("spliceType",splice),("minReadSupport", minReadSupport)))
						for mergedComb in itertools.product([("techname", comb[0]), ("techname", "byTech"), ("techname", "allSeqTechs")], [("corrLevel", corrLevel), ("corrLevel", "byCorr")] , [("capDesign", comb[1]), ("capDesign", "byCapDesign")], [("sizeFrac", comb[2]), ], [("barcodes",comb[3]), ("barcodes", "allTissues"), ("barcodes", "byTissue")], [("spliceType",splice)], [("minReadSupport",minReadSupport)]):
							AUTHORIZEDCOMBINATIONSMERGE.append((mergedComb))





AUTHORIZEDCOMBINATIONS=set(AUTHORIZEDCOMBINATIONS)
AUTHORIZEDCOMBINATIONSMERGE=set(AUTHORIZEDCOMBINATIONSMERGE)


PLOTSbyTISSUE=set(BARCODESpluSMERGED)
PLOTSbyTISSUE.add("byTissue")
PLOTSbyTISSUEnoALL=set(BARCODES)
PLOTSbyTISSUEnoALL.add("byTissue")


PLOTSbySEQTECH=set(TECHNAMES)
PLOTSbySEQTECH.add("byTech")
PLOTSbySEQTECHnoALLTECHS=set(TECHNAMES)
PLOTSbySEQTECHnoALLTECHS.add("byTech")

PLOTSbyCORRLEVEL=set(FINALCORRECTIONLEVELS)
if config["LORDEC_CORRECT"]:
	PLOTSbyCORRLEVEL.add("byCorr")
PLOTSbyCAPDESIGN=set(CAPDESIGNS)
PLOTSbyCAPDESIGN.add("byCapDesign")
PLOTSbyCAPDESIGNplusMERGED=set(PLOTSbyCAPDESIGN)
for capD in MERGEDCAPDESIGNS:
	PLOTSbyCAPDESIGNplusMERGED.add(capD)

#def sampleNameToSeqPlatform(sample_name):
#	seqPlatform=sampleAnnot.loc[sampleAnnot.sample_name == sample_name, 'seqPlatform'].values[0]
#	return seqPlatform


def returnPlotFilenames(basename):
	plotsList=[]
	for comb in itertools.product(["legendOnly"], config["PLOTFORMATS"]):
		plotsList.append(basename + "." + comb[0] + "." + comb[1])
	for comb in itertools.product(["xy", "yx"], ["wLegend", "woLegend"], config["PLOTFORMATS"]):
		plotsList.append(basename + "." + comb[0] + "." + comb[1] + "." + comb[2])
	return plotsList



def filtered_product(*args):
	found=False
	for wc_comb in itertools.product(*args):
		#print(wc_comb)
		if wc_comb in AUTHORIZEDCOMBINATIONS:
			found=True
			yield(wc_comb)
		else:
			if len(wc_comb) > 4:
				if wc_comb[2][0] in ('capDesign') and wc_comb[2][1] in MERGEDCAPDESIGNS and wc_comb[4] in (('barcodes', 'allTissues'),):
					found=True
					yield(wc_comb)
	if not found:
		pprint(" Error in function filtered_product. Args were:", file=sys.stderr)
		pprint((args), file=sys.stderr)
		quit(" Error. Could not yield any input file.", file=sys.stderr)

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
				if wc_comb[4] in (('barcodes', 'allTissues'),):
					found=True
					#print("YIELD 1")
					yield(wc_comb)

			else:
				found=True
				#print("YIELD 2")
				yield(wc_comb)
		else:
			if wc_comb[2][0] in ('capDesign') and wc_comb[2][1] in MERGEDCAPDESIGNS and wc_comb[4] in (('barcodes', 'allTissues'),):
				found=True
				#print("YIELD 3")
				yield(wc_comb)
	if not found:
		pprint(" Error in function filtered_product_merge. Args were:", file=sys.stderr)
		pprint((args), file=sys.stderr)
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
				if wc_comb[4] in (('barcodes', 'allTissues'),):
					found=True
					yield(wc_comb)

			else:
				found=True
				yield(wc_comb)
		else:
			if wc_comb[2][0] in ('capDesign') and wc_comb[2][1] in MERGEDCAPDESIGNS and wc_comb[4] in (('barcodes', 'allTissues'),):
				#print ("OKcapDesign_product_merge")
				for capD in MERGEDCAPDESIGNS[wc_comb[2][1]]:
					wc_comb2=list(wc_comb)
					wc_comb2[2]=('capDesign', capD)
					wc_comb2=tuple(wc_comb2)
					#print (wc_comb2)
					found=True
					yield(wc_comb2)
	if not found:
		pprint(" Error in function filtered_capDesign_product_merge. Args were:", file=sys.stderr)
		pprint((args), file=sys.stderr)
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
					if wc_comb[3] in (('sizeFrac', '0+'),) and wc_comb[4] in (('barcodes', 'allTissues'),):
						yieldWC=True
					else:
					#print ("NO AUTH 2")
						continue
				else:
		#		print ("AUTH 1")
					yieldWC=True
		else:
			if wc_comb[2][0] in ('capDesign') and wc_comb[2][1] in MERGEDCAPDESIGNS and wc_comb[4] in (('barcodes', 'allTissues'),):
				yieldWC=True
		if yieldWC == True:
			#print("AUTH")
			yield(wc_comb)
			found=True
	if not found:
		pprint(" Error in function filtered_merge_figures. Args were:", file=sys.stderr)
		pprint((args), file=sys.stderr)
		quit(" Error. Could not yield any input file.")


def merge_figures_params(c,bf,bt,cl, tn, splicing_status=None):
	clBool=None
	if len(FINALCORRECTIONLEVELS)>1 and cl == FINALCORRECTIONLEVELS[1]:
		clBool="Yes"
	elif cl == FINALCORRECTIONLEVELS[0]:
		clBool="No"
	else:
		clBool=cl
	seqTechFilterString=''
	removeFacetLabels=''
	if tn not in ('byTech') and c not in ('byCapDesign') and bt not in ('byTissue'):
		removeFacetLabels = """
plotFacetXy <- parse(text =paste(plotFacetXy, \\" + theme(strip.text.x = element_blank(), strip.text.y = element_blank())\\"))
plotFacetYx <- parse(text=paste(plotFacetYx, \\" + theme(strip.text.x = element_blank(), strip.text.y = element_blank())\\"))
"""		
	splicingStatusFilterString=''
	if splicing_status != None:
		if splicing_status == 'spliced':
			splicingStatusFilterString="dat <- subset(dat, spliced==1)\n"
		elif splicing_status == 'unspliced':
			splicingStatusFilterString="dat <- subset(dat, spliced==0)\n"
		elif splicing_status == 'all':
			splicingStatusFilterString=''
		else:
			quit("Error. Invalid value for splicing_status in function merge_figures_params.")
	if tn not in ('byTech'):
		seqTechFilterString="dat <- subset(dat, seqTech=='" + tn + "')\n"
	else:
		seqTechFilterString="dat <- subset(dat, seqTech!='allSeqTechs')\n"
	substSeqTechString="""
dat\$seqTech <- gsub(':$', '', dat\$seqTech)
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

horizCats <- (length(unique(dat\$correctionLevel))/2) * length(unique(dat\$capDesign)) * length(unique(dat\$tissue))
vertCats <- length(unique(dat\$seqTech))

wXyPlot = (horizCats * 0.9) +1.7
hXyPlot = (vertCats * 0.6) + 1.7

geom_textSize=1.4 # * (max(vertCats, horizCats))
themeSize = (14/5) * geom_textSize
# https://stackoverflow.com/questions/25061822/ggplot-geom-text-font-size-control/25062509
lineSize=geom_textSize/8
minorLineSize=lineSize/2

"""
	facetPlotSetup= f"""
plotFacetXy <- parse(text =paste(plotBase, \\"facet_grid( seqTech ~ capDesign + tissue, scales='free_y')\\"))
plotFacetYx <- parse(text=paste(plotBase, \\"facet_grid( capDesign + tissue ~ seqTech, scales='free_y')\\"))
{removeFacetLabels}
pXy <- eval(plotFacetXy)
pYx <- eval(plotFacetYx)

legend <- get_legend(pXy)

pXyNoLegend <- pXy + theme(legend.position='none')
pYxNoLegend <- pYx + theme(legend.position='none')

legendOnly <- grid.arrange(legend)
pXyGrob <- as.grob(pXy)
pYxGrob <- as.grob(pYx)
pXyNoLegendGrob <- as.grob(pXyNoLegend)
pYxNoLegendGrob <- as.grob(pYxNoLegend)


hLegendOnly <- convertUnit(sum(legend\$heights), 'in', valueOnly=TRUE)
wLegendOnly <- convertUnit(sum(legend\$widths), 'in', valueOnly=TRUE)


hYxPlot <- wXyPlot
wYxPlot <- hXyPlot 

hXyNoLegendPlot<- hXyPlot 
wXyNoLegendPlot<- wXyPlot - wLegendOnly

hYxNoLegendPlot<- hYxPlot
wYxNoLegendPlot<- wYxPlot - wLegendOnly
"""	

	if c not in ('byCapDesign'):
		capDesignFilterString+="dat <- subset(dat, capDesign=='" + c + "')\n"
 	if clBool not in ('byCorr'):
 		corrLevelFilterString+="dat <- subset(dat, correctionLevel=='" + clBool + "')\n"
 		xlabString=''
 		hideXaxisLabels="theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +"
 		hideYaxisLabels="theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) +"
# #	if bf in ('allFracs'):
#		sizeFracFilterString+="dat <- subset(dat, sizeFrac=='" + bf + "')\n"
#	else:
	sizeFracFilterString+="dat <- subset(dat, sizeFrac=='" + bf + "')\n"
	if bt in ('allTissues'):
		tissueFilterString+="dat <- subset(dat, tissue=='" + bt + "')\n"
	else:
		if bt in ('byTissue'):
			tissueFilterString+="dat <- subset(dat, tissue!='allTissues')\n"
		else:
			tissueFilterString+="dat <- subset(dat, tissue=='" + bt + "')\n"
	return(capDesignFilterString, corrLevelFilterString, sizeFracFilterString, tissueFilterString, substSeqTechString, substTissueString,  xlabString, hideXaxisLabels , graphDimensions, hideYaxisLabels, seqTechFilterString, splicingStatusFilterString, facetPlotSetup)


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
	return(returnSubGroup1String + " " + " ".join(sorted(returnSubGroup1StringData)) + "\n" + returnSubGroup2String + " " + " ".join(sorted(returnSubGroup2StringData)) + "\n" + returnSubGroup3String + " " + " ".join(sorted(returnSubGroup3StringData))+ "\n" + returnSubGroup4String + " " + " ".join(sorted(returnSubGroup4StringData)))

def trackChecked(t,c,cd,s,b,m):
	checked=''
	if t.find('pacBio') == 0 and (len(FINALCORRECTIONLEVELS)>1 and c == FINALCORRECTIONLEVELS[1]):
		checked='off'
	elif t.find('ont') == 0 and (len(FINALCORRECTIONLEVELS)>1 and c == FINALCORRECTIONLEVELS[1]):
		checked='off'
	else:
		if s == '0+' and b == 'allTissues':
			checked='on'
		else:
			checked='off'
	if m == '1':
		checked='off'
	return(checked)



include: "lrCorrection.smk"
if config["DEMULTIPLEX"]:
	include: "demultiplex.smk"

include: "fastqStats.smk"
include: "lrMapping.smk"
include: "srMapping.smk"
include: "polyAmapping.smk"
include: "introns.smk"
include: "processReadMappings.smk"
include: "tmClassification.smk"
include: "tmEndSupport.smk"
include: "trackHub.smk"
include: "lociEndSupport.smk"





#pseudo-rule specifying the target files we ultimately want.
rule all:
	input:
		expand(returnPlotFilenames(config["PLOTSDIR"] + "HiSS.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.stats"), filtered_merge_figures, techname=PLOTSbySEQTECHnoALLTECHS, corrLevel=PLOTSbyCORRLEVEL,  capDesign=PLOTSbyCAPDESIGN, sizeFrac=SIZEFRACS, barcodes=PLOTSbyTISSUE),
		expand("mappings/nonAnchoredMergeReads/colored/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.splicing_status:all.endSupport:{endSupport}.bed", filtered_product_merge, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNSplusMERGED, sizeFrac=SIZEFRACS, barcodes=BARCODESpluSMERGED, endSupport=ENDSUPPORTcategories, minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"], splicedStatus=TMSPLICEDSTATUScategories),
		config["STATSDATADIR"] + "all.highConfSplicedReads.stats.tsv",

		expand(FQ_CORR_PATH + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.fastq.gz", filtered_product, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS) if config["DEMULTIPLEX"] else expand(FQ_CORR_PATH + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.fastq.gz", filtered_product, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES),

		expand(returnPlotFilenames(config["PLOTSDIR"] + "readLength.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.readLength.stats"), filtered_merge_figures, techname=PLOTSbySEQTECHnoALLTECHS, corrLevel=PLOTSbyCORRLEVEL, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS) if config["DEMULTIPLEX"] else expand(returnPlotFilenames(config["PLOTSDIR"] + "readLength.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.readLength.stats"), filtered_merge_figures, techname=PLOTSbySEQTECHnoALLTECHS, corrLevel=PLOTSbyCORRLEVEL, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES), # facetted histograms of read length
		expand(FQ_CORR_PATH + "qc/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.dupl.txt", filtered_product,techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS) if config["DEMULTIPLEX"] else expand(FQ_CORR_PATH + "qc/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.{barcodes}.dupl.txt", filtered_product,techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES),

		expand(config["PLOTSDIR"] + "{techname}Corr{corrLevel}.fastq.UP.stats.{ext}", techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, ext=config["PLOTFORMATS"]) if config["DEMULTIPLEX"] else expand( DUMMY_DIR + "dummy{number}.txt", number='1'), # UP reads plots
		expand(config["PLOTSDIR"] + "{techname}Corr{corrLevel}.fastq.BC.stats.{ext}", techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, ext=config["PLOTFORMATS"]) if config["DEMULTIPLEX"]  else expand( DUMMY_DIR + "dummy{number}.txt", number='2') , # barcode reads plots
		expand(config["PLOTSDIR"] + "{techname}Corr{corrLevel}.fastq.foreignBC.stats.{ext}", techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, ext=config["PLOTFORMATS"]) if config["DEMULTIPLEX"]  else expand( DUMMY_DIR + "dummy{number}.txt", number='3') , #foreign barcode reads plots
		expand(DEMULTIPLEX_DIR + "demultiplexFastqs/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.fastq.gz", filtered_product, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES) if config["DEMULTIPLEX"]  else expand( DUMMY_DIR + "dummy{number}.txt", number='4'),
		expand("mappings/readMapping/qc/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.bam.dupl.txt",filtered_product_merge, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODESpluSMERGED),
		expand(config["PLOTSDIR"] + "{techname}Corr{corrLevel}.ambiguousBarcodes.reads.stats.{ext}", techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, ext=config["PLOTFORMATS"]) if config["DEMULTIPLEX"]  else expand( DUMMY_DIR + "dummy{number}.txt", number='6'), # ambiguous barcodes plots
		expand(config["PLOTSDIR"] + "{techname}Corr{corrLevel}_{capDesign}.adapters.location.stats.{ext}",techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, ext=config["PLOTFORMATS"]) if config["DEMULTIPLEX"]  else expand( DUMMY_DIR + "dummy{number}.txt", number='7'), #location of adapters over reads
		expand(config["PLOTSDIR"] + "{techname}Corr{corrLevel}.chimeric.reads.stats.{ext}", techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, ext=config["PLOTFORMATS"]) if config["DEMULTIPLEX"]  else expand( DUMMY_DIR + "dummy{number}.txt", number='8'), # stats on chimeric reads
		expand(config["PLOTSDIR"] + "{techname}Corr{corrLevel}.finalDemul.reads.stats.{ext}", techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, ext=config["PLOTFORMATS"]) if config["DEMULTIPLEX"]  else expand( DUMMY_DIR + "dummy{number}.txt", number='9'), #final demultiplexing stats
		expand(DEMULTIPLEX_DIR + "qc/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.demul.QC1.txt",filtered_product, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS) if config["DEMULTIPLEX"]  else expand( DUMMY_DIR + "dummy{number}.txt", number='10'), # QC on demultiplexing (checks that there is only one barcode assigned per read
		expand(DEMULTIPLEX_DIR + "qc/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.demul.QC2.txt", filtered_product, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS) if config["DEMULTIPLEX"]  else expand( DUMMY_DIR + "dummy{number}.txt", number='11'), # QC on demultiplexing (checks that there is only one barcode assigned per read
		expand(config["PLOTSDIR"] + "{techname}Corr{corrLevel}.demultiplexing.perSample.stats.{ext}", techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, ext=config["PLOTFORMATS"])  if config["DEMULTIPLEX"]  else expand( DUMMY_DIR + "dummy{number}.txt", number='12'),
		expand(returnPlotFilenames(config["PLOTSDIR"] + "lrMapping.basic.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.lrMapping.basic.stats"), filtered_merge_figures, techname=PLOTSbySEQTECHnoALLTECHS, corrLevel=PLOTSbyCORRLEVEL,  capDesign=PLOTSbyCAPDESIGN, sizeFrac=SIZEFRACS, barcodes=PLOTSbyTISSUE),
		expand(returnPlotFilenames(config["PLOTSDIR"] + "intraPriming.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.intraPriming.stats"), filtered_merge_figures, techname=PLOTSbySEQTECHnoALLTECHS, corrLevel=PLOTSbyCORRLEVEL,  capDesign=PLOTSbyCAPDESIGN, sizeFrac=SIZEFRACS, barcodes=PLOTSbyTISSUE),
		expand(returnPlotFilenames(config["PLOTSDIR"] + "lrMapping.spikeIns.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.lrMapping.spikeIns.stats"), filtered_merge_figures, techname=PLOTSbySEQTECHnoALLTECHS, corrLevel=PLOTSbyCORRLEVEL,  capDesign=PLOTSbyCAPDESIGN, sizeFrac=SIZEFRACS, barcodes=PLOTSbyTISSUE),

		expand(config["PLOTSDIR"] + "hiSeq.mapping.stats/all.hiSeq.mapping.stats.{ext}", ext=config["PLOTFORMATS"]),
		expand(config["PLOTSDIR"] + "hiSeq.SJs.stats/all.hiSeq.SJs.stats.{ext}", ext=config["PLOTFORMATS"]),
		expand ("mappings/clusterPolyAsites/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.polyAsites.clusters.bed", filtered_product_merge, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODESpluSMERGED),
		
		
		
		
		expand(returnPlotFilenames(config["PLOTSDIR"] + "merged.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.merged.stats"), filtered_merge_figures, techname=PLOTSbySEQTECHnoALLTECHS, corrLevel=PLOTSbyCORRLEVEL,  capDesign=PLOTSbyCAPDESIGN, sizeFrac=SIZEFRACS, barcodes=PLOTSbyTISSUE, ext=config["PLOTFORMATS"], minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"]),
		expand(returnPlotFilenames(config["PLOTSDIR"] + "cagePolyASupport.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.splicing_status:{splicedStatus}.cagePolyASupport.stats"), filtered_merge_figures, techname=PLOTSbySEQTECH, corrLevel=PLOTSbyCORRLEVEL,  capDesign=PLOTSbyCAPDESIGNplusMERGED, sizeFrac=SIZEFRACS, barcodes=PLOTSbyTISSUE, minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"], splicedStatus=TMSPLICEDSTATUScategories),
		expand(returnPlotFilenames(config["PLOTSDIR"] + "tmerge.vs.gencode.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.splicing_status:{splicedStatus}.endSupport:{endSupport}.vs.gencode.stats"), filtered_merge_figures, techname=PLOTSbySEQTECH, corrLevel=PLOTSbyCORRLEVEL,  capDesign=PLOTSbyCAPDESIGNplusMERGED, sizeFrac=SIZEFRACS, barcodes=PLOTSbyTISSUE, endSupport=ENDSUPPORTcategories, minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"], splicedStatus=TMSPLICEDSTATUScategories),
		expand(returnPlotFilenames(config["PLOTSDIR"] + "targetCoverage.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.targetCoverage.stats"), filtered_merge_figures, techname=PLOTSbySEQTECH, corrLevel=PLOTSbyCORRLEVEL,  capDesign=PLOTSbyCAPDESIGN, sizeFrac=SIZEFRACS, barcodes=PLOTSbyTISSUE, minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"]) if config["CAPTURE"]  else expand( DUMMY_DIR + "dummy{number}.txt", number='15'),
		#expand(returnPlotFilenames(config["PLOTSDIR"] + "spliceSites.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.spliceSites.stats"), filtered_merge_figures, techname=PLOTSbySEQTECHnoALLTECHS, corrLevel=PLOTSbyCORRLEVEL,  capDesign=PLOTSbyCAPDESIGN, sizeFrac=SIZEFRACS, barcodes=PLOTSbyTISSUE, minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"]),
		#expand(returnPlotFilenames(config["PLOTSDIR"] + "matureRNALength.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.splicing_status:{splicedStatus}.matureRNALength.box.stats"), filtered_merge_figures, techname=PLOTSbySEQTECHnoALLTECHS, corrLevel=PLOTSbyCORRLEVEL,  capDesign=PLOTSbyCAPDESIGN, sizeFrac=SIZEFRACS, barcodes=PLOTSbyTISSUEnoALL, minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"], splicedStatus=TMSPLICEDSTATUScategories),
		expand(returnPlotFilenames(config["PLOTSDIR"] + "matureRNALength.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.splicing_status:{splicedStatus}.matureRNALength.hist.stats"), filtered_merge_figures, techname=PLOTSbySEQTECHnoALLTECHS, corrLevel=PLOTSbyCORRLEVEL,  capDesign=PLOTSbyCAPDESIGN, sizeFrac=SIZEFRACS, barcodes=PLOTSbyTISSUEnoALL, minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"], splicedStatus=TMSPLICEDSTATUScategories),

		expand(returnPlotFilenames(config["PLOTSDIR"] + "polyAreads.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.polyAreads.stats"), filtered_merge_figures, techname=PLOTSbySEQTECHnoALLTECHS, corrLevel=PLOTSbyCORRLEVEL,  capDesign=PLOTSbyCAPDESIGN, sizeFrac=SIZEFRACS, barcodes=PLOTSbyTISSUE),
		expand(returnPlotFilenames(config["PLOTSDIR"] + "tmerge.novelLoci.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.endSupport:{endSupport}.novelLoci.stats"), filtered_merge_figures, techname=PLOTSbySEQTECH, corrLevel=PLOTSbyCORRLEVEL,  capDesign=PLOTSbyCAPDESIGNplusMERGED, sizeFrac=SIZEFRACS, barcodes=PLOTSbyTISSUE, endSupport=ENDSUPPORTcategories, minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"]),
		expand(returnPlotFilenames(config["PLOTSDIR"] + "tmerge.vs.Gencode.SJs.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.vs.Gencode.SJs.stats"), filtered_merge_figures, techname=PLOTSbySEQTECH, corrLevel=PLOTSbyCORRLEVEL,  capDesign=PLOTSbyCAPDESIGN, sizeFrac=SIZEFRACS, barcodes=PLOTSbyTISSUE, minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"]),
		expand(returnPlotFilenames(config["PLOTSDIR"] + "tmerge.vs.SIRVs.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.{filt}.tmerge.min{minReadSupport}reads.vs.SIRVs.stats"), filtered_merge_figures, techname=PLOTSbySEQTECH, corrLevel=PLOTSbyCORRLEVEL,  capDesign=PLOTSbyCAPDESIGN, sizeFrac=SIZEFRACS, barcodes=PLOTSbyTISSUE, minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"], filt=READFILTERS) if SIRVpresent  else expand( DUMMY_DIR + "dummy{number}.txt", number='13'),
		expand(returnPlotFilenames(config["PLOTSDIR"] + "tmerge.vs.gencode.SnPr.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.splicing_status:{splicedStatus}.endSupport:{endSupport}.vs.gencode.SnPr.stats"), filtered_merge_figures, techname=PLOTSbySEQTECH, corrLevel=PLOTSbyCORRLEVEL,  capDesign=PLOTSbyCAPDESIGN, sizeFrac=SIZEFRACS, barcodes=PLOTSbyTISSUE, endSupport=ENDSUPPORTcategories, minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"], splicedStatus=TMSPLICEDSTATUScategories),
		expand(returnPlotFilenames(config["PLOTSDIR"] + "tmerge.vs.gencode.length.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.splicing_status:bySplicingStatus.endSupport:{endSupport}.vs.gencode.length.stats"), filtered_merge_figures, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS,  capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES, endSupport=ENDSUPPORTcategories, minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"]),
		expand(returnPlotFilenames(config["PLOTSDIR"] + "tmerge.vs.gencode.length.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.splicing_status:all.endSupport:{endSupport}.vs.gencode.length.stats"), filtered_merge_figures, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS,  capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES, endSupport=ENDSUPPORTcategories, minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"]),

		expand(returnPlotFilenames(config["PLOTSDIR"] + "tmerge.vs.SIRVs.detection.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.vs.SIRVs.detection.stats"), filtered_merge_figures, techname=PLOTSbySEQTECH, corrLevel=PLOTSbyCORRLEVEL,  capDesign=PLOTSbyCAPDESIGN, sizeFrac=SIZEFRACS, barcodes=PLOTSbyTISSUE, minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"]) if SIRVpresent  else expand( DUMMY_DIR + "dummy{number}.txt", number='14'),
		expand(returnPlotFilenames(config["PLOTSDIR"] + "FLloci.gencodeOnly.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.FLloci.gencodeOnly.stats"), filtered_merge_figures, techname='allSeqTechs', corrLevel=FINALCORRECTIONLEVELS,  capDesign=PLOTSbyCAPDESIGNplusMERGED, sizeFrac=SIZEFRACS, barcodes='allTissues', minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"]) if GENOMETOPREVIOUS is not None else expand( DUMMY_DIR + "dummy{number}.txt", number='15'),
		expand(returnPlotFilenames(config["PLOTSDIR"] + "FLloci.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.FLloci.stats"), filtered_merge_figures, techname='allSeqTechs', corrLevel=FINALCORRECTIONLEVELS,  capDesign=PLOTSbyCAPDESIGNplusMERGED, sizeFrac=SIZEFRACS, barcodes='allTissues', minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"]) if GENOMETOPREVIOUS is not None else expand( DUMMY_DIR + "dummy{number}.txt", number='16'),
		expand(config["PLOTSDIR"] + "sampleComparison.stats/{capDesign}/{capDesign}_min{minReadSupport}reads.splicing_status:{splicedStatus}.endSupport:{endSupport}.heatmap.sampleComparison.simpson.png", capDesign=CAPDESIGNS, minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"], splicedStatus=TMSPLICEDSTATUScategories, endSupport=ENDSUPPORTcategories),
		expand(config["PLOTSDIR"] + "gencode.detected.length.stats/{capDesign}.tmerge.min{minReadSupport}reads.splicing_status:{splicedStatus}.endSupport:{endSupport}.gencode.detected.length.stats.{ext}", capDesign=CAPDESIGNS, minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"], splicedStatus=TMSPLICEDSTATUScategories, endSupport=ENDSUPPORTcategories, ext=config["PLOTFORMATS"]),
		expand(config["PLOTSDIR"] + "dhsVsCage5primeComparison.venn.stats/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.dhsVsCage5primeComparison.venn.stats.pdf", filtered_product, techname=TECHNAMES, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS,sizeFrac=SIZEFRACS, barcodes=BARCODES, minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"]),
		expand(returnPlotFilenames(config["PLOTSDIR"] + "geneReadCoverage.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.geneReadCoverage.stats"), filtered_merge_figures, techname=PLOTSbySEQTECHnoALLTECHS, corrLevel=PLOTSbyCORRLEVEL, capDesign=PLOTSbyCAPDESIGN, sizeFrac=SIZEFRACS, barcodes=PLOTSbyTISSUEnoALL), # facetted histograms of read length
		expand(config["PLOTSDIR"] + "readProfile/byTechCorr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.readProfile.density.png", filtered_product, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES),
		expand(config["PLOTSDIR"] + "readProfile/byTechCorr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.readProfile.heatmap.png", filtered_product, corrLevel=FINALCORRECTIONLEVELS, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES),
		expand(returnPlotFilenames(config["PLOTSDIR"] + "sequencingError.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.sequencingError.allErrors.stats"), filtered_merge_figures, techname=PLOTSbySEQTECHnoALLTECHS, corrLevel=PLOTSbyCORRLEVEL,  capDesign=PLOTSbyCAPDESIGN, sizeFrac=SIZEFRACS, barcodes=PLOTSbyTISSUE),
		expand(returnPlotFilenames(config["PLOTSDIR"] + "sequencingError.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.sequencingError.deletionsOnly.stats"), filtered_merge_figures, techname=PLOTSbySEQTECHnoALLTECHS, corrLevel=PLOTSbyCORRLEVEL,  capDesign=PLOTSbyCAPDESIGN, sizeFrac=SIZEFRACS, barcodes=PLOTSbyTISSUE),
		expand(returnPlotFilenames(config["PLOTSDIR"] + "sqantiRtsJunctions.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.sqantiRtsJunctions.stats"), filtered_merge_figures, techname=PLOTSbySEQTECHnoALLTECHS, corrLevel=PLOTSbyCORRLEVEL,  capDesign=PLOTSbyCAPDESIGN, sizeFrac=SIZEFRACS, barcodes=PLOTSbyTISSUE),
		expand(returnPlotFilenames(config["PLOTSDIR"] + "sqantiCanJunctions.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.sqantiCanJunctions.stats"), filtered_merge_figures, techname=PLOTSbySEQTECHnoALLTECHS, corrLevel=PLOTSbyCORRLEVEL,  capDesign=PLOTSbyCAPDESIGN, sizeFrac=SIZEFRACS, barcodes=PLOTSbyTISSUE),
		expand(returnPlotFilenames(config["PLOTSDIR"] + "gencode.geneDetection.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.splicing_status:{splicedStatus}.endSupport:{endSupport}.gencode.geneDetection.stats"), filtered_merge_figures, techname=PLOTSbySEQTECH, corrLevel=PLOTSbyCORRLEVEL,  capDesign=PLOTSbyCAPDESIGN, sizeFrac=SIZEFRACS, barcodes=PLOTSbyTISSUE, endSupport=ENDSUPPORTcategories, minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"], splicedStatus=TMSPLICEDSTATUScategories),
		#expand(config["STATSDATADIR"] + "all.{capDesign}.mappedReadlength.summary.tsv", capDesign=CAPDESIGNS),
		#expand(config["STATSDATADIR"] + "all.{capDesign}.readlength.summary.tsv", capDesign=CAPDESIGNS),
		expand(returnPlotFilenames(config["PLOTSDIR"] + "readToBiotypeBreakdown.stats/{techname}/Corr{corrLevel}/{capDesign}/{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}_{barcodes}.readToBiotypeBreakdown.stats"), filtered_merge_figures, techname=PLOTSbySEQTECHnoALLTECHS, corrLevel=PLOTSbyCORRLEVEL,  capDesign=PLOTSbyCAPDESIGN, sizeFrac=SIZEFRACS, barcodes=PLOTSbyTISSUE),

## temporary/intermediate
		config["STATSDATADIR"] + "all.fastq.timestamps.tsv",
		config["STATSDATADIR"] + "all.readlength.summary.tsv",
		config["STATSDATADIR"] + "all.mappedReadlength.summary.tsv",
		expand("annotations/{capDesign}.partition.gff", capDesign=CAPDESIGNS),
## temporary



		#################
		### Track hub ###
		#################
		config["TRACK_HUB_DIR"] + "hub.txt",
		config["TRACK_HUB_DIR"] + "genomes.txt",
		expand(config["TRACK_HUB_DIR"] + "{genome}/trackDb.txt", genome=GENOMES),


##
# "Dummy" rule to skip undesired targets depending on config
##
rule dummy:
	input: "/dev/null"
	output: DUMMY_DIR + "dummy{number}.txt"
	shell:
		'''
touch {output}
		'''


rule sortIndexGenome:
	input: config["GENOMESDIR"] +"{genome}.fa"
	output: 
		sorted=config["GENOMESDIR"] +"{genome}.sorted.fa",
		bioperlindex=config["GENOMESDIR"] +"{genome}.sorted.fa.index"
	shell:
		'''
uuid=$(uuidgen)
#check for duplicate sequences:
count=$(cat {input} | fgrep ">" | sort|uniq -d | wc -l)
if [ $count -gt 0 ]; then echo "$count duplicate sequence IDs found"; exit 1; fi

FastaToTbl {input} | sort -T {config[TMPDIR]} -k1,1 | TblToFasta > {config[TMPDIR]}/$uuid 
mv {config[TMPDIR]}/$uuid {output.sorted}
perl -e 'use Bio::DB::Fasta; my $chrdb = Bio::DB::Fasta->new("{output.sorted}");'

		'''

rule makeGenomeFile:
	input: config["GENOMESDIR"] +"{genome}.sorted.fa"
	output: config["GENOMESDIR"] +"{genome}.sorted.genome"
	shell:
		'''
 uuid=$(uuidgen)
FastaToTbl {input} | awk '{{print $1"\\t"length($2)}}' | sort -k1,1 > {config[TMPDIR]}/$uuid
mv {config[TMPDIR]}/$uuid {output}
		'''

rule makeGencodePartition:
	input:
		gtf=lambda wildcards: CAPDESIGNTOANNOTGTF[wildcards.capDesign],
		genome = lambda wildcards: config["GENOMESDIR"] + CAPDESIGNTOGENOME[wildcards.capDesign] + ".sorted.genome"
	output: "annotations/{capDesign}.partition.gff"
	shell:
		'''
set +eu
conda activate bedtools_env
set -eu
partitionAnnotation.sh {input.gtf} {input.genome} | sortgff > {output}

# QC:
genomeSize=$(cat {input.genome} | cut -f2|sum.sh)
testSize=$(cat {output} | awk '{print ($5-$4)+1}' | sum.sh)
if [ $testSize -ne $genomeSize ]; then
echoerr "ERROR: sum of feature sizes in output gff is not equal to genome size. The output is probably bogus."
exit 1
fi
		'''


rule makeSirvGff:
	input: lambda wildcards: CAPDESIGNTOANNOTGTF[wildcards.capDesign]
	output: "annotations/{capDesign}.SIRVs.gff"
	shell:
		'''
cat {input} | awk '$1 ~ /SIRV/' | sortgff > {output}

		'''

rule makeGencodeBed:
	input: lambda wildcards: CAPDESIGNTOANNOTGTF[wildcards.capDesign]
	output: "annotations/{capDesign}.bed"
	shell:
		'''
cat {input} | gff2bed_full.pl - |sortbed > {output}
		'''

rule simplifyGencode:
	input: lambda wildcards: CAPDESIGNTOANNOTGTF[wildcards.capDesign]
	output: "annotations/simplified/{capDesign}.gencode.simplified_biotypes.gtf"
	shell:
		'''
uuidTmpOut=$(uuidgen)
cat {input}  | simplifyGencodeGeneTypes.pl - | sort -T {config[TMPDIR]}  -k1,1 -k4,4n -k5,5n  > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
		'''

rule collapseGencode:
	input: "annotations/simplified/{capDesign}.gencode.simplified_biotypes.gtf"
	output: "annotations/simplified/{capDesign}.gencode.collapsed.simplified_biotypes.gtf"
	threads:1
	shell:
		'''
uuid=$(uuidgen)
uuidTmpOut=$(uuidgen)
cat {input} | skipcomments | sort -T {config[TMPDIR]}  -k1,1 -k4,4n -k5,5n  | tmerge --exonOverhangTolerance {config[exonOverhangTolerance]} - |sort -T {config[TMPDIR]}  -k1,1 -k4,4n -k5,5n  > {config[TMPDIR]}/$uuid
uuidL=$(uuidgen)
set +eu
conda activate bedtools_env
set -eu
bedtools intersect -s -wao -a {config[TMPDIR]}/$uuid -b {config[TMPDIR]}/$uuid | buildLoci.pl - |sort -T {config[TMPDIR]}  -k1,1 -k4,4n -k5,5n  > {config[TMPDIR]}/$uuidL
mergeToRef.pl {input} {config[TMPDIR]}/$uuidL | sort -T {config[TMPDIR]}  -k1,1 -k4,4n -k5,5n  > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}

		'''
