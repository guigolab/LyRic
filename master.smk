import glob
from collections import defaultdict
import os
import itertools
import sys
import pandas as pd
from pprint import pprint


onsuccess:
	print("Workflow finished smoothly.")

onerror:
	print("Workflow finished with errors, see log file.")

#prefix all shell scripts with the following commands:
shell.prefix('source ~/.bashrc; set +eu; conda deactivate;  set -euo pipefail; export PATH="$PWD/utils/:$PATH";')
#(it seems DRMAA does not source ~/.bashrc by itself)
include: "functions.py"

# consistent ggplot2 theme across R plots:
GGPLOT_PUB_QUALITY="theme(axis.text= element_text(size=themeSize*1.8), axis.ticks = element_line(size=lineSize), axis.line = element_line(colour = '#595959', size=lineSize), axis.title=element_text(size = themeSize*2), panel.grid.major = element_line(colour='#d9d9d9', size=lineSize),panel.grid.minor = element_line(colour='#e6e6e6', size=minorLineSize),panel.border = element_blank(),panel.background = element_blank(), strip.background = element_rect(colour='#737373',fill='white'), legend.key.size=unit(0.5,'line'), legend.title=element_text(size=themeSize*1.2), legend.text=element_text(size=themeSize), strip.text = element_text(size = themeSize))"

#list of stat plot formats to generate (only one format at a time is supported):
plotFormat=['png',]
# color sizeFracs:
sizeFrac_Rpalette="c('0+'='#b3b3b3', '0-1' ='#f765ac','1+' ='#b370f9')"

# generic long color palette for R:
long_Rpalette="c('#8dd3c7','#ffffb3','#80ff80','#ff4dff','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd','#ccebc5','#ffed6f','#c4ff4d','#ff66a3')"

# color palette for seqPlatform's and libraryPrep's
sampleAnnot_Rpalette="list(seqPlatform = c(ONT = '#8dd3c7', pacBioSI = '#ffffb3', pacBioSII = '#5d513a'), libraryPrep = c(CapTrap = '#ff4dff', SMARTer = '#bebada', SmartSeq2 ='#751aff', Teloprime = '#fb8072', directRNA = '#80b1d3', 'Rt' = '#8c8c8c', PcrOnt='#3333ff', R2C2='#ff9900', IsoSeq='#800080'))"

# color palette for gencode biotypes
simpleBiotypes_Rpalette="c('lncRNA' = '#66ccff',  'nonExonic' = '#6666ff',  'protein_coding' = '#ff8c66', miRNA = '#808000', 'misc_RNA' = '#99ff99','pseudogene' = '#d98cb3', 'rRNA' = '#d9d9d9', 'ERCC' = '#8dd3c7', 'SIRV' = '#399384')"


# which genome build corresponds to to each capture design:
CAPDESIGNTOGENOME=config["capDesignToGenome"]
# which gencode annotation GTF correspond to each capture design:
GENOMETOANNOTGTF=config["genomeToAnnotGtf"]
# non-overlapping targeted regions for each capture design
CAPDESIGNTOTARGETSGFF=config['capDesignToTargetsGff']
# mapping of capDesign to capDesign:
CAPDESIGNTOCAPDESIGN=config["capDesignToCapDesign"]
# which CAGE peak BED file corresponds to each genome build:
GENOMETOCAGEPEAKS=config["genomeToCAGEpeaks"]
# which ENCODE DHS peak BED file corresponds to each genome build:
GENOMETODHSPEAKS=config["genomeToDHSpeaks"]
# which polyA signal file corresponds to each genome build:
GENOMETOPAS=config["genomeToPAS"]
# URL for UCSC Track Hub data files:
TRACK_HUB_DATA_URL=config["TRACK_HUB_BASE_URL"] + "dataFiles/"
# 5' and 3' end support categories for transcript models (TMs):
ENDSUPPORTcategories=["all", "cagePolyASupported"] 
	# "all" = all TMs (no filter applied)
	# "cagePolyASupported" = CAGE (5') + PolyA (3') -supported TMs only
# splicing status categories of TMs:
TMSPLICEDSTATUScategories=["all", "spliced", "unspliced"]

READFILTERS=['HiSS', 'noFilt']

GENOMES=[]
GENOMETOCAPDESIGNS=defaultdict(list)
for capD in CAPDESIGNTOGENOME:
	genome=CAPDESIGNTOGENOME[capD]
	GENOMES.append(genome)
	GENOMETOCAPDESIGNS[genome].append(capD)
GENOMES=set(GENOMES)


SIRVpresent=None
try:
	config["SIRVinfo"]
except KeyError:
	SIRVpresent=False
	print ("Will NOT compare TMs to SIRVs because config['SIRVinfo'] not found.", file=sys.stderr)
else:
	SIRVpresent=True



(TECHNAMES, CAPDESIGNS, SIZEFRACS, SAMPLEREPS) = glob_wildcards("fastqs/" + "{techname, [^_/]+}_{capDesign}_{sizeFrac}_{sampleRep}.fastq.gz")
SAMPLEREPS=set(SAMPLEREPS)
CAPDESIGNS=set(CAPDESIGNS)
SIZEFRACS=set(SIZEFRACS)
TECHNAMES=set(TECHNAMES)

# no underscores allowed in wildcards, to avoid greedy matching since we use them as separators:
wildcard_constraints:
 	capDesign = "[^_/]+",
 	sizeFrac = "[^_/]+",
 	techname = "[^_/]+",
 	sampleRep = "[^_/]+",
	minReadSupport = '|'.join(config["MINIMUM_TMERGE_READ_SUPPORT"]),
	endSupport  = '|'.join(ENDSUPPORTcategories),
	splicedStatus = '|'.join(TMSPLICEDSTATUScategories),
	filt= '|'.join(READFILTERS),
	sample_name= "[^/]+",

## read sample metadata annotations into a pandas dataframe:
sampleAnnot = pd.read_table(config['SAMPLE_ANNOT'], header=0, sep='\t')

## check that samples are properly annotated:
SAMPLES, = glob_wildcards("fastqs/" + "{sample_name, [^/]+}.fastq.gz")
for sample_name in SAMPLES:
	sampleRow=sampleAnnot[sampleAnnot.sample_name == sample_name]
	assert len(sampleRow) < 2, "Duplicate found for sample " + sample_name + " in " + config['SAMPLE_ANNOT']
	assert len(sampleRow) > 0, "Sample " + sample_name + " not found in " + config['SAMPLE_ANNOT']

#drop rows of sampleAnnot whose sample_name are not in SAMPLES list:
sampleAnnot = sampleAnnot[sampleAnnot['sample_name'].isin(SAMPLES)]

##########################################################################
### Make mapping dict of sampleRep group ID (a sample, most of the times)
### to list of sampleReps belonging to this group.
### config['sampleRepGroupBy'] contains the list of attributes to group
### sampleReps on. These attributes are column names in the sample
### annotation file (config['SAMPLE_ANNOT'])
##########################################################################

# group samples contained in sampleAnnot by attributes listed in config['sampleRepGroupBy']
sampleRepGroups=sampleAnnot.groupby(config['sampleRepGroupBy'])
# initialize dict that will contain mapping of sampleRep group ID to list of sampleReps belonging to this group:
sampleRepGroupIdToSampleReps={} 
# populate sampleRepGroupIdToSampleReps:
for key, item in sampleRepGroups:
	# remove spaces from group identifiers, because they will serve as file basenames
	sampleRepGroupId='_'.join(key).replace(' ', '')
	# assign list of sample_names (dict value) to sampleRepGroupId (dict key)
 	sampleRepGroupIdToSampleReps[sampleRepGroupId]=item['sample_name'].to_list() 
################
### Done!
################

#convert sampleAnnot to dict to facilitate later access
sampleAnnot.set_index('sample_name', inplace=True)
sampleAnnotDict = sampleAnnot.to_dict('index')



########################
### make list of authorized wildcard combinations
### inspired by https://stackoverflow.com/questions/41185567/how-to-use-expand-in-snakemake-when-some-particular-combinations-of-wildcards-ar
###

AUTHORIZEDCOMBINATIONS = [] #contains combinations of wildcards corresponding to existing FASTQs or "by*" combinations when relevant
TECHNAMESplusBY=set(TECHNAMES)
TECHNAMESplusBY.add("byTech")
SAMPLEREPSplusBY=set(SAMPLEREPS)
SAMPLEREPSplusBY.add("bySampleRep")
CAPDESIGNSplusBY=set(CAPDESIGNS)
CAPDESIGNSplusBY.add("byCapDesign")

# Populate AUTHORIZEDCOMBINATIONS:
for comb in itertools.product(TECHNAMESplusBY,CAPDESIGNSplusBY,SIZEFRACS,SAMPLEREPSplusBY):
	if comb[0] == "byTech":
		filesList=glob.glob("fastqs/" + "*" + "_" + comb[1] + "_" + comb[2] + "_" + comb[3] + ".fastq.gz")
		if(len(filesList)>1): #authorize "byTech" combination only if it's relevant
			authorizeComb(comb)

	elif comb[1] == "byCapDesign":
		filesList=glob.glob("fastqs/" + comb[0] + "_" + "*" + "_" + comb[2] + "_" + comb[3] + ".fastq.gz")
		if(len(filesList)>1): #authorize "byCapDesign" combination only if it's relevant
			authorizeComb(comb)
	elif comb[3] == "bySampleRep":
		filesList=glob.glob("fastqs/" + comb[0] + "_" + comb[1] + "_" + comb[2] + "_" + "*" + ".fastq.gz")
		if(len(filesList)>1): #authorize "bySampleRep" combination only if it's relevant
			authorizeComb(comb)
	elif(os.path.isfile("fastqs/" + comb[0] + "_" + comb[1] + "_" + comb[2] + ".fastq.gz") or
	os.path.isfile("fastqs/" + comb[0] + "_" + comb[1] + "_" + comb[2]  + "_" + comb[3] + ".fastq.gz")):
		authorizeComb(comb)


AUTHORIZEDCOMBINATIONS=set(AUTHORIZEDCOMBINATIONS)

include: "fastqStats.smk"
include: "lrMapping.smk"
include: "srMapping.smk"
include: "polyAmapping.smk"
include: "introns.smk"
include: "processReadMappings.smk"
include: "tmClassification.smk"
include: "tmEndSupport.smk"
include: "trackHub.smk"





#pseudo-rule specifying the target files we ultimately want.
rule all:
	input:
		# transcriptome GTFs (per sampleRep):
		expand("output/mappings/mergedReads/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.HiSS.tmerge.min{minReadSupport}reads.endSupport:all.gff", filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, sampleRep=SAMPLEREPS, endSupport=ENDSUPPORTcategories, minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"]),
		# transcriptome GTFs (per grouped sampleReps)
		expand("output/mappings/mergedReads/groupedSampleReps/{groupedSampleRepBasename}.min{minReadSupport}reads.splicing_status:all.endSupport:all.gff.gz", groupedSampleRepBasename=sampleRepGroupIdToSampleReps.keys(), minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"]),
		# read-to-TM mapping file (per grouped sampleReps, required by LRGASP to check what each TM contains)
		expand("output/mappings/mergedReads/groupedSampleReps/{groupedSampleRepBasename}.min{minReadSupport}reads.splicing_status:all.endSupport:all.readsToTm.tsv.gz", groupedSampleRepBasename=sampleRepGroupIdToSampleReps.keys(), minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"]),		
		expand(returnPlotFilenames("output/plots/" + "HiSS.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.HiSS.stats"), filtered_product, techname=TECHNAMESplusBY,  capDesign=CAPDESIGNSplusBY, sizeFrac=SIZEFRACS, sampleRep=SAMPLEREPSplusBY) if config['produceStatPlots'] else '/dev/null',
		expand(returnPlotFilenames("output/plots/" + "readLength.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.readLength.stats"), filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, sampleRep=SAMPLEREPS) if config['produceStatPlots'] else '/dev/null', # facetted histograms of read length
		expand("output/fastqs/" + "qc/{techname}_{capDesign}_{sizeFrac}.{sampleRep}.dupl.txt", filtered_product,techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, sampleRep=SAMPLEREPS),
		expand("output/mappings/longReadMapping/qc/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.bam.dupl.txt",filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, sampleRep=SAMPLEREPS),
		expand(returnPlotFilenames("output/plots/" + "lrMapping.basic.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.lrMapping.basic.stats"), filtered_product, techname=TECHNAMESplusBY,  capDesign=CAPDESIGNSplusBY, sizeFrac=SIZEFRACS, sampleRep=SAMPLEREPSplusBY) if config['produceStatPlots'] else '/dev/null',
		expand(returnPlotFilenames("output/plots/" + "intraPriming.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.intraPriming.stats"), filtered_product, techname=TECHNAMESplusBY,  capDesign=CAPDESIGNSplusBY, sizeFrac=SIZEFRACS, sampleRep=SAMPLEREPSplusBY) if config['produceStatPlots'] else '/dev/null',
		expand(returnPlotFilenames("output/plots/" + "lrMapping.spikeIns.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.lrMapping.spikeIns.stats"), filtered_product, techname=TECHNAMESplusBY,  capDesign=CAPDESIGNSplusBY, sizeFrac=SIZEFRACS, sampleRep=SAMPLEREPSplusBY) if config['produceStatPlots'] else '/dev/null',
		expand("output/plots/" + "hiSeq.mapping.stats/all.hiSeq.mapping.stats.{ext}", ext=plotFormat) if config['produceStatPlots'] else '/dev/null',
		expand("output/plots/" + "hiSeq.SJs.stats/all.hiSeq.SJs.stats.{ext}", ext=plotFormat) if config['produceStatPlots'] else '/dev/null',
		expand(returnPlotFilenames("output/plots/" + "merged.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.min{minReadSupport}reads.merged.stats"), filtered_product, techname=TECHNAMESplusBY,  capDesign=CAPDESIGNSplusBY, sizeFrac=SIZEFRACS, sampleRep=SAMPLEREPSplusBY, ext=plotFormat, minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"]) if config['produceStatPlots'] else '/dev/null',
		expand(returnPlotFilenames("output/plots/" + "cagePolyASupport.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.min{minReadSupport}reads.splicing_status:{splicedStatus}.cagePolyASupport.stats"), filtered_product, techname=TECHNAMESplusBY,  capDesign=CAPDESIGNSplusBY, sizeFrac=SIZEFRACS, sampleRep=SAMPLEREPSplusBY, minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"], splicedStatus=TMSPLICEDSTATUScategories) if config['produceStatPlots'] else '/dev/null',
		expand(returnPlotFilenames("output/plots/" + "tmerge.vs.gencode.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.splicing_status:{splicedStatus}.endSupport:{endSupport}.vs.gencode.stats"), filtered_product, techname=TECHNAMESplusBY,  capDesign=CAPDESIGNSplusBY, sizeFrac=SIZEFRACS, sampleRep=SAMPLEREPSplusBY, endSupport=ENDSUPPORTcategories, minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"], splicedStatus=TMSPLICEDSTATUScategories) if config['produceStatPlots'] else '/dev/null',
		expand(returnPlotFilenames("output/plots/" + "targetCoverage.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.min{minReadSupport}reads.targetCoverage.stats"), filtered_product, techname=TECHNAMESplusBY,  capDesign=CAPDESIGNSplusBY, sizeFrac=SIZEFRACS, sampleRep=SAMPLEREPSplusBY, minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"]) if config["CAPTURE"] and config['produceStatPlots'] else '/dev/null',
		expand(returnPlotFilenames("output/plots/" + "matureRNALength.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.min{minReadSupport}reads.splicing_status:{splicedStatus}.matureRNALength.hist.stats"), filtered_product, techname=TECHNAMESplusBY,  capDesign=CAPDESIGNSplusBY, sizeFrac=SIZEFRACS, sampleRep=SAMPLEREPSplusBY, minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"], splicedStatus=TMSPLICEDSTATUScategories) if config['produceStatPlots'] else '/dev/null',
		expand(returnPlotFilenames("output/plots/" + "polyAreads.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.polyAreads.stats"), filtered_product, techname=TECHNAMESplusBY,  capDesign=CAPDESIGNSplusBY, sizeFrac=SIZEFRACS, sampleRep=SAMPLEREPSplusBY) if config['produceStatPlots'] else '/dev/null',
		expand(returnPlotFilenames("output/plots/" + "tmerge.novelLoci.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.endSupport:{endSupport}.novelLoci.stats"), filtered_product, techname=TECHNAMESplusBY,  capDesign=CAPDESIGNSplusBY, sizeFrac=SIZEFRACS, sampleRep=SAMPLEREPSplusBY, endSupport=ENDSUPPORTcategories, minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"]) if config['produceStatPlots'] else '/dev/null',
		expand(returnPlotFilenames("output/plots/" + "tmerge.vs.Gencode.SJs.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.vs.Gencode.SJs.stats"), filtered_product, techname=TECHNAMESplusBY,  capDesign=CAPDESIGNSplusBY, sizeFrac=SIZEFRACS, sampleRep=SAMPLEREPSplusBY, minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"]) if config['produceStatPlots'] else '/dev/null',
		expand(returnPlotFilenames("output/plots/" + "tmerge.vs.SIRVs.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.{filt}.tmerge.min{minReadSupport}reads.vs.SIRVs.stats"), filtered_product, techname=TECHNAMESplusBY,  capDesign=CAPDESIGNSplusBY, sizeFrac=SIZEFRACS, sampleRep=SAMPLEREPSplusBY, minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"], filt=READFILTERS) if SIRVpresent and config['produceStatPlots'] else '/dev/null',
		expand(returnPlotFilenames("output/plots/" + "tmerge.vs.gencode.SnPr.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.splicing_status:{splicedStatus}.endSupport:{endSupport}.vs.gencode.SnPr.stats"), filtered_product, techname=TECHNAMESplusBY,  capDesign=CAPDESIGNSplusBY, sizeFrac=SIZEFRACS, sampleRep=SAMPLEREPSplusBY, endSupport=ENDSUPPORTcategories, minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"], splicedStatus=TMSPLICEDSTATUScategories) if config['produceStatPlots'] else '/dev/null',
		expand(returnPlotFilenames("output/plots/" + "tmerge.vs.gencode.length.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.splicing_status:bySplicingStatus.endSupport:{endSupport}.vs.gencode.length.stats"), filtered_product, techname=TECHNAMES,  capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, sampleRep=SAMPLEREPS, endSupport=ENDSUPPORTcategories, minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"]) if config['produceStatPlots'] else '/dev/null',
		expand(returnPlotFilenames("output/plots/" + "tmerge.vs.gencode.length.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.splicing_status:all.endSupport:{endSupport}.vs.gencode.length.stats"), filtered_product, techname=TECHNAMES,  capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, sampleRep=SAMPLEREPS, endSupport=ENDSUPPORTcategories, minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"]) if config['produceStatPlots'] else '/dev/null',
		expand(returnPlotFilenames("output/plots/" + "tmerge.vs.SIRVs.detection.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.vs.SIRVs.detection.stats"), filtered_product, techname=TECHNAMESplusBY,  capDesign=CAPDESIGNSplusBY, sizeFrac=SIZEFRACS, sampleRep=SAMPLEREPSplusBY, minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"]) if SIRVpresent and config['produceStatPlots'] else '/dev/null',
		expand("output/plots/" + "dhsVsCage5primeComparison.venn.stats/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.dhsVsCage5primeComparison.venn.stats.pdf", filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS,sizeFrac=SIZEFRACS, sampleRep=SAMPLEREPS, minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"]) if config['produceStatPlots'] else '/dev/null',
		expand(returnPlotFilenames("output/plots/" + "geneReadCoverage.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.geneReadCoverage.stats"), filtered_product, techname=TECHNAMESplusBY, capDesign=CAPDESIGNSplusBY, sizeFrac=SIZEFRACS, sampleRep=SAMPLEREPSplusBY) if config['produceStatPlots'] else '/dev/null',
		expand("output/plots/" + "readProfile/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.readProfile.density.png", filtered_product, techname='byTech', capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, sampleRep=SAMPLEREPS) if config['produceStatPlots'] else '/dev/null',
		expand("output/plots/" + "readProfile/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.readProfile.heatmap.png", filtered_product, techname='byTech', capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, sampleRep=SAMPLEREPS) if config['produceStatPlots'] else '/dev/null',
		expand(returnPlotFilenames("output/plots/" + "sequencingError.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.sequencingError.allErrors.stats"), filtered_product, techname=TECHNAMESplusBY,  capDesign=CAPDESIGNSplusBY, sizeFrac=SIZEFRACS, sampleRep=SAMPLEREPSplusBY) if config['produceStatPlots'] else '/dev/null',
		expand(returnPlotFilenames("output/plots/" + "sequencingError.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.sequencingError.deletionsOnly.stats"), filtered_product, techname=TECHNAMESplusBY,  capDesign=CAPDESIGNSplusBY, sizeFrac=SIZEFRACS, sampleRep=SAMPLEREPSplusBY) if config['produceStatPlots'] else '/dev/null',
		expand(returnPlotFilenames("output/plots/" + "gencode.geneDetection.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.splicing_status:{splicedStatus}.endSupport:{endSupport}.gencode.geneDetection.stats"), filtered_product, techname=TECHNAMESplusBY,  capDesign=CAPDESIGNSplusBY, sizeFrac=SIZEFRACS, sampleRep=SAMPLEREPSplusBY, endSupport=ENDSUPPORTcategories, minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"], splicedStatus=TMSPLICEDSTATUScategories) if config['produceStatPlots'] else '/dev/null',
		expand(returnPlotFilenames("output/plots/" + "readToBiotypeBreakdown.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.readToBiotypeBreakdown.stats"), filtered_product, techname=TECHNAMESplusBY,  capDesign=CAPDESIGNSplusBY, sizeFrac=SIZEFRACS, sampleRep=SAMPLEREPSplusBY) if config['produceStatPlots'] else '/dev/null',
		expand(returnPlotFilenames("output/plots/" + "tmerge.ntCoverageByGenomePartition.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{sampleRep}.tmerge.min{minReadSupport}reads.endSupport:{endSupport}.ntCoverageByGenomePartition.stats"), filtered_product, techname=TECHNAMESplusBY,  capDesign=CAPDESIGNSplusBY, sizeFrac=SIZEFRACS, sampleRep=SAMPLEREPSplusBY, endSupport=ENDSUPPORTcategories, minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"]) if config['produceStatPlots'] else '/dev/null',

## temporary/intermediate
		"output/statsFiles/" + "all.highConfSplicedReads.stats.tsv" if config['produceHtmlStatsTable'] else '/dev/null',
		"output/statsFiles/" + "all.fastq.timestamps.tsv" if config['produceHtmlStatsTable'] else '/dev/null',
		"output/statsFiles/" + "all.readlength.summary.tsv" if config['produceHtmlStatsTable'] else '/dev/null',
		expand("output/statsFiles/" + "all.min{minReadSupport}reads.matureRNALengthSummary.stats.tsv", minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"]) if config['produceHtmlStatsTable'] else '/dev/null',
		expand("output/statsFiles/" + "all.tmerge.min{minReadSupport}reads.endSupport:{endSupport}.novelLoci.qc.stats.tsv", minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"], endSupport=ENDSUPPORTcategories) if config['produceHtmlStatsTable'] else '/dev/null',
## temporary



		#################
		### Track hub ###
		#################
		"output/trackHub/" + "hub.txt" if config['produceTrackHub'] else '/dev/null',
		"output/trackHub/" + "genomes.txt" if config['produceTrackHub'] else '/dev/null',
		expand("output/trackHub/" + "{genome}/trackDb.txt", genome=GENOMES) if config['produceTrackHub'] else '/dev/null',




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
		gtf=lambda wildcards: GENOMETOANNOTGTF[CAPDESIGNTOGENOME[wildcards.capDesign]],
		genome = lambda wildcards: config["GENOMESDIR"] + CAPDESIGNTOGENOME[wildcards.capDesign] + ".sorted.genome"
	conda: "envs/xtools_env.yml"
	output: "output/annotations/{capDesign}.partition.gff"
	shell:
		'''

partitionAnnotation.sh {input.gtf} {input.genome} | sortgff > {output}

# QC:
genomeSize=$(cat {input.genome} | cut -f2|sum.sh)
testSize=$(cat {output} | awk '{{print ($5-$4)+1}}' | sum.sh)
if [ $testSize -ne $genomeSize ]; then
echoerr "ERROR: sum of feature sizes in output gff is not equal to genome size. The output is probably bogus."
exit 1
fi
		'''


rule makeSirvGff:
	input: lambda wildcards: GENOMETOANNOTGTF[CAPDESIGNTOGENOME[wildcards.capDesign]]
	output: "output/annotations/{capDesign}.SIRVs.gff"
	shell:
		'''
cat {input} | awk '$1 ~ /SIRV/' | sortgff > {output}

		'''

rule makeGencodeBed:
	input: lambda wildcards: GENOMETOANNOTGTF[CAPDESIGNTOGENOME[wildcards.capDesign]]
	output: "output/annotations/{capDesign}.bed"
	shell:
		'''
cat {input} | gff2bed_full.pl - |sortbed > {output}
		'''

rule simplifyGencode:
	input: lambda wildcards: GENOMETOANNOTGTF[CAPDESIGNTOGENOME[wildcards.capDesign]]
	output: "output/annotations/simplified/{capDesign}.gencode.simplified_biotypes.gtf"
	shell:
		'''
uuidTmpOut=$(uuidgen)
cat {input}  | simplifyGencodeGeneTypes.pl - | sort -T {config[TMPDIR]}  -k1,1 -k4,4n -k5,5n  > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}
		'''

rule collapseGencode:
	input: "output/annotations/simplified/{capDesign}.gencode.simplified_biotypes.gtf"
	output: "output/annotations/simplified/{capDesign}.gencode.collapsed.simplified_biotypes.gtf"
	threads:1
	conda: "envs/xtools_env.yml"
	shell:
		'''
uuid=$(uuidgen)
uuidTmpOut=$(uuidgen)
cat {input} | skipcomments | sort -T {config[TMPDIR]}  -k1,1 -k4,4n -k5,5n  | tmerge --exonOverhangTolerance {config[exonOverhangTolerance]} - |sort -T {config[TMPDIR]}  -k1,1 -k4,4n -k5,5n  > {config[TMPDIR]}/$uuid
uuidL=$(uuidgen)

bedtools intersect -s -wao -a {config[TMPDIR]}/$uuid -b {config[TMPDIR]}/$uuid | buildLoci.pl - |sort -T {config[TMPDIR]}  -k1,1 -k4,4n -k5,5n  > {config[TMPDIR]}/$uuidL
mergeToRef.pl {input} {config[TMPDIR]}/$uuidL | sort -T {config[TMPDIR]}  -k1,1 -k4,4n -k5,5n  > {config[TMPDIR]}/$uuidTmpOut
mv {config[TMPDIR]}/$uuidTmpOut {output}

		'''
