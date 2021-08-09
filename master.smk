import glob
from collections import defaultdict
import os
import itertools
import sys
import pandas as pd
from pprint import pprint

#prefix all shell scripts with the following commands:
shell.prefix('source ~/.bashrc; set +eu; conda deactivate;  set -euo pipefail; export PATH="$PWD/utils/:$PATH";')
#(it seems DRMAA does not source ~/.bashrc by itself)

# consistent ggplot2 theme across R plots:
GGPLOT_PUB_QUALITY="theme(axis.text= element_text(size=themeSize*1.8), axis.ticks = element_line(size=lineSize), axis.line = element_line(colour = '#595959', size=lineSize), axis.title=element_text(size = themeSize*2), panel.grid.major = element_line(colour='#d9d9d9', size=lineSize),panel.grid.minor = element_line(colour='#e6e6e6', size=minorLineSize),panel.border = element_blank(),panel.background = element_blank(), strip.background = element_rect(colour='#737373',fill='white'), legend.key.size=unit(0.5,'line'), legend.title=element_text(size=themeSize*1.2), legend.text=element_text(size=themeSize), strip.text = element_text(size = themeSize))"

# color sizeFracs:
sizeFrac_Rpalette="c('0+'='#b3b3b3', '0-1' ='#f765ac','1+' ='#b370f9')"

# generic long color palette for R:
long_Rpalette="c('#8dd3c7','#ffffb3','#80ff80','#ff4dff','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd','#ccebc5','#ffed6f','#c4ff4d','#ff66a3')"

# color palette for seqPlatform's and libraryPrep's
sampleAnnot_Rpalette="list(seqPlatform = c(ONT = '#8dd3c7', pacBioSI = '#ffffb3', pacBioSII = '#5d513a'), libraryPrep = c(CapTrap = '#ff4dff', SMARTer = '#bebada', SmartSeq2 ='#751aff', Teloprime = '#fb8072', directRNA = '#80b1d3', 'Rt' = '#8c8c8c', PcrOnt='#3333ff', R2C2='#ff9900', IsoSeq='#800080'))"

# color palette for gencode biotypes
simpleBiotypes_Rpalette="c('lncRNA' = '#66ccff',  'nonExonic' = '#6666ff',  'protein_coding' = '#ff8c66', miRNA = '#808000', 'misc_RNA' = '#99ff99','pseudogene' = '#d98cb3', 'rRNA' = '#d9d9d9', 'ERCC' = '#8dd3c7', 'SIRV' = '#399384')"

# in which directory to find input FASTQ files:
LR_FASTQDIR=config["LR_FASTQDIR"]
# which genome build corresponds to to each capture design:
CAPDESIGNTOGENOME=config["capDesignToGenome"]
# which gencode annotation GTF correspond to each capture design:
CAPDESIGNTOANNOTGTF=config["capDesignToAnnotGtf"]
# which CAGE peak BED file corresponds to each genome build:
GENOMETOCAGEPEAKS=config["genomeToCAGEpeaks"]
# which ENCODE DHS peak BED file corresponds to each genome build:
GENOMETODHSPEAKS=config["genomeToDHSpeaks"]
# which polyA signal file corresponds to each genome build:
GENOMETOPAS=config["genomeToPAS"]
# where to write dummy files 
DUMMY_DIR="dummy/"

TRACK_HUB_DATA_URL=config["TRACK_HUB_BASE_URL"] + "dataFiles/"

print ("Input LR FASTQs are in: " + LR_FASTQDIR, file=sys.stderr)





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

for genome in GENOMETOCAPDESIGNS:
	mergedCapDString=''
	capDperGenome=0
	for capD in sorted(GENOMETOCAPDESIGNS[genome]):
		if capD.find('preCap') == -1:  #avoid merging pre-capture with post-capture
			mergedCapDString+=capD
			capDperGenome+=1
	if capDperGenome > 1:
		for capD in GENOMETOCAPDESIGNS[genome]:
			CAPDESIGNTOANNOTGTF[mergedCapDString]=CAPDESIGNTOANNOTGTF[capD]
		CAPDESIGNTOGENOME[mergedCapDString]=genome


for capD in CAPDESIGNTOGENOME:
	genome=CAPDESIGNTOGENOME[capD]

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



DEMULTIPLEXED_FASTQS=LR_FASTQDIR
(TECHNAMES, CAPDESIGNS, SIZEFRACS, BARCODES) = glob_wildcards(LR_FASTQDIR + "{techname, [^_/]+}_{capDesign}_{sizeFrac}_{barcodes}.fastq.gz")
BARCODES=set(BARCODES)
BARCODESUNDETER=list(BARCODES)


CAPDESIGNS=set(CAPDESIGNS)
SIZEFRACS=set(SIZEFRACS)
TECHNAMES=set(TECHNAMES)


wildcard_constraints:
 	capDesign = "[^_/]+",
 	sizeFrac = "[^_/]+",
 	techname = "[^_/]+",
 	barcodes = "[^_/]+",
	minReadSupport = '|'.join(config["MINIMUM_TMERGE_READ_SUPPORT"]),
	endSupport  = '|'.join(ENDSUPPORTcategories),
	splicedStatus = '|'.join(TMSPLICEDSTATUScategories),
	filt= '|'.join(READFILTERS),
	sample_name= "[^/]+",

## read sample metadata annotations into a pandas dataframe:
sampleAnnot = pd.read_table(config['SAMPLE_ANNOT'], header=0, sep='\t')

## check that samples are properly annotated:
SAMPLES, = glob_wildcards(LR_FASTQDIR + "{sample_name, [^/]+}.fastq.gz")
for sample_name in SAMPLES:
	sampleRow=sampleAnnot[sampleAnnot.sample_name == sample_name]
	assert len(sampleRow) < 2, "Duplicate found for sample " + sample_name + " in " + config['SAMPLE_ANNOT']
	assert len(sampleRow) > 0, "Sample " + sample_name + " not found in " + config['SAMPLE_ANNOT']

sampleAnnot.set_index('sample_name', inplace=True)
sampleAnnotDict = sampleAnnot.to_dict('index')



########################
### make list of authorized wildcard combinations
### inspired by https://stackoverflow.com/questions/41185567/how-to-use-expand-in-snakemake-when-some-particular-combinations-of-wildcards-ar
###



AUTHORIZEDCOMBINATIONS = [] #contains combinations of wildcards corresponding to existing FASTQs
TECHNAMESplusBY=set(TECHNAMES)
TECHNAMESplusBY.add("byTech")
BARCODESplusBY=set(BARCODES)
BARCODESplusBY.add("byTissue")
CAPDESIGNSplusBY=set(CAPDESIGNS)
CAPDESIGNSplusBY.add("byCapDesign")

for comb in itertools.product(TECHNAMESplusBY,CAPDESIGNSplusBY,SIZEFRACS,BARCODESplusBY):
	if(os.path.isfile(LR_FASTQDIR + comb[0] + "_" + comb[1] + "_" + comb[2] + ".fastq.gz") or
	os.path.isfile(LR_FASTQDIR + comb[0] + "_" + comb[1] + "_" + comb[2]  + "_" + comb[3] + ".fastq.gz")):
		AUTHORIZEDCOMBINATIONS.append((("techname", comb[0]),("capDesign", comb[1])))
		AUTHORIZEDCOMBINATIONS.append((("techname", comb[0]),("capDesign", comb[1]),("sizeFrac", comb[2]),("barcodes", comb[3])))
		for byComb in itertools.product([("techname", comb[0]), ("techname", "byTech")], [("capDesign", comb[1]), ("capDesign", "byCapDesign")], [("sizeFrac", comb[2]),], [("barcodes", comb[3]), ("barcodes", "byTissue")]):
			AUTHORIZEDCOMBINATIONS.append((byComb))



AUTHORIZEDCOMBINATIONS=set(AUTHORIZEDCOMBINATIONS)

def returnPlotFilenames(basename):
	plotsList=[]
	for comb in itertools.product(["legendOnly"], config["PLOTFORMATS"]):
		plotsList.append(basename + "." + comb[0] + "." + comb[1])
	for comb in itertools.product(["xy", "yx"], ["wLegend", "woLegend"], config["PLOTFORMATS"]):
		plotsList.append(basename + "." + comb[0] + "." + comb[1] + "." + comb[2])
	return plotsList



def filtered_product(*args):
#return combinations of wildcards that correspond to existing combinations in input (i.e. corresponding input FASTQ exists)
#args[0]: techname
#args[1]: capDesign
#args[2]: sizeFrac (optional)
#args[3]: barcodes (optional)
#args[>3] are irrelevant

	found=False
	if len(args) < 2:
		quit("Error in function filtered_product (wrong number of arguments, shoud be >1).")
	elif len(args) == 2:
		for wc_comb in itertools.product(*args):
			if wc_comb in AUTHORIZEDCOMBINATIONS:
				found=True
				yield(wc_comb)
	else:
		for wc_comb in itertools.product(*args):
			if wc_comb[0:4] in AUTHORIZEDCOMBINATIONS:
				if wc_comb[1] in (('capDesign', 'byCapDesign'),):
					if wc_comb[3] not in (('barcodes', 'byTissue'),): #DO NOT plot all tissues and capDesigns on the same axis (too big)
						found=True
						yield(wc_comb)
					else:
						continue
				else:
					found=True
					yield(wc_comb)
	if not found:
		pprint("Error in function filtered_product. Args were:", sys.stderr)
		pprint((args), sys.stderr)
		quit("Error. Could not yield any input file.")


def multi_figures(capDesign,sizeFrac,barcodes,techname, splicing_status=None):
	figure_settings = dict();
	figure_settings['technameFilterString']=''
	figure_settings['hideXaxisLabels']="theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +"
	removeFacetLabels=''
	if techname not in ('byTech') and capDesign not in ('byCapDesign') and barcodes not in ('byTissue'):
		removeFacetLabels = """
plotFacetXy <- parse(text =paste(plotFacetXy, \\" + theme(strip.text.x = element_blank(), strip.text.y = element_blank())\\"))
plotFacetYx <- parse(text=paste(plotFacetYx, \\" + theme(strip.text.x = element_blank(), strip.text.y = element_blank())\\"))
"""		
	figure_settings['splicingStatusFilterString']=''
	if splicing_status != None:
		if splicing_status == 'spliced':
			figure_settings['splicingStatusFilterString']="dat <- subset(dat, spliced==1)\n"
		elif splicing_status == 'unspliced':
			figure_settings['splicingStatusFilterString']="dat <- subset(dat, spliced==0)\n"
		elif splicing_status == 'all':
			figure_settings['splicingStatusFilterString']=''
		else:
			quit("Error. Invalid value for splicing_status in function multi_figures.")
	if techname not in ('byTech'):
		figure_settings['technameFilterString']="dat <- subset(dat, seqTech=='" + techname + "')\n"
	figure_settings['substSeqTechString']="""
dat\$seqTech <- gsub(':$', '', dat\$seqTech)
dat\$seqTech <- gsub(':', '\\n', dat\$seqTech)
"""
	figure_settings['substTissueString']="dat\$tissue <- gsub('" + capDesign + "_', '', dat\$tissue)\n"
	figure_settings['capDesignFilterString']=''
	figure_settings['sizeFracFilterString']=''
	figure_settings['tissueFilterString']=''
	figure_settings['graphDimensions']="""

horizCats <- length(unique(dat\$capDesign)) * length(unique(dat\$tissue))
vertCats <- length(unique(dat\$seqTech))

wXyPlot = (horizCats * 0.9) +1.7
hXyPlot = (vertCats * 0.6) + 1.7

geom_textSize=1.4 
themeSize = (14/5) * geom_textSize
# https://stackoverflow.com/questions/25061822/ggplot-geom-text-font-size-control/25062509
lineSize=geom_textSize/8
minorLineSize=lineSize/2

"""
	figure_settings['facetPlotSetup']= f"""
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

	if capDesign not in ('byCapDesign'):
		figure_settings['capDesignFilterString']+="dat <- subset(dat, capDesign=='" + capDesign + "')\n"

	figure_settings['sizeFracFilterString']+="dat <- subset(dat, sizeFrac=='" + sizeFrac + "')\n"
	if barcodes not in ('byTissue'):
		figure_settings['tissueFilterString']+="dat <- subset(dat, tissue=='" + barcodes + "')\n"
	return(figure_settings)


def trackHubSubGroupString(techname, capDesign, sizeFrac, barcodes, minReadSupport):
	techname=(("techname", techname),)
	capDesign=(("capDesign", capDesign),)
	sizeFracs=[]
	for sizeF in sizeFrac:
		sizeFracs.append(("sizeFrac", sizeF))
	barCodes=[]
	for barC in barcodes:
		barCodes.append(("barcodes", barC))
	minRS=[]
	for minReadSupport in minReadSupport:
		minRS.append(("minReadSupport", minReadSupport))
	returnSubGroup1String="subGroup1 sample Sample"
	returnSubGroup2String="subGroup2 sizeFraction Size_fraction"
	returnSubGroup3String="subGroup3 minReadSupport Min_read_support_per_TM"

	returnSubGroup1StringData=[]
	returnSubGroup2StringData=[]
	returnSubGroup3StringData=[]
	for wc_comb in itertools.product(techname,capDesign, sizeFracs, barCodes, minRS):
		if wc_comb in AUTHORIZEDCOMBINATIONS:
			returnSubGroup1StringData.append(wc_comb[3][1] + "=" + wc_comb[3][1])
			returnSubGroup2StringData.append(wc_comb[2][1] + "=" + wc_comb[2][1])
			returnSubGroup3StringData.append(wc_comb[4][1] + "=" + wc_comb[4][1])

	returnSubGroup1StringData=set(returnSubGroup1StringData)
	returnSubGroup2StringData=set(returnSubGroup2StringData)
	returnSubGroup3StringData=set(returnSubGroup3StringData)
	
	return(returnSubGroup1String + " " + " ".join(sorted(returnSubGroup1StringData)) + "\n" + returnSubGroup2String + " " + " ".join(sorted(returnSubGroup2StringData)) + "\n" + returnSubGroup3String + " " + " ".join(sorted(returnSubGroup3StringData)))


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
		expand(returnPlotFilenames(config["PLOTSDIR"] + "HiSS.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{barcodes}.HiSS.stats"), filtered_product, techname=TECHNAMESplusBY,  capDesign=CAPDESIGNSplusBY, sizeFrac=SIZEFRACS, barcodes=BARCODESplusBY),
		expand("mappings/nonAnchoredMergeReads/colored/{techname}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.splicing_status:all.endSupport:{endSupport}.bed", filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES, endSupport=ENDSUPPORTcategories, minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"]),
		config["STATSDATADIR"] + "all.highConfSplicedReads.stats.tsv",

		expand(LR_FASTQDIR + "{techname}_{capDesign}_{sizeFrac}_{barcodes}.fastq.gz", filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES),

		expand(returnPlotFilenames(config["PLOTSDIR"] + "readLength.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{barcodes}.readLength.stats"), filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES), # facetted histograms of read length
		expand(LR_FASTQDIR + "qc/{techname}_{capDesign}_{sizeFrac}.{barcodes}.dupl.txt", filtered_product,techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES),
		expand("mappings/readMapping/qc/{techname}_{capDesign}_{sizeFrac}_{barcodes}.bam.dupl.txt",filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES),
		expand(returnPlotFilenames(config["PLOTSDIR"] + "lrMapping.basic.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{barcodes}.lrMapping.basic.stats"), filtered_product, techname=TECHNAMESplusBY,  capDesign=CAPDESIGNSplusBY, sizeFrac=SIZEFRACS, barcodes=BARCODESplusBY),
		expand(returnPlotFilenames(config["PLOTSDIR"] + "intraPriming.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{barcodes}.intraPriming.stats"), filtered_product, techname=TECHNAMESplusBY,  capDesign=CAPDESIGNSplusBY, sizeFrac=SIZEFRACS, barcodes=BARCODESplusBY),
		expand(returnPlotFilenames(config["PLOTSDIR"] + "lrMapping.spikeIns.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{barcodes}.lrMapping.spikeIns.stats"), filtered_product, techname=TECHNAMESplusBY,  capDesign=CAPDESIGNSplusBY, sizeFrac=SIZEFRACS, barcodes=BARCODESplusBY),

		expand(config["PLOTSDIR"] + "hiSeq.mapping.stats/all.hiSeq.mapping.stats.{ext}", ext=config["PLOTFORMATS"]),
		expand(config["PLOTSDIR"] + "hiSeq.SJs.stats/all.hiSeq.SJs.stats.{ext}", ext=config["PLOTFORMATS"]),
		expand ("mappings/clusterPolyAsites/{techname}_{capDesign}_{sizeFrac}_{barcodes}.polyAsites.clusters.bed", filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES),
		
		
		
		
		expand(returnPlotFilenames(config["PLOTSDIR"] + "merged.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.merged.stats"), filtered_product, techname=TECHNAMESplusBY,  capDesign=CAPDESIGNSplusBY, sizeFrac=SIZEFRACS, barcodes=BARCODESplusBY, ext=config["PLOTFORMATS"], minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"]),
		expand(returnPlotFilenames(config["PLOTSDIR"] + "cagePolyASupport.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.splicing_status:{splicedStatus}.cagePolyASupport.stats"), filtered_product, techname=TECHNAMESplusBY,  capDesign=CAPDESIGNSplusBY, sizeFrac=SIZEFRACS, barcodes=BARCODESplusBY, minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"], splicedStatus=TMSPLICEDSTATUScategories),
		expand(returnPlotFilenames(config["PLOTSDIR"] + "tmerge.vs.gencode.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.splicing_status:{splicedStatus}.endSupport:{endSupport}.vs.gencode.stats"), filtered_product, techname=TECHNAMESplusBY,  capDesign=CAPDESIGNSplusBY, sizeFrac=SIZEFRACS, barcodes=BARCODESplusBY, endSupport=ENDSUPPORTcategories, minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"], splicedStatus=TMSPLICEDSTATUScategories),
		expand(returnPlotFilenames(config["PLOTSDIR"] + "targetCoverage.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.targetCoverage.stats"), filtered_product, techname=TECHNAMESplusBY,  capDesign=CAPDESIGNSplusBY, sizeFrac=SIZEFRACS, barcodes=BARCODESplusBY, minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"]) if config["CAPTURE"]  else expand( DUMMY_DIR + "dummy{number}.txt", number='15'),

		expand(returnPlotFilenames(config["PLOTSDIR"] + "matureRNALength.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{barcodes}.min{minReadSupport}reads.splicing_status:{splicedStatus}.matureRNALength.hist.stats"), filtered_product, techname=TECHNAMESplusBY,  capDesign=CAPDESIGNSplusBY, sizeFrac=SIZEFRACS, barcodes=BARCODESplusBY, minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"], splicedStatus=TMSPLICEDSTATUScategories),

		expand(returnPlotFilenames(config["PLOTSDIR"] + "polyAreads.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{barcodes}.polyAreads.stats"), filtered_product, techname=TECHNAMESplusBY,  capDesign=CAPDESIGNSplusBY, sizeFrac=SIZEFRACS, barcodes=BARCODESplusBY),
		expand(returnPlotFilenames(config["PLOTSDIR"] + "tmerge.novelLoci.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.endSupport:{endSupport}.novelLoci.stats"), filtered_product, techname=TECHNAMESplusBY,  capDesign=CAPDESIGNSplusBY, sizeFrac=SIZEFRACS, barcodes=BARCODESplusBY, endSupport=ENDSUPPORTcategories, minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"]),
		expand(returnPlotFilenames(config["PLOTSDIR"] + "tmerge.vs.Gencode.SJs.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.vs.Gencode.SJs.stats"), filtered_product, techname=TECHNAMESplusBY,  capDesign=CAPDESIGNSplusBY, sizeFrac=SIZEFRACS, barcodes=BARCODESplusBY, minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"]),
		expand(returnPlotFilenames(config["PLOTSDIR"] + "tmerge.vs.SIRVs.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{barcodes}.{filt}.tmerge.min{minReadSupport}reads.vs.SIRVs.stats"), filtered_product, techname=TECHNAMESplusBY,  capDesign=CAPDESIGNSplusBY, sizeFrac=SIZEFRACS, barcodes=BARCODESplusBY, minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"], filt=READFILTERS) if SIRVpresent  else expand( DUMMY_DIR + "dummy{number}.txt", number='13'),
		expand(returnPlotFilenames(config["PLOTSDIR"] + "tmerge.vs.gencode.SnPr.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.splicing_status:{splicedStatus}.endSupport:{endSupport}.vs.gencode.SnPr.stats"), filtered_product, techname=TECHNAMESplusBY,  capDesign=CAPDESIGNSplusBY, sizeFrac=SIZEFRACS, barcodes=BARCODESplusBY, endSupport=ENDSUPPORTcategories, minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"], splicedStatus=TMSPLICEDSTATUScategories),
		expand(returnPlotFilenames(config["PLOTSDIR"] + "tmerge.vs.gencode.length.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.splicing_status:bySplicingStatus.endSupport:{endSupport}.vs.gencode.length.stats"), filtered_product, techname=TECHNAMES,  capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES, endSupport=ENDSUPPORTcategories, minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"]),
		expand(returnPlotFilenames(config["PLOTSDIR"] + "tmerge.vs.gencode.length.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.splicing_status:all.endSupport:{endSupport}.vs.gencode.length.stats"), filtered_product, techname=TECHNAMES,  capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES, endSupport=ENDSUPPORTcategories, minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"]),

		expand(returnPlotFilenames(config["PLOTSDIR"] + "tmerge.vs.SIRVs.detection.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.vs.SIRVs.detection.stats"), filtered_product, techname=TECHNAMESplusBY,  capDesign=CAPDESIGNSplusBY, sizeFrac=SIZEFRACS, barcodes=BARCODESplusBY, minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"]) if SIRVpresent  else expand( DUMMY_DIR + "dummy{number}.txt", number='14'),
		expand(config["PLOTSDIR"] + "dhsVsCage5primeComparison.venn.stats/{techname}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.dhsVsCage5primeComparison.venn.stats.pdf", filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS,sizeFrac=SIZEFRACS, barcodes=BARCODES, minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"]),
		expand(returnPlotFilenames(config["PLOTSDIR"] + "geneReadCoverage.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{barcodes}.geneReadCoverage.stats"), filtered_product, techname=TECHNAMESplusBY, capDesign=CAPDESIGNSplusBY, sizeFrac=SIZEFRACS, barcodes=BARCODESplusBY), # facetted histograms of read length
		expand(config["PLOTSDIR"] + "readProfile/{techname}_{capDesign}_{sizeFrac}_{barcodes}.readProfile.density.png", filtered_product, techname='byTech', capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES),
		expand(config["PLOTSDIR"] + "readProfile/{techname}_{capDesign}_{sizeFrac}_{barcodes}.readProfile.heatmap.png", filtered_product, techname='byTech', capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES),
		expand(returnPlotFilenames(config["PLOTSDIR"] + "sequencingError.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{barcodes}.sequencingError.allErrors.stats"), filtered_product, techname=TECHNAMESplusBY,  capDesign=CAPDESIGNSplusBY, sizeFrac=SIZEFRACS, barcodes=BARCODESplusBY),
		expand(returnPlotFilenames(config["PLOTSDIR"] + "sequencingError.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{barcodes}.sequencingError.deletionsOnly.stats"), filtered_product, techname=TECHNAMESplusBY,  capDesign=CAPDESIGNSplusBY, sizeFrac=SIZEFRACS, barcodes=BARCODESplusBY),
		expand(returnPlotFilenames(config["PLOTSDIR"] + "gencode.geneDetection.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.splicing_status:{splicedStatus}.endSupport:{endSupport}.gencode.geneDetection.stats"), filtered_product, techname=TECHNAMESplusBY,  capDesign=CAPDESIGNSplusBY, sizeFrac=SIZEFRACS, barcodes=BARCODESplusBY, endSupport=ENDSUPPORTcategories, minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"], splicedStatus=TMSPLICEDSTATUScategories),
		expand(returnPlotFilenames(config["PLOTSDIR"] + "readToBiotypeBreakdown.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{barcodes}.readToBiotypeBreakdown.stats"), filtered_product, techname=TECHNAMESplusBY,  capDesign=CAPDESIGNSplusBY, sizeFrac=SIZEFRACS, barcodes=BARCODESplusBY),
		expand(returnPlotFilenames(config["PLOTSDIR"] + "tmerge.ntCoverageByGenomePartition.stats/{techname}/{capDesign}/{techname}_{capDesign}_{sizeFrac}_{barcodes}.tmerge.min{minReadSupport}reads.endSupport:{endSupport}.ntCoverageByGenomePartition.stats"), filtered_product, techname=TECHNAMESplusBY,  capDesign=CAPDESIGNSplusBY, sizeFrac=SIZEFRACS, barcodes=BARCODESplusBY, endSupport=ENDSUPPORTcategories, minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"]),

## temporary/intermediate
		config["STATSDATADIR"] + "all.fastq.timestamps.tsv",
		config["STATSDATADIR"] + "all.readlength.summary.tsv",
		expand("annotations/{capDesign}.partition.gff", capDesign=CAPDESIGNS),
		expand(config["STATSDATADIR"] + "all.min{minReadSupport}reads.matureRNALengthSummary.stats.tsv", minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"]),
		expand(config["STATSDATADIR"] + "all.tmerge.min{minReadSupport}reads.endSupport:{endSupport}.novelLoci.qc.stats.tsv", minReadSupport=config["MINIMUM_TMERGE_READ_SUPPORT"], endSupport=ENDSUPPORTcategories),
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
	conda: "envs/xtools_env.yml"
	output: "annotations/{capDesign}.partition.gff"
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
