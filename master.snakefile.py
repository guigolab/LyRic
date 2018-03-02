import glob
from collections import defaultdict
import os.path
from itertools import product


# # path to dropbox folder where to sync output R plots:
DROPBOX_PLOTS=config["DROPBOX_PLOTSDIR"]
GGPLOT_PUB_QUALITY=config["GGPLOT_PUB_QUALITY"]
#TECHNAME=config["TECHNAME"]
CAPDESIGNTOGENOME=config["capDesignToGenome"]
CAPDESIGNTOANNOTGTF=config["capDesignToAnnotGtf"]
#mappingDir=config["PB_MAPPINGS"]
sizeFrac_Rpalette=config["SIZEFRACTION_RPALETTE"]
#print(CAPDESIGNTOGENOME)

# no underscores allowed in wildcards, to avoid greedy matching since we use them as separators
wildcard_constraints:
 	capDesign = "[^_/]+",
 	sizeFrac = "[^_/]+",
 	techname = "[^_/]+",
# 	barcodes = "[^\.]",
# 	barcodesU ="^(merged)"

# get CAPDESIGNS (capture designs, i.e. Hv1, Hv2, Mv1) and SIZEFRACS (size fractions) variables from FASTQ file names (warning: this will generate duplicate entries):
(TECHNAMES, CAPDESIGNS, SIZEFRACS) = glob_wildcards(config["FQPATH"] + "{techname}_{capDesign}_{sizeFrac}.fastq.gz")
# remove duplicate entries:
CAPDESIGNS=set(CAPDESIGNS)
SIZEFRACS=set(SIZEFRACS)
TECHNAMES=set(TECHNAMES)


########################
### make list of forbidden combinations (missing input files)
### inspired by https://stackoverflow.com/questions/41185567/how-to-use-expand-in-snakemake-when-some-particular-combinations-of-wildcards-ar

AUTHORIZEDCOMBINATIONS = []

for comb in product(TECHNAMES,CAPDESIGNS,SIZEFRACS):
	if(os.path.isfile(config["FQPATH"] + comb[0] + "_" + comb[1] + "_" + comb[2] + ".fastq.gz")):
		tup=(("techname", comb[0]),("capDesign", comb[1]),("sizeFrac", comb[2]))
		AUTHORIZEDCOMBINATIONS.append(tup)

#print ("AUTHORIZEDCOMBINATIONS:")
#print (AUTHORIZEDCOMBINATIONS)

#print("\n\nTODO:\n  - check that read IDs appear only once in output after mappings (BAMs/BEDs)")

def filter_combinator(whitelist):
	def filtered_combinator(*args, **kwargs):
		for wc_comb in product(*args, **kwargs):
#			print ("wc_comb 0-3:")
#			print (wc_comb[0:3])
            #use only first 2/3 elements of args, which should be techname, capDesign, sizeFrac. We don't care about the rest, as they're not present in input files
			found=False
			for ac in AUTHORIZEDCOMBINATIONS:
#				print("current AUTH:")
#				print (ac)
				if(wc_comb[0:3] == ac):
					found=True
					yield(wc_comb)
					break
			if not found and wc_comb[2][0] != 'sizeFrac':
#if this point is reached it means wc_comb[0:3] was not found. Try with wc_comb[0:2]
				for ac in AUTHORIZEDCOMBINATIONS:
#				print("current AUTH:")
#				print (ac)
					if(wc_comb[0:2] == ac[0:2]):
						found=True
						yield(wc_comb)
						break
#				for currAc in ac:
#					print(currAc)

			# if (frozenset(wc_comb[0:3]) in whitelist):
			# 	print(" OK")
			# 	yield wc_comb
			# # elif (frozenset(wc_comb[0:2]) not in blacklist):
			# # 	print(" OK")
			# # 	yield wc_comb
			# else:
			# 	print(" FORBIDDEN")
	return filtered_combinator

#filtered_product = filter_combinator(FORBIDDENCOMBINATIONS)
filtered_product = filter_combinator(AUTHORIZEDCOMBINATIONS)


#check that barcode IDs and capDesign match (i.e. prefix of barcode ID matches capDesign)
def capDesignBarcodeMatch(*args):
	for wc_comb in product(*args):
#		print (wc_comb)
		capDesign=wc_comb[1][1]
		barcode=wc_comb[3][1]
#		print (capDesign, barcode)
		if (barcode).find(capDesign) == 0:
	 #or wildcards.barcodesU == 'Undeter':
#			print ("MATCH ")
			yield wc_comb


adaptersTSV = "demultiplexing/all_adapters.tsv"
f = open(adaptersTSV, 'r')
BARCODES = []
BARCODESUNDETER = []
CAPDESIGNTOBARCODES = defaultdict(list)
for line in f:
	columns = line.split("\t")
	barcodeId = columns[0].split("_")
	if columns[0] != 'UP':
		BARCODES.append(columns[0])
		BARCODESUNDETER.append(columns[0])
		CAPDESIGNTOBARCODES[barcodeId[0]].append(columns[0])

BARCODESUNDETER.append("Undeter")
BARCODES=set(BARCODES)
BARCODESUNDETER=set(BARCODESUNDETER)

include: "demultiplex.snakefile.py"
include: "fastqStats.snakefile.py"
include: "lrMapping.snakefile.py"
include: "srMapping.snakefile.py"
include: "polyAmapping.snakefile.py"

#to avoid AmbiguousRuleException: (and in that order, otherwise {capDesign}_{sizeFrac}.Undeter.fastq will be generated using rule demultiplexFastqs and not getUndeterminedReads, which will always give empty Undeter files :
#ruleorder: getUndeterminedReads > demultiplexFastqs
#ruleorder: mergeSizeFracBams > readMapping

#pseudo-rule specifying the target files we ultimately want.
rule all:
	input:
		expand(config["PLOTSDIR"] + "{techname}_{capDesign}_all.readlength.{ext}", filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, ext=config["PLOTFORMATS"]), # facetted histograms of read length
		expand(config["FQPATH"] + "qc/{techname}_{capDesign}_{sizeFrac}.dupl.txt", filtered_product,techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS),
 		expand(config["PLOTSDIR"] + "{techname}.fastq.UP.stats.{ext}", techname=TECHNAMES, ext=config["PLOTFORMATS"]), # UP reads plots
 		expand(config["PLOTSDIR"] + "{techname}.fastq.BC.stats.{ext}", techname=TECHNAMES, ext=config["PLOTFORMATS"]), # barcode reads plots
 		expand(config["PLOTSDIR"] + "{techname}.fastq.foreignBC.stats.{ext}", techname=TECHNAMES, ext=config["PLOTFORMATS"]), #foreign barcode reads plots
 		expand ("mappings/" + "hiSeq_{capDesign}.bam", capDesign=CAPDESIGNS),  # mapped short reads
 		expand ("mappings/" + "mergeSizeFracBams/{techname}_{capDesign}_{barcodes}.merged.bam", capDesignBarcodeMatch, techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES),  # mapped reads
		expand ("mappings/" + "clusterPolyAsites/{techname}_{capDesign}_{barcodes}.polyAsites.clusters.bed", capDesignBarcodeMatch, techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES),
		expand("mappings/" + "getPolyAreadsList/{techname}_{capDesign}_{barcodes}.polyAreads.list", capDesignBarcodeMatch, techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES),
		 expand("mappings/" + "makePolyABigWigs/{techname}_{capDesign}_{barcodes}.polyAsitesNoErcc.{strand}.bw", capDesignBarcodeMatch, techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES, strand=config["STRANDS"]),

# 		expand("mappings/" + "qc/{techname}_{capDesign}_{barcodes}.merged.bam.dupl.txt", techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES),
# 		expand ("mappings/" + "readBedToGff/{techname}_{capDesign}_{barcodes}.merged.gff", techname=TECHNAMES, capDesign=CAPDESIGNS, barcodes=BARCODES),
# 		expand("mappings/qualimap_reports/" + "{techname}_{capDesign}.merged2/genome_results.txt", techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS),
 		expand(config["PLOTSDIR"] + "{techname}.ambiguousBarcodes.reads.stats.{ext}", techname=TECHNAMES, ext=config["PLOTFORMATS"]), # ambiguous barcodes plots
 		expand(config["PLOTSDIR"] + "{techname}_{capDesign}.adapters.location.stats.{ext}",techname=TECHNAMES, capDesign=CAPDESIGNS, ext=config["PLOTFORMATS"]), #location of adapters over reads
 		expand(config["PLOTSDIR"] + "{techname}.chimeric.reads.stats.{ext}", techname=TECHNAMES, ext=config["PLOTFORMATS"]), # stats on chimeric reads
 		expand(config["PLOTSDIR"] + "{techname}.finalDemul.reads.stats.{ext}", techname=TECHNAMES, ext=config["PLOTFORMATS"]), #final demultiplexing stats
 		expand(config["DEMULTIPLEX_DIR"] + "qc/{techname}_{capDesign}_{sizeFrac}.demul.QC1.txt",filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS), # QC on demultiplexing (checks that there is only one barcode assigned per read
 		expand(config["DEMULTIPLEX_DIR"] + "qc/{techname}_{capDesign}_{sizeFrac}.demul.QC2.txt", filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS), # QC on demultiplexing (checks that there is only one barcode assigned per read
 		expand(config["DEMULTIPLEX_DIR"] + "qc/{techname}_{capDesign}_{sizeFrac}.{barcodes}.demul.QC3.txt",filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS, barcodes=BARCODES),
 		expand(config["PLOTSDIR"] + "{techname}.demultiplexing.perSample.stats.{ext}", techname=TECHNAMES, ext=config["PLOTFORMATS"]),
# 		expand(config["PLOTSDIR"] + "{techname}.mapping.perSample.perFraction.stats.{ext}", techname=TECHNAMES, ext=config["PLOTFORMATS"])


