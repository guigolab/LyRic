
def authorizeComb(comb):
	global AUTHORIZEDCOMBINATIONS
	AUTHORIZEDCOMBINATIONS.append(
		(("techname", comb[0]), ("capDesign", comb[1])))
	AUTHORIZEDCOMBINATIONS.append(
		(("techname", comb[0]), ("capDesign", comb[1]), ("sizeFrac", comb[2]), ("barcodes", comb[3])))


def returnPlotFilenames(basename):
	plotsList = []
	for comb in itertools.product(["legendOnly"], config["PLOTFORMATS"]):
		plotsList.append(basename + "." + comb[0] + "." + comb[1])
	for comb in itertools.product(["xy", "yx"], ["wLegend", "woLegend"], config["PLOTFORMATS"]):
		plotsList.append(basename + "." +
						 comb[0] + "." + comb[1] + "." + comb[2])
	return plotsList


def filtered_product(*args):
	# return combinations of wildcards that correspond to existing combinations in input (i.e. corresponding input FASTQ exists)
	# args[0]: techname
	# args[1]: capDesign
	# args[2]: sizeFrac (optional)
	# args[3]: barcodes (optional)
	# args[>3] are ignored by filter but returned if args[0:3] pass

	found = False
	#print("NEW ARGS: ")
	# pprint(args)
	if len(args) < 2:
		quit("Error in function filtered_product (wrong number of arguments, shoud be >1).")
	elif len(args) == 2:
		for wc_comb in itertools.product(*args):
			if wc_comb in AUTHORIZEDCOMBINATIONS:
				found = True
				yield(wc_comb)
	else:
		for wc_comb in itertools.product(*args):
			if wc_comb[0:4] in AUTHORIZEDCOMBINATIONS:
				found = True
				yield(wc_comb)
			# else:
			# 	print ("NOT FOUND")
			# 	pprint(wc_comb)

#	if not found:
#		return None
		#pprint("Error in function filtered_product. Args were:", sys.stderr)
		#pprint((args), sys.stderr)
		#quit("Error. Could not yield any input file.")


def multi_figures(capDesign, sizeFrac, barcodes, techname, splicing_status=None):
	figure_settings = dict()
	figure_settings['technameFilterString'] = ''
	figure_settings['hideXaxisLabels'] = "theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +"
	removeFacetLabels = ''
	if techname not in ('byTech') and capDesign not in ('byCapDesign') and barcodes not in ('byTissue'):
		removeFacetLabels = """
plotFacetXy <- parse(text =paste(plotFacetXy, \\" + theme(strip.text.x = element_blank(), strip.text.y = element_blank())\\"))
plotFacetYx <- parse(text=paste(plotFacetYx, \\" + theme(strip.text.x = element_blank(), strip.text.y = element_blank())\\"))
"""
	figure_settings['splicingStatusFilterString'] = ''
	if splicing_status != None:
		if splicing_status == 'spliced':
			figure_settings['splicingStatusFilterString'] = "dat <- subset(dat, spliced==1)\n"
		elif splicing_status == 'unspliced':
			figure_settings['splicingStatusFilterString'] = "dat <- subset(dat, spliced==0)\n"
		elif splicing_status == 'all':
			figure_settings['splicingStatusFilterString'] = ''
		else:
			quit("Error. Invalid value for splicing_status in function multi_figures.")
	if techname not in ('byTech'):
		figure_settings['technameFilterString'] = "dat <- subset(dat, seqTech=='" + \
			techname + "')\n"
	figure_settings['substSeqTechString'] = """
dat\$seqTech <- gsub(':$', '', dat\$seqTech)
dat\$seqTech <- gsub(':', '\\n', dat\$seqTech)
"""
	figure_settings['substTissueString'] = "dat\$tissue <- gsub('" + \
		capDesign + "_', '', dat\$tissue)\n"
	figure_settings['capDesignFilterString'] = ''
	figure_settings['sizeFracFilterString'] = ''
	figure_settings['tissueFilterString'] = ''
	figure_settings['graphDimensions'] = """

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
	figure_settings['facetPlotSetup'] = f"""
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
		figure_settings['capDesignFilterString'] += "dat <- subset(dat, capDesign=='" + \
			capDesign + "')\n"

	figure_settings['sizeFracFilterString'] += "dat <- subset(dat, sizeFrac=='" + \
		sizeFrac + "')\n"
	if barcodes not in ('byTissue'):
		figure_settings['tissueFilterString'] += "dat <- subset(dat, tissue=='" + \
			barcodes + "')\n"
	return(figure_settings)


def trackHubSubGroupString(techname, capDesign, sizeFrac, barcodes, minReadSupport):
	techname = (("techname", techname),)
	capDesign = (("capDesign", capDesign),)
	sizeFracs = []
	for sizeF in sizeFrac:
		sizeFracs.append(("sizeFrac", sizeF))
	barCodes = []
	for barC in barcodes:
		barCodes.append(("barcodes", barC))
	minRS = []
	for minReadSupport in minReadSupport:
		minRS.append(("minReadSupport", minReadSupport))
	returnSubGroup1String = "subGroup1 sample Sample"
	returnSubGroup2String = "subGroup2 sizeFraction Size_fraction"
	returnSubGroup3String = "subGroup3 minReadSupport Min_read_support_per_TM"

	returnSubGroup1StringData = []
	returnSubGroup2StringData = []
	returnSubGroup3StringData = []
	for wc_comb in itertools.product(techname, capDesign, sizeFracs, barCodes, minRS):
		if wc_comb in AUTHORIZEDCOMBINATIONS:
			returnSubGroup1StringData.append(
				wc_comb[3][1] + "=" + wc_comb[3][1])
			returnSubGroup2StringData.append(
				wc_comb[2][1] + "=" + wc_comb[2][1])
			returnSubGroup3StringData.append(
				wc_comb[4][1] + "=" + wc_comb[4][1])

	returnSubGroup1StringData = set(returnSubGroup1StringData)
	returnSubGroup2StringData = set(returnSubGroup2StringData)
	returnSubGroup3StringData = set(returnSubGroup3StringData)

	return(returnSubGroup1String + " " + " ".join(sorted(returnSubGroup1StringData)) + "\n" + returnSubGroup2String + " " + " ".join(sorted(returnSubGroup2StringData)) + "\n" + returnSubGroup3String + " " + " ".join(sorted(returnSubGroup3StringData)))