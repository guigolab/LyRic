#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

# test if there is at least two arguments: if not, return an error
if (length(args) != 2) {
  stop("Exactly two arguments must be given\n\t1: input sample annotation TSV;\n\t2: output HTML", call.=FALSE)
}

inputSampleAnnotation=args[1]
outputHtml=args[2]

print(paste("Input: ", inputSampleAnnotation))
print(paste("Output: ", outputHtml))

library(formattable)
library(dplyr)
library(tidyr)
library(DT)
#library(kableExtra)

concatMetadata <- function(df){
	return(mutate(df, sample_name=paste(seqTech, capDesign, sizeFrac, tissue, sep='_'), seqTech=NULL, capDesign=NULL, sizeFrac=NULL, tissue=NULL, correctionLevel=NULL))
}

dat <- read.table("Rinput/all.readlength.summary.tsv", header=TRUE, sep='\t')
dat <- concatMetadata(dat)

sampleAnnot <- read.table(inputSampleAnnotation, header=TRUE, sep='\t')

fastqTimestamps <- read.table("Rinput/all.fastq.timestamps.tsv", header=TRUE, sep='\t')
mappingStats <- read.table("Rinput/all.basic.mapping.stats.tsv", header=TRUE, sep='\t')
mappingStats <- concatMetadata(mappingStats)
mappingStats <- select(mappingStats, -totalReads)

# mappedLength <- read.table("Rinput/all.mappedReadlength.summary.tsv", header=TRUE, sep='\t')
# mappedLength <- concatMetadata(mappedLength)
# mappedLength <- select(mappedLength, -n,-median,  -mean)
# mappedLength <- rename(mappedLength, maxM = max)

mergedStats <- read.table("Rinput/all.min2reads.merged.stats.tsv", header=TRUE, sep='\t')
mergedStats <- spread(mergedStats, category, count)
mergedStats <- concatMetadata(mergedStats)

ntCoverageStats <- read.table("Rinput/all.tmerge.min2reads.endSupport:all.vs.ntCoverageByGenomePartition.stats.tsv", header=TRUE, sep='\t')
ntCoverageStats <- concatMetadata(ntCoverageStats)
ntCoverageStats <- filter(ntCoverageStats, splicingStatus=="all")
ntCoverageStats <- select(ntCoverageStats, -splicingStatus)
ntCoverageSummary <- ntCoverageStats %>% group_by(sample_name) %>% summarise(sum=sum(nt_coverage_count))
ntCoverageSummary <- rename(ntCoverageSummary, 'genomic nts\ncovered'= sum)
ntCoverageStatsCounts <- select(ntCoverageStats, -nt_coverage_percent)
ntCoverageStatsCounts <- spread(ntCoverageStatsCounts, partition, nt_coverage_count)
ntCoverageStatsPercent <- select(ntCoverageStats, -nt_coverage_count)
ntCoverageStatsPercent <- spread(ntCoverageStatsPercent, partition, nt_coverage_percent)

ntCoverageStats <- inner_join(ntCoverageStatsCounts, ntCoverageStatsPercent, by='sample_name')
ntCoverageStats <- relocate(ntCoverageStats, intergenic.x, .after=sample_name)
ntCoverageStats <- relocate(ntCoverageStats, intergenic.y, .after=intergenic.x)
ntCoverageStats <- relocate(ntCoverageStats, CDS.x, .after=intergenic.y)
ntCoverageStats <- relocate(ntCoverageStats, CDS.y, .after=CDS.x)
ntCoverageStats <- relocate(ntCoverageStats, UTR.x, .after=CDS.y)
ntCoverageStats <- relocate(ntCoverageStats, UTR.y, .after=UTR.x)
ntCoverageStats <- relocate(ntCoverageStats, exonOfNCT.x, .after=UTR.y)
ntCoverageStats <- relocate(ntCoverageStats, exonOfNCT.y, .after=exonOfNCT.x)
ntCoverageStats <- relocate(ntCoverageStats, exonOfPseudo.x, .after=exonOfNCT.y)
ntCoverageStats <- relocate(ntCoverageStats, exonOfPseudo.y, .after=exonOfPseudo.x)
ntCoverageStats <- relocate(ntCoverageStats, intron.x, .after=exonOfPseudo.y)
ntCoverageStats <- relocate(ntCoverageStats, intron.y, .after=intron.x)



dat <- inner_join(sampleAnnot, dat, by='sample_name')
dat <- inner_join(fastqTimestamps, dat, by='sample_name')
#dat <- inner_join(dat, mappedLength, by='sample_name')
dat <- inner_join(dat, mappingStats, by='sample_name')
dat <- inner_join(dat, mergedStats, by='sample_name')
dat <- inner_join(dat, ntCoverageSummary, by='sample_name')
dat <- inner_join(dat, ntCoverageStats, by='sample_name')

dat <- mutate(dat, '%\nHCGM\nreads' =  HCGMreads / mappedReads)
dat <- relocate(dat, '%\nHCGM\nreads', .after=HCGMreads)
dat <- relocate(dat, date_sequenced, .after=sample_name)
dat <- relocate(dat, seqPlatform, .after=seqCenter)
dat <- relocate(dat, species, .after=seqPlatform)
dat <- relocate(dat, tissue, .after=species)
dat <- relocate(dat, biosample_id, .after=sample_name)
dat <- relocate(dat, reverse_transcriptase, .after=libraryPrep)
dat <- mutate(dat, 'merge\nrate' = mergedTMs / HCGMreads)
dat <- relocate(dat, 'merge\nrate', .after=mergedTMs)
dat$n <- comma(dat$n, digits=0)
dat$median <- comma(dat$median, digits=0)
dat$mean <- comma(dat$mean, digits=0)
dat$max <- comma(dat$max, digits=0)
#dat$maxM <- comma(dat$maxM, digits=0)
dat$FASTQ_modified <- as.Date(dat$FASTQ_modified, format="%Y-%m-%d")
dat$date_sequenced <- as.Date(as.character(dat$date_sequenced), format="%Y%m%d")
dat$percentMappedReads <- percent(dat$percentMappedReads, digits=0)
dat$'%\nHCGM\nreads' <- percent(dat$'%\nHCGM\nreads', digits=0)
dat$HCGMreads <- comma(dat$HCGMreads, digits=0)
dat$mergedTMs <- comma(dat$mergedTMs, digits=0)
dat$'merge\nrate' <- percent(dat$'merge\nrate', digits=2)
dat$'genomic nts\ncovered' <- comma(dat$'genomic nts\ncovered',digits=0)

dat$intergenic.x <- comma(dat$intergenic.x, digits=0)
dat$CDS.x <- comma(dat$CDS.x, digits=0)
dat$UTR.x <- comma(dat$UTR.x, digits=0)
dat$exonOfNCT.x <- comma(dat$exonOfNCT.x, digits=0)
dat$exonOfPseudo.x <- comma(dat$exonOfPseudo.x, digits=0)
dat$intron.x <- comma(dat$intron.x, digits=0)

dat$intergenic.y <- percent(dat$intergenic.y, digits=0)
dat$CDS.y <- percent(dat$CDS.y, digits=0)
dat$UTR.y <- percent(dat$UTR.y, digits=0)
dat$exonOfNCT.y <- percent(dat$exonOfNCT.y, digits=0)
dat$exonOfPseudo.y <- percent(dat$exonOfPseudo.y, digits=0)
dat$intron.y <- percent(dat$intron.y, digits=0)


#dat <- select(dat, -mappedReads)
dat <- rename(dat, 'capped\nspike-ins?' = cappedSpikeIns)
dat <- rename(dat, 'FASTQ\nlast\nmodified' = FASTQ_modified)
dat <- rename(dat, 'date\nsequenced' = date_sequenced)
dat <- rename(dat, 'reverse\ntranscriptase' = reverse_transcriptase)
dat <- rename(dat, '# reads' = n)
dat <- rename(dat, 'median\nread\nlength' = median)
dat <- rename(dat, 'mean\nread\nlength' = mean)
dat <- rename(dat, 'max\nread\nlength\n(FASTQ)' = max)
#dat <- rename(dat, 'max\nmapped\nread\nlength\n(BAM)' = maxM)
dat <- rename(dat, '%\nmapped\nreads'=percentMappedReads)
dat <- rename(dat, 'seq\nplatform' = seqPlatform)
dat <- rename(dat, 'seq\ncenter' = seqCenter)
dat <- rename(dat, 'library\nprep' = libraryPrep)
dat <- rename(dat, 'size\nfraction' = sizeFrac)
dat <- rename(dat, 'cell\nfraction' = cellFrac)
dat <- rename(dat, 'bio\nrep' = bioRep)
dat <- rename(dat, 'tech\nrep' = techRep)
dat <- rename(dat, 'capture\ndesign' = captureDesign)
dat <- rename(dat, '#\nHCGM\nreads' = HCGMreads)
dat <- rename(dat, '#\nmerged\nTMs\n(min.\n2 reads)' =mergedTMs)
dat <- rename(dat, '# TM nts over\nintergenic\nregions' = intergenic.x, '% TM nts over\nintergenic\nregions' = intergenic.y, '# TM nts over\nCDS\nregions'=CDS.x, '% TM nts over\nCDS\nregions'=CDS.y, '# TM nts over\nUTR\nregions'=UTR.x, '% TM nts over\nUTR\nregions'=UTR.y, '# TM nts over\nexons of\nnoncoding\ntranscripts' =exonOfNCT.x, '% TM nts over\nexons of\nnoncoding\ntranscripts' =exonOfNCT.y, '# TM nts over\nexons of\npseudogenes' = exonOfPseudo.x, '% TM nts over\nexons of\npseudogenes' = exonOfPseudo.y, '# TM nts over\nintrons' = intron.x, '% TM nts over\nintrons' = intron.y)


tb <- formattable(dat, 
#	align = c("l",rep("r", NCOL(dat) -1)), 
	list('# reads' = color_tile("#def7e9", "#45ba78"),
	'median\nread\nlength' = color_tile("#e6f7ff", "#0088cc"),
	'mean\nread\nlength' = color_tile("#e6f7ff", "#0088cc"),
	'max\nread\nlength\n(FASTQ)' = color_tile("#e6f7ff", "#0088cc"),
	'%\nmapped\nreads' = color_tile("#def7e9", "#45ba78"),
#	'max\nmapped\nread\nlength\n(BAM)' = color_tile("#e6f7ff", "#0088cc"),
	'#\nHCGM\nreads' = color_tile("#def7e9", "#45ba78"),
	'%\nHCGM\nreads' = color_tile("#def7e9", "#45ba78"),
	'#\nmerged\nTMs\n(min.\n2 reads)' = color_tile('#ffe6cc','#ff8c1a'),
	'merge\nrate' = color_tile('#ffe6cc','#ff8c1a'),
	'genomic nts\ncovered' = color_tile("#f9ecf2", "#c6538c"),
	'# TM nts over\nintergenic\nregions' = color_tile("#f9ecf2", "#c6538c"),
	'% TM nts over\nintergenic\nregions' = color_tile("#f9ecf2", "#c6538c"),
	'# TM nts over\nCDS\nregions'=color_tile("#f9ecf2", "#c6538c"),
	'% TM nts over\nCDS\nregions'=color_tile("#f9ecf2", "#c6538c"),
	'# TM nts over\nUTR\nregions'=color_tile("#f9ecf2", "#c6538c"),
	'% TM nts over\nUTR\nregions'=color_tile("#f9ecf2", "#c6538c"),
	'# TM nts over\nexons of\nnoncoding\ntranscripts' =color_tile("#f9ecf2", "#c6538c"),
	'% TM nts over\nexons of\nnoncoding\ntranscripts' =color_tile("#f9ecf2", "#c6538c"),
	'# TM nts over\nexons of\npseudogenes' = color_tile("#f9ecf2", "#c6538c"),
	'% TM nts over\nexons of\npseudogenes' = color_tile("#f9ecf2", "#c6538c"),
	'# TM nts over\nintrons' = color_tile("#f9ecf2", "#c6538c"),
	'% TM nts over\nintrons' = color_tile("#f9ecf2", "#c6538c")


#	'species' = formatter("span",
#		style = x ~ style(display = "block", color = "black", "border-radius" = "0px", 'background-color' =  ifelse(x == 'Human', '#ffcce0',  ifelse( x== 'Mouse', '#ffff80', NA)))),
#	'cell\nfraction' = formatter("span",
#		style = x ~ style(display = "block", color = "black", "border-radius" = "0px", 'background-color' =  ifelse(x == 'T', '#d9d9d9', ifelse( x=='C', '#ffb3ff', ifelse( x== 'N', '#ccebff', NA))))),
#	'seq\nplatform' = formatter("span",
#		style = x ~ style(display = "block", color = "black", "border-radius" = "0px", 'background-color' =  ifelse(x == 'ONT', '#8dd3c7', ifelse( x=='pacBioSI', '#ffffb3', ifelse( x== 'pacBioSII', '#baac91', NA))))),
#	'size\nfraction' = formatter("span",
#		style = x ~ style(display = "block", color = "black", "border-radius" = "0px", 'background-color' =  ifelse(x == '1+', '#b370f9', NA))),
#	'library\nprep' = formatter("span",
#		style = x ~ style(display = "block", color = "black", "border-radius" = "0px", 'background-color' =  ifelse(x == 'CapTrap', '#ff99ff', ifelse( x=='SMARTer', '#bebada', ifelse( x== 'Teloprime', '#fdbcb5', ifelse(x=='directRNA', '#b3d0e5', ifelse(x== 'Rt', '#d9d9d9', ifelse(x=='PcrOnt', '#9999ff', NA))))))))

)
) 
#tb %>% as.htmlwidget() %>% htmlwidgets::saveWidget(file="./html/summary_table_test.html")
#kbl(tb, escape=TRUE,format='latex') %>%  kable_styling("hover", full_width = F)

#htmlString <- tb %>% select('# reads') %>% rename (reads = '# reads') %>% mutate(reads = color_tile("white", "orange")(reads)) %>% kable('html', escape=F)
#write(htmlString, file='html/test.html')


tbDt <- as.datatable(tb[order(tb$tissue), ], extensions = c('Buttons','ColReorder','FixedColumns', 'FixedHeader', 'KeyTable','RowGroup', 'Select', 'SearchPanes'), 
filter='top', 
rownames=FALSE, 
selection='none',
options= list(
	dom = 'PBlfrtip',
	buttons = c('csv','colvis'),
	buttons = list(list(extend='colvis', column = c(2,3,4))),
	searchPanes = list(threshold = 0),
	columnDefs = list(list(searchPanes = list(show=TRUE), targets = c(0,5,6,7,8,9,10,11,12,13,16,17))),
	colReorder=TRUE, 
#	fixedColumns = TRUE,
	fixedHeader = TRUE, 
	pageLength =400,
	#scrollX = TRUE,
	keys = TRUE,
	rowId =0,
#	rowGroup = list(dataSrc = c(7)),
	initComplete = JS(
    "function(settings, json) {",
    "$('body').css({'font-family': 'Helvetica'});",
    "}"
  )
	)
	)

#tbDt <- tbDt %>% formatStyle('species', backgroundColor = styleEqual (c('Human', 'Mouse'), c('#ffcce0', '#ffff80')))
#tbDt <- tbDt %>% formatStyle('seq_platform', backgroundColor = styleEqual (c('ONT', 'pacBioSI', 'pacBioSII'), c('#8dd3c7', '#ffffb3', '#baac91')))
#tbDt <- tbDt %>% formatStyle('library\nprep', backgroundColor = styleEqual (c('CapTrap', 'SMARTer', 'Teloprime', 'directRNA', 'Rt', 'PcrOnt'), c('#ff99ff', '#bebada', '#fdbcb5', '#b3d0e5', '#d9d9d9', '#9999ff')))


# unique(dat$seqPlatform)
DT::saveWidget(tbDt, paste0(getwd(), "/", outputHtml), title='GENCODE - CRG Dashboard of Experiments')
