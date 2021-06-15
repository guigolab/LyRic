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

hissStats <- read.table("Rinput/all.HiSS.stats.tsv", header=TRUE, sep='\t')
hissStats <- concatMetadata(hissStats)
hissStats <- spread(hissStats, category, count)
hissStats <- rename(hissStats, HCGMspliced='HCGM-spliced', HCGMmono='HCGM-mono')
hissStats <- mutate(hissStats, percentHcgmSpliced=HCGMspliced/(HCGMmono+HCGMspliced))
hissStats <- select(hissStats, sample_name, percentHcgmSpliced)
hissStats <- rename(hissStats, '%\nspliced\nHCGMs' = percentHcgmSpliced)

mergedStats <- read.table("Rinput/all.min2reads.merged.stats.tsv", header=TRUE, sep='\t')
mergedStats <- spread(mergedStats, category, count)
mergedStats <- concatMetadata(mergedStats)

tmLengthStats <- read.table('Rinput/all.min2reads.matureRNALengthSummary.stats.tsv', header=TRUE, sep='\t')
tmLengthStats <- concatMetadata(tmLengthStats)
tmLengthStats <- filter(tmLengthStats, category=='CLS_TMs')
tmMedLengthStats <- spread(select(tmLengthStats, -max), category, med)
tmMedLengthStats <- rename(tmMedLengthStats, 'TM\nmedian\nlength'=CLS_TMs)
tmMaxLengthStats <- spread(select(tmLengthStats, -med), category, max)
tmMaxLengthStats <- rename(tmMaxLengthStats, 'TM\nmax\nlength'=CLS_TMs)
tmLengthStats <- inner_join(tmMedLengthStats, tmMaxLengthStats, by='sample_name')

sirvAccuracyStats <- read.table('Rinput/all.HiSS.tmerge.min2reads.vs.SIRVs.stats.tsv', header=TRUE, sep='\t')
sirvAccuracyStats <- concatMetadata(sirvAccuracyStats)
sirvAccuracyStats <- filter(sirvAccuracyStats, level=='Transcript')
sirvAccuracyStats <- spread(select(sirvAccuracyStats, -level), metric, value)
sirvAccuracyStats <- rename(sirvAccuracyStats, 'SIRV\nTx\nSn' = 'Sn')
sirvAccuracyStats <- rename(sirvAccuracyStats, 'SIRV\nTx\nPr' = 'Pr')

cagePolyASupportStats <- read.table('Rinput/all.min2reads.splicing_status:all.cagePolyASupport.stats.tsv', header=TRUE, sep='\t')
cagePolyASupportStats <- select(cagePolyASupportStats, -count)
cagePolyASupportStats <- spread(cagePolyASupportStats, category, percent)
cagePolyASupportStats <- select(cagePolyASupportStats, -cageOnly, -noCageNoPolyA, -polyAOnly)
cagePolyASupportStats <- concatMetadata(cagePolyASupportStats)
cagePolyASupportStats <- rename(cagePolyASupportStats, '%\nFL TMs' = cageAndPolyA)

novelLociStats <- read.table('Rinput/all.tmerge.min2reads.endSupport:all.novelLoci.stats.tsv', header=TRUE, sep='\t')
novelLociStats <- concatMetadata(novelLociStats)
novelLociStats <- select(novelLociStats, -percent)
novelLociStats <- spread(novelLociStats, category, count)
novelLociStats <- rename(novelLociStats, 'novelIntergenicLoci'=intergenic, 'novelIntronicLoci'=intronic)
novelLociStats <- mutate(novelLociStats, 'novelLoci'=novelIntergenicLoci+novelIntronicLoci)
novelLociStats <- mutate(novelLociStats, percentIntergenicNovelLoci=novelIntergenicLoci/novelLoci)
novelLociStats <- select(novelLociStats, -novelIntergenicLoci, -novelIntronicLoci)

novelFlLociStats <- read.table('Rinput/all.tmerge.min2reads.endSupport:cagePolyASupported.novelLoci.stats.tsv', header=TRUE, sep='\t')
novelFlLociStats <- concatMetadata(novelFlLociStats)
novelFlLociStats <- select(novelFlLociStats, -percent)
novelFlLociStats <- spread(novelFlLociStats, category, count)
novelFlLociStats <- rename(novelFlLociStats, 'novelIntergenicFlLoci'=intergenic, 'novelIntronicFlLoci'=intronic)
novelFlLociStats <- mutate(novelFlLociStats, 'novelFlLoci'=novelIntergenicFlLoci+novelIntronicFlLoci)
novelFlLociStats <- select(novelFlLociStats, -novelIntergenicFlLoci, -novelIntronicFlLoci)

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
dat <- inner_join(dat, hissStats, by='sample_name')
dat <- inner_join(dat, mergedStats, by='sample_name')
dat <- inner_join(dat, tmLengthStats, by='sample_name')
dat <- inner_join(dat, sirvAccuracyStats, by='sample_name')
dat <- inner_join(dat, cagePolyASupportStats, by='sample_name')
dat <- inner_join(dat, ntCoverageSummary, by='sample_name')
dat <- inner_join(dat, ntCoverageStats, by='sample_name')
dat <- inner_join(dat, novelLociStats, by='sample_name')
dat <- inner_join(dat, novelFlLociStats, by='sample_name')

dat <- mutate(dat, '%\nHCGM\nreads' =  HCGMreads / mappedReads)
dat <- relocate(dat, '%\nHCGM\nreads', .after=HCGMreads)
dat <- relocate(dat, '%\nspliced\nHCGMs', .after='%\nHCGM\nreads')
dat <- relocate(dat, date_sequenced, .after=sample_name)
dat <- relocate(dat, seqPlatform, .after=seqCenter)
dat <- relocate(dat, species, .after=seqPlatform)
dat <- relocate(dat, tissue, .after=species)
dat <- relocate(dat, biosample_id, .after=sample_name)
dat <- relocate(dat, reverse_transcriptase, .after=libraryPrep)
dat <- mutate(dat, 'merge\nrate' = mergedTMs / HCGMreads)
dat <- mutate(dat, '# novel\nloci\nper\nmillion\nmapped\nreads'= novelLoci / (mappedReads/1000000) )
dat <- relocate(dat, '# novel\nloci\nper\nmillion\nmapped\nreads', .after=novelLoci)
dat <- mutate(dat, '% novel\nloci\nFL'= novelFlLoci/novelLoci)
dat <- rename(dat, '# novel\nloci'=novelLoci, '% novel\nloci\nintergenic' = percentIntergenicNovelLoci)
dat <- select(dat, -mappedReads, -novelFlLoci)
dat$'# novel\nloci' <- comma(dat$'# novel\nloci', digits=0)
dat$'# novel\nloci\nper\nmillion\nmapped\nreads' <- comma(dat$'# novel\nloci\nper\nmillion\nmapped\nreads', digits=0)
dat$'% novel\nloci\nintergenic' <- percent(dat$'% novel\nloci\nintergenic', digits=0)
dat$'% novel\nloci\nFL' <- percent(dat$'% novel\nloci\nFL', digits=0)

dat <- relocate(dat, 'merge\nrate', .after=mergedTMs)
dat$n <- comma(dat$n, digits=0)
dat$median <- comma(dat$median, digits=0)
dat$mean <- comma(dat$mean, digits=0)
dat$max <- comma(dat$max, digits=0)
dat$FASTQ_modified <- as.Date(dat$FASTQ_modified, format="%Y-%m-%d")
dat$date_sequenced <- as.Date(as.character(dat$date_sequenced), format="%Y%m%d")
dat$percentMappedReads <- percent(dat$percentMappedReads, digits=0)
dat$'%\nHCGM\nreads' <- percent(dat$'%\nHCGM\nreads', digits=0)
dat$'%\nspliced\nHCGMs' <- percent(dat$'%\nspliced\nHCGMs', digits=0)
dat$HCGMreads <- comma(dat$HCGMreads, digits=0)
dat$mergedTMs <- comma(dat$mergedTMs, digits=0)
dat$'merge\nrate' <- percent(dat$'merge\nrate', digits=2)
dat$'TM\nmedian\nlength' <- comma(dat$'TM\nmedian\nlength', digits=0)
dat$'TM\nmax\nlength' <- comma(dat$'TM\nmax\nlength', digits=0)
dat$'SIRV\nTx\nSn' <- percent(dat$'SIRV\nTx\nSn'/100, digits=0)
dat$'SIRV\nTx\nPr' <- percent(dat$'SIRV\nTx\nPr'/100, digits=0)
dat$'%\nFL TMs' <- percent(dat$'%\nFL TMs', digits=0)
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


dat <- rename(dat, 'capped\nspike-ins?' = cappedSpikeIns)
dat <- rename(dat, 'FASTQ\nlast\nmodified' = FASTQ_modified)
dat <- rename(dat, 'date\nsequenced' = date_sequenced)
dat <- rename(dat, 'reverse\ntranscriptase' = reverse_transcriptase)
dat <- rename(dat, '# reads' = n)
dat <- rename(dat, 'median\nread\nlength' = median)
dat <- rename(dat, 'mean\nread\nlength' = mean)
dat <- rename(dat, 'max\nread\nlength\n(FASTQ)' = max)
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
	list('# reads' = color_tile("#def7e9", "#45ba78"),
	'median\nread\nlength' = color_tile("#e6f7ff", "#0088cc"),
	'mean\nread\nlength' = color_tile("#e6f7ff", "#0088cc"),
	'max\nread\nlength\n(FASTQ)' = color_tile("#e6f7ff", "#0088cc"),
	'%\nmapped\nreads' = color_tile("#def7e9", "#45ba78"),
	'#\nHCGM\nreads' = color_tile("#def7e9", "#45ba78"),
	'%\nHCGM\nreads' = color_tile("#def7e9", "#45ba78"),
	'%\nspliced\nHCGMs' = color_tile("#def7e9", "#45ba78"), 
	'#\nmerged\nTMs\n(min.\n2 reads)' = color_tile('#ffe6cc','#ff8c1a'),
	'merge\nrate' = color_tile('#ffe6cc','#ff8c1a'),
	'TM\nmedian\nlength'= color_tile("#e6f7ff", "#0088cc"),
	'TM\nmax\nlength'= color_tile("#e6f7ff", "#0088cc"),
	'%\nFL TMs' = color_tile('#f8ecf8', '#d17ad1'),
	'SIRV\nTx\nSn'=color_tile('#ffeee6','#ff7733'),
	'SIRV\nTx\nPr'=color_tile('#ecf9f2', '#339966'),
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
	'% TM nts over\nintrons' = color_tile("#f9ecf2", "#c6538c"),
	'# novel\nloci' = color_tile('#e6ffff', '#00cccc'),
	'# novel\nloci\nper\nmillion\nmapped\nreads' = color_tile('#e6ffff', '#00cccc'),
	'% novel\nloci\nintergenic' = color_tile('#e6ffff', '#00cccc'),
	'% novel\nloci\nFL' = color_tile('#e6ffff', '#00cccc')


)
) 
#tb %>% as.htmlwidget() %>% htmlwidgets::saveWidget(file="./html/summary_table_test.html")
#kbl(tb, escape=TRUE,format='latex') %>%  kable_styling("hover", full_width = F)

#htmlString <- tb %>% select('# reads') %>% rename (reads = '# reads') %>% mutate(reads = color_tile("white", "orange")(reads)) %>% kable('html', escape=F)
#write(htmlString, file='html/test.html')


tbDt <- as.datatable(tb[order(tb$tissue), ], extensions = c('Buttons','ColReorder','FixedColumns', 'FixedHeader', 'KeyTable','RowGroup', 'Select', 'SearchPanes'), 
#filter='top', 
rownames=FALSE, 
selection='none',
options= list(
	dom = 'PBlfrtip',
	buttons = c('csv','colvis'),
	buttons = list(list(extend='colvis', column = c(2,3,4))),
	searchPanes = list(threshold = 0),
	columnDefs = list(list(searchPanes = list(show=TRUE), targets = c(0,5,6,7,8,9,10,11,12,13,16,17))),
	colReorder=TRUE, 
	fixedColumns = list(leftColumns=1),
	fixedHeader = TRUE, 
	pageLength =400,
#	scrollX = TRUE,
	autowidth=TRUE,
	keys = TRUE,
	rowId =0,
	searching = FALSE,
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
