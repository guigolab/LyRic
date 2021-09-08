#!/usr/bin/env Rscript

library(formattable)
library(tidyverse)
library(DT)

args = commandArgs(trailingOnly=TRUE)

# test if there is at least two arguments: if not, return an error
if (length(args) != 15) {
  stop("Exactly 15 arguments must be given\n", call.=FALSE)
}

outputHtml=args[1]
inputSampleAnnotation=args[2]
allFastqTimeStamps=args[3]
allReadLengths=args[4]
allBasicMappingStats=args[5]
allHissStats=args[6]
allMergedStats=args[7]
allMatureRnaLengthStats=args[8]
allTmergeVsSirvStats=args[9]
allCagePolyASupportStats=args[10]
allNovelLociStats=args[11]
allNovelFlLociStats=args[12]
allNovelLociQcStats=args[13]
allNovelFlLociQcStats=args[14]
allNtCoverageStats=args[15]

print(paste("Output: ", outputHtml))

concatMetadata <- function(df){
	head(df)
	return(mutate(df, sample_name=paste(seqTech, capDesign, sizeFrac, sampleRep, sep='_'), seqTech=NULL, capDesign=NULL, sizeFrac=NULL, sampleRep=NULL))
}

dat <- read.table(allReadLengths, header=TRUE, sep='\t')
dat <- concatMetadata(dat)

sampleAnnot <- read.table(inputSampleAnnotation, header=TRUE, sep='\t')

fastqTimestamps <- read.table(allFastqTimeStamps, header=TRUE, sep='\t')
mappingStats <- read.table(allBasicMappingStats, header=TRUE, sep='\t')
mappingStats <- concatMetadata(mappingStats)
mappingStats <- select(mappingStats, -totalReads)

hissStats <- read.table(allHissStats, header=TRUE, sep='\t')
hissStats <- concatMetadata(hissStats)
hissStats <- spread(hissStats, category, count)
hissStats <- rename(hissStats, HCGMspliced='HCGM-spliced', HCGMmono='HCGM-mono')
hissStats <- mutate(hissStats, percentHcgmSpliced=HCGMspliced/(HCGMmono+HCGMspliced))
hissStats <- select(hissStats, sample_name, percentHcgmSpliced)
hissStats <- rename(hissStats, '% spliced HCGMs' = percentHcgmSpliced)

mergedStats <- read.table(allMergedStats, header=TRUE, sep='\t')
mergedStats <- spread(mergedStats, category, count)
mergedStats <- concatMetadata(mergedStats)

tmLengthStats <- read.table(allMatureRnaLengthStats, header=TRUE, sep='\t')
tmLengthStats <- concatMetadata(tmLengthStats)
tmLengthStats <- filter(tmLengthStats, category=='CLS_TMs')
tmMedLengthStats <- spread(select(tmLengthStats, -max), category, med)
tmMedLengthStats <- rename(tmMedLengthStats, 'TM median length'=CLS_TMs)
tmMaxLengthStats <- spread(select(tmLengthStats, -med), category, max)
tmMaxLengthStats <- rename(tmMaxLengthStats, 'TM max length'=CLS_TMs)
tmLengthStats <- inner_join(tmMedLengthStats, tmMaxLengthStats, by='sample_name')

sirvAccuracyStats <- read.table(allTmergeVsSirvStats, header=TRUE, sep='\t')
sirvAccuracyStats <- concatMetadata(sirvAccuracyStats)
sirvAccuracyStats <- filter(sirvAccuracyStats, level=='Transcript')
sirvAccuracyStats <- spread(select(sirvAccuracyStats, -level), metric, value)
sirvAccuracyStats <- rename(sirvAccuracyStats, 'SIRV Tx Sn' = 'Sn')
sirvAccuracyStats <- rename(sirvAccuracyStats, 'SIRV Tx Pr' = 'Pr')

cagePolyASupportStats <- read.table(allCagePolyASupportStats, header=TRUE, sep='\t')
cagePolyASupportStats <- select(cagePolyASupportStats, -count)
cagePolyASupportStats <- spread(cagePolyASupportStats, category, percent)
cagePolyASupportStats <- select(cagePolyASupportStats, -cageOnly, -noCageNoPolyA, -polyAOnly)
cagePolyASupportStats <- concatMetadata(cagePolyASupportStats)
cagePolyASupportStats <- rename(cagePolyASupportStats, '% FL TMs' = cageAndPolyA)

novelLociStats <- read.table(allNovelLociStats, header=TRUE, sep='\t')
novelLociStats <- concatMetadata(novelLociStats)
novelLociStats <- select(novelLociStats, -percent)
novelLociStats <- spread(novelLociStats, category, count)
novelLociStats <- rename(novelLociStats, 'novelIntergenicLoci'=intergenic, 'novelIntronicLoci'=intronic)
novelLociStats <- mutate(novelLociStats, 'novelLoci'=novelIntergenicLoci+novelIntronicLoci)
novelLociStats <- mutate(novelLociStats, percentIntergenicNovelLoci=novelIntergenicLoci/novelLoci)
novelLociStats <- select(novelLociStats, -novelIntergenicLoci, -novelIntronicLoci)

novelFlLociStats <- read.table(allNovelFlLociStats, header=TRUE, sep='\t')
novelFlLociStats <- concatMetadata(novelFlLociStats)
novelFlLociStats <- select(novelFlLociStats, -percent)
novelFlLociStats <- spread(novelFlLociStats, category, count)
novelFlLociStats <- rename(novelFlLociStats, 'novelIntergenicFlLoci'=intergenic, 'novelIntronicFlLoci'=intronic)
novelFlLociStats <- mutate(novelFlLociStats, 'novelFlLoci'=novelIntergenicFlLoci+novelIntronicFlLoci)
novelFlLociStats <- select(novelFlLociStats, -novelIntergenicFlLoci, -novelIntronicFlLoci)

novelLociQcStats <- read.table(allNovelLociQcStats, header=TRUE, sep='\t')
novelLociQcStats <- concatMetadata(novelLociQcStats)
novelLociQcStats <- mutate(novelLociQcStats, countSpliced=total-countMono, percentSpliced=1-percentMono)
novelLociQcStats <- select(novelLociQcStats, -total, -minLengthAllTms, -minLengthMonoTms, -minLengthSplicedTms, -countMono, -percentMono)
novelLociQcStats <- relocate(novelLociQcStats, percentSpliced, .after=countSpliced)
novelLociQcStats <- relocate(novelLociQcStats, percentRepeats, .after=countRepeats)

novelFlLociQcStats <- read.table(allNovelFlLociQcStats, header=TRUE, sep='\t')
novelFlLociQcStats <- concatMetadata(novelFlLociQcStats)
novelFlLociQcStats <- mutate(novelFlLociQcStats, countSpliced=total-countMono, percentSpliced=1-percentMono)
novelFlLociQcStats <- select(novelFlLociQcStats, -total, -minLengthAllTms, -minLengthMonoTms, -minLengthSplicedTms, -countMono, -percentMono)
novelFlLociQcStats <- relocate(novelFlLociQcStats, percentSpliced, .after=countSpliced)
novelFlLociQcStats <- relocate(novelFlLociQcStats, percentRepeats, .after=countRepeats)
novelFlLociQcStats <- rename(novelFlLociQcStats, countFlLociSpliced=countSpliced, '% FL spliced novel loci'=percentSpliced, '# FL novel loci on repeats'=countRepeats, '% FL novel loci on repeats'=percentRepeats, 'median length all TMs (FL novel loci)'=medianLengthAllTms, 'max length all TMs (FL novel loci)'=maxLengthAllTms, 'median length monoexonic TMs (FL novel loci)'=medianLengthMonoTms, 'max length monoexonic TMs (FL novel loci)'=maxLengthMonoTms, 'median length spliced TMs (FL novel loci)'=medianLengthSplicedTms, 'max length spliced TMs (FL novel loci)'=maxLengthSplicedTms)


ntCoverageStats <- read.table(allNtCoverageStats, header=TRUE, sep='\t')
ntCoverageStats <- concatMetadata(ntCoverageStats)
ntCoverageStats <- filter(ntCoverageStats, splicingStatus=="all")
ntCoverageStats <- select(ntCoverageStats, -splicingStatus)
ntCoverageSummary <- ntCoverageStats %>% group_by(sample_name) %>% summarise(sum=sum(nt_coverage_count))
ntCoverageSummary <- rename(ntCoverageSummary, 'genomic nts covered'= sum)
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
dat <- inner_join(dat, novelLociQcStats, by='sample_name')
dat <- inner_join(dat, novelFlLociQcStats, by='sample_name')

dat <- mutate(dat, '% HCGM reads' =  HCGMreads / mappedReads)
dat <- relocate(dat, '% HCGM reads', .after=HCGMreads)
dat <- relocate(dat, '% spliced HCGMs', .after='% HCGM reads')
dat <- relocate(dat, date_sequenced, .after=sample_name)
dat <- relocate(dat, seqPlatform, .after=seqCenter)
dat <- relocate(dat, species, .after=seqPlatform)
dat <- relocate(dat, tissue, .after=species)
dat <- relocate(dat, biosample_id, .after=sample_name)
dat <- relocate(dat, reverse_transcriptase, .after=libraryPrep)
dat <- mutate(dat, 'merge rate' = mergedTMs / HCGMreads)
dat <- mutate(dat, '# novel loci per million mapped reads'= novelLoci / (mappedReads/1000000) )
dat <- mutate(dat, '# FL spliced novel loci per million mapped reads' = countFlLociSpliced / (mappedReads/1000000) )
dat <- rename(dat, '# FL spliced novel loci'=countFlLociSpliced)
dat <- relocate(dat, '# novel loci per million mapped reads', .after=novelLoci)
dat <- relocate(dat, '# FL spliced novel loci per million mapped reads', .after='# FL spliced novel loci')

dat <- mutate(dat, '% novel loci FL'= novelFlLoci/novelLoci)
dat <- relocate(dat, '% novel loci FL', .after=novelFlLoci)
dat <- rename(dat, '# novel loci'=novelLoci, '% novel loci intergenic' = percentIntergenicNovelLoci, '# novel loci FL'=novelFlLoci, '# spliced novel loci'=countSpliced, '% spliced novel loci'=percentSpliced, '# novel loci on repeats'=countRepeats, '% novel loci on repeats'=percentRepeats, 'median length all TMs (novel loci)'=medianLengthAllTms, 'max length all TMs (novel loci)'=maxLengthAllTms, 'median length monoexonic TMs (novel loci)'=medianLengthMonoTms, 'max length monoexonic TMs (novel loci)'=maxLengthMonoTms, 'median length spliced TMs (novel loci)'=medianLengthSplicedTms, 'max length spliced TMs (novel loci)'=maxLengthSplicedTms)
dat <- relocate(dat, '# spliced novel loci', '% spliced novel loci', .after='% novel loci intergenic')
dat <- select(dat, -mappedReads)
dat$'# novel loci' <- comma(dat$'# novel loci', digits=0)
dat$'# novel loci per million mapped reads' <- comma(dat$'# novel loci per million mapped reads', digits=0)
dat$'% novel loci intergenic' <- percent(dat$'% novel loci intergenic', digits=0)
dat$'% novel loci FL' <- percent(dat$'% novel loci FL', digits=0)
dat$'# novel loci FL' <- comma(dat$'# novel loci FL', digits=0)
dat$'# spliced novel loci' <- comma(dat$'# spliced novel loci', digits=0)
dat$'% spliced novel loci' <- percent(dat$'% spliced novel loci', digits=0)
dat$'# novel loci on repeats' <- comma(dat$'# novel loci on repeats', digits=0)
dat$'% novel loci on repeats' <- percent(dat$'% novel loci on repeats', digits=0)
dat$'# FL spliced novel loci per million mapped reads'<- comma(dat$'# FL spliced novel loci per million mapped reads', digits=0)

dat$'median length all TMs (novel loci)' <- comma(dat$'median length all TMs (novel loci)', digits=0)
dat$'max length all TMs (novel loci)' <- comma(dat$'max length all TMs (novel loci)', digits=0)
dat$'median length monoexonic TMs (novel loci)' <- comma(dat$'median length monoexonic TMs (novel loci)', digits=0)
dat$'max length monoexonic TMs (novel loci)' <- comma(dat$'max length monoexonic TMs (novel loci)', digits=0)
dat$'median length spliced TMs (novel loci)' <- comma(dat$'median length spliced TMs (novel loci)', digits=0)
dat$'max length spliced TMs (novel loci)' <- comma(dat$'max length spliced TMs (novel loci)', digits=0)
dat$'# FL spliced novel loci' <- comma(dat$'# FL spliced novel loci', digits=0)
dat$'% FL spliced novel loci' <- percent(dat$ '% FL spliced novel loci', digits=0)
dat$'# FL novel loci on repeats' <- comma(dat$'# FL novel loci on repeats', digits=0)
dat$'% FL novel loci on repeats'  <- percent(dat$'% FL novel loci on repeats', digits=0)
dat$'median length all TMs (FL novel loci)' <- comma(dat$'median length all TMs (FL novel loci)', digits=0)
dat$'max length all TMs (FL novel loci)' <- comma(dat$'max length all TMs (FL novel loci)', digits=0)
dat$'median length monoexonic TMs (FL novel loci)' <- comma(dat$'median length monoexonic TMs (FL novel loci)', digits=0)
dat$'max length monoexonic TMs (FL novel loci)' <- comma(dat$'max length monoexonic TMs (FL novel loci)', digits=0)
dat$'median length spliced TMs (FL novel loci)' <- comma(dat$'median length spliced TMs (FL novel loci)', digits=0)
dat$'max length spliced TMs (FL novel loci)' <- comma(dat$'max length spliced TMs (FL novel loci)', digits=0)

dat <- relocate(dat, 'merge rate', .after=mergedTMs)
dat$n <- comma(dat$n, digits=0)
dat$median <- comma(dat$median, digits=0)
dat$mean <- comma(dat$mean, digits=0)
dat$max <- comma(dat$max, digits=0)
dat$FASTQ_modified <- as.Date(dat$FASTQ_modified, format="%Y-%m-%d")
dat$date_sequenced <- as.Date(as.character(dat$date_sequenced), format="%Y%m%d")
dat$percentMappedReads <- percent(dat$percentMappedReads, digits=0)
dat$'% HCGM reads' <- percent(dat$'% HCGM reads', digits=0)
dat$'% spliced HCGMs' <- percent(dat$'% spliced HCGMs', digits=0)
dat$HCGMreads <- comma(dat$HCGMreads, digits=0)
dat$mergedTMs <- comma(dat$mergedTMs, digits=0)
dat$'merge rate' <- percent(dat$'merge rate', digits=2)
dat$'TM median length' <- comma(dat$'TM median length', digits=0)
dat$'TM max length' <- comma(dat$'TM max length', digits=0)
dat$'SIRV Tx Sn' <- percent(dat$'SIRV Tx Sn'/100, digits=0)
dat$'SIRV Tx Pr' <- percent(dat$'SIRV Tx Pr'/100, digits=0)
dat$'% FL TMs' <- percent(dat$'% FL TMs', digits=0)
dat$'genomic nts covered' <- comma(dat$'genomic nts covered',digits=0)

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


dat <- rename(dat, 'capped spike-ins?' = cappedSpikeIns)
dat <- rename(dat, 'FASTQ last modified' = FASTQ_modified)
dat <- rename(dat, 'date sequenced' = date_sequenced)
dat <- rename(dat, 'reverse transcriptase' = reverse_transcriptase)
dat <- rename(dat, '# reads' = n)
dat <- rename(dat, 'median read length' = median)
dat <- rename(dat, 'mean read length' = mean)
dat <- rename(dat, 'max read length (FASTQ)' = max)
dat <- rename(dat, '% mapped reads'=percentMappedReads)
dat <- rename(dat, 'seq platform' = seqPlatform)
dat <- rename(dat, 'seq center' = seqCenter)
dat <- rename(dat, 'library prep' = libraryPrep)
dat <- rename(dat, 'size fraction' = sizeFrac)
dat <- rename(dat, 'cell fraction' = cellFrac)
dat <- rename(dat, 'bio rep' = bioRep)
dat <- rename(dat, 'tech rep' = techRep)
dat <- rename(dat, 'capture design' = captureDesign)
dat <- rename(dat, '# HCGM reads' = HCGMreads)
dat <- rename(dat, '# merged TMs (min. 2 reads)' =mergedTMs)
dat <- rename(dat, '# TM nts over intergenic regions' = intergenic.x, '% TM nts over intergenic regions' = intergenic.y, '# TM nts over CDS regions'=CDS.x, '% TM nts over CDS regions'=CDS.y, '# TM nts over UTR regions'=UTR.x, '% TM nts over UTR regions'=UTR.y, '# TM nts over exons of noncoding transcripts' =exonOfNCT.x, '% TM nts over exons of noncoding transcripts' =exonOfNCT.y, '# TM nts over exons of pseudogenes' = exonOfPseudo.x, '% TM nts over exons of pseudogenes' = exonOfPseudo.y, '# TM nts over introns' = intron.x, '% TM nts over introns' = intron.y)


tb <- formattable(dat, 
	list('# reads' = color_tile("#def7e9", "#45ba78"),
	'median read length' = color_tile("#e6f7ff", "#0088cc"),
	'mean read length' = color_tile("#e6f7ff", "#0088cc"),
	'max read length (FASTQ)' = color_tile("#e6f7ff", "#0088cc"),
	'% mapped reads' = color_tile("#def7e9", "#45ba78"),
	'# HCGM reads' = color_tile("#def7e9", "#45ba78"),
	'% HCGM reads' = color_tile("#def7e9", "#45ba78"),
	'% spliced HCGMs' = color_tile("#def7e9", "#45ba78"), 
	'# merged TMs (min. 2 reads)' = color_tile('#ffe6cc','#ff8c1a'),
	'merge rate' = color_tile('#ffe6cc','#ff8c1a'),
	'TM median length'= color_tile("#e6f7ff", "#0088cc"),
	'TM max length'= color_tile("#e6f7ff", "#0088cc"),
	'% FL TMs' = color_tile('#f8ecf8', '#d17ad1'),
	'SIRV Tx Sn'=color_tile('#ffeee6','#ff7733'),
	'SIRV Tx Pr'=color_tile('#ffeee6','#ff7733'),
	'genomic nts covered' = color_tile("#f9ecf2", "#c6538c"),
	'# TM nts over intergenic regions' = color_tile("#f9ecf2", "#c6538c"),
	'% TM nts over intergenic regions' = color_tile("#f9ecf2", "#c6538c"),
	'# TM nts over CDS regions'=color_tile("#f9ecf2", "#c6538c"),
	'% TM nts over CDS regions'=color_tile("#f9ecf2", "#c6538c"),
	'# TM nts over UTR regions'=color_tile("#f9ecf2", "#c6538c"),
	'% TM nts over UTR regions'=color_tile("#f9ecf2", "#c6538c"),
	'# TM nts over exons of noncoding transcripts' =color_tile("#f9ecf2", "#c6538c"),
	'% TM nts over exons of noncoding transcripts' =color_tile("#f9ecf2", "#c6538c"),
	'# TM nts over exons of pseudogenes' = color_tile("#f9ecf2", "#c6538c"),
	'% TM nts over exons of pseudogenes' = color_tile("#f9ecf2", "#c6538c"),
	'# TM nts over introns' = color_tile("#f9ecf2", "#c6538c"),
	'% TM nts over introns' = color_tile("#f9ecf2", "#c6538c"),
	'# novel loci' = color_tile('#e6ffff', '#00cccc'),
	'# novel loci per million mapped reads' = color_tile('#e6ffff', '#00cccc'),
	'% novel loci intergenic' = color_tile('#e6ffff', '#00cccc'),
	'% novel loci FL' = color_tile('#e6ffff', '#00cccc'),
	'# novel loci FL' = color_tile('#e6ffff', '#00cccc'),
	'# spliced novel loci' = color_tile('#e6ffff', '#00cccc'),
	'% spliced novel loci' = color_tile('#e6ffff', '#00cccc'),
	'# novel loci on repeats' = color_tile('#e6ffff', '#00cccc'),
	'% novel loci on repeats' = color_tile('#e6ffff', '#00cccc')



)
) 

tbDt <- as.datatable(tb[order(tb$tissue), ], extensions = c('Buttons','ColReorder','FixedColumns', 'FixedHeader', 'KeyTable','RowGroup', 'Select', 'SearchPanes'), 
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
	autowidth=TRUE,
	keys = TRUE,
	rowId =0,
	initComplete = JS(
    "function(settings, json) {",
    "$('body').css({'font-family': 'Helvetica'});",
    "}"
  )
	)
	)
#dir.create('./html/')
DT::saveWidget(tbDt, paste0(getwd(), "/", outputHtml), title='Dashboard of Experiments')
write_tsv(dat, paste0(tools::file_path_sans_ext(outputHtml), ".tsv", sep=''))
