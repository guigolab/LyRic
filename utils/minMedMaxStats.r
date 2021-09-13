#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

input=args[1]
# input should be 1 column with numeric values, no header

library(tidyverse)
library(data.table)
dat<-fread(input, header=F)

#summarise function doesn't like empty data frames, so we need to populate it if empty:
if(length(dat)==0){
dat <-data.frame(NA)
}
names(dat) <- c('values')
#dat
dat %>%
  summarise(med=median(values), max=max(values), min=min(values)) -> datSumm

#datSumm
write(datSumm$min, '')
write(datSumm$med, '')
write(datSumm$max, '')
