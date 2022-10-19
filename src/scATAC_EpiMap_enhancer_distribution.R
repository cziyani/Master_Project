# Script to explore EpiMap enhancer as well as scATAC_seq enhancer distribution

library(data.table)
library(ggplot2)

## Load and process Epimap-enhancer data
df = fread("atac_peak_file.bed", header = F, sep = "\t")
colnames(df) = c("chr", "start", "end")

## Plot the distribution
df$size = abs(df$end - df$start)
dp = ggplot(df  , aes(x=size)) + geom_histogram(bins = 100) 
summary(df$size)

## Load and process SHARE-seq scATAC-seq enhancer 
dt = fread("enhancers_processing_task1_sorted_merge.bed", header = F, sep = "\t")
colnames(dt) = c("chr", "start", "end")

## Plot the distribution
dt$size = abs(dt$end - dt$start)
dp1 = ggplot(dt  , aes(x=size)) + geom_histogram(bins = 100) 
