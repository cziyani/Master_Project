##########explore EpiMap enhancer after merging as well as scATAC_seq enhancer distribution

##add the packages
library(data.table)
library(ggplot2)

##epimap-enhancer
df = fread("atac_peak_file.bed", header = F, sep = "\t")
colnames(df) = c("chr", "start", "end")

###plot
df$size = abs(df$end - df$start)
dp = ggplot(df  , aes(x=size)) + geom_histogram(bins = 100) 
summary(df$size)

##atac-enhancer
dt = fread("enhancers_processing_task1_sorted_merge.bed", header = F, sep = "\t")
colnames(dt) = c("chr", "start", "end")

###plot
dt$size = abs(dt$end - dt$start)
dp1 = ggplot(dt  , aes(x=size)) + geom_histogram(bins = 100) 
