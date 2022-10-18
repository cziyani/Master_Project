# Script to explore SHARE-seq coex significant associations.

library(data.table)
library(ggplot2)
library(plyr)
library(dplyr)
library(splitstackshape)
library(qdap)

## Load and process the SHARE-seq coex data
df = fread("coex_peak3.tsv", header = T, sep = "\t")
df[corrSign == "-"]$corr = -df[corrSign == "-"]$corr
df$gene = data.table(unlist(lapply(df$pairID,function(x)unlist(strsplit(x,"[|]"))[1])))$V1
df$enh = data.table(unlist(lapply(df$pairID,function(x)unlist(strsplit(x,"[|]"))[2])))$V1
df$chr = data.table(unlist(lapply(df$enh,function(x)unlist(strsplit(x,"[_]"))[1])))$V1
df$start=data.table(unlist(lapply(df$enh,function(x)unlist(strsplit(x,"[_]"))[2])))$V1
df$end=data.table(unlist(lapply(df$enh,function(x)unlist(strsplit(x,"[_]"))[3])))$V1
df = df[!is.na(corr)]


## Apply a filter: consider only the corr > 0.05 & FDR 5%
dt = df[corr > 0.05]
dt$p.adj = p.adjust(dt$corrPval, method = "BH")
dt = dt[dt$p.adj < 0.05] 


## Count the number of enhancer per gene
enh_per_gene1 = data.table(table(dt$gene)) 
summary(enh_per_gene1)
colnames(enh_per_gene1) = c("gene" , "NBenh")
## Count the number of gene that have more than 1 enh
enh_per_gene1[NBenh > 1] 
## Merge "dt" with the list of the number of enhancer per gene
dt1 = merge(dt, enh_per_gene1, by = "gene") 
## keep only genes that have more than 1 enh 
dt1 = dt1[NBenh > 1] 
summary(dt1)
dt1= dt1[,-c(2,3,4,5,7,12)]
colnames(coexDT) = c("gene", "chr", "start", "end")
## Add the char "chr" and remove spaces
dt1$chr = sub("^", "chr", dt1$chr)
dt1$chr= gsub( " ", "", dt1$chr) 
dt1$start= gsub( " ", "", dt1$start) 
dt1$end= gsub( " ", "", dt1$end) 
## Create a tag column
dt1$tag = apply(dt1,1,function(x)paste(x[4],x[5],x[6],sep ="_"))
dt1 = dt1[,-c(4,5,6)]

## Save the file
write.table(dt1, "coex-peak3-corr-0.05-pval-cor-col.bed",col.names = F , row.names = F , quote = F, sep = "\t" )
