# Script to create simplified files to compare SHARE-seq with other datasets (EPIMAP, ABCmodel)

library(data.table)
library(ggplot2)
library(plyr)
library(gprofiler2)


## Load and process SHARE-seq coex data
df = fread("coex_peak3.tsv", header = T, sep = "\t")
df[corrSign == "-"]$corr = -df[corrSign == "-"]$corr
df$gene = data.table(unlist(lapply(df$pairID, function(x) unlist(strsplit(x,"[|]"))[1])))$V1
df$enh = data.table(unlist(lapply(df$pairID, function(x) unlist(strsplit(x,"[|]"))[2])))$V1
df$chr = data.table(unlist(lapply(df$enh, function(x) unlist(strsplit(x,"[_]"))[1])))$V1
df$start=data.table(unlist(lapply(df$enh, function(x) unlist(strsplit(x,"[_]"))[2])))$V1
df$end=data.table(unlist(lapply(df$enh, function(x) unlist(strsplit(x,"[_]"))[3])))$V1
df = df[!is.na(corr)]
share_enh = df[,.(chr, start, end, gene, corr)]
share_enh$chr = sub("^", "chr", share_enh$chr )
colnames(share_enh)[5] = "score"
share_enh$tag = apply(share_enh, 1, function(x) paste(x[1],x[2],x[3],x[4],sep = "_"))
write.table(share_enh, "share_enhancer.bed",col.names = F , row.names = F , quote = F, sep = "\t" )



## Load and process EPIMAP data
de = fread("links_by_group.lymphoblastoid.tsv", header = T, sep = "\t")
epimap_enh = de[,.(chr, start, end, gene, score)]
write.table(epimap_enh, "epimap-enhancer.bed",col.names = F , row.names = F , quote = F, sep = "\t" )



## Load and process ABCmodel data
de = fread("GM12878-Roadmap", header = T, sep = "\t")
colnames(de)[7] <- "gene"
abc_model = de[,.(chr, start, end, gene, ABC.Score)]
##convert gene name to gene ID
abc_model_n = gconvert(query = abc_model$gene, organism = "hsapiens", target="ENSG")
colnames(abc_model_n)[2] <- "gene"
ABC_final <-join(abc_model, abc_model_n, type="left")
ABC_final = ABC_final[!is.na(input_number)]
ABC_model_final = ABC_final[,.(chr, start, end, target, ABC.Score)]
colnames(ABC_model_final)[4] = "gene"
write.table(ABC_model_final,"ABCmodel_enhancer.bed",col.names = F , row.names = F , quote = F, sep = "\t" )

