

library(data.table)
library(ggplot2)
library(plyr)
library(gprofiler2)


##read SHARE-data
df = fread("coex_peak3.tsv", header = T, sep = "\t")
df[corrSign == "-"]$corr = -df[corrSign == "-"]$corr
df$gene = data.table(unlist(lapply(df$pairID, function(x) unlist(strsplit(x,"[|]"))[1])))$V1
df$enh = data.table(unlist(lapply(df$pairID, function(x) unlist(strsplit(x,"[|]"))[2])))$V1
df$chr = data.table(unlist(lapply(df$enh, function(x) unlist(strsplit(x,"[_]"))[1])))$V1
df$start=data.table(unlist(lapply(df$enh, function(x) unlist(strsplit(x,"[_]"))[2])))$V1
df$end=data.table(unlist(lapply(df$enh, function(x) unlist(strsplit(x,"[_]"))[3])))$V1
## remove the NaN
df = df[!is.na(corr)]
####create another data.table
share_enh = df[,.(chr, start, end, gene, corr)]
share_enh$chr = sub("^", "chr", share_enh$chr )
##change the col name
colnames(share_enh)[5] = "score"
##create a new col called tag
share_enh$tag = apply(share_enh, 1, function(x) paste(x[1],x[2],x[3],x[4],sep = "_"))
##save the file
write.table(share_enh, "share_enhancer.bed",col.names = F , row.names = F , quote = F, sep = "\t" )



##read EPIMAP-data
de = fread("links_by_group.lymphoblastoid.tsv", header = T, sep = "\t")
####create another data.table
epimap_enh = de[,.(chr, start, end, gene, score)]
##save the file
write.table(epimap_enh, "epimap-enhancer.bed",col.names = F , row.names = F , quote = F, sep = "\t" )



##read the abc-model data
de = fread("GM12878-Roadmap", header = T, sep = "\t")
##rename the 7th col
colnames(de)[7] <- "gene"
##create another data.table 
abc_model = de[,.(chr, start, end, gene, ABC.Score)]
##convert gene name to gene ID : it will as the same time remove nan
abc_model_n = gconvert(query = abc_model$gene, organism = "hsapiens", target="ENSG")
##rename the second col
colnames(abc_model_n)[2] <- "gene"
#join 2 data table, if we don t give "by" it will by default search 
#for the common col if not it will just join both table
ABC_final <-join(abc_model, abc_model_n, type="left")
##remove the na from the table
ABC_final = ABC_final[!is.na(input_number)]
####create another data.table
ABC_model_final = ABC_final[,.(chr, start, end, target, ABC.Score)]
##rename the 4th col
colnames(ABC_model_final)[4] = "gene"
##save the file
write.table(ABC_model_final,"ABCmodel_enhancer.bed",col.names = F , row.names = F , quote = F, sep = "\t" )

