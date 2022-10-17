########## Explore all enhancer within gene-enhancer associations list
########## The produced file will be used while exploring and plotting 
########## Enhancers co-activity


library(data.table)
library(ggplot2)

#####load and process data

data = fread("coex_peak3.tsv",header=T, stringsAsFactors = FALSE)

##keep only the first col containing 
####gene-enhancer associations 
data = data[,.(pairID)]

##split enhancer and gene coordiantes
data$gene = data.table(unlist(lapply(data$pairID, function(x) unlist(strsplit(x,"[|]"))[1])))$V1
data$enh = data.table(unlist(lapply(data$pairID, function(x) unlist(strsplit(x,"[|]"))[2])))$V1

##remove the col
data$pairID = NULL

####save file
write.table(data,"coex_gene_enhancer.tsv",quote=F,sep="\t",col.names=T,row.names=F)
