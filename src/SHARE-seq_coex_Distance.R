# Script to Explore SHARE-seq coex gene-enhancer Distance 

library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)


## Load and process SHARE-seq coex data
df = fread("coex_peak3.tsv", header = T, sep = "\t")
df$gene = data.table(unlist(lapply(df$pairID, function(x) unlist(strsplit(x,"[|]"))[1])))$V1
df$enh = data.table(unlist(lapply(df$pairID, function(x) unlist(strsplit(x,"[|]"))[2])))$V1


## Load genes coordiantes
dg = fread("genes.bed", header = F, sep = "\t")
dg$gene = data.table(unlist(lapply(dg$V4, function(x) unlist(strsplit(x,"[.]"))[1])))$V1
colnames(dg) = c("chr","start1", "end1", "gene-info" , "info", "strand", "gene" )


## Merge based on gene coordinates (SHARE data) 
gene_share = merge(df, dg, by="gene")
gene_share$start = data.table(unlist(lapply(gene_share$enh, function(x) unlist(strsplit(x,"[_]"))[2])))$V1
gene_share$end= data.table(unlist(lapply(gene_share$enh, function(x) unlist(strsplit(x,"[_]"))[3])))$V1


## Midpoint calculation 
gene_share$mean = (as.numeric(gene_share$start )+ as.numeric(gene_share$end)) / 2

## Distance calculation
gene_share$distance = 0
gene_share$tss = gene_share$start1
gene_share[strand == "-"]$tss = gene_share[strand == "-"]$end1
gene_share$distance = abs(gene_share$tss - gene_share$mean)


## Filter : apply a cutoff
gene_share1 = gene_share[ corr > 0.05]
gene_share1$p.adj = p.adjust(gene_share1$corrPval, method = "BH")
gene_share1= gene_share1[gene_share1$p.adj < 0.05]
colnames(gene_share1)[20] <- "distance1"


## Plot the SHARE-seq coex gene-enhancer Distance based on midpoint calculation
ggplot() +
  geom_density(aes(gene_share$distance, fill = "All-tested"), alpha = .5, data = gene_share) +
  geom_density(aes(gene_share1$distance1, fill = "Associated"), alpha = .6, data = gene_share1)+
  theme_bw()+
  ggtitle("Associated: Corr>0.05, FDR 5%")+
  theme(plot.title=element_text( hjust=1, vjust=-8, size = 8, face="bold"))+
  labs(x="Gene-enhancer distance", y = "Density")+
  theme(legend.title=element_blank())
