# Script to plot enhancers co-activity across genes and cells

library(data.table)
library(ggplot2)
library(dplyr)

## Accept arguments from bash script
args <- commandArgs(trailingOnly = TRUE)

## Produce a "png" plot
png(paste(args,".png",sep=""), res = 300, width = 12, height = 8,units = "in")

## Load and process SHARE-seq coex the data : 
data = fread("gene_enhancer_cell.out",header=T, stringsAsFactors = FALSE)
d = data.table(table(data$gene))
## read the gene-tss file from Genecode
tssData = fread("genes_tss.bed",sep="\t",header=F, stringsAsFactors = FALSE)
tssData$gene = data.table(unlist(lapply(tssData$V4, function(x) unlist(strsplit(x,"[.]"))[1])))$V1
tssData1 = unique(tssData[gene == args])
## The gene-tss is the start of the gene when the strand-sign is +
tssData1$tss = tssData1$V2
## The gene-tss is the end of the gene when the strand-sign is -
tssData1[V6 == "-"]$tss = tssData1[V6 == "-"]$end
## The gene-tss is the start of the gene when the strand-sign is +
tssData$tss = tssData$V2
## The gene-tss is the end of the gene when the strand-sign is -
tssData[V6 == "-"]$tss = tssData[V6 == "-"]$end
##create a data table which contain all the other genes-tss
tssData2 = tssData[!(tssData$tss %in% tssData1$tss),]
## whithin 1MB from the chosen gene tss
#tssData2 = tssData2[tss <= (tssData1$tss+1000000) & tss >= (tssData1$tss-1000000)]
##match chr between the chosen gene and the other gene that are next to it 
#tssData2 = tssData2[(tssData2$V1 %in% tssData1$V1),]
example = data[gene == args]

## Find all the combinations by cell and order it by enhancers
combinationDT = data.table()
for (c in unique(example$cell)){
  tags = unique(example[cell == c][order(tag)]$tag)
  text = ""
  for (tag in tags){text = paste(text, tag, sep = "|")}
  
  combinationDT = rbind(combinationDT, data.table(set = text, cell = c, num = length(tags)))
}
countDT = data.table(table(combinationDT$set))
combinationCountDT = merge(combinationDT,countDT, by.x = "set", by.y = "V1")
example$chr = data.table(unlist(lapply(example$tag, function(x) unlist(strsplit(x,"_"))[1])))$V1
example$start = as.numeric(data.table(unlist(lapply(example$tag, function(x) unlist(strsplit(x,"_"))[2])))$V1)
example$end = as.numeric(data.table(unlist(lapply(example$tag, function(x) unlist(strsplit(x,"_"))[3])))$V1)
example = merge(example, combinationCountDT, by = "cell")

## Read the coex_enhancer file
all_enhancer = fread("coex_gene_enhancer.tsv",header=T, stringsAsFactors = FALSE)
## keep unique enhancers
all_enhancer= unique(all_enhancer)
all_enhancer = all_enhancer[gene %in% unique(example$gene)]
all_enhancer$chr1 = data.table(unlist(lapply(all_enhancer$enh, function(x) unlist(strsplit(x,"_"))[1])))$V1
all_enhancer$start1 = as.numeric(data.table(unlist(lapply(all_enhancer$enh, function(x) unlist(strsplit(x,"_"))[2])))$V1)
all_enhancer$end1 = as.numeric(data.table(unlist(lapply(all_enhancer$enh, function(x) unlist(strsplit(x,"_"))[3])))$V1)
all_enhancer$chr1 = sub("^", "chr", all_enhancer$chr1)
## Order by the num of enhancer and the total num of cells
example = example[order(num, N, -set),]
## rank data based on the cells
rankdf = data.table(cell=unique(example$cell))
rankdf$rank = seq(1,nrow(rankdf))
example = merge(example, rankdf, by="cell")

## Enhancer coordinate 
enhancers = unique(example[,.(start,end)])
## The non-associated enhancer coordinate 
all_enh = unique(all_enhancer[,.(start1, end1)])
all_enh=all_enh[!(all_enh$start1 %in% enhancers$start),]
## The max of the rank
maxcoord = max(example$rank)

## Plot
ggplot() +
  #geom_rect(data = example, aes(xmin = start, xmax = end, ymin = rank, ymax = rank +1, fill = tag, alpha=5), size = 4)+
  geom_linerange(data = example, aes(x=(start+end)/2, ymin=rank, ymax=rank+1, color=tag), size =2)+
  ##gene localization
  geom_rect(data = tssData1, aes(xmin = V2, xmax = V3, ymin = maxcoord+ maxcoord/50, ymax =maxcoord+ maxcoord/ 25), size = 2, color = "grey") +
  ##enhancers that are regulating these gene
  geom_rect(data = enhancers, aes(xmin = start, xmax = end, ymin = maxcoord+ maxcoord/50, ymax = maxcoord+ maxcoord/ 25), size = 2, color = "black") +
  ##all the other enhnacers
  geom_rect(data = all_enh, aes(xmin = start1, xmax = end1, ymin = maxcoord+ maxcoord/50, ymax = maxcoord+ maxcoord/ 25), size = 0.5, color = "green") +
  ##chosen gene tss
  geom_point(data = tssData1, aes(x = V2, y = max(example$rank)+maxcoord/50), size = 3, shape = 1) +
  
  ##gene that are next to the chosen gene withen 1MB
  #geom_point(data = tssData2, aes(x = tss, y = max(example$rank)+maxcoord/50), size = 0.5, shape = 2)+
  
  ggtitle(args) +
  ylab("Cell rank") +
  xlab("Genomic coordinate") +
  theme_minimal() +
  geom_path() +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
                text = element_text(size=10), 
                panel.background = element_rect(colour = "black", fill = "white", size = 1),
                axis.text.x  = element_text(hjust=0.8,face = "bold", colour="black")
  )

dev.off()
