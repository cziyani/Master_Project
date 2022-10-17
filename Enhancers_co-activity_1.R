# Script to create a list of enhancers that are active simultaneously in the same gene(s) and cell(s) where the gene expressed.

library(data.table)
library(ggplot2)


## Load and process gene expression data
geneDT= fread("GSM4156603_GM12878.rep3.rna.counts.ensg.final.binary.sparse.bed", header = F, sep = "\t")
data.table(table(geneDT$V2))
colnames(geneDT) = c("gene", "cell")
## remove sequencing identifiers code from Cell ID cordinate
geneDT$cell = gsub(",P1.50","",as.character(geneDT$cell))
geneDT$cell = gsub(",P1.51","",as.character(geneDT$cell))


## Load and process Peaks data
enhDT= fread("peak_mergenhancer_overlap_F_0.5_sorted.out", header = F, sep = "\t")
colnames(enhDT) = c("chr", "start", "end", "cell")
## Remove code identifiers from the cell ID coordinate
enhDT$cell = gsub(",P1.02","",as.character(enhDT$cell))
enhDT$cell = gsub(",P1.03","",as.character(enhDT$cell))
enhDT$chr= gsub( " ", "", enhDT$chr) 
enhDT$start= gsub( " ", "", enhDT$start) 
enhDT$end= gsub( " ", "", enhDT$end) 
enhDT$tag = apply(enhDT,1,function(x)paste(x[1],x[2],x[3],sep ="_"))
enhDT = enhDT[,-c(1,2,3)]


## Load and process Share-seq coex data
coexDT=fread("coex-peak3-corr-0.05.bed", header = F, sep = "\t")
data.table(table(coexDT$gene)) ##5215 gene
coexDT = coexDT[,-c(2,3,4,5,6,7,8,12)]
colnames(coexDT) = c("gene", "chr", "start", "end")
coexDT$chr = sub("^", "chr", coexDT$chr)
coexDT$chr= gsub( " ", "", coexDT$chr) 
coexDT$start= gsub( " ", "", coexDT$start) 
coexDT$end= gsub( " ", "", coexDT$end) 
coexDT$tag = apply(coexDT,1,function(x)paste(x[2],x[3],x[4],sep ="_"))
coexDT = coexDT[,-c(2,3,4)]

                   
## Filter the previous datas
geneDT = geneDT[gene %in% unique(coexDT$gene)]
geneDT = geneDT[cell %in% unique(enhDT$cell)]
enhDT = enhDT[tag %in% unique(coexDT$tag)]
enhDT = enhDT[cell %in% unique(geneDT$cell)]

                   
## Merge data
d1 = unique(merge(unique(coexDT), unique(geneDT), by = "gene", allow.cartesian = T))
d2 = unique(merge(unique(coexDT), unique(enhDT), by = "tag", allow.cartesian = T))
mergedData = merge(d1, d2, by = c("gene","tag","cell")) 

                   
# Save file
write.table(mergedData,"gene_enhancer_cell.out",quote=F,sep="\t",col.names=T,row.names=F)
