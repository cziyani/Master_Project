###########create a list of enhancer that are active simultaneously in gene and cell where gene expressed


library(data.table)
library(ggplot2)


##########Load/process data

geneDT= fread("GSM4156603_GM12878.rep3.rna.counts.ensg.final.binary.sparse.bed", header = F, sep = "\t")
data.table(table(geneDT$V2))
colnames(geneDT) = c("gene", "cell")

##remove sequencing identifiers code from Cell ID cordinate
geneDT$cell = gsub(",P1.50","",as.character(geneDT$cell))
geneDT$cell = gsub(",P1.51","",as.character(geneDT$cell))

##read the peaks data----
enhDT= fread("peak_mergenhancer_overlap_F_0.5_sorted.out", header = F, sep = "\t")

##rename col
colnames(enhDT) = c("chr", "start", "end", "cell")

##remove code identifiers from the cell ID coordiante
enhDT$cell = gsub(",P1.02","",as.character(enhDT$cell))
enhDT$cell = gsub(",P1.03","",as.character(enhDT$cell))

### remove space
enhDT$chr= gsub( " ", "", enhDT$chr) 
enhDT$start= gsub( " ", "", enhDT$start) 
enhDT$end= gsub( " ", "", enhDT$end) 

##create enhancers coordiante
enhDT$tag = apply(enhDT,1,function(x)paste(x[1],x[2],x[3],sep ="_"))

###remove col
enhDT = enhDT[,-c(1,2,3)]


##read coex data----
coexDT=fread("coex-peak3-corr-0.05.bed", header = F, sep = "\t")
data.table(table(coexDT$gene)) ##5215 gene

##remove col
coexDT = coexDT[,-c(2,3,4,5,6,7,8,12)]

##rename col
colnames(coexDT) = c("gene", "chr", "start", "end")

##add the char "chr
coexDT$chr = sub("^", "chr", coexDT$chr)
coexDT$chr= gsub( " ", "", coexDT$chr) 
coexDT$start= gsub( " ", "", coexDT$start) 
coexDT$end= gsub( " ", "", coexDT$end) 

##create a tag
coexDT$tag = apply(coexDT,1,function(x)paste(x[2],x[3],x[4],sep ="_"))
coexDT = coexDT[,-c(2,3,4)]

##filter data
geneDT = geneDT[gene %in% unique(coexDT$gene)]
geneDT = geneDT[cell %in% unique(enhDT$cell)]
enhDT = enhDT[tag %in% unique(coexDT$tag)]
enhDT = enhDT[cell %in% unique(geneDT$cell)]

##merge data
d1 = unique(merge(unique(coexDT), unique(geneDT), by = "gene", allow.cartesian = T))
d2 = unique(merge(unique(coexDT), unique(enhDT), by = "tag", allow.cartesian = T))

###merge data
mergedData = merge(d1, d2, by = c("gene","tag","cell")) 

####Save file
write.table(mergedData,"gene_enhancer_cell.out",quote=F,sep="\t",col.names=T,row.names=F)
