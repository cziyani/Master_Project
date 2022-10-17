# Script to explore SHARE-seq coex and EpiMap datasets  

library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)
library(ggpubr)


## Load and process SHARE-seq coex data
df = fread("coex_peak3.tsv", header = T, sep = "\t")
df$gene = data.table(unlist(lapply(df$pairID, function(x) unlist(strsplit(x,"[|]"))[1])))$V1
df$enh = data.table(unlist(lapply(df$pairID, function(x) unlist(strsplit(x,"[|]"))[2])))$V1
df = df[!is.na(corr)] 
## Count how many enhancer(s) per gene in the coex data
enh_per_gene = data.table(table(df$gene))
summary(enh_per_gene) 
## Plot the Number of Enhancer per gene
ggplot(enh_per_gene , aes(x=N)) + geom_histogram(binwidth = 1, color = "lightblue")+
  labs(title= "Number of enhancer per gene",
       y="gene", x = "enhancer")          
## Filter : Consider 5% FDR and corr > 0.05
dt = df[corr > 0.05]
dt$p.adj = p.adjust(dt$corrPval, method = "BH")
dt = dt[dt$p.adj < 0.05] 
## Count how many enhancer(s) per gene after the cutoff
enh_per_gene_SHARE = data.table(table(dt$gene))
## Plot
P1 = ggplot(enh_per_gene_SHARE , aes(x=N)) + geom_histogram(binwidth = 1, color = "lightblue", alpha = 1)+
  labs(y="Number of gene", x = "Number of enhancer") +
  ggtitle("            Min: 1
                       Median: 3
                       Mean: 4,6")  +
  theme(plot.title=element_text( hjust=1, vjust=-17, face='bold', size = 8))





## Load and process EPIMAP data
de = fread("links_by_group.lymphoblastoid.tsv", header = T, sep = "\t")
de = subset(de, chr != "chrX" & chr != "chrY")
de$tag = apply(de, 1, function(x) paste(x[5],x[6],x[7],sep = "_"))
unique(data.table(table(de$tag)))
## Count the number of enhancer per gene
enh_per_gene_epimap = data.table(table(de$gene))
summary(enh_per_gene_epimap)
## Plot
P2 = ggplot(enh_per_gene_epimap, aes(x=N)) + geom_histogram(binwidth = 10,color = "orange", alpha = 0.4) +
  labs( y="Number of gene", x = "Number of enhancer") +
  ggtitle("                  Min: 1
                             Median: 28
                             Mean: 47") +
  theme(plot.title=element_text( hjust=1, vjust=-17, face = "bold", size = 8))

               

## Combine both Plots on the same page
ggarrange(P1, P2 , 
          labels = c("(A) SHARE-seq", "(B) EpiMap"),
          ncol = 2, nrow = 1)
