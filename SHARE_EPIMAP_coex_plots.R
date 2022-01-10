################ Explore SHARE-seq and EpiMap datasets  


library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)
library(ggpubr)

#############load and process data

df = fread("coex_peak3.tsv", header = T, sep = "\t")

##create 2 new col : gene / enh
df$gene = data.table(unlist(lapply(df$pairID, function(x) unlist(strsplit(x,"[|]"))[1])))$V1
df$enh = data.table(unlist(lapply(df$pairID, function(x) unlist(strsplit(x,"[|]"))[2])))$V1

## remove the NaN
df = df[!is.na(corr)] 

##count how many enhancer ??er gene i our data
enh_per_gene = data.table(table(df$gene))
summary(enh_per_gene) 

########## plot
ggplot(enh_per_gene , aes(x=N)) + geom_histogram(binwidth = 1, color = "lightblue")+
  labs(title= "Number of enhancer per gene",
       y="gene", x = "enhancer")



## consider 5% FDR and corr > 0.05
dt = df[corr > 0.05]
dt$p.adj = p.adjust(dt$corrPval, method = "BH")
dt = dt[dt$p.adj < 0.05] 

##number of enh per gene after the cutoff
enh_per_gene_SHARE = data.table(table(dt$gene))

#######plot
P1 = ggplot(enh_per_gene_SHARE , aes(x=N)) + geom_histogram(binwidth = 1, color = "lightblue", alpha = 1)+
  labs(y="Number of gene", x = "Number of enhancer") +
  ggtitle("            Min: 1
                       Median: 3
                       Mean: 4,6")  +
  theme(plot.title=element_text( hjust=1, vjust=-17, face='bold', size = 8))




#1)gene-enh epimap association data----
de = fread("links_by_group.lymphoblastoid.tsv", header = T, sep = "\t")
de = subset(de, chr != "chrX" & chr != "chrY")

###enhancer coordiante
de$tag = apply(de, 1, function(x) paste(x[5],x[6],x[7],sep = "_"))

unique(data.table(table(de$tag)))

##plot the number of enhancer per gene
enh_per_gene_epimap = data.table(table(de$gene))
summary(enh_per_gene_epimap)

####Plot
P2 = ggplot(enh_per_gene_epimap, aes(x=N)) + geom_histogram(binwidth = 10,color = "orange", alpha = 0.4) +
  labs( y="Number of gene", x = "Number of enhancer") +
  ggtitle("                  Min: 1
                             Median: 28
                             Mean: 47") +
  theme(plot.title=element_text( hjust=1, vjust=-17, face = "bold", size = 8))


#### combine Plots
ggarrange(P1, P2 , 
          labels = c("(A) SHARE-seq", "(B) EpiMap"),
          ncol = 2, nrow = 1)



