######## create a list of enhnacer per gene and cell

library(data.table)
library(ggplot2)
library(plyr)
library(dplyr)
library(splitstackshape)
library(qdap)
library(magrittr)

######load and process data

final_data= fread("gene_enhancer_cell.out", header = T, sep = "\t")

####create a list of enhancer
final_data1 = final_data %>% group_by(gene, cell) %>%mutate(enhancerList = toString(tag)) %>%as.data.frame()
#final_data1 =  setDT(data1)[, list(Values=paste(V2, collapse=",")) ,c("V1","V3","V4")]

final_data1 = final_data1[,-c(2)]
final_data1 = unique(final_data1) 

###count the number of cells where gene expressed and enhnacer are active 
s = data.table(table(final_data1$gene, final_data1$enhancerList)) 
s = s[N > 0]

####renqme col
colnames(s) = c("gene", "enhnacerList", "NBcell")

###save file
write.table(s, "gene_enhanerList.out",col.names = T , row.names = F , quote = F, sep = "\t" )
