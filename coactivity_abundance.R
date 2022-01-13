
library(data.table)
library(ggplot2)
library(dplyr)
library(plyr)
library(tidyverse)
library(ggthemes)


##read the data: list of enh that are active in the same gene and cell
data = fread("gene_enhanerList.out",header=T, stringsAsFactors = FALSE)

##count how many enh exit in the list
data$count =lengths(strsplit(data$enhnacerList, split = ','))

##keep only cell that have min 2 co-active enh
data = data[count>=2]

##read the data
df = fread("gene_enhanerList.out",header=T, stringsAsFactors = FALSE)

##remove the tag col
df$enhnacerList = NULL

##create data.table on gene
df1 = data.table(table(df$gene))

##rename col
colnames(df1) = c("gene", "total_NBcell")

##merge data
final_result = merge(data, df1, by="gene")

##find the percentage of cells
final_result$cell = final_result$NBcell/final_result$total_NBcell

##keep only unique associations
final_result= unique(final_result)

##give the same genes the same max number of enh
dt = final_result[, countMax := max(count), by = .(gene)][, enhnacerList := NULL][]

##create a col called evidence
dt$evidence = 0

##create evidence col based on cell percentage
dt[cell >= 0.05]$evidence = 1 

##create a new data.table
dt1 = dt[,.(gene, countMax, evidence)]


dt1 =dt1 %>%
  arrange(gene, desc(evidence)) %>% 
  group_by(gene, countMax) %>% 
  slice(1L) 

##data.table
dt1 = data.table(dt1)

##data table on evidence and count
dt2 = data.table(table(dt1$evidence, dt1$countMax))

##rename col
colnames(dt2) = c("evidence", "NBenh", "NBgene")
dt2$NBenh = as.numeric(dt2$NBenh)

##plot
ggplot(dt2,aes(dt2$NBenh, dt2$NBgene, fill = evidence)) +
  geom_bar(stat="identity")+
  xlab("Number of enhancer") + 
  ylab("Number of gene")+
  theme_bw()+
  ggtitle("Number of enhancer per gene > 2")+
  theme(plot.title=element_text( hjust=0.95, vjust=-10, size = 8, face="bold"))







