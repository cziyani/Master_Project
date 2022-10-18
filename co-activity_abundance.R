# Script to explore enhancers co-activity abundance / Frequency

library(data.table)
library(ggplot2)
library(dplyr)
library(plyr)
library(tidyverse)
library(ggthemes)


## Load and process the data : list of enhancers that are active in the same gene and cell
data = fread("gene_enhanerList.out",header=T, stringsAsFactors = FALSE)
## Count how many enhancers exit in the list
data$count =lengths(strsplit(data$enhnacerList, split = ','))
## keep only cell that have min 2 co-active enh
data = data[count>=2]

## Load and process the data for the second time
df = fread("gene_enhanerList.out",header=T, stringsAsFactors = FALSE)
df$enhnacerList = NULL
df1 = data.table(table(df$gene))
colnames(df1) = c("gene", "total_NBcell")
final_result = merge(data, df1, by="gene")
final_result$cell = final_result$NBcell/final_result$total_NBcell
final_result= unique(final_result)
## Give the same gene(s) the same max number of enh
dt = final_result[, countMax := max(count), by = .(gene)][, enhnacerList := NULL][]
## Create a col called evidence
dt$evidence = 0
## Fill the evidence column based on cell percentage
dt[cell >= 0.05]$evidence = 1 ##1756
dt1 = dt[,.(gene, countMax, evidence)]
## Remove redundancy
dt1 =dt1 %>%
  arrange(gene, desc(evidence)) %>% 
  group_by(gene, countMax) %>% 
  slice(1L) 
dt1 = data.table(dt1)
dt2 = data.table(table(dt1$evidence, dt1$countMax))
colnames(dt2) = c("evidence", "NBenh", "NBgene")
dt2$NBenh = as.numeric(dt2$NBenh)
summary(dt2$evidence)

## Plot
ggplot(dt2,aes(dt2$NBenh, dt2$NBgene, fill = evidence)) +
  geom_bar(stat="identity")+
  xlab("Number of enhancer") + 
  ylab("Number of gene")+
  theme_bw()+
  ggtitle("Number of enhancer per gene > 2")+
  theme(plot.title=element_text( hjust=0.95, vjust=-10, size = 8, face="bold"))
