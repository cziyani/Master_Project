# Script to explore distance between pairs of co-active enhancers

library(data.table)
library(ggplot2)
library(dplyr)
library(plyr)
library(ggthemes)

## Load and process the SHARE-seq coex data : list of enhancer per gene 
data = fread("gene_enhanerList.out",header=T, stringsAsFactors = FALSE) 
data=unique(data)
data$gene = NULL
data$count =lengths(strsplit(data$enhnacerList, split = ','))
## keep only genes that have 2 enh
data = data[count == 2]
## Create a column for each enhancer in the list
data$enh1 = data.table(unlist(lapply(data$enhnacerList, function(x) unlist(strsplit(x,"[,]"))[1])))$V1
data$enh2 = data.table(unlist(lapply(data$enhnacerList, function(x) unlist(strsplit(x,"[,]"))[2])))$V1
data$enhnacerList=NULL                                   
data$start1= data.table(unlist(lapply(data$enh1, function(x) unlist(strsplit(x,"[_]"))[2])))$V1
data$end1 = data.table(unlist(lapply(data$enh1, function(x) unlist(strsplit(x,"[_]"))[3])))$V1
data$start2 = data.table(unlist(lapply(data$enh2, function(x) unlist(strsplit(x,"[_]"))[2])))$V1
data$end2 = data.table(unlist(lapply(data$enh2, function(x) unlist(strsplit(x,"[_]"))[3])))$V1
distanceData =  data[,.(start1,end1, start2,end2, NBcell,count )]

## Calculate the midpoint between pairs of enhancers for Distnace calculation
distanceData$mean1 = (as.numeric(distanceData$start1)+ as.numeric(distanceData$end1)) / 2
distanceData$mean2 = (as.numeric(distanceData$start2)+ as.numeric(distanceData$end2)) / 2
distanceData = distanceData[,-c(1:4)]
## Sum up cells
dt2 = distanceData[, sum(NBcell), by=c("mean1", "mean2", "count")]
## Calculate the distance between pairs of enhnacers regulating the same gene
dt2$distance = abs(as.numeric(dt2$mean2) - as.numeric(dt2$mean1))

## Box-plot
dt2$bins<- cut(as.numeric(dt2$distance), breaks = c(0, 50000, 100000, 200000, 500000, 1000000, 1500000,2000000),
               labels = c("0-50KB", "50KB-100KB", "100KB-200KB", "200KB-500KB", "500KB-1000KB", "1000KB-1500KB", "1500KB-2000KB" ),
               include.lowest = TRUE)
colnames(dt2) = c("mean1", "mean2","count" ,"NBcell", "distance", "bins")


df <- ddply(dt2,
            .(bins),
            summarise,
            min = min(NBcell),
            q1 = quantile(NBcell,0.25),
            med = median(NBcell),
            q3 = quantile(NBcell,0.75),
            max= max(NBcell),
            lab = length(bins))

colnames(dt2)[6] = "Distance"

ggplot(dt2, aes(dt2$Distance,dt2$NBcell, fill=Distance)) +
  geom_boxplot(alpha = 0.4)+
  scale_y_continuous(limits=c(0,100))+
  theme(legend.position="none") +
  #geom_text(data = df,aes(bins, min, label = lab), vjust = 0.5)
  stat_summary(fun.y = median, fun.ymax = length,
               geom = "text", aes(label = ..ymax..), vjust = -0.5)+
  labs(y="Number of cells", x = "Distance between co-active enhancers")+
  theme_bw()+
  theme(plot.title=element_text( hjust=1, vjust=-17, size = 10, face="bold"),
        axis.title.x = element_text(colour = "black", face="bold"),
        axis.title.y = element_text(colour = "black", face="bold"))
