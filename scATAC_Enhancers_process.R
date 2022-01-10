# Script to explore scATAC-seq data and LCL enhancer annotation

library(ggplot2)
library(data.table)
library(ggthemes)
library(grid)
library(gridExtra)
library(ggpubr)

########
##load and process data
df = fread("peak_mergenhancer_overlap_F_0.5_sorted.out.gz", header = F, sep = "\t")

##remove chr x and y
df = subset(df, V1 != "chrX" & V1 != "chrY")

##create a new column containing enhancer coordinate
df$tag = apply(df, 1, function(x) paste(x[1],x[2],x[3],sep = "_"))

##rename col
colnames(df) = c("chr", "start", "end", "cell","tag")

##count how many enhancer(s) per cell in our data
enh_per_cell = data.table(table(df$cell)) 

##count how many cell(s) per enhancer in our data
cell_per_enh = data.table(table(df$tag)) 

##summary
summary(enh_per_cell$N)
summary(cell_per_enh)


##plot the number of enhancers per cells
P1 = ggplot(enh_per_cell , aes(x=N)) + 
  geom_histogram(binwidth = 50, color="black", fill="lightblue") + 
  xlim(c(-50,2000))+
  labs(x="Number of Enhancer", y="Number of cells")+ theme_bw()+
  ggtitle("Min : 1
           Median : 11
           Mean : 88") +
  theme(plot.title=element_text( hjust=1, vjust=-17, size = 10, face = "bold")) 



##plot the number of cells per enh
P2=ggplot(cell_per_enh , aes(x=N)) +
  geom_histogram(binwidth = 50, color="black", fill="orange", alpha = 0.4) +
  xlim(c(-50,3000))+ 
  labs(x="Number of Cell" , y="Number of enhnacers")+ theme_bw()+
  ggtitle("Min : 1
           Median : 193
           Mean : 344") +
  theme(plot.title=element_text( hjust=1, vjust=-17, size = 10, face = "bold"))



###put both plots in the same page
ggarrange(P1, P2 , 
          labels = c("A", "B"),
          ncol = 2, nrow = 1)

