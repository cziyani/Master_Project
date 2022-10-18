# Script to explore scATAC-seq data and LCL enhancer annotation

library(ggplot2)
library(data.table)
library(ggthemes)
library(grid)
library(gridExtra)
library(ggpubr)
 

## Load and process the data
df = fread("peak_mergenhancer_overlap_F_0.5_sorted.out.gz", header = F, sep = "\t")
## Remove chr x and y
df = subset(df, V1 != "chrX" & V1 != "chrY")
df$tag = apply(df, 1, function(x) paste(x[1],x[2],x[3],sep = "_"))
colnames(df) = c("chr", "start", "end", "cell","tag")
enh_per_cell = data.table(table(df$cell)) 
cell_per_enh = data.table(table(df$tag)) 

                        
## Plot the number of enhancer(s) per cell
P1 = ggplot(enh_per_cell , aes(x=N)) + 
  geom_histogram(binwidth = 50, color="black", fill="lightblue") + 
  xlim(c(-50,2000))+
  labs(x="Number of Enhancer", y="Number of cells")+ theme_bw()+
  ggtitle("Min : 1
           Median : 11
           Mean : 88") +
  theme(plot.title=element_text( hjust=1, vjust=-17, size = 10, face = "bold")) 


## Plot the number of cell(s) per enhancer
P2=ggplot(cell_per_enh , aes(x=N)) +
  geom_histogram(binwidth = 50, color="black", fill="orange", alpha = 0.4) +
  xlim(c(-50,3000))+ 
  labs(x="Number of Cell" , y="Number of enhnacers")+ theme_bw()+
  ggtitle("Min : 1
           Median : 193
           Mean : 344") +
  theme(plot.title=element_text( hjust=1, vjust=-17, size = 10, face = "bold"))


## Arrange both plots on the same page
ggarrange(P1, P2 , 
          labels = c("A", "B"),
          ncol = 2, nrow = 1)
