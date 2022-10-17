# Script to calculate the correlation between our approach SHARE-seq coex with other studies approch : EpiMap and ABCmodel(Activity By Contact model)

library(data.table)
library(ggplot2)
require(plyr)
library(gridExtra)
library(grid)
library(ggpubr)



# Correlation between SHARE-seq and EPIMAP:
final_result = fread("share_epimap_overlap_matched.out", header = F, sep = "\t")
colnames(final_result) = c("chr", "start", "end", "gene", "cor", "tag", "chr1", "start1", "end1", "gene1", "score")
final_result[cor > 0.5] 
## Correlation test:
cor.test(final_result$cor, final_result$score, method ="spearman") 
## Plot x/y:
p1 = ggplot(final_result, aes(x=cor, y=score) ) +
  geom_bin2d() +
  geom_smooth(method = "lm")+
  theme_bw()+
  labs(title = "Spearman R:0.18, P-val<2.2e-16", x="SHARE-seq correlation", y="Epimap Score" )+
  theme(plot.title=element_text( hjust=1, vjust=-7, size = 8, face="bold"),
        aspect.ratio = 1,
        legend.key.size = unit(0.25, 'cm'),
        axis.title.x = element_text(colour = "black", size = 8),
        axis.title.y = element_text(colour = "black", size = 8)
  )



# Correlation between SHARE-seq and ABCmodel:
final_result = fread("share_abcmodel_overlap_matched.out", header = F, sep = "\t")
colnames(final_result) = c("chr", "start", "end", "gene", "cor", "tag", "chr1", "start1", "end1", "gene1", "ABC.Score")
## Correlation test:
cor.test(final_result$cor, final_result$ABC.Score, method ="spearman") ##0.04858158
## Plot x/y:
p2= ggplot(final_result, aes(as.numeric(cor), ABC.Score)) +
  geom_bin2d() +
  geom_smooth(method = "lm")+
  theme_bw()+
  labs(title = "Spearman R:0.048, P-val<2.2e-16", x="SHARE-seq correlation", y="ABCmodel Score" ) +
  theme(plot.title=element_text( hjust=1, vjust=-7, size = 8, face="bold"),
        aspect.ratio = 1,
        legend.key.size = unit(0.25, 'cm'),
        axis.title.x = element_text(colour = "black", size = 8),
        axis.title.y = element_text(colour = "black", size = 8))




# Correlation between EPIMAP and ABCmodel:
final_result = fread("epimap_abcmodel_overlap_matched.out", header = F, sep = "\t")
final_result$tag = apply(final_result, 1, function(x) paste(x[6],x[7],x[8],x[9],sep = "_"))
colnames(final_result) = c("chr", "start", "end", "gene", "score", "chr1", "start1", "end1", "gene1", "ABC.score", "tag1")
## Correlation test:
cor.test(final_result$score, final_result$ABC.score, method ="spearman") ##0.01823965 
## Plot x/y:
p3 = ggplot(final_result, aes(score,ABC.score)) + 
  geom_bin2d() +
  geom_smooth(method = "lm")+
  labs(title = "Spearman R:0.018, P-val<9.974e-11" ,x="Epimap Score", y="ABCmodel Score")+
  theme_bw()+
  theme(plot.title=element_text( hjust=1, vjust=-7, size = 8, face="bold"),
        aspect.ratio = 1,
        legend.key.size = unit(0.25, 'cm'),
        axis.title.x = element_text(colour = "black", size = 8),
        axis.title.y = element_text(colour = "black", size = 8),
        
  )

## Arrange multiple plots on the same page
ggarrange(p1,p2,p3, ncol = 3, nrow = 1)



