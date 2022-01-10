######## explore observed and expected enhancers combinations


library(data.table)
library(ggplot2)
library(plyr)
library(dplyr)
library(splitstackshape)
library(qdap)
library(magrittr)
library(purrr)
library(ggpubr)


######load and process data

df= fread("gene_enhanerList.out", header = T, sep = "\t")
df1 = data.table(table(df$gene)) 
colnames(df1) = c("gene", "obs_comb")
df$count = lengths(strsplit(df$enhnacerList, split = ','))
max(df$count) ##23


dt= fread("gene_enhancer_cell.out", header = T, sep = "\t") ##2224692 associations
##remove the cell column
dt$cell = NULL
##keep unique associations
dt = unique(dt) 
##data table on the gene
dt1 = data.table(table(dt$gene)) 
## keep gene that have more than 1 enhancer
###because we are interested by stydyin g co-active enhancer within the same gene
dt1 = dt1[N > 1]

##factorial : expected combinations
dt1$fact = as.integer(factorial(dt1$N))+1

##remove the NA
dt1 = dt1[!is.na(dt1$fact)] ##remove NA because we have some gene that have 20 enhancer
colnames(dt1) = c("gene","nb_enh","exp_comb")

##merge
dt_final = merge(df1, dt1, by = "gene")

##check if the oberved combinations are not bigger than the expected ones 
dt_final[obs_comb > exp_comb]

##density plot
d1 = data.table(dt_final$gene, stack(dt_final[,2:4]))
summary(d1)

d1 = d1[ind != "nb_enh"]

colnames(d1)[3] = "Group"

##apply a wilcoxon test
wilcox.test(df1$obs_comb, dt1$exp_comb)

###plot
ggplot(d1, aes(x=log10(values), fill = Group)) + 
    geom_density(stat = "density", alpha = 0.5, adjust = 5) +  
  scale_fill_brewer(palette = "Set2")+
  labs(x="Gene-enhancer combinations", y = "Density")+
  theme(aspect.ratio = 1)+
  theme(legend.title=element_blank())+
  ggtitle("Observed mean: 92
           Expected mean: 13454982
           Wilcoxon p-value: 2.2e-16")+
  theme(plot.title=element_text( hjust=1, vjust=-17, size = 10, face="bold"),
        axis.title.x = element_text(colour = "black", face="bold"),
        axis.title.y = element_text(colour = "black", face="bold"))
  


#####explore and validate previous results

dt_final$obs_log = log10(dt_final$obs_comb)
dt_final$exp_log = log10(dt_final$exp_comb)

unique(dt_final$exp_comb)

dp2 = ggplot( dt_final, aes(x = obs_comb, y = exp_comb ))  +
  geom_point( position = position_jitter(width = 0.01, height = 0.01), alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  scale_x_log10(limits = c(1,max(dt_final$exp_comb))) +
  scale_y_log10(limits = c(1,max(dt_final$exp_comb))) +
  theme_linedraw() + 
  labs(x="observed combinations", y="expected combinations")+
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=24), 
        legend.text=element_text(size=20), legend.title=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black", size = 1),
        panel.border = element_rect(colour = "black", fill=NA, size=1))



ggplot(dt_final, aes(x=obs_comb, y = exp_comb)) + 
  geom_point(position = "jitter")

##create a data table on exp
c = data.table(table(dt_final$exp_comb))

##kepp genes that have 3 or 4 enh
d = dt_final[nb_enh == 4] ##600 gene that have 4 enhancers

##keep only the cases when the obs_comb is too diff from the exp_comb
d2 = d[obs_comb != exp_comb] ##214 gene that have 3/4 enhancers whre the exp combination
#are diff from the observed ones

##keep only the gene column
d2 = d2[,.(gene)]

##save it
write.table(d2, "gene_list",col.names = T , row.names = F , quote = F, sep = "\t" )
