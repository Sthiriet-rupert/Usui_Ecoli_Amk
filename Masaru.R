

## Heatmap Mut vs MIC
# Biofilm (same for Planktonic data)

data<-read.table("R_data/Mutation-MIC_Biofilm.txt",header = T,sep = "\t")

names<-data$Mut
names<-c("Mut",names)
colnames(data)<-names

library(reshape2)

dat<-melt(data=data,id.vars = c("Mut"),variable.name = "Mut2",value.name = "Freq")

library(ggplot2)

namesx<-data$Mut
namesy<-rev(namesx)

hp<-ggplot(data = dat, aes(x=Mut, y=Mut2, fill=Freq)) + 
    geom_tile(color = "black")+
    scale_fill_gradientn(values=c(1, .4, .1,0), colours=c("red", "orange", "green", "white"), limit = c(0,128),
                         name="Frequency") + theme_minimal()+ 
    theme(axis.text.x = element_text(vjust = 1, 
                                     size = 12, hjust = 1))+
    coord_fixed()


hp + theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size=12),
    axis.text.x = element_text(angle = 45,hjust=-0.01,vjust=0.2),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position="bottom",
    legend.direction = "horizontal")+
    guides(fill = guide_colorbar(barwidth = 12, barheight = 1,
                                 title.vjust = 0.9,title.hjust = 0))   + scale_x_discrete(position = "top",limits=namesx) +  scale_y_discrete(limits=namesy) + theme(plot.margin = margin(.2,3.5,.5,.5, "cm"))



ggsave("FIgures/Mut_vs_MIC_Biofilm.tiff",width=11,height = 10,device = "tiff",dpi = 700)


## Boxplot + Anova on clones MICs
# Color points by selection plate

data<-read.table("Anova_MIC_V2.txt",header = T,sep = "\t")

library(ggpubr)

comp<-list(c("Biofilm-5MIC","Biofilm-80MIC"),c("Biofilm-5MIC","Planktonic-5MIC"),c("Biofilm-5MIC","Planktonic-80MIC"),c("Biofilm-80MIC","Planktonic-5MIC"),c("Biofilm-80MIC","Planktonic-80MIC"),c("Planktonic-5MIC","Planktonic-80MIC"))

p <- ggboxplot(data, x = "Condition", y = "MIC",outlier.shape = NA) +geom_point(position=position_jitterdodge(jitter.width=2, dodge.width = 0), pch=21, aes(fill=factor(SelectionPlate),color=factor(SelectionPlate)))+ scale_color_manual(values=c("#f18d7b", "#b8ce55", "#56B4E9","#bb7de0"))
p + stat_compare_means(method = "kruskal",label.y = 200) + stat_compare_means(comparisons = comp,label="p.signif",label.y=c(140,150,160,170,180,190)) + xlab("") + ylab("Clone MIC (µg / mL)")+ scale_y_continuous(breaks = c(25,50,75,100,125))
ggsave("Figures/Boxplot_Anova_CloneMICs_ColorPoints.tiff",width=7,height = 7,device = "tiff",dpi = 700)



## Heatmaps for mutated genes in the final populations
# (same for the heatmaps with all SNPs)


data<-read.table("R_data/Mut_Genes_Cycle10_Biofilm.txt",header = T,sep = "\t")

library(reshape2)
library(ggplot2)

names<-data$Mutation
namesy<-rev(names)


dat<-melt(data=data,id.vars = c("Mutation"),measure.vars = c("X1","X2","X3","X4","X5","X6","X7","X8","X9"),variable.name = "Sample",value.name = "Freq")

hp<-ggplot(data = dat, aes(x=Sample, y=Mutation, fill=Freq)) + 
    geom_tile(color = "black")+
    scale_fill_gradientn(values=c(1, .5, 0), colours=c("red", "#f08080", "white"), limit = c(0,100),
                         name="Frequency") + theme_minimal()+ 
    theme(axis.text.x = element_text(vjust = 1, 
                                     size = 12, hjust = 1))+
    coord_fixed()


hp + theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size=12),
    axis.text.x = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position="bottom",
    legend.direction = "horizontal")+
    guides(fill = guide_colorbar(barwidth = 12, barheight = 1,
                                 title.vjust = 0.9,title.hjust = 0))   + scale_x_discrete(position = "top")+  scale_y_discrete(limits=namesy)

ggsave("Figures/Heatmap_Genes_Biof.tiff",width=7,height = 10,device = "tiff",dpi = 700)



## Boxplot entropy

library(ggpubr)

data<-read.table("R_data/Entropy.txt",header = T,sep = "\t")

comp<-list(c("Domain I","FusA"),c("Domain II","FusA"),c("Domain III","FusA"),c("Domain IV","FusA"),c("Domain V","FusA"))

p <- ggboxplot(data, x = "Domain", y = "Entropy",color="Domain",outlier.shape = NA,add = "jitter") 
p + stat_compare_means(method = "anova",label.y = 0.45) + stat_compare_means(comparisons = comp,label="p.signif",label.y=c(0.26,0.3,0.35,0.4,0.45)) + scale_y_continuous(breaks = c(0.05,0.10,0.15,0.20,0.25))  + xlab("") + ylab("Entropy")+theme(axis.title.y = element_text(size=18),axis.text.x = element_text(angle = 45,vjust=1,hjust = 1),legend.position="none")

ggsave("Figures/Boxplot_Entropy.tiff",width=7,height = 7,device = "tiff",dpi = 700)



