## Produce the data to be used in R
## Use the lolopop python package
# python3.9 -m pip install lolipop
# lolipop lineage --genotypes --no-outline -o ./Documents/Stan/Projects/EvolTolAB/Masaru/MullerPlots/Biofilm_5 -i ./Documents/Stan/SandBox/Biofilm_5.tsv

# ./Documents/Stan/SandBox/Biofilm_5.tsv looks like this:

#	Gene	Genotype	1	2	3	5	7	9	10
#	SbmA H161P	SbmA H161P	0	0	0	0	0	0	46
#	SbmA L376P	SbmA L376P	0	0	0	44	33	57	14
#	SbmA R167C	SbmA R167C	0	0	0	0	0	0	19
#	SbmA W111*	SbmA W111*	0	0	0	0	0	0	7
#	SbmA W99*	SbmA W99*	0	0	0	0	0	0	6

## Produce the MullerPlots using R from the data previously generated

library(ggplot2)
library(ggmuller)

setwd("~/Documents/Stan/Projects/EvolTolAB/Masaru/MullerPlots/Tables_For_R")
population <- read.table("Biofilm_1.populations.tsv",sep="\t", header=TRUE)
edges <- read.table("Biofilm_1.edges.tsv", sep="\t",header=TRUE)

Muller_df <- get_Muller_df(edges, population)
palette <- c("#F4A460","#EE82EE",	"white")

ggplot(Muller_df, aes_string(x = "Generation", y = "Frequency", group = "Group_id", fill = "Identity", colour = "Identity")) +
geom_area() + theme_classic() +
theme(legend.position = "right",legend.key = element_rect(colour = "black",linetype="solid")) +
guides(linetype = FALSE, color = FALSE) +
scale_y_continuous(labels = 25 * (0:4), name = "Frequency", expand=c(0,0)) +
scale_x_continuous(expand=c(0,0),breaks = seq(0,10,2),name="Cycles") +
scale_fill_manual(name = "Genotype", values = palette) +
scale_color_manual(values = palette) + theme(text = element_text(size=20,face="bold"))

ggsave("MullerPlot_Biofilm-1.tiff",width=6,height = 5,device = "tiff",dpi = 700)




Plk2
palette <- c("#FF7F00","white")

Plk1
palette <- c("#32CD32","#1E90FF","white")

Biofilm6
palette <- c("#E31A1C","#EE82EE","#FF7F50","#1F78B4","white")

Biofilm5
palette <- c("#E31A1C","#DAA520","#EE82EE","#1F78B4","#CD853F","white")

Biofilm4
palette <- c("#FDBF6F","#33A02C","#E31A1C","#1F78B4","white")

Biofilm3
palette <- c("#E31A1C","#1F78B4","#FB9A99","white")

Biofilm2
palette <- c("#FDBF6F",	"#4169E1","#228B22","white")

Biofilm1
palette <- c("#F4A460","#EE82EE",	"white")
