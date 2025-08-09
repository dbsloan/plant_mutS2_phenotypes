###R code adapted from Warren et al. 2023 MBE. doi: 10.1093/molbev/msad163
###Modified by Shady Kuster

## Generate trees of aaRSs with targeting predictions ## 

library (tidyverse)
library(ggtree)
library(aplot)
library(treeio)


#read in tree, root it, and make ggtree object
mutS2_tree = read.tree ("../RAxML_bipartitions.nwk")
outgroup <- c("Bacillus_subtilis_subsp_CAB14818.1", "Helicobacter_pylori_AAD07685.1",
             "Borreliella_burgdorferi_AAC66481.1","Aquifex_aeolicus_AAC07247.1")
mutS2_tree <- root(ape::unroot(mutS2_tree), outgroup)

tree = ggtree(mutS2_tree) + 
  geom_tiplab() + 
  geom_treescale(offset = 0.1, x = 0, y=6) + 
  geom_text2(aes(subset = !isTip, label=label))

#read in prediction data and extract the mt and plastid data
target_dat = read.table("targeting_predictions.txt", header=TRUE, sep = "\t")

mito = target_dat %>% filter(Organelle=="Mitochondria")
mito[,4] <- round(mito[,4], 2)

chloro = target_dat %>% filter(Organelle=="Chloroplast")
chloro[,4] <- round(chloro[,4], 2)

#make the mt plot
mplot = ggplot(mito, aes(x=Algorithm, y=ID, fill=value)) +
  geom_tile(width=.95, height=.95, color="gray", size = 0.5) + 
  theme_tree2() +
  scale_fill_gradient(limits=c(0,1), low="white", high="chocolate4", name="Mitochondrial Targeting") +
  theme(legend.position = "right", legend.direction = "horizontal") +
  geom_text(aes(label=value), size = 2)
mplot 

#make the plastid plot
cplot = ggplot(chloro, aes(x=Algorithm, y=ID, fill=value)) +
  geom_tile(width=.95, height=.95, color="gray", size = 0.25) + 
  theme_tree2() + 
  scale_fill_gradient(limits=c(0,1), low="white", high="darkolivegreen4", name = "Plastid Targeting") +
  theme(legend.position = "right", legend.direction = "horizontal") +
  geom_text(aes(label=value), size = 2)
cplot

#combine the plots and make pdf
plot<-mplot %>% insert_left(tree) %>% insert_right(cplot)
plot
pdf(file="Targeting.pdf", width=10, height=6)
  print(plot)
graphics.off()
