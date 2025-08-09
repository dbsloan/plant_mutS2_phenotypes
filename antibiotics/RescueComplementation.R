library(tidyverse)
library(patchwork)

X2023_Sept_rescue <- read.csv("RescueComplementation.csv")

#This code replaces the placeholders values for drug concentration with their actual values
X2023_Sept_rescue$DC = factor(X2023_Sept_rescue$DC, labels = c("Control", "Ethanol", "Chl 2.0", "Ery 2.0", "Spec 2.0"))
#This code replaces the BA value for genotype with AB
X2023_Sept_rescue$Genotype = factor(X2023_Sept_rescue$Genotype, levels = c("A", "B", "AB", "W", "AC1", "AC2", "BC1", "BC2"))


FvFm.plot <- ggplot(X2023_Sept_rescue, aes(x=Genotype, y=FvFm, group=DC)) +
  geom_point(cex=1.5,pch=1.0,position=position_jitter(w=0.1,h=0), alpha=0.5) +
  stat_summary(fun.data = 'mean_se', geom = 'errorbar', width=0.1) +
  stat_summary(fun.data = 'mean_se', geom = 'pointrange') +
  #  geom_point(data=mean.FvFm.data, aes(x=Genotype, y=FvFm)) +
  facet_grid(~DC) +
  ylab("Fv/Fm")


Root.plot <- ggplot(X2023_Sept_rescue, aes(x=Genotype, y=root_mm, group=DC)) +
  geom_point(cex=1.5,pch=1.0,position=position_jitter(w=0.1,h=0), alpha=0.5) +
  stat_summary(fun.data = 'mean_se', geom = 'errorbar', width=0.1) +
  stat_summary(fun.data = 'mean_se', geom = 'pointrange') +
  #  geom_point(data=mean.Root.data, aes(x=Genotype, y=root_mm)) +
  facet_grid(~DC) +
  ylab("Root Length (mm)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

#Ths uses the patchwork library to combine the two saved plots in vertical orientation (Fv/Fm on top of root)
combined.plot <- FvFm.plot / Root.plot

#This code reformats the table in a way that will allow for making FvFm and root plots with a single ggplot call
X2023_Sept_rescue.long <- X2023_Sept_rescue %>%
  pivot_longer(
    cols = c(FvFm, root_mm),
    names_to = "Phenotype",
    values_to = "PhenotypeValue"
  )

#Rename Root to Root Length (mm) and FvFm to Fv/Fm
X2023_Sept_rescue.long$Phenotype = factor(X2023_Sept_rescue.long$Phenotype, labels = c("Fv/Fm", "Root Length (mm)"))

#calculating averages but this has been commented out below
mean.phenotype.data <- X2023_Sept_rescue.long %>% 
  filter(!is.na(PhenotypeValue)) %>% 
  group_by(Genotype, DC, Phenotype) %>% 
  summarise(PhenotypeValue = mean(PhenotypeValue))


#This is a hack to allow for both manually specifying the y axis limits and have them different for the different rows in the grid (ggplot lets you do one of those but not both)
#it is implemented with teh geom_blank lines in the ggplot call below
X2023_Sept_rescue.long <- X2023_Sept_rescue.long %>%
  mutate(ymin = case_when(
    Phenotype == "Fv/Fm" ~ 0,
    Phenotype == "Root Length (mm)" ~ 0,
    TRUE ~ NA_real_
  ),
  ymax = case_when(
    Phenotype == "Fv/Fm" ~ 1,
    Phenotype == "Root Length (mm)" ~ 35,
    TRUE ~ NA_real_
  ))


combined.plot.grid <- ggplot(X2023_Sept_rescue.long, aes(x=Genotype, y=PhenotypeValue)) +
  geom_point(cex=0.8,pch=1,position=position_jitter(w=0.1,h=0), alpha=0.25) +
  stat_summary(fun.data = 'mean_se', geom = 'errorbar', width=0.4, linewidth=0.25) +
  stat_summary(fun.data = 'mean_se', geom = 'pointrange', size=0.125, linewidth=0.25) +
  #  geom_point(data=mean.phenotype.data, aes(x=Genotype, y=PhenotypeValue)) +
  facet_grid(Phenotype~DC, scales="free_y") +
  geom_blank(aes(y = ymin)) +  #part of the hack to manually scale y axis      
  geom_blank(aes(y = ymax)) +  #part of the hack to manually scale y axis   
  labs(y = NULL) +
  theme_bw() +
  theme(
  axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
  axis.text = element_text(size=6),
  axis.title = element_text(size=7, face="bold"),
  strip.text = element_text(size=7, margin = margin(t = 2, r = 2, b = 2, l = 2))
)

combined.plot.grid

ggsave("RescueComplementation.pdf", width=6, height=2.5)