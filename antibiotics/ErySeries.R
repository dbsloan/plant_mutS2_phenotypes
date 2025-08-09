library(tidyverse)
library(patchwork)

Feb2023_ERY <- read.csv("ErySeries.csv")

#This code replaces the placeholders values for drug concentration with their actual values
Feb2023_ERY$DC = factor(Feb2023_ERY$DC, labels = c("Control", "Ethanol", "Ery 0.4", "Ery 2", "Ery 10"))
#This code replaces the BA value for genotype with AB
Feb2023_ERY$Genotype = factor(Feb2023_ERY$Genotype, labels = c("A", "B", "AB", "W"))


mean.FvFm.data <- Feb2023_ERY %>% 
  filter(!is.na(FvFm)) %>% 
  group_by(Genotype, DC) %>% 
  summarise(FvFm = mean(FvFm))

mean.Root.data <- Feb2023_ERY %>% 
  filter(!is.na(Root)) %>% 
  group_by(Genotype, DC) %>% 
  summarise(Root = mean(Root))



FvFm.plot <- ggplot(Feb2023_ERY, aes(x=Genotype, y=FvFm, group=DC)) +
  geom_point(cex=1.5,pch=1.0,position=position_jitter(w=0.1,h=0), alpha=0.5) +
  stat_summary(fun.data = 'mean_se', geom = 'errorbar', width=0.1) +
  stat_summary(fun.data = 'mean_se', geom = 'pointrange') +
#  geom_point(data=mean.FvFm.data, aes(x=Genotype, y=FvFm)) +
  facet_grid(~DC) +
  coord_cartesian(ylim = c(0, 0.9)) +
  ylab("Fv/Fm")

Root.plot <- ggplot(Feb2023_ERY, aes(x=Genotype, y=Root, group=DC)) +
  geom_point(cex=1.5,pch=1.0,position=position_jitter(w=0.1,h=0), alpha=0.5) +
  stat_summary(fun.data = 'mean_se', geom = 'errorbar', width=0.1) +
  stat_summary(fun.data = 'mean_se', geom = 'pointrange') +
#  geom_point(data=mean.Root.data, aes(x=Genotype, y=Root)) +
  facet_grid(~DC) +
  ylab("Root Length (mm)")

#Ths uses the patchwork library to combine the two saved plots in vertical orientation (Fv/Fm on top of root)
combined.plot <- FvFm.plot / Root.plot


#This code reformats the table in a way that will allow for making FvFm and root plots with a single ggplot call
Feb2023_ERY.long <- Feb2023_ERY %>%
  pivot_longer(
    cols = c(FvFm, Root),
    names_to = "Phenotype",
    values_to = "PhenotypeValue"
  )

#Rename Root to Root Length (mm) and FvFm to Fv/Fm
Feb2023_ERY.long$Phenotype = factor(Feb2023_ERY.long$Phenotype, labels = c("Fv/Fm", "Root Length (mm)"))

mean.phenotype.data <- Feb2023_ERY.long %>% 
  filter(!is.na(PhenotypeValue)) %>% 
  group_by(Genotype, DC, Phenotype) %>% 
  summarise(PhenotypeValue = mean(PhenotypeValue))


#This is a hack to allow for both manually specifying the y axis limits and have them different for the different rows in the grid (ggplot lets you do one of those but not both)
#it is implemented with teh geom_blank lines in the ggplot call below
Feb2023_ERY.long <- Feb2023_ERY.long %>%
  mutate(ymin = case_when(
    Phenotype == "Fv/Fm" ~ 0,
    Phenotype == "Root Length (mm)" ~ 0,
    TRUE ~ NA_real_
  ),
  ymax = case_when(
    Phenotype == "Fv/Fm" ~ 1,
    Phenotype == "Root Length (mm)" ~ 60,
    TRUE ~ NA_real_
  ))


combined.plot.grid <- ggplot(Feb2023_ERY.long, aes(x=Genotype, y=PhenotypeValue)) +
  geom_point(cex=1,pch=1,position=position_jitter(w=0.1,h=0), alpha=0.25) +
  stat_summary(fun.data = 'mean_se', geom = 'errorbar', width=0.4, linewidth=0.25) +
  stat_summary(fun.data = 'mean_se', geom = 'pointrange', size=0.2, linewidth=0.25) +
#  geom_point(data=mean.phenotype.data, aes(x=Genotype, y=PhenotypeValue)) +
  facet_grid(Phenotype~DC, scales="free_y") +
  geom_blank(aes(y = ymin)) +  #part of the hack to manually scale y axis      
  geom_blank(aes(y = ymax)) +  #part of the hack to manually scale y axis   
  labs(y = NULL) +
  theme_bw() +
  theme(
    axis.text = element_text(size=6),
    axis.title = element_text(size=7, face="bold"),
    strip.text = element_text(size=7, margin = margin(t = 1.5, r = 1.5, b = 1.5, l = 1.5))
  )

combined.plot.grid

ggsave("ErySeries.pdf", width=4, height=2.5)
