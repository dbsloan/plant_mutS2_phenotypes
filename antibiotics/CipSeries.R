library(tidyverse)


CIP_alldata3 <- read.csv("CipSeries.csv")

mean.aveLeaf.data <- CIP_alldata3 %>% group_by(Genotype, DC) %>% summarise(aveLeaf = mean(aveLeaf))

#This code replaces the treatment names so that concentration comes after CIP to match order in other plots
CIP_alldata3$DC = factor(CIP_alldata3$DC, labels = c("Control", "CIP 0.25", "CIP 0.5", "CIP 0.75", "CIP 1.0"))


ggplot(CIP_alldata3, aes(x=Genotype, y=aveLeaf, group=DC)) +
  geom_point(cex=0.8,pch=1,position=position_jitter(w=0.1,h=0), alpha=0.25) +
  stat_summary(fun.data = 'mean_se', geom = 'errorbar', width=0.4, linewidth=0.25) +
  stat_summary(fun.data = 'mean_se', geom = 'pointrange', size=0.2, linewidth=0.25) +
#  geom_point(data=mean.aveLeaf.data, aes(x=Genotype, y=aveLeaf)) +
  facet_grid(~DC) +
  scale_x_discrete(limits = c("A","B","AB", "W")) +
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100), limits = c(0, 115)) +
  ylab("Percent with True Leaves") +
  theme_bw() +
  theme(
    axis.text = element_text(size=6),
    axis.title = element_text(size=7, face="bold"),
    strip.text = element_text(size=7, margin = margin(t = 1.5, r = 1.5, b = 1.5, l = 1.5))
  )

ggsave("CipSeries.pdf", width=4, height=1.6)
