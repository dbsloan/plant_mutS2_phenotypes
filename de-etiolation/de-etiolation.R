library(tidyverse)

greening = read.csv("de-etiolation.csv")

greening$Genotype = factor(greening$Genotype, levels=c("A", "B", "AB", "W"))

greening %>%
  ggplot(aes(x=Genotype, y=PercentGreen)) +
  geom_point (cex=0.75, pch=1.0,position=position_jitter(w=0.1,h=0), alpha=0.5) +
  stat_summary(fun.data = 'mean_se', geom = 'errorbar', width=0.2, linewidth=0.5) +
  stat_summary(fun.data = 'mean_se', geom = 'pointrange', cex=0.25) +
  facet_wrap(~Experiment, nrow=1) +
  theme_bw() +
  theme(legend.position = "none", 
        axis.title = element_text(size=7, face="bold"),
        axis.text = element_text(size=6),
        strip.text = element_text(size=6, margin = margin(1.5, 1.5, 1.5, 1.5))
        ) +
  ylab("Green Plants (%)")

ggsave("de-etiotation.pdf", width=4.5, height=1.8)