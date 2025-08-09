library(tidyverse)

overlap = read.table("overlap_freq.txt", header=TRUE)

overlap$Genotype = factor(overlap$Genotype, levels = c("WT", "mutS2A/B"))
overlap$Treatment = factor(overlap$Treatment, levels = c("Control", "Cipro"))

ggplot(data=overlap, aes(x=Overlap, y=Frequency)) +
  geom_bar(stat="identity", fill="gray50") +
  geom_line(data=overlap,aes(y=Simulation), color="firebrick3") +
  facet_grid(Treatment ~ Genotype) +
  theme_bw() +
  xlim(c(-1,31)) +
  ylim(c(0,.42)) +
  xlab ("Sequence Overlap at Breakpoints (bp)") +
  theme(panel.grid.minor = element_blank(), axis.title = element_text(size=7, face="bold"), axis.text = element_text(size=6), strip.text = element_text(size=7, margin = margin(t = 1.5, b = 1.5, l = 1.5, r = 1.5)))

ggsave ("Overlap_plot.pdf", width=3.25, height=2.5)