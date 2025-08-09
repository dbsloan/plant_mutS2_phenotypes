library(tidyverse)
library(scales)

pos_data_mutS2 = read.table ("hotspots.mutS2.txt", header=TRUE)
pos_data_WT = read.table ("hotspots.WT.txt", header=TRUE)

pos_data <- bind_rows(
  pos_data_mutS2 %>% mutate(Genotype = "mutS2A/B"),
  pos_data_WT %>% mutate(Genotype = "WT")
)

ggplot(data=pos_data, aes(x=Position, y=Read_Count)) +
  geom_bar(stat="identity") +
  facet_wrap(~Genotype, nrow=2, scales="free_y") +
  theme_bw() +
  xlab ("Position in Plastid Genome (bp)") +
  ylab ("Discordant Read Count") +
  scale_x_continuous(labels = comma) +
  theme(panel.grid = element_blank(), axis.title = element_text(size=7, face="bold"), axis.text = element_text(size=6), strip.text = element_text(size=7, margin = margin(t = 1.5, b = 1.5, l = 1.5, r = 1.5)))

ggsave ("Hotspot_plot.pdf", width=3.25, height=3.25)