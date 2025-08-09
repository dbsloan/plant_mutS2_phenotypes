library(tidyverse)
library(car)

mutS2_data <- read.csv("discord.summary.mito.mod.csv")

mutS2_data$Genotype <- factor(mutS2_data$Genotype, levels = c("WT", "mutS2A/B"))
mutS2_data$Treatment <- factor(mutS2_data$Treatment, levels = c("Control", "Cipro"))


mutS2_data <- mutS2_data %>% mutate(DiscordantRatePerMM = 1e6*(same_neg_R1_lower_cut+same_neg_R1_higher_cut+same_pos_R1_lower_cut+same_pos_R1_higher_cut+outward_cut)/MappedPairs)

ggplot(data=mutS2_data, aes(y=DiscordantRatePerMM, x=Genotype)) +
  facet_wrap(~Treatment) +
  stat_summary(aes(group = Rep), geom = "line", fun = mean, color="gray") +
  geom_point(size=1.5, stroke = 0) +
  stat_summary(geom = "point", fun = "mean", size = 12, shape = 95, alpha=0.5, position=position_dodge(width = 0.8), mapping=aes(y=DiscordantRatePerMM, x=Genotype)) + theme_bw() +
#  ylim(0,1150) +
  ylab ("Mitogenome Discordant Read Pairs") +
  theme(panel.grid.minor = element_blank(), axis.title = element_text(size=7, face="bold"), axis.text = element_text(size=6), strip.text = element_text(size=7, margin = margin(t = 1.5, b = 1.5, l = 1.5, r = 1.5)))

ggsave ("Discord_plot.mito.pdf", width=3.25, height=2.5)

mutS2_data_cipro = mutS2_data %>%
  filter(Treatment == "Cipro") %>% 
  group_by(Rep, Genotype) %>%
  summarize(Grouped_DiscordantRatePerMM = sum(DiscordantRatePerMM), .groups = "drop") %>%
  pivot_wider(names_from = Genotype, values_from = Grouped_DiscordantRatePerMM)

t.test(data=mutS2_data_cipro, x=mutS2_data_cipro$WT, y=mutS2_data_cipro$`mutS2A/B`, paired=TRUE)


mutS2_data_control = mutS2_data %>%
  filter(Treatment == "Control") %>% 
  group_by(Rep, Genotype) %>%
  summarize(Grouped_DiscordantRatePerMM = sum(DiscordantRatePerMM), .groups = "drop") %>%
  pivot_wider(names_from = Genotype, values_from = Grouped_DiscordantRatePerMM)

t.test(data=mutS2_data_control, x=mutS2_data_control$WT, y=mutS2_data_control$`mutS2A/B`, paired=TRUE)
