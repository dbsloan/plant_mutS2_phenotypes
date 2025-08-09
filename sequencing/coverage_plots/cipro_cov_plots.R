library (tidyverse)
library (scales)

cov_data_mito = read.table("depth.mito.sw.mod.txt", header=TRUE)
cov_data_plastid = read.table("depth.plastid.sw.mod.txt", header=TRUE)

cov_data_mito_long = cov_data_mito %>%
  pivot_longer(
    cols = c(C10_Mut_Control_3, C11_WT_Cipro_3, C12_Mut_Cipro_3, 
             C17_WT_Control_4, C18_Mut_Control_4, C19_WT_Cipro_4, 
             C01_WT_Control_1, C20_Mut_Cipro_4, C21_WT_Control_5, 
             C22_Mut_Control_5, C23_WT_Cipro_5, C24_Mut_Cipro_5, 
             C02_Mut_Control_1, C03_WT_Cipro_1, C04_Mut_Cipro_1, 
             C05_WT_Control_2, C06_Mut_Control_2, C07_WT_Cipro_2, 
             C08_Mut_Cipro_2, C09_WT_Control_3), 
    names_to = "Library",
    values_to = "Coverage"
  )


cov_data_plastid_long = cov_data_plastid %>%
  pivot_longer(
    cols = c(C10_Mut_Control_3, C11_WT_Cipro_3, C12_Mut_Cipro_3, 
             C17_WT_Control_4, C18_Mut_Control_4, C19_WT_Cipro_4, 
             C01_WT_Control_1, C20_Mut_Cipro_4, C21_WT_Control_5, 
             C22_Mut_Control_5, C23_WT_Cipro_5, C24_Mut_Cipro_5, 
             C02_Mut_Control_1, C03_WT_Cipro_1, C04_Mut_Cipro_1, 
             C05_WT_Control_2, C06_Mut_Control_2, C07_WT_Cipro_2, 
             C08_Mut_Cipro_2, C09_WT_Control_3), 
    names_to = "Library",
    values_to = "Coverage"
  )


ggplot (data = cov_data_mito_long, aes(x=Position/1000, y=Coverage)) + 
  geom_bar(stat="identity") +
  facet_wrap(~Library, ncol=4) +
  theme_bw() +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma) +
  coord_cartesian(ylim = c(0, 5000)) +
  xlab ("Position (kb)") +
  ylab ("Read Depth") +
  theme(legend.position = "none", panel.grid = element_blank())

ggsave ("mito_cipro_coverage.pdf", height=6.7, width=6.7)


ggplot (data = cov_data_plastid_long, aes(x=Position/1000, y=Coverage)) + 
  geom_bar(stat="identity") +
  facet_wrap(~Library, ncol=4) +
  theme_bw() +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma) +
  coord_cartesian(ylim = c(0, 50000)) +
  xlab ("Position (kb)") +
  ylab ("Read Depth") +
  theme(legend.position = "none", panel.grid = element_blank())

ggsave ("plastid_cipro_coverage.pdf", height=6.7, width=6.7)

