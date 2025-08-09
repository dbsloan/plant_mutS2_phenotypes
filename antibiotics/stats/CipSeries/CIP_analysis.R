#spreadsheet CIP_alldata3
library(readxl)
CIP_alldata3 <- read_excel("Sloan lab/NIH grant/MutS2/Manuscript/Code_spreadsheets/CIP/CIP_alldata3.xlsx", 
                           col_types = c("text", "text", "text", 
                                         "text", "text", "numeric", "numeric", 
                                         "text", "numeric"))
View(CIP_alldata3)


library(lme4)
library(lmerTest)
library(pbkrtest)
library(emmeans)
library(multcomp)
library(multcompView)

Model <- lmer(aveLeaf ~ Genotype*DC + (1|Exp) + (1|ExpRep), data = CIP_alldata3)
anova(Model)
emmeans(Model, pairwise ~ Genotype|DC)

#use emmeans with the tukey test letters
Model <- lmer(aveLeaf ~ Genotype*DC + (1|Exp) + (1|ExpRep), data = CIP_alldata3)
Leaf_ee <- emmeans(Model, pairwise ~ Genotype|DC, adjust = "tukey")
Leaf_cld <- cld(Leaf_ee, alpha = 0.05, Letters = letters)
Leaf_cld