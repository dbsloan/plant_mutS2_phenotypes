#spreadsheet is X2023_April_all_drugs

library(lme4)
library(lmerTest)
library(pbkrtest)
library(emmeans)
library(multcomp)
library(multcompView)

#This code replaces the placeholders values for drug concentration with their actual values
X2023_April_all_drugs$DC = factor(X2023_April_all_drugs$DC, labels = c("Control", "Ethanol", "Chl 0.4", "Chl 2.0", "Ery 0.4", "Ery 2.0", "Hyg 0.4","Hyg 2.0", "Spec 0.4", "Spec 2.0"))
#This code replaces the BA value for genotype with AB
X2023_April_all_drugs$Genotype = factor(X2023_April_all_drugs$Genotype, labels = c("A", "B", "AB", "W"))


#FvFm analysis
Model <- lmer(FvFm ~ Genotype*DC + (1|Plate) + (1|Block) + (1|F2), data = X2023_April_all_drugs)
anova(Model)
emmeans(Model, pairwise ~ Genotype|DC)

#Tukey significance groups
Model <- lmer(FvFm ~ Genotype*DC + (1|Plate) + (1|Block) + (1|F2), data = X2023_April_all_drugs)
FvFm_em <- emmeans(Model, pairwise ~ Genotype|DC, adjust = "tukey")
FvFm_cld <- cld(FvFm_em, alpha = 0.05, Letters = letters)
FvFm_cld


#Root analysis
Model <- lmer(Root ~ Genotype*DC + (1|Plate) + (1|Block) + (1|F2), data = X2023_April_all_drugs)
anova(Model)
emmeans(Model, pairwise ~ Genotype|DC)

#Tukey significance groups
Model <- lmer(Root ~ Genotype*DC + (1|Plate) + (1|Block) + (1|F2), data = X2023_April_all_drugs)
FvFm_em <- emmeans(Model, pairwise ~ Genotype|DC, adjust = "tukey")
FvFm_cld <- cld(FvFm_em, alpha = 0.05, Letters = letters)
FvFm_cld