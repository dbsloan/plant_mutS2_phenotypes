#spreadsheet (2023_Sept_rescue)
library(readxl)
X2023_Sept_rescue <- read_excel("Sloan lab/NIH grant/MutS2/Manuscript/Code_spreadsheets/Rescue_Sept 2023/2023_Sept_rescue.xlsx", 
                                col_types = c("numeric", "text", "text", 
                                              "text", "text", "text", "text", "text", "text", "numeric", "numeric"))
View(X2023_Sept_rescue)


library(lme4)
library(lmerTest)
library(pbkrtest)
library(emmeans)
library(multcomp)
library(multcompView)


#This code replaces the placeholders values for drug concentration with their actual values
X2023_Sept_rescue$DC = factor(X2023_Sept_rescue$DC, labels = c("Control", "Ethanol", "Chl 2.0", "Ery 2.0", "Spec 2.0"))
#This code replaces the BA value for genotype with AB
X2023_Sept_rescue$Genotype = factor(X2023_Sept_rescue$Genotype, levels = c("A", "B", "AB", "W", "AC1", "AC2", "BC1", "BC2"))


#FvFm analysis
Model <- lmer(FvFm ~ Genotype*DC + (1|plate) + (1|block), data = X2023_Sept_rescue)
anova(Model)
emmeans(Model, pairwise ~ Genotype|DC)

#Tukey significance groups
Model <- lmer(FvFm ~ Genotype*DC + (1|plate) + (1|block), data = X2023_Sept_rescue)
FvFm_ee <- emmeans(Model, pairwise ~ Genotype|DC, adjust = "tukey")
FvFm_cld <- cld(FvFm_ee, alpha = 0.05, Letters = letters)
FvFm_cld


#Root analysis
Model <- lmer(root_mm ~ Genotype*DC + (1|plate) + (1|block), data = X2023_Sept_rescue)
anova(Model)
emmeans(Model, pairwise ~ Genotype|DC)

#Tukey significance groups
Model <- lmer(root_mm ~ Genotype*DC + (1|plate) + (1|block), data = X2023_Sept_rescue)
root_ee <- emmeans(Model, pairwise ~ Genotype|DC, adjust = "tukey")
root_cld <- cld(root_ee, alpha = 0.05, Letters = letters)
root_cld