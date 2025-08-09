#(spreadsheet Feb2023_ERY)


library(lme4)
library(lmerTest)
library(pbkrtest)
library(emmeans)
library(multcomp)
library(multcompView)

#This code replaces the placeholders values for drug concentration with their actual values
Feb2023_ERY$DC = factor(Feb2023_ERY$DC, labels = c("Control", "Ethanol", "Ery 0.4", "Ery 2", "Ery 10"))
#This code replaces the BA value for genotype with AB
Feb2023_ERY$Genotype = factor(Feb2023_ERY$Genotype, labels = c("A", "B", "AB", "W"))

#analysis for FvFm
Model <- lmer(FvFm ~ Genotype*DC + (1|plate) + (1|F2), data = Feb2023_ERY)
anova(Model)
emmeans(Model, pairwise ~ Genotype|DC)

#Tukey significance groups
Model <- lmer(FvFm ~ Genotype*DC + (1|plate) + (1|F2), data = Feb2023_ERY)
FvFm_ee <- emmeans(Model, pairwise ~ Genotype|DC, adjust = "tukey")
FvFm_cld <- cld(FvFm_ee, alpha = 0.05, Letters = letters)
FvFm_cld



#analysis for Root length
Model <- lmer(Root ~ Genotype*DC + (1|plate) + (1|F2), data = Feb2023_ERY)
anova(Model)
emmeans(Model, pairwise ~ Genotype|DC)


#Tukey significance groups
Model <- lmer(Root ~ Genotype*DC + (1|plate) + (1|F2), data = Feb2023_ERY)
Root_ee <- emmeans(Model, pairwise ~ Genotype|DC, adjust = "tukey")
Root_cld <- cld(Root_ee, alpha = 0.05, Letters = letters)
Root_cld