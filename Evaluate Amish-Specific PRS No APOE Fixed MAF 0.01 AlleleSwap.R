### Evaluate Amish-Specific PRS split by site (IN and OH)

library(ROCR)
library(pROC)
library(DescTools)
library(ggpubr)
library(cowplot)
library(tidyverse)

setwd("D:/Research/AD PRS/PRS")

OH_pheno = read.table("D:/Research/AD PRS/GWAS/GENESIS/PV_Data_AD_OH.txt", header=T)
IN_pheno = read.table("D:/Research/AD PRS/GWAS/GENESIS/PV_Data_AD_IN.txt", header=T)

OH_PRS = read.table("IN_Train_OH_Test_MAF01_No_APOE_AlleleSwap.best", header=T)
IN_PRS = read.table("OH_Train_IN_Test_MAF01_No_APOE_AlleleSwap.best", header=T)

OH_PRS = OH_PRS[match(OH_pheno$IID, OH_PRS$IID),]
IN_PRS = IN_PRS[match(IN_pheno$IID, IN_PRS$IID),]

OH_pheno$PRS = OH_PRS$PRS
IN_pheno$PRS = IN_PRS$PRS

table(OH_pheno$dx)
table(IN_pheno$dx)

## Making APOE Count Variables

# OH
for (i in 1:nrow(OH_pheno)){
  if (is.na(OH_pheno$APOE[i])) {OH_pheno$apoe2[i] = NA}
  else if (OH_pheno$APOE[i] == "2|2") {OH_pheno$apoe2[i] = 2}
  else if (OH_pheno$APOE[i] == "2|3") {OH_pheno$apoe2[i] = 1}
  else if (OH_pheno$APOE[i] == "2|4") {OH_pheno$apoe2[i] = 1}
  else if (OH_pheno$APOE[i] == "3|3") {OH_pheno$apoe2[i] = 0}
  else if (OH_pheno$APOE[i] == "3|4") {OH_pheno$apoe2[i] = 0}
  else if (OH_pheno$APOE[i] == "4|4") {OH_pheno$apoe2[i] = 0}
}

for (i in 1:nrow(OH_pheno)){
  if (is.na(OH_pheno$APOE[i])) {OH_pheno$apoe4[i] = NA}
  else if (OH_pheno$APOE[i] == "2|2") {OH_pheno$apoe4[i] = 0}
  else if (OH_pheno$APOE[i] == "2|3") {OH_pheno$apoe4[i] = 0}
  else if (OH_pheno$APOE[i] == "2|4") {OH_pheno$apoe4[i] = 1}
  else if (OH_pheno$APOE[i] == "3|3") {OH_pheno$apoe4[i] = 0}
  else if (OH_pheno$APOE[i] == "3|4") {OH_pheno$apoe4[i] = 1}
  else if (OH_pheno$APOE[i] == "4|4") {OH_pheno$apoe4[i] = 2}
}

table(OH_pheno$APOE, OH_pheno$apoe2)
table(OH_pheno$APOE, OH_pheno$apoe4)

# IN

for (i in 1:nrow(IN_pheno)){
  if (is.na(IN_pheno$APOE[i])) {IN_pheno$apoe2[i] = NA}
  else if (IN_pheno$APOE[i] == "2|2") {IN_pheno$apoe2[i] = 2}
  else if (IN_pheno$APOE[i] == "2|3") {IN_pheno$apoe2[i] = 1}
  else if (IN_pheno$APOE[i] == "2|4") {IN_pheno$apoe2[i] = 1}
  else if (IN_pheno$APOE[i] == "3|3") {IN_pheno$apoe2[i] = 0}
  else if (IN_pheno$APOE[i] == "3|4") {IN_pheno$apoe2[i] = 0}
  else if (IN_pheno$APOE[i] == "4|4") {IN_pheno$apoe2[i] = 0}
}

for (i in 1:nrow(IN_pheno)){
  if (is.na(IN_pheno$APOE[i])) {IN_pheno$apoe4[i] = NA}
  else if (IN_pheno$APOE[i] == "2|2") {IN_pheno$apoe4[i] = 0}
  else if (IN_pheno$APOE[i] == "2|3") {IN_pheno$apoe4[i] = 0}
  else if (IN_pheno$APOE[i] == "2|4") {IN_pheno$apoe4[i] = 1}
  else if (IN_pheno$APOE[i] == "3|3") {IN_pheno$apoe4[i] = 0}
  else if (IN_pheno$APOE[i] == "3|4") {IN_pheno$apoe4[i] = 1}
  else if (IN_pheno$APOE[i] == "4|4") {IN_pheno$apoe4[i] = 2}
}

table(IN_pheno$APOE, IN_pheno$apoe2)
table(IN_pheno$APOE, IN_pheno$apoe4)

## OH Models

# Covariate only model

cov_mod_OH = glm(dx ~ sex + Age, data=OH_pheno)
summary(cov_mod_OH)
Cstat(cov_mod_OH)

# APOE only model

apoe_mod_OH = glm(dx ~ APOE, data=OH_pheno)
summary(apoe_mod_OH)
Cstat(apoe_mod_OH)

apoe_count_mod_OH  = glm(dx ~ apoe2 + apoe4, data=OH_pheno)
summary(apoe_count_mod_OH)
Cstat(apoe_count_mod_OH)

# Cov + APOE only model

cov_apoe_mod_OH  = glm(dx ~ sex + Age + apoe2 + apoe4, data=OH_pheno)
summary(cov_apoe_mod_OH)
Cstat(cov_apoe_mod_OH)

# PRS models

prs_mod_OH = glm(dx ~ PRS, data=OH_pheno)
summary(prs_mod_OH)
Cstat(prs_mod_OH)

prs_apoe_mod_OH = glm(dx ~ PRS + apoe2 + apoe4, data=OH_pheno)
summary(prs_apoe_mod_OH)
Cstat(prs_apoe_mod_OH)

sex_age_prs_mod_OH = glm(dx ~ sex + Age + PRS, data=OH_pheno)
summary(sex_age_prs_mod_OH)
Cstat(sex_age_prs_mod_OH)

full_mod_OH = glm(dx ~ sex + Age + PRS + apoe2 + apoe4, data=OH_pheno)
summary(full_mod_OH)
Cstat(full_mod_OH)

## IN Models

# Covariate only model

cov_mod_IN = glm(dx ~ sex + Age, data=IN_pheno)
summary(cov_mod_IN)
Cstat(cov_mod_IN)

# APOE only model

apoe_mod_IN = glm(dx ~ APOE, data=IN_pheno)
summary(apoe_mod_IN)
Cstat(apoe_mod_IN)

apoe_count_mod_IN  = glm(dx ~ apoe2 + apoe4, data=IN_pheno)
summary(apoe_count_mod_IN)
Cstat(apoe_count_mod_IN)

# Cov + APOE only model

cov_apoe_mod_IN  = glm(dx ~ sex + Age + apoe2 + apoe4, data=IN_pheno)
summary(cov_apoe_mod_IN)
Cstat(cov_apoe_mod_IN)

# PRS models

prs_mod_IN = glm(dx ~ PRS, data=IN_pheno)
summary(prs_mod_IN)
Cstat(prs_mod_IN)

prs_apoe_mod_IN = glm(dx ~ PRS + apoe2 + apoe4, data=IN_pheno)
summary(prs_apoe_mod_IN)
Cstat(prs_apoe_mod_IN)

sex_age_prs_mod_IN = glm(dx ~ sex + Age + PRS, data=IN_pheno)
summary(sex_age_prs_mod_IN)
Cstat(sex_age_prs_mod_IN)

full_mod_IN = glm(dx ~ sex + Age + PRS + apoe2 + apoe4, data=IN_pheno)
summary(full_mod_IN)
Cstat(full_mod_IN)


## Other Follow-up Information

boxplot(OH_pheno$PRS ~ OH_pheno$dx)
boxplot(IN_pheno$PRS ~ IN_pheno$dx)

OH_pheno_cases = subset(OH_pheno, dx == 1)
OH_pheno_controls = subset(OH_pheno, dx == 0)
summary(OH_pheno_cases$PRS)
summary(OH_pheno_controls$PRS)
t.test(OH_pheno_cases$PRS, OH_pheno_controls$PRS)

IN_pheno_cases = subset(IN_pheno, dx == 1)
IN_pheno_controls = subset(IN_pheno, dx == 0)
summary(IN_pheno_cases$PRS)
summary(IN_pheno_controls$PRS)
t.test(IN_pheno_cases$PRS, IN_pheno_controls$PRS)


### Rerun analysis with 1:1 design

## OH

set.seed(12345)
OH_pheno_controls_small = OH_pheno_controls[sample(nrow(OH_pheno_controls), 78),]
OH_pheno_even = rbind(OH_pheno_controls_small, OH_pheno_cases)

# Covariate only model

cov_mod_OH_even = glm(dx ~ sex + Age, data=OH_pheno_even)
summary(cov_mod_OH_even)
Cstat(cov_mod_OH_even)

# APOE only model

apoe_mod_OH_even = glm(dx ~ APOE, data=OH_pheno_even)
summary(apoe_mod_OH_even)
Cstat(apoe_mod_OH_even)

apoe_count_mod_OH_even  = glm(dx ~ apoe2 + apoe4, data=OH_pheno_even)
summary(apoe_count_mod_OH_even)
Cstat(apoe_count_mod_OH_even)

# Cov + APOE model

cov_apoe_mod_OH_even  = glm(dx ~ sex + Age + apoe2 + apoe4, data=OH_pheno_even)
summary(cov_apoe_mod_OH_even)
Cstat(cov_apoe_mod_OH_even)

# PRS models

prs_mod_OH_even = glm(dx ~ PRS, data=OH_pheno_even)
summary(prs_mod_OH_even)
Cstat(prs_mod_OH_even)

prs_apoe_mod_OH_even = glm(dx ~ PRS + apoe2 + apoe4, data=OH_pheno_even)
summary(prs_apoe_mod_OH_even)
Cstat(prs_apoe_mod_OH_even)

sex_age_prs_mod_OH_even = glm(dx ~ sex + Age + PRS, data=OH_pheno_even)
summary(sex_age_prs_mod_OH_even)
Cstat(sex_age_prs_mod_OH_even)

full_mod_OH_even = glm(dx ~ sex + Age + PRS + apoe2 + apoe4, data=OH_pheno_even)
summary(full_mod_OH_even)
Cstat(full_mod_OH_even)

## IN

set.seed(12345)
IN_pheno_controls_small = IN_pheno_controls[sample(nrow(IN_pheno_controls), 78),]
IN_pheno_even = rbind(IN_pheno_controls_small, IN_pheno_cases)

# Covariate only model

cov_mod_IN_even = glm(dx ~ sex + Age, data=IN_pheno_even)
summary(cov_mod_IN_even)
Cstat(cov_mod_IN_even)

# APOE only model

apoe_mod_IN_even = glm(dx ~ APOE, data=IN_pheno_even)
summary(apoe_mod_IN_even)
Cstat(apoe_mod_IN_even)

apoe_count_mod_IN_even  = glm(dx ~ apoe2 + apoe4, data=IN_pheno_even)
summary(apoe_count_mod_IN_even)
Cstat(apoe_count_mod_IN_even)

# Cov + APOE model

cov_apoe_mod_IN_even  = glm(dx ~ sex + Age + apoe2 + apoe4, data=IN_pheno_even)
summary(cov_apoe_mod_IN_even)
Cstat(cov_apoe_mod_IN_even)

# PRS models

prs_mod_IN_even = glm(dx ~ PRS, data=IN_pheno_even)
summary(prs_mod_IN_even)
Cstat(prs_mod_IN_even)

prs_apoe_mod_IN_even = glm(dx ~ PRS + apoe2 + apoe4, data=IN_pheno_even)
summary(prs_apoe_mod_IN_even)
Cstat(prs_apoe_mod_IN_even)

sex_age_prs_mod_IN_even = glm(dx ~ sex + Age + PRS, data=IN_pheno_even)
summary(sex_age_prs_mod_IN_even)
Cstat(sex_age_prs_mod_IN_even)

full_mod_IN_even = glm(dx ~ sex + Age + PRS + apoe2 + apoe4, data=IN_pheno_even)
summary(full_mod_IN_even)
Cstat(full_mod_IN_even)

### Rerun analysis with lower age cutoff design

## OH

OH_pheno_age75 = subset(OH_pheno, Age >= 75)

# Covariate only model

cov_mod_OH_age75 = glm(dx ~ sex + Age, data=OH_pheno_age75)
summary(cov_mod_OH_age75)
Cstat(cov_mod_OH_age75)

# APOE only model

apoe_mod_OH_age75 = glm(dx ~ APOE, data=OH_pheno_age75)
summary(apoe_mod_OH_age75)
Cstat(apoe_mod_OH_age75)

apoe_count_mod_OH_age75  = glm(dx ~ apoe2 + apoe4, data=OH_pheno_age75)
summary(apoe_count_mod_OH_age75)
Cstat(apoe_count_mod_OH_age75)

# Cov + APOE only model

cov_apoe_mod_OH_age75  = glm(dx ~ apoe2 + apoe4 + sex + Age, data=OH_pheno_age75)
summary(cov_apoe_mod_OH_age75)
Cstat(cov_apoe_mod_OH_age75)

# PRS models

prs_mod_OH_age75 = glm(dx ~ PRS, data=OH_pheno_age75)
summary(prs_mod_OH_age75)
Cstat(prs_mod_OH_age75)

prs_apoe_mod_OH_age75 = glm(dx ~ PRS + apoe2 + apoe4, data=OH_pheno_age75)
summary(prs_apoe_mod_OH_age75)
Cstat(prs_apoe_mod_OH_age75)

sex_age_prs_mod_OH_age75 = glm(dx ~ sex + Age + PRS, data=OH_pheno_age75)
summary(sex_age_prs_mod_OH_age75)
Cstat(sex_age_prs_mod_OH_age75)

full_mod_OH_age75 = glm(dx ~ sex + Age + PRS + apoe2 + apoe4, data=OH_pheno_age75)
summary(full_mod_OH_age75)
Cstat(full_mod_OH_age75)

## IN


IN_pheno_age75 = subset(IN_pheno, Age >= 75)

# Covariate only model

cov_mod_IN_age75 = glm(dx ~ sex + Age, data=IN_pheno_age75)
summary(cov_mod_IN_age75)
Cstat(cov_mod_IN_age75)

# APOE only model

apoe_mod_IN_age75 = glm(dx ~ APOE, data=IN_pheno_age75)
summary(apoe_mod_IN_age75)
Cstat(apoe_mod_IN_age75)

apoe_count_mod_IN_age75  = glm(dx ~ apoe2 + apoe4, data=IN_pheno_age75)
summary(apoe_count_mod_IN_age75)
Cstat(apoe_count_mod_IN_age75)

# Cov + APOE model

cov_apoe_mod_IN_age75  = glm(dx ~ apoe2 + apoe4 + sex + Age, data=IN_pheno_age75)
summary(cov_apoe_mod_IN_age75)
Cstat(cov_apoe_mod_IN_age75)

# PRS models

prs_mod_IN_age75 = glm(dx ~ PRS, data=IN_pheno_age75)
summary(prs_mod_IN_age75)
Cstat(prs_mod_IN_age75)

prs_apoe_mod_IN_age75 = glm(dx ~ PRS + apoe2 + apoe4, data=IN_pheno_age75)
summary(prs_apoe_mod_IN_age75)
Cstat(prs_apoe_mod_IN_age75)

sex_age_prs_mod_IN_age75 = glm(dx ~ sex + Age + PRS, data=IN_pheno_age75)
summary(sex_age_prs_mod_IN_age75)
Cstat(sex_age_prs_mod_IN_age75)

full_mod_IN_age75 = glm(dx ~ sex + Age + PRS + apoe2 + apoe4, data=IN_pheno_age75)
summary(full_mod_IN_age75)
Cstat(full_mod_IN_age75)

#### Plots for poster
boxplot(OH_pheno$PRS ~ OH_pheno$dx)
boxplot(IN_pheno$PRS ~ IN_pheno$dx)

OH_pheno$PRS2 = OH_pheno$PRS - mean(OH_pheno$PRS)
summary(OH_pheno$PRS2)
OH_pheno$PRS2 = OH_pheno$PRS2 * 1/sd(OH_pheno$PRS2)
sd(OH_pheno$PRS2)
summary(OH_pheno$PRS2)

IN_pheno$PRS2 = IN_pheno$PRS - mean(IN_pheno$PRS)
summary(IN_pheno$PRS2)
IN_pheno$PRS2 = IN_pheno$PRS2 * 1/sd(IN_pheno$PRS2)
sd(IN_pheno$PRS2)
summary(IN_pheno$PRS2)

p1=ggplot(OH_pheno, aes(x=as.factor(dx), y=PRS2, fill=as.factor(dx))) +
  geom_violin() +
  geom_boxplot(width=.2) +
  theme(legend.position = "none") +
  scale_fill_brewer(palette = "Spectral") +
  ggtitle("IN-Trained PRS Distribution in OH") +
  theme(plot.title = element_text(hjust=0.5, size=22), axis.text = element_text(size=16), axis.title=element_text(size=16, face="bold")) +
  xlab("Dementia Status") +
  ylab("Polygenic Risk Score") +
  scale_x_discrete(labels=c("Unaffected", "Affected"))

p2=ggplot(IN_pheno, aes(x=as.factor(dx), y=PRS2, fill=as.factor(dx))) +
  geom_violin() +
  geom_boxplot(width=.2) +
  theme(legend.position = "none") +
  scale_fill_brewer(palette = "Spectral") +
  ggtitle("OH-Trained PRS Distribution in IN") +
  theme(plot.title = element_text(hjust=0.5, size = 22), axis.text = element_text(size=16), axis.title=element_text(size=16, face="bold")) +
  xlab("Dementia Status") +
  ylab("Polygenic Risk Score") +
  scale_x_discrete(labels=c("Unaffected", "Affected"))

ggarrange(p1, p2) # Use 1500 x 500


### Plots for BGSS poster and updated paper

p4=ggplot(OH_pheno, aes(x=as.factor(dx), y=PRS2, fill=as.factor(dx))) +
  geom_violin() +
  geom_boxplot(width=.2) +
  theme(legend.position = "none") +
  scale_fill_brewer(palette = "Spectral") +
  ggtitle("IN-Trained PRS Distribution in OH") +
  theme(plot.title = element_text(hjust=0.5, size=20), axis.text = element_text(size=20), axis.title=element_text(size=20, face="bold")) +
  xlab("Dementia Status") +
  ylab("Polygenic Risk Score") +
  scale_x_discrete(labels=c("Unaffected", "Affected"))

p3=ggplot(IN_pheno, aes(x=as.factor(dx), y=PRS2, fill=as.factor(dx))) +
  geom_violin() +
  geom_boxplot(width=.2) +
  theme(legend.position = "none") +
  scale_fill_brewer(palette = "Spectral") +
  ggtitle("OH-Trained PRS Distribution in IN") +
  theme(plot.title = element_text(hjust=0.5, size = 20), axis.text = element_text(size=20), axis.title=element_text(size=20, face="bold")) +
  xlab("Dementia Status") +
  ylab("Polygenic Risk Score") +
  scale_x_discrete(labels=c("Unaffected", "Affected"))

ggarrange(p3, p4) # Use 1000 x 380

Fig2=ggarrange(p3, p4) # Use 1500 x 500

tiff("Figure2.tiff", units = "in", width = 12, height = 6, res = 300)
Fig2
dev.off()
