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

EU_PRS_in_OH = read.table("Jansen_Weights_OH_Test_No_APOE_MAF05.best", header=T)
IN_PRS_in_OH = read.table("IN_Train_OH_Test_MAF05_No_APOE_AlleleSwap.best", header=T)
OH_PRS_in_OH = read.table("Posthoc PRSs/OH_Train_OH_Test_MAF05_No_APOE.best", header=T)

EU_PRS_in_IN = read.table("Jansen_Weights_IN_Test_No_APOE_MAF05.best", header=T)
OH_PRS_in_IN = read.table("OH_Train_IN_Test_MAF05_No_APOE_AlleleSwap.best", header=T)
IN_PRS_in_IN = read.table("Posthoc PRSs/IN_Train_IN_Test_MAF05_No_APOE.best", header=T)

EU_PRS_in_OH = EU_PRS_in_OH[match(OH_pheno$IID, EU_PRS_in_OH$IID),]
IN_PRS_in_OH = IN_PRS_in_OH[match(OH_pheno$IID, IN_PRS_in_OH$IID),]
OH_PRS_in_OH = OH_PRS_in_OH[match(OH_pheno$IID, OH_PRS_in_OH$IID),]

EU_PRS_in_IN = EU_PRS_in_IN[match(IN_pheno$IID, EU_PRS_in_IN$IID),]
OH_PRS_in_IN = OH_PRS_in_IN[match(IN_pheno$IID, OH_PRS_in_IN$IID),]
IN_PRS_in_IN = IN_PRS_in_IN[match(IN_pheno$IID, IN_PRS_in_IN$IID),]

OH_pheno$EU_PRS = EU_PRS_in_OH$PRS
OH_pheno$IN_PRS = IN_PRS_in_OH$PRS
OH_pheno$OH_PRS = OH_PRS_in_OH$PRS

IN_pheno$EU_PRS = EU_PRS_in_IN$PRS
IN_pheno$OH_PRS = OH_PRS_in_IN$PRS
IN_pheno$IN_PRS = IN_PRS_in_IN$PRS

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

### Check correlations
cor(OH_pheno$EU_PRS, OH_pheno$IN_PRS)
cor(OH_pheno$EU_PRS, OH_pheno$OH_PRS)
cor(OH_pheno$OH_PRS, OH_pheno$IN_PRS)

cor(IN_pheno$EU_PRS, IN_pheno$IN_PRS)
cor(IN_pheno$EU_PRS, IN_pheno$OH_PRS)
cor(IN_pheno$OH_PRS, IN_pheno$IN_PRS)

