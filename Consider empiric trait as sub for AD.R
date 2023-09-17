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

OH_PRS = read.table("IN_Train_OH_Test_MAF05_No_APOE_AlleleSwap.best", header=T)
IN_PRS = read.table("OH_Train_IN_Test_MAF05_No_APOE_AlleleSwap.best", header=T)

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

# PRS models

prs_mod_OH = glm(dx ~ PRS, data=OH_pheno)
summary(prs_mod_OH)
Cstat(prs_mod_OH)

## IN Models

# PRS models

prs_mod_IN = glm(dx ~ PRS, data=IN_pheno)
summary(prs_mod_IN)
Cstat(prs_mod_IN)

## Make simulated traits
# OH: 730 total, 78 cases, 652 controls
# IN: 789 total, 74 cases, 715 controls
#OHt1 = sample(0:1, 730, replace=T, prob = c(0.95, 0.05))
OH_pheno$t1 = 1:730
OH_pheno$t2 = 1:730
OH_pheno$t3 = 1:730
OH_pheno$t4 = 1:730
OH_pheno$t5 = 1:730
OH_pheno$t6 = 1:730
OH_pheno$t7 = 1:730
OH_pheno$t8 = 1:730
OH_pheno$t9 = 1:730
OH_pheno$t10 = 1:730
OHt1 = sample.int(n = 730, size = 78, replace=F)
OHt2 = sample.int(n = 730, size = 78, replace=F)
OHt3 = sample.int(n = 730, size = 78, replace=F)
OHt4 = sample.int(n = 730, size = 78, replace=F)
OHt5 = sample.int(n = 730, size = 78, replace=F)
OHt6 = sample.int(n = 730, size = 78, replace=F)
OHt7 = sample.int(n = 730, size = 78, replace=F)
OHt8 = sample.int(n = 730, size = 78, replace=F)
OHt9 = sample.int(n = 730, size = 78, replace=F)
OHt10 = sample.int(n = 730, size = 78, replace=F)

IN_pheno$t1 = 1:789
IN_pheno$t2 = 1:789
IN_pheno$t3 = 1:789
IN_pheno$t4 = 1:789
IN_pheno$t5 = 1:789
IN_pheno$t6 = 1:789
IN_pheno$t7 = 1:789
IN_pheno$t8 = 1:789
IN_pheno$t9 = 1:789
IN_pheno$t10 = 1:789
INt1 = sample.int(n = 789, size = 74, replace=F)
INt2 = sample.int(n = 789, size = 74, replace=F)
INt3 = sample.int(n = 789, size = 74, replace=F)
INt4 = sample.int(n = 789, size = 74, replace=F)
INt5 = sample.int(n = 789, size = 74, replace=F)
INt6 = sample.int(n = 789, size = 74, replace=F)
INt7 = sample.int(n = 789, size = 74, replace=F)
INt8 = sample.int(n = 789, size = 74, replace=F)
INt9 = sample.int(n = 789, size = 74, replace=F)
INt10 = sample.int(n = 789, size = 74, replace=F)

for (i in 1:nrow(OH_pheno)){
  if (OH_pheno$t1[i] %in% OHt1) {OH_pheno$t1[i] = 1}
  else {
    OH_pheno$t1[i] = 0
  }
}

for (i in 1:nrow(OH_pheno)){
  if (OH_pheno$t2[i] %in% OHt2) {OH_pheno$t2[i] = 1}
  else {
    OH_pheno$t2[i] = 0
  }
}

for (i in 1:nrow(OH_pheno)){
  if (OH_pheno$t3[i] %in% OHt1) {OH_pheno$t3[i] = 1}
  else {
    OH_pheno$t3[i] = 0
  }
}

for (i in 1:nrow(OH_pheno)){
  if (OH_pheno$t4[i] %in% OHt1) {OH_pheno$t4[i] = 1}
  else {
    OH_pheno$t4[i] = 0
  }
}

for (i in 1:nrow(OH_pheno)){
  if (OH_pheno$t5[i] %in% OHt1) {OH_pheno$t5[i] = 1}
  else {
    OH_pheno$t5[i] = 0
  }
}

for (i in 1:nrow(OH_pheno)){
  if (OH_pheno$t6[i] %in% OHt6) {OH_pheno$t6[i] = 1}
  else {
    OH_pheno$t1[i] = 0
  }
}

for (i in 1:nrow(OH_pheno)){
  if (OH_pheno$t7[i] %in% OHt1) {OH_pheno$t7[i] = 1}
  else {
    OH_pheno$t7[i] = 0
  }
}

for (i in 1:nrow(OH_pheno)){
  if (OH_pheno$t8[i] %in% OHt1) {OH_pheno$t8[i] = 1}
  else {
    OH_pheno$t8[i] = 0
  }
}

for (i in 1:nrow(OH_pheno)){
  if (OH_pheno$t9[i] %in% OHt1) {OH_pheno$t9[i] = 1}
  else {
    OH_pheno$t9[i] = 0
  }
}

for (i in 1:nrow(OH_pheno)){
  if (OH_pheno$t10[i] %in% OHt10) {OH_pheno$t10[i] = 1}
  else {
    OH_pheno$t10[i] = 0
  }
}

for (i in 1:nrow(IN_pheno)){
  if (IN_pheno$t1[i] %in% INt1) {IN_pheno$t1[i] = 1}
  else {
    IN_pheno$t1[i] = 0
  }
}

for (i in 1:nrow(IN_pheno)){
  if (IN_pheno$t2[i] %in% INt2) {IN_pheno$t2[i] = 1}
  else {
    IN_pheno$t2[i] = 0
  }
}

for (i in 1:nrow(IN_pheno)){
  if (IN_pheno$t3[i] %in% INt1) {IN_pheno$t3[i] = 1}
  else {
    IN_pheno$t3[i] = 0
  }
}

for (i in 1:nrow(IN_pheno)){
  if (IN_pheno$t4[i] %in% INt1) {IN_pheno$t4[i] = 1}
  else {
    IN_pheno$t4[i] = 0
  }
}

for (i in 1:nrow(IN_pheno)){
  if (IN_pheno$t5[i] %in% INt1) {IN_pheno$t5[i] = 1}
  else {
    IN_pheno$t5[i] = 0
  }
}

for (i in 1:nrow(IN_pheno)){
  if (IN_pheno$t6[i] %in% INt6) {IN_pheno$t6[i] = 1}
  else {
    IN_pheno$t1[i] = 0
  }
}

for (i in 1:nrow(IN_pheno)){
  if (IN_pheno$t7[i] %in% INt1) {IN_pheno$t7[i] = 1}
  else {
    IN_pheno$t7[i] = 0
  }
}

for (i in 1:nrow(IN_pheno)){
  if (IN_pheno$t8[i] %in% INt1) {IN_pheno$t8[i] = 1}
  else {
    IN_pheno$t8[i] = 0
  }
}

for (i in 1:nrow(IN_pheno)){
  if (IN_pheno$t9[i] %in% INt1) {IN_pheno$t9[i] = 1}
  else {
    IN_pheno$t9[i] = 0
  }
}

for (i in 1:nrow(IN_pheno)){
  if (IN_pheno$t10[i] %in% INt10) {IN_pheno$t10[i] = 1}
  else {
    IN_pheno$t10[i] = 0
  }
}

### Get AUCs

# OH
t1_mod_OH = glm(t1 ~ PRS, data=OH_pheno)
Cstat(t1_mod_OH)

t2_mod_OH = glm(t2 ~ PRS, data=OH_pheno)
Cstat(t2_mod_OH)

t3_mod_OH = glm(t3 ~ PRS, data=OH_pheno)
Cstat(t3_mod_OH)

t4_mod_OH = glm(t4 ~ PRS, data=OH_pheno)
Cstat(t4_mod_OH)

t5_mod_OH = glm(t5 ~ PRS, data=OH_pheno)
Cstat(t5_mod_OH)

t6_mod_OH = glm(t6 ~ PRS, data=OH_pheno)
Cstat(t6_mod_OH)

t7_mod_OH = glm(t7 ~ PRS, data=OH_pheno)
Cstat(t7_mod_OH)

t8_mod_OH = glm(t8 ~ PRS, data=OH_pheno)
Cstat(t8_mod_OH)

t9_mod_OH = glm(t9 ~ PRS, data=OH_pheno)
Cstat(t9_mod_OH)

t10_mod_OH = glm(t10 ~ PRS, data=OH_pheno)
Cstat(t10_mod_OH)


# IN

t1_mod_IN = glm(t1 ~ PRS, data=IN_pheno)
Cstat(t1_mod_IN)

t2_mod_IN = glm(t2 ~ PRS, data=IN_pheno)
Cstat(t2_mod_IN)

t3_mod_IN = glm(t3 ~ PRS, data=IN_pheno)
Cstat(t3_mod_IN)

t4_mod_IN = glm(t4 ~ PRS, data=IN_pheno)
Cstat(t4_mod_IN)

t5_mod_IN = glm(t5 ~ PRS, data=IN_pheno)
Cstat(t5_mod_IN)

t6_mod_IN = glm(t6 ~ PRS, data=IN_pheno)
Cstat(t6_mod_IN)

t7_mod_IN = glm(t7 ~ PRS, data=IN_pheno)
Cstat(t7_mod_IN)

t8_mod_IN = glm(t8 ~ PRS, data=IN_pheno)
Cstat(t8_mod_IN)

t9_mod_IN = glm(t9 ~ PRS, data=IN_pheno)
Cstat(t9_mod_IN)

t10_mod_IN = glm(t10 ~ PRS, data=IN_pheno)
Cstat(t10_mod_IN)

