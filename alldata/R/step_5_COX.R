##
##        step 5 COX分析：建立模型
##              2021.9.15
##


rm(list = ls())
options(stringsAsFactors = F)

library(dplyr)
library(survival)
library(survminer)
library(glmnet)
source('./source/Staticalf.R')
source('./source/forstf.R')


## 1.加载数据
load("../rdata/step_4_Pyro_DAT.Rdata")
load('../rdata/step_4_TCGA_clinc.Rdata')
COX_dat <- Pyro_DAT
COX_clinc <- TCGA_clinc

## 2.整理数据
COX_tdat <- as.data.frame(t(COX_dat))
COX_tdat$submitter <- rownames(COX_tdat)
COX_mdat <- inner_join(COX_clinc,COX_tdat,by = "submitter")
rownames(COX_mdat) <- COX_mdat$submitter
# COX_mdat <- COX_mdat[!is.na(COX_mdat$OS.time),]

## 3.单因素cox
uncox1 <- uncoxf(indat = COX_mdat,ingene = rownames(COX_dat))
uncox <- uncox1[uncox1$pvalue < 0.05,]
uncox

## 4.lasso回归
lasso <- lassof(indat = COX_mdat,ingene = rownames(uncox),insep = 5)
lasso

## 5.多因素cox回归
mucox <- mucoxf(indat = COX_mdat,ingene = rownames(uncox))
mucox

## 6.保存数据
write.csv(uncox1,file = '../output/step_5_cox_uncox1.csv',quote = F)
write.csv(uncox,file = '../output/step_5_cox_uncox.csv',quote = F)
write.csv(lasso,file = '../output/step_5_cox_lasso.csv',quote = F)
write.csv(mucox,file = '../output/step_5_cox_mucox.csv',quote = F)
save(COX_mdat,lasso,file = '../rdata/step_5_COX.Rdata')

## 7.绘制森林图
uncox1 <- read.csv('../output/step_5_cox_uncox1.csv',row.names = 1)
uncox1 <- read.csv('../output/step_5_cox_uncox1.csv',row.names = 1)
Forstf(indat = uncox1,instep = 5)
Forstf(indat = uncox,instep = '5_9cox')

