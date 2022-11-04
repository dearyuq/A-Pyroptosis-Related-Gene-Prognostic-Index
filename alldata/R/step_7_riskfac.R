##
##  step_8 风险因子关联图
##  2021.9.16
##  5个基因的风险因子图
##  

# 空间设置
rm(list = ls())
options(stringsAsFactors = F)
options(scipen=200)

source('./source/factorf.R')

## 1,加载数据
load('../rdata/step_5_COX.Rdata')
load("../rdata/step_6_Risk_mdat.Rdata")
load('../rdata/step_6_2_ICGC_DATclinc.Rdata')

fac_mix <- ICGC_clinc
factGene <- lasso$gene

factorf(indat = fac_mix,inGene = factGene,insep = '7_icgc2')




