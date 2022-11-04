##
##    step10:  高低风险组 免疫活性比较
##            2021.9.14
##

## 工作空间
rm(list = ls())
options(stringsAsFactors = F)
suppressMessages(library(GSEABase))
source('./source/immunef.R')

## 1.加载数据
load("./rdata/step_4_TCGA_DAT.Rdata")
load('./rdata/step_6_Risk_mdat.Rdata')
load('./rdata/step_6_2_ICGC_DATclinc.Rdata')
geneSets <- getGmt("../input/addinfo/Hhimmune.symbol.gmt")

## 2.整理数据
ICGC_clinc$submitter <- ICGC_clinc$icgc_donor_id


## 3.免疫ssGSEA，图片输出到本地fig/
ssgst <- Imm1316f(indat = ICGC_DAT,inclinc = ICGC_clinc,inSet = geneSets,ingroup = 'risk',insep = '12_ICGC')


# 4.保存数据
write.csv(ssgst,file = "../output/step_10_ssGSEA_result.csv")


