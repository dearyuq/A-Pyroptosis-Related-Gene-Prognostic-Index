## 
##    step6.2:使用ICGC 进行外部验证
##            2021.9.14
## 
## 

## 空间设置
rm(list = ls())
options(stringsAsFactors = F,scipen=200)
library(tidyr)
library(dplyr)
library(survival)
library(survminer)
source('./source/rocf.R')


#  前期准备：简单准备数据
# ICGC_DAT <- data.table::fread("F:/DataBase/ICGC/LIHC/exp_seq.LIRI-JP.tsv",data.table = F)
# ICGC_DAT <- ICGC_DAT[,c(1,5,8,9)]
# ICGC_DAT = separate(ICGC_DAT, submitted_sample_id, sep = "_",into = c("ID", "type"))
# ICGC_DAT = ICGC_DAT[ICGC_DAT$type == 'Cancer',]
# dat_meger <- data.frame(icgc_donor_id = unique(ICGC_DAT$icgc_donor_id))
# ICGC_Gene <- unique(ICGC_DAT$gene_id)
# 
# for (i in ICGC_Gene) {
#   data_gene <- ICGC_DAT[ICGC_DAT$gene_id == i,]
#   data_gene <- data_gene[!duplicated(data_gene$icgc_donor_id),]
#   data_gene <- data_gene[,c(1,5)]
#   names(data_gene)[2] <- i
#   dat_meger <- merge(dat_meger,data_gene,by = 'icgc_donor_id',all = T)
# }
# dat_meger <- as.data.frame(t(dat_meger))
# colnames(dat_meger) <- dat_meger[1,]
# dat_meger <- dat_meger[-1,]
# #character->numric
# ICGC_DAT <- apply(as.matrix(dat_meger),2,function(x) as.numeric(x))
# rownames(ICGC_DAT) <- rownames(dat_meger)
# ICGC_DAT <- as.data.frame(t(ICGC_DAT))
# # 替换每一列的NA值为该列的平均值
# x[is.na(x)]=mean(x,na.rm = T)
# ICGC_DAT2 = apply(ICGC_DAT,2,function(x){
#   return(x)
# })
# ICGC_DAT <- as.data.frame(t(ICGC_DAT2))
# 
# save(ICGC_DAT,file = 'rdata/step_6_2_ICGC_DAT.Rdata')

# 1.加载需要的数据
load('../rdata/step_5_COX.Rdata')
load('../rdata/step_6_2_ICGC_DAT.Rdata')

ICGC_clinc <- data.table::fread("F:/DataBase/ICGC/LIHC/donor.LIRI-JP.tsv",data.table = F)
rownames(ICGC_clinc) <- ICGC_clinc$icgc_donor_id
ICGC_clinc <- ICGC_clinc[ICGC_clinc$icgc_donor_id %in% colnames(ICGC_DAT),]
ICGC_clinc <- ICGC_clinc[,c(1,5,6,9,15,17,19,21)]
ICGC_clinc$OS.time <- signif(ICGC_clinc$donor_survival_time/31,3)
ICGC_clinc$OS <- ifelse(ICGC_clinc$donor_vital_status == 'alive',0,1)
ICGC_clinc$Age <- ICGC_clinc$donor_age_at_diagnosis
ICGC_clinc$Gender <- ICGC_clinc$donor_sex
ICGC_clinc$Stage <- ICGC_clinc$donor_tumour_stage_at_diagnosis


# 2.准备COX数据
COX_dat <- log2(ICGC_DAT + 1)
COX_dat <- COX_dat[rownames(COX_dat) %in% lasso$gene,]
COX_clinc <- ICGC_clinc

COX_tdat <- as.data.frame(t(COX_dat))
COX_tdat$icgc_donor_id <- rownames(COX_tdat)
COX_mdat <- inner_join(COX_clinc,COX_tdat,by = "icgc_donor_id")
rownames(COX_mdat) <- COX_mdat$icgc_donor_id
COX_mdat$OS.time <- as.numeric(COX_mdat$OS.time)
COX_mdat$OS <- as.numeric(COX_mdat$OS)


## 2.计算分线得分
Risk_mdat <- theriskf(indat = COX_mdat,ingene = lasso$gene,incoeff = lasso$coeff)
Risk_mdat$RiskScore
median(Risk_mdat$RiskScore)
Risk_mdat$risk = ifelse(Risk_mdat$RiskScore > median(Risk_mdat$RiskScore),'highrisk','lowrisk')
table(Risk_mdat$risk)


## 3.ROC绘图  4x4.5  file = step_6_2_ICGC_ROC
rocf(Risk_mdat)

## 4.高低风险组KM
KM_dat <- Risk_mdat

diff=survdiff(Surv(OS.time, OS) ~ risk,data = KM_dat)
pValue=1-pchisq(diff$chisq,df=1)
fit <- survfit(Surv(OS.time, OS) ~ risk, data = KM_dat)
plot6_2 <- ggsurvplot(fit,palette = c("#E7B800", "#2E9FDF"),
                    risk.table =F,pval =TRUE,pval.size= 4.3,
                    conf.int =F,xlab ="Time (mouth)",ylab = 'Survival rate',
                    ggtheme =theme_light(), 
                    xlim = c(0,60),
                    surv.median.line = "hv",
                    legend = c(0.82,0.30),
                    legend.title = "",legend.labs = c("High risk", "Low risk"),
                    pval.coord=c(1,0.03))
plot6_2
pdf(file = '../fig/step_6_2_ICGC_KM.pdf',width = 4,height = 4)
print(plot6_2)
dev.off()

## 5.重定向数据
ICGC_DAT <- ICGC_DAT
ICGC_clinc <- Risk_mdat
ICGC_clinc$RiskScore <- ICGC_clinc$score
ICGC_clinc$Stage0 <- ifelse(ICGC_clinc$Stage %in% c(1,2),'i/ii','iii/iv')

## 6.保存数据
save(ICGC_DAT,ICGC_clinc,file = '../rdata/step_6_2_ICGC_DATclinc.Rdata')


