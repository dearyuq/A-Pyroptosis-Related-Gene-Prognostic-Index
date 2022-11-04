##
##        step 6: ROC分析
##              2020.7.13
##

rm(list = ls())
options(stringsAsFactors = F)
library(survival)
library(survminer)

source('./source/rocf.R')

## 1.加载cox结果
load('../rdata/step_5_COX.Rdata')
load('../rdata/step_1_Pyro_DAT.Rdata')

intersect(lasso$gene,rownames(Pyro_DAT))

## 2.计算分线得分
Risk_mdat <- theriskf(indat = COX_mdat,ingene = lasso$gene,incoeff = lasso$coeff)
median(Risk_mdat$RiskScore)   # 1.654612

## 3.ROC绘图   4x4.5  step_6_riskscore_ROC
rocf(Risk_mdat)

## 4.高低风险组KM
KM_dat <- Risk_mdat

diff <- survdiff(Surv(OS.time, OS) ~ RiskGroup,data = KM_dat)
pValue <- 1-pchisq(diff$chisq,df=1)
fit <- survfit(Surv(OS.time, OS) ~ RiskGroup, data = KM_dat)
plot6 <- ggsurvplot(fit,palette = c("#E7B800", "#2E9FDF"),
                    risk.table =F,pval =TRUE,pval.size= 4.3,
                    conf.int = F,xlab ="Time (month)",ylab = 'Survival rate',
                    ggtheme =theme_light(),
                    xlim = c(0,60),
                    surv.median.line = "hv",
                    legend = c(0.82,0.83),
                    legend.title = "",legend.labs = c("High risk", "Low risk"),
                    pval.coord=c(1,0.03))
plot6
pdf(file = '../fig/step_6_riskscore_KM.pdf',width = 4,height = 4)
print(plot6)
dev.off()

## 5.保存数据
save(Risk_mdat,file = '../rdata/step_6_Risk_mdat.Rdata')

