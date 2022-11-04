##
##  step_8 独立预后
##        2021.9.16
##        能否作为独立的预后因素
##

## 空间设置
rm(list = ls())
options(stringsAsFactors = F)
library(ComplexHeatmap)
library(circlize)
source('./source/Staticalf.R')
source('./source/rocf.R')


## 1,加载数据
load("./rdata/step_6_Risk_mdat.Rdata")
load("./rdata/step_6_2_ICGC_DATclinc.Rdata")

names(Risk_mdat)
Risk_mdat$RiskScore <- Risk_mdat$score
Risk_mdat$Stage0 <- ifelse(Risk_mdat$Stage %in% c('i','ii'),'i/ii','iii/iv')
TCGA_clinc <- Risk_mdat

## 2.重定向变量
coxdat = 'TCGA_clinc'
invarlist <- c("RiskScore",'Age','Gender','Stage0')

## 3,临床COX
uncox <- uncoxf(indat = get(coxdat),ingene = invarlist)
uncox
mucox <- mucoxf(indat = get(coxdat),ingene = invarlist)
mucox

## 4.绘制森林图
library(forestplot)
coxlist <- c('uncox','mucox')
collist <- c('green','red')
for (i in 1:2) {
  indat <- get(coxlist[i])
  hrtext <- cbind(c(' ',rownames(indat)),
                  c('P.value',indat$pvalue),
                  c('Hazard ratio',indat$HR_95_CI))
  hrdat <- data.frame('mean'  = c(NA,indat$HR),
                      'lower' = c(NA,indat$HR95L),
                      'upper' = c(NA,indat$HR95H))
  hrdat <- apply(as.matrix(hrdat),2,function(x) as.numeric(x))
  plot1 <- forestplot(hrtext, hrdat,
                      boxsize = 0.2,
                      zero = 1,lwd.zero = 3,lwd.ci = 3,
                      xlab="Hazard ratio",
                      graphwidth = unit(4,'cm'),
                      lineheight = unit(1.8,'cm'),
                      colgap=unit(4,"mm"),
                      txt_gp=fpTxtGp(label = gpar(cex = 1),
                                     ticks = gpar(cex = 1),
                                     xlab  = gpar(cex = 1),
                                     title = gpar(cex = 1)),
                      col=fpColors(box=collist[i],line="darkblue", summary="royalblue")
                      )
  print(plot1)
  pdf(file = paste0('./fig/step_8_',coxdat,'_',coxlist[i],'_forst.pdf'),width = 6,height = 6)
  print(plot1)
  dev.off()
}

## 6.模型基因的表达热图
load('../rdata/step_5_COX.Rdata')
load('../rdata/step_2_TCGA_DEG.Rdata')
Risk_mdat <- Risk_mdat[order(Risk_mdat$risk),]
htdat <- Risk_mdat[,c(rownames(lasso))]
htdat <- as.data.frame(t(htdat))
TCGA_DEG_INFO2 <- TCGA_DEG_INFO[rownames(TCGA_DEG_INFO) %in% rownames(htdat),]

x = TCGA_DEG_INFO2$logFC
names(x) = rownames(TCGA_DEG_INFO2)
FCgene = names(sort(x))
ht_dat = htdat[FCgene,]

ann_col = data.frame(Age = ifelse(Risk_mdat$Age > 60,'>60','<=60'),
                     Stage = Risk_mdat$Stage,
                     Gender = Risk_mdat$Gender,
                     # OS = Risk_mdat$Status,
                     Risk = Risk_mdat$risk)
rownames(ann_col) = colnames(ht_dat)
ann_colors = list(Age    = c('>60'="#FFCC99",'<=60'="#CCFF99"),
                  Stage  = c("i"="#99CCCC","ii"="#FFCC99","iii"="#FFCCCC","iv"="#FF9999"),
                  Gender = c("female"="#FFFF99","male"="#CCCCCC"),
                  OS     = c('Alive'='#99CC66','Dead'='#FF6666'),
                  Risk   = c('lowrisk' = '#00CCCC','highrisk' = '#FF9999'))
main_col = colorRamp2(c(-2, 0, 2), c("green", "white", "red"))
plot6 <- pheatmap(as.matrix(ht_dat),name = ' ',
                  show_rownames = T,show_colnames = F,
                  cluster_rows=T,cluster_cols = F,
                  # main="Heatmap of top 200 lncRNA",
                  annotation_row = NA,
                  annotation_col = ann_col,
                  annotation_colors =ann_colors,
                  scale="row",color = main_col,
                  treeheight_row = 15)
plot6
pdf(file="../fig/step_8_independ_7gene_heatmap.pdf",width = 6,height = 5.5)
plot6
dev.off()


## 7.诺莫图
indat = Risk_mdat; invar = invarlist; instep = 8; inmm = 100
{
  library(rms)
  library(survival)
  Nmdat <- indat
  Nmdat <- Nmdat[,c(invar,'OS','OS.time')]
  Nmdat$OS.time <- as.numeric(Nmdat$OS.time)*30
  dd <- datadist(Nmdat)
  options(datadist="dd")
  multivarl <- as.formula(paste0('Surv(OS.time,OS)~', 
                                 paste(invar, sep = '', collapse = '+')))
  
  coxm_1 <- cph(formula = multivarl,data=Nmdat,surv=T,x=T,y=T,time.inc = 365)
  surv <- Survival(coxm_1)
  surv1 <- function(x) surv(1*365,x)
  surv3 <- function(x) surv(3*365,x)
  surv5 <- function(x) surv(5*365,x)
  nomo <- nomogram(coxm_1,fun = list(surv1,surv3,surv5),lp = T,
                   funlabel = c('1-year survival Probability','3-year survival Probability','5-year survival Probability'),
                   maxscale = 100,fun.at = c(0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1))
  pdf(file = paste0('fig/step_',instep,'_nomodule.pdf'),width = 4,height = 4)
  plot(nomo,lplabel = 'Linear Preadictor',
       xfrac = .35,varname.label = T,varname.label.sep = '=',ia.space = .2,
       tck = NA,tcl = 0.2,lmgp = 0.3,
       points.label = 'Points',total.points.label = 'Total Points',
       total.sep.page = F,
       cap.labels = F,cex.var = 0.53,cex.axis = 0.53,lwd = 0.53,
       label.every = 1,col.grid = gray(c(0.8,0.95)))
  dev.off()
  ## 校准曲线
  f1 <- cph(formula = multivarl, x=T, y=T, surv=T, data=Nmdat, time.inc=365)
  cal1 <- calibrate(f1, cmethod="KM", method="boot", u=1*365, m=inmm, B=1000)
  f3 <- cph(formula = multivarl, x=T, y=T, surv=T, data=Nmdat, time.inc=3*365)
  cal3 <- calibrate(f3, cmethod="KM", method="boot", u=3*365, m=inmm, B=1000)
  f5 <- cph(formula = multivarl, x=T, y=T, surv=T, data=Nmdat, time.inc=5*365)
  cal5 <- calibrate(f5, cmethod="KM", method="boot", u=5*365, m=inmm, B=1000)
  
  pdf(file = paste0('fig/step_',instep,'_nomoadjust.pdf'),width = 4,height = 4)
  plot(cal1,lwd=1,lty=1, cex.axis = 1,cex.lab=1,
       errbar.col = '#666699',
       xlab='Nomogram-Predicted Probability',
       ylab='Actual',
       col = '#666699',subtitles = F,
       xlim = c(0,1),ylim = c(0,1))
  plot(cal3,lwd=1,lty=1, cex.axis = 1,cex.lab=1,add=T,
       errbar.col = '#339933',
       col = '#339933',subtitles = F,
       xlim = c(0,1),ylim = c(0,1))
  plot(cal5,lwd=1,lty=1, cex.axis = 1,cex.lab=1,add=T,
       errbar.col = '#FF0033',
       col = '#FF0033',subtitles = F,
       xlim = c(0,1),ylim = c(0,1))
  abline(0,1,lty=1,lwd=1)
  legend("bottomright",legend=c("1 - year","3 - year","5 - year"), 
         col=c("#666699","#339933","#FF0033"),
         # x.intersp=3, y.intersp=0.8,
         lty= 1 ,lwd= 4,
         bty = "n",
         seg.len=1,cex=1)
  dev.off()
  
  c_in <- coxph(multivarl,data=Nmdat)
  sum.surv <- summary(c_in)
  print(sum.surv$concordance)
}











