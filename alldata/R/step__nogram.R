##
## step 9 诺莫图
## 2020.11.30
##

# 空间设置
{
  rm(list = ls())
  options(stringsAsFactors = F)
  
  library(rms)
  library(survival)
}

# 加载数据
load("tmp/step_8_ROC.Rdata")
LIHC <- ROC_dat

names(LIHC)
LIHC <- LIHC[,c('score',"OS","OS.time")]
LIHC$OS.time <- LIHC$OS.time*30

head(LIHC)
dd=datadist(LIHC)
options(datadist="dd")


# ****** 构建logk模型，绘制诺莫图 ******
{
  f <- lrm(OS ~ score, data =  LIHC,x = T,y = T)
  f
  nom <- nomogram(f, fun=plogis, lp=F, funlabel="Risk")
  plot(nom)
  
  cail <- calibrate(f,method = 'boot',B = 100)
  plot(cail,xlim = c(0,1.0),ylim = c(0,1.0))
}

# ****** 绘制COX回归中位生存时间的Nomogram图 ******
## 2.1构建COX比例风险模型  4x4
{
  coxm_1 <- cph(formula = Surv(OS.time,OS)~ score,data=LIHC,surv=T,x=T,y=T,time.inc = 365)
  surv <- Survival(coxm_1)
  surv1 <- function(x) surv(1*365,x)
  surv3 <- function(x) surv(3*365,x)
  surv5 <- function(x) surv(5*365,x)
  x <- nomogram(coxm_1,fun = list(surv1,surv3,surv5),lp = T,
                funlabel = c('1-year survival Probability','3-year survival Probability','5-year survival Probability'),
                maxscale = 100,fun.at = c(0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1))
  plot(x,lplabel = 'Linear Preadictor',
       xfrac = .35,varname.label = T,varname.label.sep = '=',ia.space = .2,
       tck = NA,tcl = 0.2,lmgp = 0.3,
       points.label = 'Points',total.points.label = 'Total Points',
       total.sep.page = F,
       cap.labels = F,cex.var = 0.53,cex.axis = 0.53,lwd = 0.53,
       label.every = 1,col.grid = gray(c(0.8,0.95)))
}

## 2.2构建校准曲线   4x4
{
  #1年校准曲线
  f1 <- cph(Surv(OS.time, OS) ~ score, x=T, y=T, surv=T, data=LIHC, time.inc=365)
  cal1 <- calibrate(f1, cmethod="KM", method="boot", u=1*365, m=15, B=1000)
  plot(cal1,lwd=1,lty=1, cex.axis = 1,cex.lab=1,
       errbar.col = '#666699',
       xlab='Nomogram-Predicted Probability',
       ylab='Actual',
       col = '#666699',subtitles = F,#副标题
       xlim = c(0,1),ylim = c(0,1))
  
  #3年校准曲线
  f3 <- cph(Surv(OS.time, OS) ~ score, x=T, y=T, surv=T, data=LIHC, time.inc=3*365)
  cal3 <- calibrate(f3, cmethod="KM", method="boot", u=3*365, m=15, B=1000)
  plot(cal3,lwd=1,lty=1, cex.axis = 1,cex.lab=1,add=T,
       errbar.col = '#339933',
       col = '#339933',subtitles = F,
       xlim = c(0,1),ylim = c(0,1))
  
  #5年校准曲线
  f5 <- cph(Surv(OS.time, OS) ~ score, x=T, y=T, surv=T, data=LIHC, time.inc=5*365)
  cal5 <- calibrate(f5, cmethod="KM", method="boot", u=5*365, m=15, B=1000)
  plot(cal5,lwd=1,lty=1, cex.axis = 1,cex.lab=1,add=T,
       errbar.col = '#FF0033',
       col = '#FF0033',subtitles = F,
       xlim = c(0,1),ylim = c(0,1))
  
  #加对角线
  abline(0,1,lty=1,lwd=1)#,col=c(rgb(0,0,255,maxColorValue= 255)))
  #加图例
  legend("bottomright",legend=c("1 - year","3 - year","5 - year"), 
         col=c("#666699","#339933","#FF0033"),
         # x.intersp=3, y.intersp=0.8,
         lty= 1 ,lwd= 4,
         bty = "n",#bty框的类型
         seg.len=1,cex=1)
}

## 2.3模型验证,采用c-index
{
  c_in <- coxph(Surv(OS.time,OS)~ score,data=LIHC)
  sum.surv <- summary(c_in)
  c_index <- sum.surv$concordance
  c_index
}




