##
##    step3:  PPI分析与焦亡基因相关性分析
##            
##

################################################
## 第一部分：PPI分析
# 在String中完成
################################################




################################################
## 第二部分：焦亡基因相关性分析

rm(list = ls())
options(stringsAsFactors = F)
library(ComplexHeatmap)

## 1.加载数据
load('../rdata/step_1_Pyro_DAT.Rdata')

## 2.准备相关性分析数据
cor_dat <- as.data.frame(t(Pyro_DAT))

# 3.计算相关性
matcar.cor <- cor(cor_dat)
matcar.cor

# 4.绘图
library(corrplot)
col<- colorRampPalette(c("blue", "white", "red"))(20)
heatmap(x = matcar.cor, col = col, symm = TRUE)

pdf(file = 'fig/step_3_pyro_corr.pdf',width = 6,height = 6)
heatmap(x = matcar.cor, col = col, symm = TRUE)
dev.off()

# 5.导出到文本
outtab <- data.frame()
for (i in 1:ncol(matcar.cor)) {
  dat <- data.frame(va1 = rownames(matcar.cor),
                    va2 = matcar.cor[,i],
                    va3 = colnames(matcar.cor)[i])
  rownames(dat) <- NULL
  outtab <- rbind(outtab,dat)
}
write.csv(outtab,file = 'output/step_3_pyro_coor.csv')



