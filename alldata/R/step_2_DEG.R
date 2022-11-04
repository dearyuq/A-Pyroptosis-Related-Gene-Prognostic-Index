##
##    step2:  差异分析
##            2021.9.13
##            使用limma差异分析

## 工作空间
rm(list = ls())
options(stringsAsFactors = F)

library(dplyr)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
suppressMessages(library(VennDiagram))
source('./source/DEGf.R')

## 1.加载数据
load("../rdata/step_1_TCGA_DAT.Rdata")
load('../rdata/step_1_TCGA_clinc.Rdata')
load('../rdata/step_1_Pyro_DAT.Rdata')

## 2.差异分析
grouplist <- ifelse(as.numeric(substr(colnames(TCGA_fpkm),14,15)) > 9,'Normal','Tumor')
table((grouplist))
batoh = 'Normal-Tumor'
TCGA_DEG_INFO <- limmaf(indat = TCGA_fpkm, ingroup = grouplist,
                        inp = 0.05,inlogfc = 1,batoh = batoh)
Pyro_DEG_Info <- TCGA_DEG_INFO[TCGA_DEG_INFO$SYMBOL %in% rownames(Pyro_DAT),]
TCGA_DEG_DAT <- TCGA_fpkm[TCGA_DEG_INFO[TCGA_DEG_INFO$State != 'not',]$SYMBOL,]
table(TCGA_DEG_INFO$State)

## 3.差异焦亡基因矩阵
Pyro_DEG_gene <- intersect(rownames(TCGA_DEG_DAT),rownames(Pyro_DAT))
Pyro_DEG_Dat <-  TCGA_DEG_DAT[rownames(TCGA_DEG_DAT) %in% Pyro_DEG_gene,]


## 4.保存数据
# save(TCGA_DEG_INFO,TCGA_DEG_DAT,grouplist,file = '../rdata/step_2_TCGA_DEG.Rdata')
# save(Pyro_DEG_Info,Pyro_DEG_Dat,grouplist,file = '../rdata/step_2_Pyro_DEG.Rdata')


## 5.1火山图
vol_dat = TCGA_DEG_INFO
vol_dat$P.value = -log10(vol_dat$P.Value)
vol_dat$logFC = vol_dat$logFC
plot1 <- ggplot(vol_dat, aes(x = logFC, y = P.value, colour=State)) +
  geom_point(alpha=0.7, size=1) + 
  # xlim(-12,12) + ylim(0,40) +
  scale_color_manual(values=c("#546de5", "black","#ff4757"))+
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  labs(x="logFC",y="-log10(Pvalue)")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank(),
        legend.key = element_rect(colour ="black"))
plot1
ggsave(plot = plot1,filename = '../fig/step_2_DEGall.pdf',width = 4,height = 4)

## 5.2热图:绘制全部焦亡基因
ht_info <- TCGA_DEG_INFO
x = ht_info$logFC
names(x) = rownames(ht_info)
FCgene = c(names(head(sort(x),100)),names(tail(sort(x),100)))
ht_dat = TCGA_fpkm[FCgene,]
# ht_dat = log2(ht_dat)
ann_col = data.frame(Groups=grouplist)
rownames(ann_col) = colnames(ht_dat)
ann_colors = list(Groups = c("Normal"="#00CCCC","Tumor"="#FF9999"))
main_col = colorRamp2(c(-2, 0, 2), c("green", "white", "red"))
plot2 <- pheatmap(as.matrix(ht_dat),name = ' ',
                  show_rownames = F,show_colnames = F,
                  cluster_rows= T,cluster_cols = F,
                  # main="Heatmap of top 200 lncRNA",
                  annotation_row = NA,
                  annotation_col = ann_col,
                  annotation_colors =ann_colors,
                  scale="row",color = main_col,
                  treeheight_row = 15)
plot2
pdf(file='../fig/step_2_heatmap.pdf',width = 6,height =6)
plot2
dev.off()

## 6.韦恩图
air1 <- length(rownames(TCGA_DEG_DAT))
air2 <- length(rownames(Pyro_DAT))
air3 <- length(rownames(Pyro_DEG_Info[Pyro_DEG_Info$State != 'not',]))
venn.plot <- draw.pairwise.venn(area1 = air1,
                                area2 = air2,
                                cross.area = air3,
                                category = c("DEG","Pyroptosis-related gene"),
                                fill = c("blue", "red"),
                                lty = "blank",
                                cex = 2,                        # 1 2 区域内部数字的字体大小 
                                cat.cex = 1,                    # 分类名称的字体大小 
                                cat.dist = 0.01,                # 分类名称距离边的距离 实际调整 
                                cat.just = rep(list(c(0.34, -3)), 2), #分类名称的位置,圈内或者圈外
                                ext.pos = 30,                   # 线的角度 默认是正上方12点位置 
                                ext.dist = -0.05,               # 外部线的距离  跟根据圆圈的大小适当调整
                                ext.length = 0.85,              # 外部线长度 
                                ext.line.lwd = 2,               # 外部线的宽度 
                                ext.line.lty = "dashed" ,       # 外部线为虚线
                                scaled = F)
pdf("../fig/step_2_pyroptos_DEG_veen.pdf",width = 4,height = 4)
grid.draw(venn.plot)
dev.off()
