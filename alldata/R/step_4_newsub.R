## 
##    根据差异焦亡基因结果分亚群
##    无监督一致性聚类
##    2021.9.14
##    剔除正常组织数据？ - 需要剔除，预后模型不是标志物，会影响模型的准确度
## 

## 空间设置
rm(list = ls())
options(stringsAsFactors = F)
library(pacman)
pacman::p_load(ConsensusClusterPlus,dplyr)
source('./source/subclusterf.R')


## 1. 加载数据
load('../rdata/step_2_Pyro_DEG.Rdata')
load('../rdata/step_1_TCGA_clinc.Rdata')
load('../rdata/step_1_TCGA_DAT.Rdata')
load('../rdata/step_1_Pyro_DAT.Rdata')

## 2去除正常样本
TCGA_clinc <- TCGA_clinc[TCGA_clinc$type == 'Tumor',]
TCGA_fpkm <- TCGA_fpkm[,TCGA_clinc$submitter]
Pyro_DAT <- Pyro_DAT[,TCGA_clinc$submitter]
Pyro_DEG_Dat <- Pyro_DEG_Dat[,TCGA_clinc$submitter]

# 3.聚类:调用自己的函数
newcluster(indat = Pyro_DEG_Dat,ink = 6)

# 4.添加分簇到临床信息，并保存到原始数据
cluster <- read.csv('cluster/cluster.k=2.consensusClass.csv',header = F)
names(cluster) <- c('submitter','cluster')

TCGA_clinc <- inner_join(TCGA_clinc,cluster,by = 'submitter')
TCGA_clinc$newsub <- ifelse(TCGA_clinc$cluster == 1,'cluster1','cluster2')

# 5.保存数据。（注意：这部分已剔除掉正常组织）
save(TCGA_fpkm,file = 'rdata/step_4_TCGA_DAT.Rdata')
save(TCGA_clinc,file = 'rdata/step_4_TCGA_clinc.Rdata')
save(Pyro_DAT,file = 'rdata/step_4_Pyro_DAT.Rdata')
save(Pyro_DEG_Dat,Pyro_DEG_Info,file = 'rdata/step_4_Pyro_DEG_Dat.Rdata')


########################################################
## 第二部分    新分簇评估：临床关系，生存
##    2021.9.14
########################################################

## 空间设置
rm(list = ls())
options(stringsAsFactors = F)
library(pacman)
pacman::p_load(ggplot2,ggpubr,ggsci,ComplexHeatmap,circlize)
library(survival)
library(survminer)

## 1. 加载数据
load('../rdata/step_4_TCGA_clinc.Rdata')
load('../rdata/step_4_Pyro_DEG_Dat.Rdata')


## 2. 亚组与临床关系: 绘制箱型图
names(TCGA_clinc)
# color2 = pal_lancet("lanonc")(9)[c(4,6)]
color2 = c('#ff9900','#146eb4')
comparison = list(c("cluster1","cluster2"))

boxp5 = ggboxplot(TCGA_clinc, size = 1, bxp.errorbar = F,
                  x="newsub", y="OS.time", color = "newsub",
                  palette = color2, add = "jitter",
                  order = c("cluster1","cluster2")) + 
  stat_compare_means(comparisons = comparison, method = "t.test") +
  theme(text = element_text(size = 25),legend.text = element_text(size = 14),
        legend.position = 'none')
boxp5
pdf(file=paste0("../fig/step_4_subtype_",'age',".pdf"),width = 4,height = 4)
boxp5
dev.off()


## 3.新亚组与临床特征热图
TCGA_clinc <- TCGA_clinc[order(TCGA_clinc$cluster),]
Pyro_DEG_Dat <- Pyro_DEG_Dat[,TCGA_clinc$submitter]

x = Pyro_DEG_Info$logFC
names(x) = rownames(Pyro_DEG_Info)
FCgene = names(sort(x))
ht_dat = Pyro_DEG_Dat[FCgene,]

ann_col = data.frame(Cluster=TCGA_clinc$newsub,
                     Stage = TCGA_clinc$Stage,
                     Gender = TCGA_clinc$Gender,
                     OS = TCGA_clinc$vital_status.demographic)
rownames(ann_col) = colnames(ht_dat)
ann_colors = list(Cluster = c("cluster1"="#99CC33","cluster2"="#FFCC00"),
                  Age = c("#66CCFF","#CC3399"),
                  Stage = c("i"="#99CCCC","ii"="#FFCC99","iii"="#FFCCCC","iv"="#FF9999","NA"="white"),
                  Gender = c("female"="#00CCCC","male"="#FF9999"),
                  OS = c('Not Reported'='white','Alive'='#99CC66','Dead'='#FF6666'))
main_col = colorRamp2(c(-4, 0, 4), c("green", "white", "red"))
plot6 <- pheatmap(as.matrix(ht_dat),name = ' ',
                  show_rownames = T,show_colnames = F,
                  cluster_rows=T,cluster_cols = F,
                  # main="Heatmap of top 200 lncRNA",
                  annotation_row = NA,
                  annotation_col = ann_col,
                  annotation_colors =ann_colors,
                  scale="row",color = main_col)
plot6
pdf(file="../fig/step_4_Newsub_clust12_heatmap.pdf",width = 6,height = 4)
plot6
dev.off()

## 4.新亚组生存分析
KM_dat <- TCGA_clinc
diff=survdiff(Surv(OS.time, OS) ~ cluster,data = KM_dat)
pValue=1-pchisq(diff$chisq,df=1)
fit <- survfit(Surv(OS.time, OS) ~ cluster, data = KM_dat)
plot7 <- ggsurvplot(fit,data = KM_dat,
                    palette = c("#2E9FDF", "#E7B800"),
                    risk.table =F,pval =TRUE,pval.size= 4.3,
                    conf.int = F,xlab ="Time (month)",ylab = 'Survival rate',
                    ggtheme =theme_light(),
                    xlim = c(0,60),
                    surv.median.line = "hv",
                    legend = c(0.82,0.83),
                    legend.title = "",legend.labs = c("cluster1", "cluster2"),
                    pval.coord=c(1,0.03))
plot7
pdf(file = '../fig/step_4_cluster_KM.pdf',width = 4,height = 4)
plot7
dev.off()

