##
## step10:高低分组,差异基因,通路富集分析
##
## 2021.3.30
##


## 工作空间
rm(list = ls())
options(stringsAsFactors = F)

library(dplyr)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)
source('./source/DEGf.R')

## 1.加载数据
load("../rdata/step_4_TCGA_DAT.Rdata")
load('../rdata/step_6_Risk_mdat.Rdata')
load("../rdata/step_5_COX.Rdata")

TCGA_fpkm <- TCGA_fpkm[,colnames(TCGA_fpkm) %in% Risk_mdat$submitter]

## 2.差异分析
grouplist <- Risk_mdat$risk
batoh = 'highrisk-lowrisk'
Risk_limma_INFO <- limmaf(indat = TCGA_fpkm,ingroup = grouplist,
                          inp = 0.05,inlogfc = 1,batoh = batoh)
table(Risk_limma_INFO$State)
modgeneFC <- Risk_limma_INFO[Risk_limma_INFO$SYMBOL %in% lasso$gene,]

## 3.火山图
vol_dat = Risk_limma_INFO
vol_dat$P.value = -log10(vol_dat$P.Value)
vol_dat$logFC = vol_dat$logFC

plot9 <- ggplot(vol_dat, aes(x = logFC, y = P.value, colour=stand)) +
  geom_point(alpha=0.7, size=1) +
  scale_color_manual(values=c("#546de5", "black","#ff4757"))+
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  labs(x="logFC",y="-log10(Pvalue)")+
  # xlim(-5,5) + ylim(0,30) +
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank(),
        legend.key = element_rect(colour ="black"))
plot9
ggsave(plot = plot9,filename = '../fig/step_9_RiskGroup_vol.pdf',width = 4,height = 4)

## 4.保存数据
save(Risk_limma_INFO,file = '../rdata/step_9_RiskGroup_DEG.Rdata')

##  5.通路富集分析
names(Risk_limma_INFO)
nrDEG = subset(Risk_limma_INFO,P.Value < 0.05 & abs(logFC) > 1)

nrDEG_2 <- bitr(unique(nrDEG$SYMBOL), fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = org.Hs.eg.db)
nrDEG_2 <- merge(nrDEG,nrDEG_2,by.x = 'SYMBOL',by.y = 'SYMBOL')

gene_up <- nrDEG_2[nrDEG_2$stand == 'up','ENTREZID']
gene_down <- nrDEG_2[nrDEG_2$stand == 'down','ENTREZID']
gene_all <- as.character(nrDEG_2[,'ENTREZID'])

## 5.2绘图
#KEGG
kk.all <- enrichKEGG(gene         = gene_all,
                     organism     = 'hsa',
                     pvalueCutoff = 0.5,
                     qvalueCutoff = 0.9)
pdf(file = '../fig/step_9_RiskGroup_KEGG.pdf',width = 8,height = 8)
barplot(kk.all,showCategory=30, title ="Enrichment KEGG")
dev.off()

#GO
go.all <- enrichGO(gene          = gene_all, 
                   OrgDb         = org.Hs.eg.db,
                   ont           = 'all', 
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.5, 
                   qvalueCutoff  = 0.9, 
                   readable      = TRUE)
pdf(file = '../fig/step_9_RiskGroup_GO.pdf',width = 8,height = 8)
barplot(go.all,split="ONTOLOGY",title = 'Enrichment GO') + 
  facet_grid(ONTOLOGY~., scale="free") + 
  scale_y_discrete(labels=function(x) str_wrap(x, width=60))
dev.off()


## 6.挑选指定通路展示
# 将data.fram转为可以画通路图的格式: plodat会报错，但可以运行
allsult <- go.all@result
allsult <- allsult[allsult$p.adjust < 0.05,]
allsult$Description
write.csv(allsult,file = '../output/step_9_allsult.csv')

# allsult <- allsult[c(1,2),]
plodat <- new("enrichResult",result = allsult)

barplot(plodat, showCategory=30, title ="Enrichment KEGG")

barplot(plodat, split="ONTOLOGY", title = 'Enrichment GO') + 
  facet_grid(ONTOLOGY~., scale="free") + 
  scale_y_discrete(labels=function(x) str_wrap(x, width=60))


