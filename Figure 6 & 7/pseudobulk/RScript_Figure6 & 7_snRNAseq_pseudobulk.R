# =======================================================================================================================
# Alpha cell inflammation during human pancreas aging and type 2 diabetes and its reversal by calorie restriction in mice
# -----------------------------------------------------------------------------------------------------------------------
# Author:      Michael Schleh
# Repository:  https://github.com/ArrojoDrigoLab/AlphaCellAging
# Journal:     Science Advances

# Description:
#   Generates Figure 6F-H & Supplementary Figure S12F-H, and Supplementary Figure S13.
#   80 week old mouse islet alpha and beta cells subjected to 8 weeks HFD or CR 
#     - A) DESeq on either i) Beta, or ii) Alpha cell normalized counts files
#     - B) Pathway Analysis - From KEGG files
#     - C) Differential Expression figure generation

# Usage:
#   1. Download packages and set the working directory.
#   2. Ensure the required input files (raw gene counts and pathways) exist in `Figure 6/' (see README).
#   3. Run: RScript_Figure6 & 7_snRNAseq_pseudobulk.R
#
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Setup: Working 
# ----------------------------------------------------------------------------
#1A. Working Directory 
#1. Setup library and call in dataframes and create metadata file.----
#Library
#PACKAGES
library(tidyverse)
library(GEOquery)
library(DESeq2)
library(apeglm)
library(ggforce)
library(ggrepel)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(reshape2)
library(pheatmap)
library(ggforce)
library(msigdbr)
library(clusterProfiler)

#1 Set your working directory, build CountMatrix, and DDS----
setwd("...")

################################################
########## Choose Alpha or beta cells ########## 
################################################
CellType <- "pseudobulk_countdata_Betacells.csv"
CellType <- "pseudobulk_countdata_Alphacells.csv"

### Create counts matrix
countdata <- read.csv(CellType, header=TRUE, row.names=1)
countdata <- as.matrix(countdata)

sampleDiet <- factor(c(rep("AL",2), rep("CR",2)))
coldata <- data.frame(row.names=colnames(countdata),sampleDiet)

#DDS file
dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design=~ sampleDiet)
dds$group <- factor(paste0(dds$sampleDiet))
design(dds) <- ~ group
dds$group <- relevel(dds$group, ref="AL") #Default = set to M0WT, change accordingly
dds <- DESeq(dds)
dds
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
normalized_counts <- counts(dds, normalized=TRUE)

#Differential expression analysis
# resultsNames(dds)
res_CR_vs_AL <- lfcShrink(dds,coef="group_CR_vs_AL",type="apeglm")
table(res_CR_vs_AL$padj<0.05)
res_CR_vs_AL <- res_CR_vs_AL[order(res_CR_vs_AL$padj), ]

#DEGs and normalized counts data
#Beta Cells
Beta_DEG<-data.frame(res_CR_vs_AL) %>% filter(padj<.05)
Beta_DEG_normalized <- data.frame(normalized_counts %>% filter(rownames(normalized_counts) %in% rownames(Beta_DEG)))

#Alpha Cells
Alpha_DEG<-data.frame(res_CR_vs_AL) %>% filter(padj<.05)
Alpha_DEG_normalized <- data.frame(normalized_counts) %>% filter(rownames(data.frame(normalized_counts)) %in% rownames(Alpha_DEG))





#Write results to an excel file
#BetaCells----
write.csv(normalized_counts, file=".../Beta_NormalizedCounts.csv")
write.csv(res_CR_vs_AL, file=".../Beta_DEGs.csv")
#AlphaCells
write.csv(normalized_counts, file=".../Alpha_NormalizedCounts.csv")
write.csv(res_CR_vs_AL, file=".../Alpha_DESeq.csv")


############# Figure 6G #########################
#David Analysis by Cluster
#AlphaCells----
Pathways <- rbind.data.frame(data.frame(read.csv("KEGG_Down_CR_Alpha_2026.csv", header = TRUE, sep = ","),
                                     Dir="Down"),
                          data.frame(read.csv("KEGG_Down_CR_Alpha_2026.csv", header = TRUE, sep = ","),
                                     Dir="Up"))
Pathways <- Pathways %>% filter(FDR<0.05)
Pathways <- Pathways %>% dplyr::filter(Term %in% c("Cell adhesion molecules",
                                      "Antigen processing and presentation",
                                      "Allograft rejection",
                                      "Type I diabetes mellitus",
                                      "Cellular senescence",
                                      "ECM-receptor interaction",
                                      "Glucagon signaling pathway",
                                      "AMPK signaling pathway",
                                      "Insulin signaling pathway"))
Pathways$FDRNew <- ifelse(Pathways$Dir=="Up",-log(Pathways$FDR),-log(Pathways$FDR)*-1)

                          
ggplot(Pathways, aes(x=FDRNew,y=reorder(Term,FDRNew))) +
  geom_segment(aes(yend = reorder(Term,FDRNew), xend = 0),color ="black",size=0.8,alpha=0.5) +
  geom_point(aes(size=Count,fill=Dir),shape=21,color="black") +
  labs(x=bquote(~-log[10]~ 'FDR')) +
  geom_vline(xintercept = 0,color="black",size=1) +
  scale_fill_manual(values=c("white","#F71480")) +
  theme(plot.title = element_blank(),
        axis.text = element_text(size=12,color="black"),
        axis.title = element_text(size=12,color="black"),
        axis.title.y = element_blank(),
        axis.line = element_line(size=1),
        legend.position = "right",
        panel.background = element_blank())


############# Figure 6G #########################
#David Analysis by Cluster
#AlphaCells----
Pathways <- rbind.data.frame(data.frame(read.csv("KEGG_Down_CR_Alpha_2026.csv", header = TRUE, sep = ","),
                                        Dir="Down"),
                             data.frame(read.csv("KEGG_Down_CR_Alpha_2026.csv", header = TRUE, sep = ","),
                                        Dir="Up"))
Pathways <- Pathways %>% filter(FDR<0.05)
Pathways <- Pathways %>% dplyr::filter(Term %in% c("Cell adhesion molecules",
                                                   "Antigen processing and presentation",
                                                   "Allograft rejection",
                                                   "Type I diabetes mellitus",
                                                   "Cellular senescence",
                                                   "ECM-receptor interaction",
                                                   "Glucagon signaling pathway",
                                                   "AMPK signaling pathway",
                                                   "Insulin signaling pathway"))
Pathways$FDRNew <- ifelse(Pathways$Dir=="Up",-log(Pathways$FDR),-log(Pathways$FDR)*-1)


ggplot(Pathways, aes(x=FDRNew,y=reorder(Term,FDRNew))) +
  geom_segment(aes(yend = reorder(Term,FDRNew), xend = 0),color ="black",size=0.8,alpha=0.5) +
  geom_point(aes(size=Count,fill=Dir),shape=21,color="black") +
  labs(x=bquote(~-log[10]~ 'FDR')) +
  geom_vline(xintercept = 0,color="black",size=1) +
  scale_fill_manual(values=c("white","#F71480")) +
  theme(plot.title = element_blank(),
        axis.text = element_text(size=12,color="black"),
        axis.title = element_text(size=12,color="black"),
        axis.title.y = element_blank(),
        axis.line = element_line(size=1),
        legend.position = "right",
        panel.background = element_blank())





#PCA Plots (Not shown in manuscript)
#Figure 4A:  PCA PLOT----
vst <- vst(dds)
intgroup=c("sampleGroup")

# Principal components analysis, Default is set up to M0 polarization
PCAdata1 <- plotPCA(vst,intgroup=c("sampleDiet"),returnData = TRUE)
PCAdata1$sampleDiet <- factor(PCAdata1$sampleDiet, levels = c("AL:CR"))
percentVar <- round(100 * attr(PCAdata1, "percentVar")) 

ggplot(PCAdata1, aes(x = PC1, y = PC2, fill=group, color = group,)) + 
   geom_mark_ellipse(aes(fill = group,
                         color = group,
                         label=NULL)) +
  geom_point(aes(fill=group,shape=group),size =6,shape=21) +
  scale_fill_manual(values=c("white","#F71480")) +
  scale_color_manual(values=c("black","black")) +
  scale_x_continuous(limits=c(-20,35)) +
  scale_y_continuous(limits=c(-20,15)) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) + 
  ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
  ggtitle("CellType") +
  theme(plot.title = element_text(hjust=0.5,size=20,face="bold"),
        axis.text = element_text(size=20,color="black"),
        axis.title = element_text(size=20,color="black"),
        axis.line = element_line(size=1),
        legend.position = "none",
        panel.background = element_blank())

#Heatmap Figures ----

#END
