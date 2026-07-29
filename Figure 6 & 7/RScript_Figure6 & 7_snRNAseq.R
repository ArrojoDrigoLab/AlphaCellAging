# =======================================================================================================================
# Alpha cell inflammation during human pancreas aging and type 2 diabetes and its reversal by calorie restriction in mice
# -----------------------------------------------------------------------------------------------------------------------
# Author:      Michael Schleh
# Repository:  https://github.com/ArrojoDrigoLab/AlphaCellAging
# Journal:     Science Advances

# Description:
#   Generates Figure 6 and 7 panels for mouse snRNA-seq in old mice fed AL vs CR diets:
#     - A) RDS setup, and cleaning
#     - B) Single cell clustering
#     - C) Population analysis
#     - D) pseudobulk analysis
#     - E) CellChat

# Usage:
#   1. Download packages and set the working directory.
#   2. Ensure the required input files exist in `Figure 6 and 7/' (see README).
#   3. Run: Rscript_Figure1_Islet_Proteomics.R
#
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Setup: Working 
# ----------------------------------------------------------------------------
#1A. Working Directory 
setwd("") <---
#setwd("C:/Users/schlehmw/Box/Mike_OldMice_snRNseq_Bender")
setwd("C:/Users/schlehmw/Box/Mike - Nature Aging/Source Data/Figure 6 & 7")

#1B. Download necessary packages
install.packages(c("tidyverse", "Seurat", "hdf5r", "patchwork", "sctransform","viridis", "enrichR", "cowplot", "Matrix", "reshape2","pheatmap", "png",
                   "RColorBrewer", "data.table","ggalluvial", "ggridges", "ggpubr", "scales","ggh4x", "rstatix"))

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("S4Vectors","SingleCellExperiment","apeglm","DESeq2"))

if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("sqjin/CellChat")

#1B. Library packages
library(tidyverse)
library(Seurat)
library(hdf5r)
library(ggplot2)
library(patchwork)
library(sctransform)
library(dplyr)
library(viridis)
library(viridisLite)
library(enrichR)
library(cowplot)
#library(Matrix.utils)
library(Matrix)
library(reshape2)
library(S4Vectors)
library(SingleCellExperiment)
library(pheatmap)
library(apeglm)
library(png)
library(DESeq2)
library(RColorBrewer)
library(data.table)
library(ggalluvial)
library(ggridges)
library(CellChat)
library(ggpubr)
library(scales)
library(ggh4x)
library(rstatix)
library(plyr)
library(patchwork)



#1C - Set parallel workers
options(future.globals.maxSize = 1e9)




##### 1D Load CellBender Pre-processed Cells #####
#### Skip to line 193 to use pre-processed object
#Load CellBender H5 Files:----
MS01 <- Read10X_h5("C:/Users/schlehmw/Box/Mike_OldMice_snRNseq_Bender/al1_bender_filtered_seurat.h5", use.names = TRUE, unique.features = TRUE)
MS02 <- Read10X_h5("C:/Users/schlehmw/Box/Mike_OldMice_snRNseq_Bender/al2_bender_filtered_seurat.h5", use.names = TRUE, unique.features = TRUE)
MS03 <- Read10X_h5("C:/Users/schlehmw/Box/Mike_OldMice_snRNseq_Bender/cr1_bender_filtered_seurat.h5", use.names = TRUE, unique.features = TRUE)
MS04 <- Read10X_h5("C:/Users/schlehmw/Box/Mike_OldMice_snRNseq_Bender/cr2_bender_filtered_seurat.h5", use.names = TRUE, unique.features = TRUE)

AL1 <- CreateSeuratObject(counts = MS01, project = "MS05", min.cells = 3, min.features = 200)
AL2 <- CreateSeuratObject(counts = MS02, project = "MS05", min.cells = 3, min.features = 200)
CR1 <- CreateSeuratObject(counts = MS03, project = "MS05", min.cells = 3, min.features = 200)
CR2 <- CreateSeuratObject(counts = MS04, project = "MS05", min.cells = 3, min.features = 200)

#Add group metadata
AL1$condition <- "AL_old"
AL2$condition <- "AL_old"
CR1$condition <- "CR_old"
CR2$condition <- "CR_old"

#Add individual metadata
AL1$identifier <- "AL1"
AL2$identifier <- "AL2"
CR1$identifier <- "CR1"
CR2$identifier <- "CR2"

# Retrieve the counts matrix from the "RNA" assay----
counts <- GetAssayData(AL1, assay = "RNA", slot = "counts")
counts[counts < 0] <- 0
AL1[["RNA"]]@layers$counts <- counts
# Option 1: Check the minimum value in the counts matrix
min_val <- min(GetAssayData(AL1, assay = "RNA", slot = "counts"))
print(min_val)  # Should be 0 or greater

counts <- GetAssayData(AL2, assay = "RNA", slot = "counts")
counts[counts < 0] <- 0
AL2[["RNA"]]@layers$counts <- counts
min_val <- min(GetAssayData(AL2, assay = "RNA", slot = "counts"))
print(min_val)  # Should be 0 or greater

counts <- GetAssayData(CR1, assay = "RNA", slot = "counts")
counts[counts < 0] <- 0
CR1[["RNA"]]@layers$counts <- counts
min_val <- min(GetAssayData(CR1, assay = "RNA", slot = "counts"))
print(min_val)  # Should be 0 or greater

counts <- GetAssayData(CR2, assay = "RNA", slot = "counts")
counts[counts < 0] <- 0
CR2[["RNA"]]@layers$counts <- counts
min_val <- min(GetAssayData(CR2, assay = "RNA", slot = "counts"))
print(min_val)  # Should be 0 or greater

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats----
AL1[["percent.mt"]] <- PercentageFeatureSet(AL1, pattern = "^MT-")
AL2[["percent.mt"]] <- PercentageFeatureSet(AL2, pattern = "^MT-")
CR1[["percent.mt"]] <- PercentageFeatureSet(CR1, pattern = "^MT-")
CR2[["percent.mt"]] <- PercentageFeatureSet(CR2, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(AL1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(AL2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(CR1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(CR2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# run standard anlaysis workflow----
AL1 <- NormalizeData(AL1)
AL1 <- FindVariableFeatures(AL1)
AL1 <- ScaleData(AL1)
AL1 <- RunPCA(AL1)

AL2 <- NormalizeData(AL2)
AL2 <- FindVariableFeatures(AL2)
AL2 <- ScaleData(AL2)
AL2 <- RunPCA(AL2)

CR1 <- NormalizeData(CR1)
CR1 <- FindVariableFeatures(CR1)
CR1 <- ScaleData(CR1)
CR1 <- RunPCA(CR1)

CR2 <- NormalizeData(CR2)
CR2 <- FindVariableFeatures(CR2)
CR2 <- ScaleData(CR2)
CR2 <- RunPCA(CR2)

# Find integration anchors and cluster----
anchors <- FindIntegrationAnchors(object.list = list(AL1, AL2, CR1, CR2), anchor.features = 4000, reduction = "rpca")

# Integrate data
integrated_data <- IntegrateData(anchors = anchors, dims = 1:30)
integrated_data

integrated_data <- ScaleData(integrated_data)
integrated_data <- RunPCA(integrated_data)

#Clustering
integrated_data <- FindNeighbors(integrated_data, dims = 1:30)
integrated_data <- FindClusters(integrated_data, resolution = 0.5)
integrated_data <- RunUMAP(integrated_data, dims = 1:30, reduction = "pca")

#Get cluster marker genes
islet.markers <- FindAllMarkers(integrated_data, only.pos = TRUE)
islet.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)



#####2 Read pre-processed RDS File
#Cellranger files pre-processed by CellBender (Python app) to remove technical background and noise before downstream analysis in Seurat (See Zenodo for pre-processed CellRanger samples)
integrated_mouse_scRNAseq <- readRDS("integrated_mouse_scRNAseq.rds")

###########################################################
########### Figure 6 & Supplementary Figure S12 ########### 
###########################################################

###Supplementary Figure S12F and G - See "pseudobulk" script

#Figure 6E - Dimplot
DimPlot(integrated_mouse_scRNAseq, reduction = "umap", group.by = "cell_type", label = TRUE, pt.size = 1.5)


# Supplementary Figure S12A
DimPlot(integrated_mouse_scRNAseq, reduction = "umap", 
        split.by = "condition", label = TRUE, pt.size = 1.5)

# Supplementary Figure S12B
#Dotplot for BCell markers
MarkerIdentifiers <- c("Ins1","Ins2","Ero1b","Gcg","Sst","Pecam1","Ppy")
EndocrineMarkers <- betacells <- subset(integrated_mouse_scRNAseq, idents = c("Beta","Alpha","Delta","EC","Ppy"), invert = FALSE)


DotPlot(EndocrineMarkers, features = MarkerIdentifiers, dot.min = 0, col.min = 0, dot.scale=10) +
  #coord_flip() +
  RotatedAxis() +
  # scale_y_discrete(
  # limits = c("B naive", "CD4 naive", "Treg", "Macrophage","B memory","CD8 effector","B proliferating","B activated",
  #            "DC","T proliferating","NK","B plasma")) +
  geom_point(aes(size = pct.exp), shape = 21, colour = "black", stroke = 0.5) +
  scale_fill_viridis_c() +
  labs(x = NULL, y = NULL) +
  scale_colour_viridis(option = "viridis") +
  theme(plot.margin = unit(c(1, 1, 1, 2), "cm"),
        axis.text.y=element_text(face="italic",),
        #legend.position = 'none'
        )


# Supplementary Figure S12C
metadata <- integrated_mouse_scRNAseq@meta.data
endocrine_cells <- metadata %>% filter(cell_type %in% c("Beta", "Alpha", "Delta", "Ppy"))
ggplot(endocrine_cells, aes(x = condition, fill = cell_type)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values=c("#F8766D","#DB8E00","#AEA100","#00C1A7"))+
  scale_y_continuous(labels = percent_format()) +
  scale_x_discrete(labels=c("AL","CR"))+
  labs(
    x = "",
    y = "Proportion",
    fill = "Cell Type"
  ) +
  theme_classic() 


#Supplementary Figure S12 E-J:
#Subset beta cells for pseudobulk/comparative analysis----
betacells <- subset(integrated_mouse_scRNAseq, idents = "Beta", invert = FALSE)
DimPlot(betacells, reduction = "umap", group.by = "cell_type", label = TRUE, split.by = "condition", pt.size = 2)

#Supplementary Figure S12I: VLNPLOT FOR SIGNIFICANT BETA CELL HSP & PROTEIN PROCESSING GENES
feature_used = c("Hsph1","Dnajb1","Hsp90ab1",
                 "Hspa8","Calr","Txnip")

VlnPlot(betacells, features="Calr",
        split.by = 'condition',pt.size=0,
        cols = c("white","#F71480")) +
  scale_y_continuous(limits = c(0.1,5)) +
  labs(title=feature_used) +
  theme(legend.position = 'none',
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())


#Re-Clustering data
betacells <- FindNeighbors(betacells, dims = 1:30)
betacells <- FindClusters(betacells, resolution = 0.3)
betacells <- RunUMAP(betacells, dims = 1:30, reduction = "pca")

#Supplementary Figure S12E - Visualize Beta Cell Clusters X Group
#i)
DimPlot(betacells, reduction = "umap", label = TRUE, split.by = "condition", pt.size = 2)
#ii)
#barplot
betacell_clusters <- betacells@meta.data
ggplot(betacell_clusters, aes(x = condition, fill = seurat_clusters)) +
  geom_bar(position = "fill") +
  #scale_fill_manual(values=c("#F8766D","#DB8E00","#AEA100","#00C1A7"))+
  scale_y_continuous(labels = percent_format()) +
  scale_x_discrete(labels=c("AL","CR"))+
  labs(
    x = "",
    y = "Proportion",
    fill = "Cell Type"
  ) +
  theme_classic() 

#### Supplementary Figure S12F and G ####
#Export beta cells for pseudobulk analyses ----
#See "RScript_Figure6 &_pseudobulk.R" for analyses of fig. S12F and G.

#filter low quality transcripts
betacells_filtered <- subset(betacells, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & nCount_RNA >200) #Ignore standard workflow because cells are already normalized/scaled for seurat
betacells_filtered$samples <- paste0(betacells_filtered$condition, betacells_filtered$identifier)

#Aggregates counts across steps
DefaultAssay(betacells_filtered)
pseudo_beta <- AggregateExpression(betacells_filtered,
                              group.by = c("seurat_clusters","samples"),
                              assays = 'RNA',
                              slot = "counts",
                              return.seurat = FALSE)
pseudo_beta <- pseudo_beta$RNA

#Data wrangling
pseudo_beta.t <- t(pseudo_beta) #transpose
pseudo_beta.t <- as.data.frame(pseudo_beta.t) #convert to data.frame
splitRows <- gsub('_.*', '', rownames(pseudo_beta.t)) # get values where to split
pseudo_beta.split <- split.data.frame(pseudo_beta.t, f = factor(splitRows)) # split data.frame

#fix colnames and transpose
pseudo_beta.split.modified <- lapply(pseudo_beta.split, function(x){
  rownames(x) <- gsub('.*_(.*)', '\\1', rownames(x))
  t(x)
})

counts_beta0 <- data.frame(pseudo_beta.split.modified$g0)
counts_beta1 <- data.frame(pseudo_beta.split.modified$g1)
counts_beta2 <- data.frame(pseudo_beta.split.modified$g2)
counts_beta3 <- data.frame(pseudo_beta.split.modified$g3)
counts_beta4 <- data.frame(pseudo_beta.split.modified$g4)
counts_beta5 <- data.frame(pseudo_beta.split.modified$g5)
counts_beta6 <- data.frame(pseudo_beta.split.modified$g6)


#### Supplementary Figure S12J-K ####
#Nebulosa Plot for Beta Cell Expression Intensity
#Samples must be subset into old "AL" and "CR" groups, and run independently, and scaled consistently.  Nebulosa does not factor for 

#features = "Hspa1a" & "Dnajb1"
Nebulosa::plot_density(betacells, features = "Hspa1a", size = 0.4) +
  scale_fill_gradient(low = "grey90", high = "red") +
  theme_void() +
  labs(title = "Hspa1a") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 12),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line = element_line(linewidth = 1),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "none"
  )

#Figure 6I-J
#Subset alpha cells----
alphacells <- subset(integrated_mouse_scRNAseq, idents = "Alpha", invert = FALSE)
FeaturePlot(alphacells, features = c("Gcg"), split.by = "condition")
DimPlot(alphacells, reduction = "umap", group.by = "cell_type", label = TRUE, split.by = "condition", pt.size = 2)

#Scale out unnecessary clusters/false positive
AlphaCell_Level2 <- ScaleData(alphacells, verbose = FALSE)
AlphaCell_Level2 <- RunPCA(AlphaCell_Level2, features = rownames(AlphaCell_Level2), npcs = 30, verbose = FALSE)
ElbowPlot(AlphaCell_Level2, ndims = 30)

AlphaCell_Level2 <- FindNeighbors(AlphaCell_Level2, dims = 1:20, verbose = TRUE) #15 minutes
AlphaCell_Level2 <- FindClusters(alphacells, resolution = c(0.1)) #45 minutes
AlphaCell_Level2 <- RunUMAP(AlphaCell_Level2, dims = 1:20) #40 minutes

alphacells <- AlphaCell_Level2 

DimPlot(alphacells,
        split.by = "condition",
        reduction = "umap", label = TRUE) + ggtitle("Resolution 0.6")
alphacells <- subset(alphacells, idents = c("0", "1"))


#Figure 6G & H
#Export alpha cells for pseudobulk analyses ----
#Check consistent counts per cell fraction
view(alphacells@meta.data)
table(Idents(alphacells), alphacells$condition)

#filter low quality transcripts
alphacells_filtered <- subset(alphacells, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & nCount_RNA >200) #Ignore standard workflow because cells are already normalized/scaled for seurat
alphacells_filtered$samples <- paste0(alphacells_filtered$condition, alphacells_filtered$identifier)
#view(alphacells_filtered@meta.data)

#Aggregates counts across steps
DefaultAssay(alphacells_filtered)
pseudo_alpha <- AggregateExpression(alphacells_filtered,
                                    group.by = c("cell_type","samples"),
                                    assays = 'RNA',
                                    slot = "counts",
                                    return.seurat = FALSE)
pseudo_alpha <- pseudo_alpha$RNA
# transpose
pseudo_alpha.t <- t(pseudo_alpha)
# convert to data.frame
pseudo_alpha.t <- as.data.frame(pseudo_alpha.t)
# get values where to split
splitRows <- gsub('_.*', '', rownames(pseudo_alpha.t))
# split data.frame
pseudo_alpha.split <- split.data.frame(pseudo_alpha.t,
                                      f = factor(splitRows))
# fix colnames and transpose
pseudo_alpha.split.modified <- lapply(pseudo_alpha.split, function(x){
  rownames(x) <- gsub('.*_(.*)', '\\1', rownames(x))
  t(x)
  
})
counts_alpha <- data.frame(pseudo_alpha.split.modified)
write.csv(counts_alpha, file='C:/Users/schlehmw/Box/Mike_OldMice_snRNseq_Bender/pseudobulk/AlphaCells/Alpha_FilteredV2.csv')

#Figure 6I: Alpha Cell MHC-I genes
feature_used = c("H2-Q7","H2-Q6","H2-K1")

VlnPlot(alphacells, features=feature_used,
        split.by = 'condition',pt.size=0,
        cols = c("white","#F71480")) +
  scale_y_continuous(limits = c(0.1,5)) +
  labs(title=feature_used) +
  theme(legend.position = 'none',
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

#### Figure 6J ####
#Alpha Cell MHC-I Module Score
MHCI_genes <- c("H2-K1","H2-D1","H2-Q2","H2-Q5","H2-Q6","H2-Q7","H2-Q12")
Alpha_Module <- AddModuleScore(
  alphacells,
  features = list(MHCI_genes),         #######<----- Change to Phenotyping
  name = "MHCIScore"
)

#Subset the AL for figure purposes----
Alpha_Module_AL <- subset(Alpha_Module, subset = condition == "AL_old")
Alpha_Module_AL <- FeaturePlot(
  Alpha_Module_AL,
  features = "MHCIScore1",
  reduction = "umap",
  pt.size = 0.5,
) +
  scale_y_continuous(limits=c(-5,-15),
                     labels = c(-5, -2.5, -0, 2.5, 5)) +
  scale_x_continuous(limits=c(-8,0)) +
  scale_color_viridis(option = "viridis",limits = c(-0.35, 0.9)) +
  ggtitle("AL") +
  theme(legend.position = 'none')

#Subset the CR for figure purposes
Alpha_Module_CR <- subset(Alpha_Module, subset = condition == "CR_old")

Alpha_Module_CR <- FeaturePlot(
  Alpha_Module_CR,
  features = "MHCIScore1",
  reduction = "umap",
  pt.size = 0.5,
) +
  scale_y_continuous(limits=c(-5,-15)) +
  scale_x_continuous(limits=c(-8,0)) +
  
  scale_color_viridis(option = "viridis",limits = c(-0.46, 2.2)) +
  ggtitle("CR") +
  theme(legend.position = 'none',
        axis.title.y = element_blank(),
        axis.text.y=element_blank())

Alpha_Module_AL | Alpha_Module_CR

#Module Scoring
Module <- data.frame(ID = Alpha_Module@meta.data$identifier,
                     Diet = Alpha_Module@meta.data$condition,
                     MHCI = Alpha_Module@meta.data$MHCIScore1)
Module$scale <- scale(Alpha_Module$MHCIScore1)

Module_ID <- Module %>% 
  group_by(Diet,ID) %>% 
  summarise(MHCI_mean=mean(scale),
            sem_group  = sd(scale) / sqrt(n()),)

MHCI_Module <- ggplot(Module, aes(x = Diet, y = scale)) +
  geom_boxplot(aes(fill = Diet),width = 0.6,color = "black",linewidth = 0.8,coef = 1.5)+  # Whiskers extend to 1.5 * IQR
  scale_fill_manual(values = c("white", "#F71480")) +
  geom_jitter(width = 0.12, size = 2, alpha = 0.01, color = "black") +
  # Wilcoxon p-value, centered above groups
  stat_compare_means(method = "wilcox.test", label = "p.format", hjust = 0,label.y = max(Module$scale) + 0.1) +
  labs(title = NULL, x = NULL) +
  theme_classic() +
  scale_x_discrete(labels=c("AL","CR")) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    axis.title.y = element_blank(),
    axis.text = element_text(size = 12),
    legend.position = "none"
  )
MHCI_Module



#################################################################################
########## Figure 7 - Immune Cell Phenotype & Supplementary Figure S12 ########## 
#################################################################################
#Subset immune cells----
immunecells <- subset(integrated_mouse_scRNAseq, idents = c("B-cell", "T-cell", "Macrophages", "prol_B-cell", "APC"))

#Figure 7A
DimPlot(immunecells, reduction = "umap", group.by = "cell_type", label = TRUE, pt.size = 2) +
  scale_color_brewer(palette = "Dark2")


#Figure 7B
#Calculate Total Proportion of Immune Cells----
#Proportion immune cells vs total islet cells
IC_Proportion <- read.csv("Proportion_ImmuneCells_per_islet.csv",header=T,sep=",")  #<--- Read csv file from Github
IC_Proportion$Category <- factor(IC_Proportion$Category, levels = c("B-cell", "T-cell", "Macrophages", "prol_B-cell","APC"))

ggplot(IC_Proportion, aes(x = Group, y = Value, fill = Category)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_brewer(palette = "Dark2")+
  scale_x_discrete(labels=c("AL","CR"))+
  labs(
    x = "",
    y = "Proportion of total islet cells",
    fill = "Cell type")+
  theme(plot.title = element_blank()) +
  theme_classic()

#ii) Proportion immune cells vs total immune cells
#Data not shown in manuscript
immunemetadata <- immunecells@meta.data
table(immunecells$cell_type,immunecells$condition)

ggplot(immunemetadata, aes(x = condition, fill = cell_type)) +
  scale_fill_brewer(palette = "Dark2")+
  geom_bar(position = "fill") +
  scale_y_continuous(labels = percent_format()) +
  scale_x_discrete(labels=c("AL","CR"))+
  labs(
    x = "",
    y = "Proportion",
    fill = "Cell Type"
  ) +
  theme_classic()

#### Figure 7C-ES
#### Subset T cells----
Tcells <- subset(integrated_mouse_scRNAseq, idents = "T-cell", invert = FALSE)
#Sanity check for T cells
DimPlot(Tcells, reduction = "umap", 
        group.by = "cell_type", 
        label = TRUE, 
        split.by = "condition", 
        pt.size = 2)

#Re-normalize/scale data data
Tcells <- NormalizeData(Tcells)      # Normalize expression
Tcells <- FindVariableFeatures(Tcells)  # Identify highly variable genes
Tcells <- ScaleData(Tcells)          # Scale data
Tcells <- RunPCA(Tcells)             # Perform PCA
Tcells <- FindNeighbors(Tcells, dims = 1:30)
Tcells <- FindClusters(Tcells, resolution = 0.3)  # Cluster T cells
Tcells <- RunUMAP(Tcells, dims = 1:30)  # Run UMAP for visualization

#Clustering
Tcells <- RunUMAP(Tcells, dims = 1:30, reduction = "pca")

#Figure 7C
DimPlot(Tcells,reduction = "umap",label =,cols="OrRd",
        #split.by = 'condition'
        ) +
  labs(title="T cells") +
  theme(plot.title = element_text(hjust = 0.5,face = "bold",size = 12),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = 'none',
    panel.spacing = unit(0.8, "lines")
  )
#Fig 7Cii) Split by group (not shown in manuscript)
DimPlot(Tcells, reduction = "umap", label = TRUE, split.by = "condition", pt.size = 2,cols = "OrRd")

#Identify new clusters
T_markers <- FindAllMarkers(Tcells, only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25)
#Sanity Check
FeaturePlot(Tcells, features = c("Cd8a"),split.by = "condition")

#Figure 6D
Tcellmetadata <- Tcells@meta.data
ggplot(Tcellmetadata, aes(x = condition, fill = seurat_clusters)) +
  scale_fill_brewer(palette = "OrRd")+
  geom_bar(position = "fill") +
  scale_y_continuous(labels = percent_format()) +
  scale_x_discrete(labels=c("AL","CR"))+
  labs(
    x = "",
    y = "Proportion",
    fill = "Cell Type"
  ) +
  theme_classic()

#Figure 6E
#FeaturesPlot for Tcell Fractions
FeaturePlot(
  Tcells,
  features = c("Gzmk","Ccr7","Il2ra","Txk"),
  split.by = "condition",
  cols = c("lightgray","red"),
  alpha = 0.5
) +
  scale_colour_gradient(
    limits = c(0, 1),
    low = "lightgray",
    high = "red",
    name = "Expression (0–1)"
  ) +
  scale_fill_gradient(
    limits = c(0, 1),
    low = "lightgray",
    high = "red",
    name = "Expression (0–1)"
  )


########Supplementary Figure S14F-K
#Macrophages
Macrophages <- subset(integrated_mouse_scRNAseq, idents = "Macrophages", invert = FALSE)
DimPlot(Macrophages, reduction = "umap", group.by = "cell_type", label = TRUE, split.by = "condition", pt.size = 2)

#Rescale Data
Macrophages <- ScaleData(Macrophages)          # Scale data
Macrophages <- RunPCA(Macrophages, features = VariableFeatures(Macrophages))
ElbowPlot(Macrophages,ndims=50)
Macrophages <- FindNeighbors(Macrophages, dims = 1:20)
Macrophages <- FindClusters(Macrophages, resolution = 0.6)
Macrophages <- RunUMAP(Macrophages, dims = 1:20)

DimPlot(Macrophages, reduction = "umap", label = F,cols="PuRd") +
  labs(title="Marcrophages") +
  theme(plot.title = element_text(hjust=0.5))

Mac_markers <- FindAllMarkers(Macrophages, only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25)

new.cluster.ids <- c(
  "M2-like",             #Cluster 0
  "M1-like",             #Cluster 1
  "M2-like",           #Cluster 2
  "M1-like")                #Cluster 3

table(Idents(Macrophages))
names(new.cluster.ids) <- levels(Macrophages)
Macrophages <- RenameIdents(Macrophages, new.cluster.ids)

#Supplementary Figure S12G
#Whole islet macrophage marker
Nebulosa::plot_density(integrated_mouse_scRNAseq,features = "Adgre1",size = 0.4) +
  #scale_color_viridis_c(option = "viridis") +
  scale_fill_viridis_c(option = "viridis") +
  # scale_x_continuous(limits = c(-12, 13)) +
  # scale_y_continuous(limits = c(-15, 8)) +
  theme_void() +
  labs(title = "F4/80 / Adgre1") +
  theme(plot.title = element_text(hjust = 0.5, face='bold.italic',size=12),
        panel.border = element_rect(color = "black",fill = NA,linewidth = 1),
        axis.line = element_line(linewidth = 1),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")

#Macrophage-specific M1- vs M2-like markers
DimPlot(Macrophages, reduction = "umap", 
        #split.by = "c",
        label = TRUE, 
        repel = F,
        cols = c("#4E79A7","#F28E2B"))

#Supplementary Figure 13H
#Dotplot for macrophage markers
MarkerIdentifiers <- c("Nos1", #M2
                       "Cxcl12", #M2
                       "Mrc1", #M2
                       "Il10", #M2
                       "Ccl5",
                       "Il6", #M1
                       "Ly6a", #M1 
                       "Il12rb2") #M1

DotPlot(Macrophages, features = MarkerIdentifiers, dot.min = 0, col.min = 0, dot.scale=10) +
  coord_flip() +
  RotatedAxis() +
  # scale_y_discrete(
    # limits = c("B naive", "CD4 naive", "Treg", "Macrophage","B memory","CD8 effector","B proliferating","B activated",
    #            "DC","T proliferating","NK","B plasma")) +
  geom_point(aes(size = pct.exp), shape = 21, colour = "black", stroke = 0.5) +
  scale_fill_viridis_c() +
  labs(x = NULL, y = NULL) +
  scale_colour_viridis(option = "viridis") +
  theme(plot.margin = unit(c(1, 1, 1, 2), "cm"),
        axis.text.y=element_text(face="italic"))

#Supplementary Figure 13I
#Macrophage populations by polarization
Macrophage <- Macrophages_metadata %>% filter(Cluster %in% c("M1-like","M2-like"))
Macrophages_metadata$Cluster <- Idents(Macrophages)

my_colors <- c("M2-like" = "#4E79A7",
               "M1-like" = "#F28E2B")

ggplot(Macrophages_metadata, aes(x = condition, fill = Cluster)) +
  geom_bar(position = "fill",color='black') +
  scale_fill_manual(values = my_colors) +
  scale_y_continuous(labels = percent_format()) +
  scale_x_discrete(labels=c("AL","CR"))+
  labs(
    title = "Macrophages",
    x = "",
    y = "% Macrophages",
    fill = "population") +
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5),
        legend.position='none')

#Supplementary Figure 13J
#Bargraph Per Group
Macrophage_Percent_Population <- read.csv("Macrophages/Macrophage_Percent_Population.csv")  #<--- Found in Github
Summary_Group <- Macrophage_Percent_Population %>% 
  group_by(Group,Cluster) %>%
  summarise(mean_prop=mean(Percent),
            sem = sd(Percent)/1.44)

ggplot() +
  geom_bar(data = Summary_Group,aes(x = Group, y = mean_prop, fill = Group),
           stat = "identity",color = "black",width = 0.65,alpha = 0.75) +
  geom_errorbar(data = Summary_Group,aes(x = Group,ymin = mean_prop - sem,ymax = mean_prop + sem),
                width = 0.2) +
  geom_jitter(data = Macrophage_Percent_Population,aes(x = Group, y = Percent),
              width = 0.12, size = 2, alpha = 0.5) +
  facet_wrap2(~ Cluster,nrow = 1,scales = "free_y",
              strip = strip_themed(background_x = elem_list_rect(fill = my_colors,color = "black",linewidth = 0.5),
                                   text_x = elem_list_text(face = "bold", color = "black"))) +
  scale_fill_manual(values=c("white","#F71480")) +
  labs(y = "% Macrophage", x = NULL) +
  theme_classic() +
  theme(
    legend.position = "none",
    strip.background = element_blank())



#M1 to M2 Ratio
M1_M2 <- read.csv("Macrophage_M1_M2 Ratio.csv")

M1_M2_Bar <- M1_M2 %>% 
  group_by(Group) %>%
  summarise(mean_prop=mean(M1_M2_Ratio),
            sem = sd(M1_M2_Ratio)/1.44)

ggplot() +
  geom_bar(data = M1_M2_Bar,aes(x = Group, y = mean_prop, fill = Group),
           stat = "identity",color = "black",width = 0.65,alpha = 0.75) +
  geom_errorbar(data = M1_M2_Bar,aes(x = Group,ymin = mean_prop - sem,ymax = mean_prop + sem),
                width = 0.2) +
  geom_jitter(data = M1_M2,aes(x = Group, y = M1_M2_Ratio),
              width = 0.12, size = 2, alpha = 0.5) +
  scale_fill_manual(values=c("white","#F71480")) +
  labs(y = "ratio", x = NULL) +
  theme_classic() +
  theme(
    legend.position = "none",
    strip.background = element_blank())

ggsave(filename = "Immune Cells/Macrophages/M1_M2_ratio.png",
       dpi = 300,width = 2,height = 2.5,units = "in",bg = "white")

##### Not Used in Manuscript - B Cells
##### BCells -----
BCells <- subset(integrated_mouse_scRNAseq, idents = c("B-cell","prol_B-cell"), invert = FALSE)
DimPlot(BCells, reduction = "umap", group.by = "cell_type", label = TRUE, split.by = "condition", pt.size = 2)

BCells <- ScaleData(BCells)                                            #Scale data
BCells <- RunPCA(BCells, features = VariableFeatures(BCells))          #RunPCA
ElbowPlot(BCells,ndims=50)                                             #Elbow (30dims)
BCells <- FindNeighbors(BCells, dims = 1:30)
BCells <- FindClusters(BCells, resolution = 0.2)
BCells <- RunUMAP(BCells, dims = 1:30)

DimPlot(BCells, reduction = "umap", label = F,cols="Spectral") +
  labs(title="B Cells") +
  theme(plot.title = element_text(hjust=0.5))

BCell_markers <- FindAllMarkers(BCells, only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25)
write.csv(BCell_markers, "Immune Cells/B Cells/BCell_ClusterMarkers.csv")

new.cluster.ids <- c(
  "B Follicular",    #Cluster 0
  "ABC",             #Cluster 1
  "B Proliferating",             #Cluster 2
  "B Proliferating",             #Cluster 3
  "ABC",        #Cluster 4
  "B Plasma",        #Cluster 5
  "B Plasma",             #Cluster 6
  "B Plasma")                  #Cluster 7
  

table(Idents(BCells))
names(new.cluster.ids) <- levels(BCells)
BCells <- RenameIdents(BCells, new.cluster.ids)



DimPlot(BCells, reduction = "umap", 
        #split.by = "condition",
        label = F, 
        repel = F,
        cols = "Spectral")


#Nebulosa Plot for Macrophages
Nebulosa::plot_density(integrated_mouse_scRNAseq,features = "Ki67",size = 0.4) +  #<--- Ki67 is example marker for proliferating B Cells
  #scale_color_viridis_c(option = "viridis") +
  scale_fill_viridis_c(option = "viridis") +
  # scale_x_continuous(limits = c(-12, 13)) +
  # scale_y_continuous(limits = c(-15, 8)) +
  theme_void() +
  labs(title = "Ms4a1") +
  theme(plot.title = element_text(hjust = 0.5, face='bold.italic',size=12),
        panel.border = element_rect(color = "black",fill = NA,linewidth = 1),
        axis.line = element_line(linewidth = 1),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")



#Dotplot for BCell markers
MarkerIdentifiers <- c("Ighd", #B Follicular
                       "Cd74", #B Follicular
                       "H2-Aa", #B Follicular
                       "Itgax", #ABC
                       "Tbx21",#ABC
                       "Psmb8",#
                       "Mki67", #B Proliferating
                       "Mcm6", #B Proliferating
                       "Prdm1", #Plasma B
                       "Sdc1", #Plasma B
                       "Igkc")  #Plasma B

DotPlot(BCells, features = MarkerIdentifiers, dot.min = 0, col.min = 0, dot.scale=10) +
  coord_flip() +
  RotatedAxis() +
  # scale_y_discrete(
  # limits = c("B naive", "CD4 naive", "Treg", "Macrophage","B memory","CD8 effector","B proliferating","B activated",
  #            "DC","T proliferating","NK","B plasma")) +
  geom_point(aes(size = pct.exp), shape = 21, colour = "black", stroke = 0.5) +
  scale_fill_viridis_c() +
  labs(x = NULL, y = NULL) +
  scale_colour_viridis(option = "viridis") +
  theme(plot.margin = unit(c(1, 1, 1, 2), "cm"),
        axis.text.y=element_text(face="italic",),
        legend.position = 'none')

#BCell subpopulations subpopulations of total B Cells:
BCell_Metadata <- data.frame(BCells@meta.data)
BCell_Metadata$Cluster <- Idents(BCells)

ggplot(BCell_Metadata, aes(x = condition, fill = Cluster)) +
  geom_bar(position = "fill",color='black') +
  scale_fill_brewer(palette = "Spectral") +
scale_y_continuous(labels = percent_format()) +
  scale_x_discrete(labels=c("AL","CR"))+
  labs(
    title = NULL,
    x = "",
    y = "% B Cells",
    fill = "population") +
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5),
        legend.position='none'
        )

#B Cell Bargraph Per Group
#3B_Proportion of all clusters by individual donor ----
BCell_Metadata <- as.data.frame(BCells@meta.data)
BCell_Metadata$Cluster <- Idents(BCells)

prop_subject <- BCell_Metadata %>%
  dplyr::count(identifier, condition, Cluster) %>%
  group_by(identifier) %>%
  mutate(total_cells = sum(n), prop = n / total_cells) %>%
  ungroup()

prop_summary <- prop_subject %>%
  group_by(condition, Cluster) %>%
  summarise(mean_prop = mean(prop),
            sem_prop  = sd(prop) / sqrt(n()),
            .groups = "drop")

n_facets <- length(unique(prop_summary$Cluster))
mycolors <- c("#D7191C","#FDAE61","#ABDDA4","#2B83BA")

ggplot() +
  geom_bar(data = prop_summary, aes(x = condition, y = mean_prop,fill=condition), stat = "identity", color = "black", width = 0.65, alpha = 0.75) +
  geom_errorbar(data = prop_summary, aes(x = condition, ymin = mean_prop - sem_prop, ymax = mean_prop + sem_prop), width = 0.2) +
  geom_jitter(data = prop_subject, aes(x = condition, y = prop), width = 0.12, size = 2, alpha = 0.5) +
  facet_wrap2(~ Cluster, nrow = 1, scales = "free_y", strip = strip_themed(background_x = elem_list_rect(fill = mycolors, color = "black",
                                                                                                          linewidth = 0.5), text_x = elem_list_text(face = "bold",color = "black"))) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values=c("white","#F71480")) +
  scale_x_discrete(labels=c("AL","CR")) +
  labs(y = "% B cells", x = NULL) +
  theme_classic() +
  theme(
    legend.position = "none",
    strip.background = element_blank()
  )
ggsave(
  filename = "Immune Cells/B Cells/BCell_Marker_XCluster.png",  # Change path if needed
  width = 6, height = 2.5, dpi = 300
)

#Statistics
ttest_results <- prop_subject %>%
  mutate(prop_log1p = log1p(prop)) %>%
  group_by(Cluster) %>%
  t_test(prop_log1p ~ condition)   # Welch t-test by default


######################################################
##################### CELL CHAT ###################### 
####### Figure 7I-J & Supplementary Figure S15 ####### 
######################################################

#Open CellChat Files
options(stringsAsFactors = FALSE)
seurat_object <- readRDS("integrated_mouse_scRNAseq.rds") 
class(seurat_object)

counts <- GetAssayData(seurat_object, assay = "integrated", slot = "data")
counts[counts < 0] <- 0
seurat_object[["integrated"]]@data <- counts
min_val <- min(GetAssayData(seurat_object, assay = "integrated", slot = "data"))
print(min_val)  # Should be 0 or greater

##SUBSETTING OBJECT
data <- SetIdent(seurat_object, value = seurat_object@meta.data$condition)
ND_data <- subset(x = data, idents = "AL_old")
table(ND_data@meta.data$cell_type)

CR_data <- subset(x = data, idents = "CR_old")
table(CR_data@meta.data$cell_type)

#Sanity Checks
DimPlot(data, reduction = "umap", label = FALSE, pt.size = 1, split.by = "condition")
DimPlot(ND_data, reduction = "umap", label = TRUE, pt.size = 3) + NoLegend()
DimPlot(CR_data, reduction = "umap", label = TRUE, pt.size = 3) + NoLegend()

##CREATING CELLCHAT OBJECT
cellChat_ND <- createCellChat(object = ND_data, group.by = "cell_type", assay = "integrated")
cellChat_ND <- updateCellChat(cellChat_ND)
cellChat_ND

cellChat_CR <- createCellChat(object = CR_data, group.by = "cell_type", assay = "integrated")
cellChat_CR <- updateCellChat(cellChat_CR)
cellChat_CR

## Set the ligand-receptor interaction database
#Before users can employ CellChat to infer cell-cell communication, they need to set the ligand-receptor interaction database and identify over-expressed ligands or receptors. 
#Our database CellChatDB is a manually curated database of literature-supported ligand-receptor interactions in both human and mouse.
#CellChatDB v2 contains ~3,300 validated molecular interactions, including ~40% of secrete autocrine/paracrine signaling interactions, ~17% of extracellular matrix (ECM)-receptor interactions, ~13% of cell-cell contact interactions and ~30% non-protein signaling.
#Compared to CellChatDB v1, CellChatDB v2 adds more than 1000 protein and non-protein interactions such as metabolic and synaptic signaling. It should be noted that for molecules that are not directly related to genes measured in scRNA-seq, CellChat v2 estimates the expression of ligands and receptors using those molecules’ key mediators or enzymes for potential communication mediated by non-proteins. 
#CellChatDB v2 also adds additional functional annotations of ligand-receptor pairs, such as UniProtKB keywords (including biological process, molecular function, functional class, disease, etc), subcellular location and relevance to neurotransmitter. 
#Users can update CellChatDB by adding their own curated ligand-receptor pairs. Please check the [tutorial on updating the ligand-receptor interaction database CellChatDB](https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/Update-CellChatDB.html). 
#When analyzing human samples, use the database **`CellChatDB.human`**; when analyzing mouse samples, use the database **`CellChatDB.mouse`**. CellChatDB categorizes ligand-receptor pairs into different types, including “Secreted Signaling”, “ECM-Receptor”, “Cell-Cell Contact” and “Non-protein Signaling”. By default, the “Non-protein Signaling” are not used. 

# use CellChatDB.mouse if running on mouse data
CellChatDB <- CellChatDB.mouse 
showDatabaseCategory(CellChatDB)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling

# Only uses the Secreted Signaling from CellChatDB v1
#  CellChatDB.use <- subsetDB(CellChatDB, search = list(c("Secreted Signaling"), c("CellChatDB v1")), key = c("annotation", "version"))

# use all CellChatDB except for "Non-protein Signaling" for cell-cell communication analysis
# CellChatDB.use <- subsetDB(CellChatDB)

# use all CellChatDB for cell-cell communication analysis
# simply use the default CellChatDB. We do not suggest to use it in this way because CellChatDB v2 includes "Non-protein Signaling" (i.e., metabolic and synaptic signaling). 
CellChatDB.use <- CellChatDB 

# set the used database in the object
cellChat_ND@DB <- CellChatDB.use
cellChat_CR@DB <- CellChatDB.use

##### Preprocessing the expression data for cell-cell communication analysis
#To infer the cell state-specific communications, CellChat identifies over-expressed ligands or receptors in one cell group and then identifies over-expressed ligand-receptor interactions if either ligand or receptor are over-expressed. 
#We also provide a function to project gene expression data onto protein-protein interaction (PPI) network. Specifically, a diffusion process is used to smooth genes’ expression values based on their neighbors’ defined in a high-confidence experimentally validated protein-protein network. 
#This function is useful when analyzing single-cell data with shallow sequencing depth because the projection reduces the dropout effects of signaling genes, in particular for possible zero expression of subunits of ligands/receptors. One might be concerned about the possible artifact introduced by this diffusion process, however, it will only introduce very weak communications. By default CellChat uses the raw data (i.e., `object@data.signaling`) instead of the projected data. To use the projected data, users should run the function `projectData` before running `computeCommunProb`, and then set `raw.use = FALSE` when running `computeCommunProb`. 

# subset the expression data of signaling genes for saving computation cost
# This step is necessary even if using the whole database
cellChat_ND <- subsetData(cellChat_ND) 
cellChat_CR <- subsetData(cellChat_CR)
future::plan("multisession", workers = 12) # do parallel
cellChat_ND <- identifyOverExpressedGenes(cellChat_ND)
cellChat_CR <- identifyOverExpressedGenes(cellChat_CR)
cellChat_ND <- identifyOverExpressedInteractions(cellChat_ND)
cellChat_CR <- identifyOverExpressedInteractions(cellChat_CR)

# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
#cellChat_ND <- projectData(cellChat_ND, PPI.mouse)
#cellChat_CR <- projectData(cellChat_CR, PPI.mouse)

# Part II: Inference of cell-cell communication network
#CellChat infers the biologically significant cell-cell communication by assigning each interaction with a probability value and peforming a permutation test. CellChat models the probability of cell-cell communication by integrating gene expression with prior known knowledge of the interactions between signaling ligands, receptors and their cofactors using the law of mass action.
#CAUTION: The number of inferred ligand-receptor pairs clearly depends on the **method for calculating the average gene expression per cell group**. By default, CellChat uses a statistically robust mean method called 'trimean', which produces fewer interactions than other methods. However, we find that CellChat performs well at predicting stronger interactions, which is very helpful for narrowing down on interactions for further experimental validations. In `computeCommunProb`, we provide an option for using other methods, such as 5% and 10% truncated mean, to calculating the average gene expression. Of note, 'trimean' approximates 25% truncated mean, implying that the average gene expression is zero if the percent of expressed cells in one group is less than 25%. To use 10% truncated mean, USER can set `type = "truncatedMean"` and `trim = 0.1`. To determine a proper value of trim, CellChat provides a function `computeAveExpr`, which can help to check the average expression of signaling genes of interest, e.g, `computeAveExpr(cellchat, features = c("CXCL12","CXCR4"), type =  "truncatedMean", trim = 0.1)`. Therefore, if well-known signaling pathways in the studied biological process are not predicted, users can try `truncatedMean` with lower values of `trim` to change the method for calculating the average gene expression per cell group. 
#When analyzing unsorted single-cell transcriptomes, under the assumption that abundant cell populations tend to send collectively stronger signals than the rare cell populations, CellChat can also consider the effect of cell proportion in each cell group in the probability calculation. USER can set `population.size = TRUE`. 

## Compute the communication probability and infer cellular communication network

ptm = Sys.time()
options(future.globals.maxSize= 891289600)
future::plan("multisession", workers = 12) # do parallel

cellChat_ND@idents <- droplevels(cellChat_ND@idents)
cellChat_ND@meta$ident <- droplevels(cellChat_ND@meta$ident)
cellChat_ND <- computeCommunProb(cellChat_ND, type = "triMean", population.size = TRUE, raw.use = TRUE)

cellChat_CR@idents <- droplevels(cellChat_CR@idents)
cellChat_CR@meta$ident <- droplevels(cellChat_CR@meta$ident)
cellChat_CR <- computeCommunProb(cellChat_CR, type = "triMean", population.size = TRUE, raw.use = TRUE)

#The key parameter for this analysis is `type`, the method for computing the average gene expression per cell group. By default `type = "triMean"`, producing fewer but stronger interactions. When setting `type = "truncatedMean"`, a value should be assigned to `trim`, producing more interactions. Please check above in detail on the **method for calculating the average gene expression per cell group**.
#Users can filter out the cell-cell communication if there are only few cells in certain cell groups. By default, the minimum number of cells required in each cell group for cell-cell communication is 10. 

cellChat_ND <- filterCommunication(cellChat_ND, min.cells = 50)
cellChat_CR <- filterCommunication(cellChat_CR, min.cells = 50)


## Infer the cell-cell communication at a signaling pathway level
#CellChat computes the communication probability on signaling pathway level by summarizing the communication probabilities of all ligands-receptors interactions associated with each signaling pathway.  
#NB: The inferred intercellular communication network of each ligand-receptor pair and each signaling pathway is stored in the slot 'net' and 'netP', respectively.

cellChat_ND <- computeCommunProbPathway(cellChat_ND)
cellChat_CR <- computeCommunProbPathway(cellChat_CR)

## Calculate the aggregated cell-cell communication network 
#CellChat calculates the aggregated cell-cell communication network by counting the number of links or summarizing the communication probability. Users can also calculate the aggregated network among a subset of cell groups by setting `sources.use` and `targets.use`.
cellChat_ND <- aggregateNet(cellChat_ND)
cellChat_CR <- aggregateNet(cellChat_CR)

### (A)	Compute and visualize the network centrality scores
# Compute the network centrality scores
future::plan("multisession", workers = 1) # do parallel
options(future.seed=TRUE)
cellChat_ND <- netAnalysis_computeCentrality(cellChat_ND, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
cellChat_CR <- netAnalysis_computeCentrality(cellChat_CR, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways

#Supplementary Figure S15A
#1 - Interactions and Interaction Strengths Graph
ptm = Sys.time()
groupSize_ND <- as.numeric(table(cellChat_ND@idents))
groupSize_CR <- as.numeric(table(cellChat_CR@idents))
par(mfrow = c(1,2), xpd=TRUE)
AL_circle <- netVisual_circle(cellChat_ND@net$count, vertex.weight = groupSize_ND, weight.scale = T, label.edge= F, title.name = "AL - Number of interactions")
AL_circle <- netVisual_circle(cellChat_ND@net$weight, vertex.weight = groupSize_ND, weight.scale = T, label.edge= F, title.name = "AL - Interaction weights/strength")

netVisual_circle(cellChat_CR@net$count, vertex.weight = groupSize_CR, weight.scale = T, label.edge= F, title.name = "CR - Number of interactions")
netVisual_circle(cellChat_CR@net$weight, vertex.weight = groupSize_CR, weight.scale = T, label.edge= F, title.name = "CR - Interaction weights/strength")

##Pathway Breakdown - AL
mat_ND <- cellChat_ND@net$weight
par(mfrow = c(3,8), xpd=TRUE)
for (i in 1:nrow(mat_ND)) {
  mat2 <- matrix(0, nrow = nrow(mat_ND), ncol = ncol(mat_ND), dimnames = dimnames(mat_ND))
  mat2[i, ] <- mat_ND[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize_ND, weight.scale = T, edge.weight.max = max(mat_ND), title.name = rownames(mat_ND)[i])
}

##Pathway Breakdown - CR
mat_CR <- cellChat_CR@net$weight
par(mfrow = c(3,8), xpd=TRUE)
for (i in 1:nrow(mat_CR)) {
  mat2 <- matrix(0, nrow = nrow(mat_CR), ncol = ncol(mat_CR), dimnames = dimnames(mat_CR))
  mat2[i, ] <- mat_CR[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize_CR, weight.scale = T, edge.weight.max = max(mat_CR), title.name = rownames(mat_CR)[i])
}

##Part III: Visualization of cell-cell communication network

pathways.show <- c("CD45")
cellChat_ND@netP$pathways
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
#vertex.receiver = seq(1,4) # a numeric vector. 
#netVisual_aggregate(cellChat_ND, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellChat_ND, signaling = pathways.show, layout = "chord")
netVisual_aggregate(cellChat_CR, signaling = pathways.show, layout = "chord")

ggsave(filename = ".png",
  width = 3,
  height = 3,
  dpi = 300)



#################### END #################### 