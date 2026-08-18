# =======================================================================================================================
# Alpha cell inflammation during human pancreas aging and type 2 diabetes and its reversal by calorie restriction in mice
# -----------------------------------------------------------------------------------------------------------------------
# Author:      Michael Schleh
# Repository:  https://github.com/ArrojoDrigoLab/AlphaCellAging
# Journal:     Science Advances

# Description:
#   Generates Figure 4C-D,G,H-I, and Supplementary Figure 6 and 7
#     - A) RDS setup, cleaning, and single cell clustering
#     - B) Second-level single-cell clustering
#     - C) Population analysis

# Usage:
#   1. Download packages and set the working directory.
#   2. Ensure the required input files exist in `Figure 4/' (see README).
#   3. Run: RScript_Figure 4_CODEX (Young vs Old).R
#
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Setup: Working 
# ----------------------------------------------------------------------------
#1A. Working Directory 
setwd("")    #<--- Set to source file location

#1B. Download necessary packages
install.packages(c(
  "tidyverse", "stringr", "reshape2", "patchwork",
  "readxl", "pheatmap", "RColorBrewer", "broom",
  "Matrix", "Seurat", "viridis", "future",
  "tibble", "devtools", "ggh4x"
))

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("Nebulosa"))

#1B. Library packages
library(tidyverse)
library(stringr)
library(reshape2)
library(patchwork)
library(readxl)
library(pheatmap)
library(RColorBrewer)
library(broom)
library(Matrix)
library(Seurat)
library(viridis)
library(future)
library(tibble)
library(devtools)
library(ggh4x)
library(Nebulosa)

#1C. Future settings
# Set future memory limit to 10 GB
options(future.globals.maxSize = 10 * 1024^3)
# Use 10 parallel workers
plan("multisession", workers = 10)

#Source Metadata
PancDB_Metadata <- read.csv("MetaData_AgingCohort.csv",header=T,sep=",")
Images_used <- c("HPAP-013","HPAP-016","HPAP-017","HPAP-027","HPAP-040","HPAP-045",
                 "HPAP-053","HPAP-056","HPAP-066","HPAP-093","HPAP-094","HPAP-095",
                 "HPAP-112","HPAP-155","HPAP-169","HPAP-174","HPAP-175")

Metadata_measure <- c("age_years","bmi","hba1c","c_peptide_ng_ml")

  PancDB_Metadata <- PancDB_Metadata %>% select(c("donor_ID","sex","bmi",
                                                "age_years","hba1c","c_peptide_ng_ml",
                                                "allocation_via","cause_of_death","liver_pancreas_removal_time_minutes",
                                                "pancreas_condition","prs_score"))
  
PancDB_Metadata <- PancDB_Metadata %>% filter(donor_ID %in% Images_used)
PancDB_Metadata$Group <- ifelse(PancDB_Metadata$age_years>40,"old","young")
PancDB_Metadata$Group <- factor(PancDB_Metadata$Group, levels = c("young","old"))

  
Metadata_Summary <- PancDB_Metadata %>%
  group_by(Group) %>%
  summarise(across(all_of(Metadata_measure),
    list(mean=~round(mean(.x, na.rm=T), 2),
         sd=~round(sd(.x,na.rm=T), 2),
         sem=~round(sd(.x,na.rm=T) / sqrt(sum(!is.na(.x))),2))))

#Supplementary Figure 6B
#i) Age
ggplot() +
  geom_col(data=Metadata_Summary, aes(x=Group,y=age_years_mean,fill=Group),color='black') +
  geom_errorbar(data=Metadata_Summary,mapping = aes(x=Group, ymin = age_years_mean-age_years_sem, ymax = age_years_mean+age_years_sem,),width=0.5) +
  geom_jitter(data=PancDB_Metadata, aes(x=Group,y=age_years),width=0.2,alpha=0.5) +
  scale_fill_manual(values=c("#E0DCB8","#607C83"))+
  labs(title="Age",
       y="years")+
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5),
        legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=12))

#Supplementary Figure 6B
#ii) BMI
ggplot() +
  geom_col(data=Metadata_Summary, aes(x=Group,y=bmi_mean,fill=Group),color='black') +
  geom_errorbar(data=Metadata_Summary,mapping = aes(x=Group, ymin = bmi_mean-bmi_sem, ymax = bmi_mean+bmi_sem,),width=0.5) +
  geom_jitter(data=PancDB_Metadata, aes(x=Group,y=bmi),width=0.2,alpha=0.5) +
  scale_fill_manual(values=c("#E0DCB8","#607C83"))+
  labs(title="BMI",
       y=bquote("kg/m"^2))+
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5),
        legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=12))

#Supplementary Figure 6B
#iii) C-Peptide
ggplot() +
  geom_col(data=Metadata_Summary, aes(x=Group,y=c_peptide_ng_ml_mean,fill=Group),color='black') +
  geom_errorbar(data=Metadata_Summary,mapping = aes(x=Group, ymin = c_peptide_ng_ml_mean-c_peptide_ng_ml_sem, ymax = c_peptide_ng_ml_mean+c_peptide_ng_ml_sem,),width=0.5) +
  geom_jitter(data=PancDB_Metadata, aes(x=Group,y=c_peptide_ng_ml),width=0.2,alpha=0.5) +
  scale_fill_manual(values=c("#E0DCB8","#607C83"))+
  labs(title="C-Peptide",
       y="ng/mL")+
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5),
        legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=12))
  
#Supplementary Figure 6B
#iii) HbA1c
ggplot() +
  geom_col(data=Metadata_Summary, aes(x=Group,y=hba1c_mean,fill=Group),color='black') +
  geom_errorbar(data=Metadata_Summary,mapping = aes(x=Group, ymin = hba1c_mean-hba1c_sem, ymax = hba1c_mean+hba1c_sem,),width=0.5) +
  geom_jitter(data=PancDB_Metadata, aes(x=Group,y=hba1c),width=0.2,alpha=0.5) +
  scale_fill_manual(values=c("#E0DCB8","#607C83"))+
  labs(title="HbA1c",
       y="%")+
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5),
        legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=12))

  
  
#Upload per-cell intensity data from Zenodo (https://doi.org/10.5281/zenodo.20496994)
CODEX <- read.table(file="Figure 4_CODEX_Cell_measurements_merged")



#Metadata and pre-processing
CODEX_Metadata <- CODEX[-1] %>%
    mutate(ID = substr(Image, 1, 8)) %>%
    tibble::rownames_to_column("CellID")  # ensure each row has a unique ID
  
CODEX_Metadata <- data.frame(CellID=rownames(CODEX_Metadata),
                             Image=CODEX_Metadata$Image,
                             CentroidX=CODEX_Metadata$Centroid.X.µm,
                             Centroidy=CODEX_Metadata$Centroid.Y.µm,
                             DonorID=CODEX_Metadata$ID,
                             Group=CODEX_Metadata$group_name,
                             row.names = rownames(CODEX_Metadata))

raw_markers <- c("GCG..Cell..Mean",             #Alpha
                 "CD45..Cell..Mean",            #Immune
                 "CD19..Cell..Mean",            #B cells
                 "CD39L3..Cell..Mean",          #ENTPD
                 "GP2..Cell..Mean",             #Pancreatic Progenitor
                 "NKX2.2..Nucleus..Mean",       #Pancreas
                 "CHGA..Cell..Mean",            #Endocrine
                 "PD.L1..Cell..Mean",           #Immune Cell Exhaustion
                 "ATP1A1..Cell..Mean",          #Na+/K+-ATPase
                 "SELP..Cell..Mean",            #Endothelial
                 "CD90..Cell..Mean",            #MSC 
                 "CD141..Cell..Mean",           #cD1Cs
                 "COL4A1..Cell..Mean",          #ECM
                 "COL6..Cell..Mean",            #ECM
                 "CD3..Cell..Mean",             #T cells
                 "CD8..Cell..Mean",             #CD8 T cells
                 "SST..Cell..Mean",             #Delta cells
                 "CD4..Cell..Mean",             #CD4T cells
                 "ACTA2..Cell..Mean",           #Smooth Muscle
                 "COL1A1..Cell..Mean",          #Connective Tissue 
                 "PAX6..Nucleus..Mean",         #Islet
                 "PPY..Cell..Mean",             #PP cells
                 "SOX9..Nucleus..Mean",         #Pancreas Progenitor
                 "CD11c..Cell..Mean",           #Dendritic Cells
                 "NKX6.1..Nucleus..Mean",       #Beta Cells
                 "KRT..Cell..Mean",             #Epithelial
                 "EPCAM..Cell..Mean",           #Epithelial
                 "CD31..Cell..Mean",            #Endothelial
                 "PDX1..Nucleus..Mean",         #Pancreas
                 "GHRL..Cell..Mean",            #Ghrelin
                 "CD163..Cell..Mean",           #Macrophage
                 "HLA.DR..Cell..Mean",          #MHC-II
                 "CD68..Cell..Mean",            #Macrophage
                 "CD66b..Cell..Mean",           #Granulocytes
                 "PNLIP..Cell..Mean",           #Acinar cells
                 "CPEP..Cell..Mean",            #Beta Cells
                 "TUBB3..Cell..Mean",           #Neuron
                 "Ki67..Nucleus..Mean",         #Proliferation Status
                 "LYVE1..Cell..Mean")           #Lymphatic

CODEX_marker <- CODEX %>% select(matches(raw_markers))
CODEX_Final <- data.frame(CODEX_Metadata,
                          CODEX_marker)

#HPAP-013_Body
{
HPAP013_Body <- CODEX_Final %>% filter(Image=="HPAP-013_CODEX_Body-of-pancreas_OCT.ome.tif")
HPAP013_Body <- CreateSeuratObject(counts = t(HPAP013_Body[,-c(1:6)]), assay = "CODEX",meta.data = HPAP013_Body[,c(1:6)])
HPAP013_Body <- NormalizeData(HPAP013_Body)
HPAP013_Body <- FindVariableFeatures(HPAP013_Body)
HPAP013_Body <- ScaleData(HPAP013_Body)
HPAP013_Body <- RunPCA(HPAP013_Body)
}

#HPAP-016_Body
{
HPAP016_Body <- CODEX_Final %>% filter(Image=="HPAP-016_CODEX_Body-of-pancreas_OCT.ome.tif")
HPAP016_Body <- CreateSeuratObject(counts = t(HPAP016_Body[,-c(1:6)]), assay = "CODEX",meta.data = HPAP016_Body[,c(1:6)])
HPAP016_Body <- NormalizeData(HPAP016_Body)
HPAP016_Body <- FindVariableFeatures(HPAP016_Body)
HPAP016_Body <- ScaleData(HPAP016_Body)
HPAP016_Body <- RunPCA(HPAP016_Body)
}

#HPAP-017_Body
{
  HPAP017_Body <- CODEX_Final %>% filter(Image=="HPAP-017_CODEX_Body-of-pancreas_OCT.ome.tif")
  HPAP017_Body <- CreateSeuratObject(counts = t(HPAP017_Body[,-c(1:6)]), assay = "CODEX",meta.data = HPAP017_Body[,c(1:6)])
  HPAP017_Body <- NormalizeData(HPAP017_Body)
  HPAP017_Body <- FindVariableFeatures(HPAP017_Body)
  HPAP017_Body <- ScaleData(HPAP017_Body)
  HPAP017_Body <- RunPCA(HPAP017_Body)
}

#HPAP-027_Body
{
  HPAP027_Body <- CODEX_Final %>% filter(Image=="HPAP-027_CODEX_Head-of-pancreas_OCT.ome.tif")
  HPAP027_Body <- CreateSeuratObject(counts = t(HPAP027_Body[,-c(1:6)]), assay = "CODEX",meta.data = HPAP027_Body[,c(1:6)])
  HPAP027_Body <- NormalizeData(HPAP027_Body)
  HPAP027_Body <- FindVariableFeatures(HPAP027_Body)
  HPAP027_Body <- ScaleData(HPAP027_Body)
  HPAP027_Body <- RunPCA(HPAP027_Body)
}

#HPAP-040_Head
{
  HPAP040_Head <- CODEX_Final %>% filter(Image=="HPAP-040_CODEX_Head-of-pancreas_OCT.ome.tif")
  HPAP040_Head <- CreateSeuratObject(counts = t(HPAP040_Head[,-c(1:6)]), assay = "CODEX",meta.data = HPAP040_Head[,c(1:6)])
  HPAP040_Head <- NormalizeData(HPAP040_Head)
  HPAP040_Head <- FindVariableFeatures(HPAP040_Head)
  HPAP040_Head <- ScaleData(HPAP040_Head)
  HPAP040_Head <- RunPCA(HPAP040_Head)
}

#HPAP-045_Tail
{
  HPAP045_Tail <- CODEX_Final %>% filter(Image=="HPAP-045_CODEX_Tail-of-pancreas_OCT.ome.tif")
  HPAP045_Tail <- CreateSeuratObject(counts = t(HPAP045_Tail[,-c(1:6)]), assay = "CODEX",meta.data = HPAP045_Tail[,c(1:6)])
  HPAP045_Tail <- NormalizeData(HPAP045_Tail)
  HPAP045_Tail <- FindVariableFeatures(HPAP045_Tail)
  HPAP045_Tail <- ScaleData(HPAP045_Tail)
  HPAP045_Tail <- RunPCA(HPAP045_Tail)
}

#HPAP-053_Head
{
  HPAP053_Head <- CODEX_Final %>% filter(Image=="HPAP-053_CODEX_Head-of-pancreas_OCT.ome.tif")
  HPAP053_Head <- CreateSeuratObject(counts = t(HPAP053_Head[,-c(1:6)]), assay = "CODEX",meta.data = HPAP053_Head[,c(1:6)])
  HPAP053_Head <- NormalizeData(HPAP053_Head)
  HPAP053_Head <- FindVariableFeatures(HPAP053_Head)
  HPAP053_Head <- ScaleData(HPAP053_Head)
  HPAP053_Head <- RunPCA(HPAP053_Head)
}

#HPAP-056_Head
{
  HPAP056_Head <- CODEX_Final %>% filter(Image=="HPAP-056_CODEX_Head-of-pancreas_OCT.ome.tif")
  HPAP056_Head <- CreateSeuratObject(counts = t(HPAP056_Head[,-c(1:6)]), assay = "CODEX",meta.data = HPAP056_Head[,c(1:6)])
  HPAP056_Head <- NormalizeData(HPAP056_Head)
  HPAP056_Head <- FindVariableFeatures(HPAP056_Head)
  HPAP056_Head <- ScaleData(HPAP056_Head)
  HPAP056_Head <- RunPCA(HPAP056_Head)
}

#HPAP-066_Tail
{
  HPAP066_Tail <- CODEX_Final %>% filter(Image=="HPAP-066_CODEX_Tail-of-pancreas_OCT.ome.tif")
  HPAP066_Tail <- CreateSeuratObject(counts = t(HPAP066_Tail[,-c(1:6)]), assay = "CODEX",meta.data = HPAP066_Tail[,c(1:6)])
  HPAP066_Tail <- NormalizeData(HPAP066_Tail)
  HPAP066_Tail <- FindVariableFeatures(HPAP066_Tail)
  HPAP066_Tail <- ScaleData(HPAP066_Tail)
  HPAP066_Tail <- RunPCA(HPAP066_Tail)
}

#HPAP-093_Body
{
  HPAP093_Body <- CODEX_Final %>% filter(Image=="HPAP-093_CODEX_Body-of-pancreas_OCT.ome.tif")
  HPAP093_Body <- CreateSeuratObject(counts = t(HPAP093_Body[,-c(1:6)]), assay = "CODEX",meta.data = HPAP093_Body[,c(1:6)])
  HPAP093_Body <- NormalizeData(HPAP093_Body)
  HPAP093_Body <- FindVariableFeatures(HPAP093_Body)
  HPAP093_Body <- ScaleData(HPAP093_Body)
  HPAP093_Body <- RunPCA(HPAP093_Body)
}

#HPAP-094_Body
{
  HPAP094_Body <- CODEX_Final %>% filter(Image=="HPAP-094_CODEX_Body-of-pancreas_OCT.ome.tif")
  HPAP094_Body <- CreateSeuratObject(counts = t(HPAP094_Body[,-c(1:6)]), assay = "CODEX",meta.data = HPAP094_Body[,c(1:6)])
  HPAP094_Body <- NormalizeData(HPAP094_Body)
  HPAP094_Body <- FindVariableFeatures(HPAP094_Body)
  HPAP094_Body <- ScaleData(HPAP094_Body)
  HPAP094_Body <- RunPCA(HPAP094_Body)
}

#HPAP-095_Head
{
  HPAP095_Head <- CODEX_Final %>% filter(Image=="HPAP-095_CODEX_Head-of-pancreas_OCT.ome.tif")
  HPAP095_Head <- CreateSeuratObject(counts = t(HPAP095_Head[,-c(1:6)]), assay = "CODEX",meta.data = HPAP095_Head[,c(1:6)])
  HPAP095_Head <- NormalizeData(HPAP095_Head)
  HPAP095_Head <- FindVariableFeatures(HPAP095_Head)
  HPAP095_Head <- ScaleData(HPAP095_Head)
  HPAP095_Head <- RunPCA(HPAP095_Head)
}

#HPAP-112_Head
{
  HPAP112_Head <- CODEX_Final %>% filter(Image=="HPAP-112_CODEX_Head-of-pancreas_OCT.ome.tif")
  HPAP112_Head <- CreateSeuratObject(counts = t(HPAP112_Head[,-c(1:6)]), assay = "CODEX",meta.data = HPAP112_Head[,c(1:6)])
  HPAP112_Head <- NormalizeData(HPAP112_Head)
  HPAP112_Head <- FindVariableFeatures(HPAP112_Head)
  HPAP112_Head <- ScaleData(HPAP112_Head)
  HPAP112_Head <- RunPCA(HPAP112_Head)
}

#HPAP-155_Body
{
  HPAP155_Body <- CODEX_Final %>% filter(Image=="HPAP-155_CODEX_Body-of-pancreas_OCT.ome.tif")
  HPAP155_Body <- CreateSeuratObject(counts = t(HPAP155_Body[,-c(1:6)]), assay = "CODEX",meta.data = HPAP155_Body[,c(1:6)])
  HPAP155_Body <- NormalizeData(HPAP155_Body)
  HPAP155_Body <- FindVariableFeatures(HPAP155_Body)
  HPAP155_Body <- ScaleData(HPAP155_Body)
  HPAP155_Body <- RunPCA(HPAP155_Body)
}

#HPAP-169_Tail
{
  HPAP169_Tail <- CODEX_Final %>% filter(Image=="HPAP-169_CODEX_Tail-of-pancreas_OCT.ome.tif")
  HPAP169_Tail <- CreateSeuratObject(counts = t(HPAP169_Tail[,-c(1:6)]), assay = "CODEX",meta.data = HPAP169_Tail[,c(1:6)])
  HPAP169_Tail <- NormalizeData(HPAP169_Tail)
  HPAP169_Tail <- FindVariableFeatures(HPAP169_Tail)
  HPAP169_Tail <- ScaleData(HPAP169_Tail)
  HPAP169_Tail <- RunPCA(HPAP169_Tail)
}

#HPAP-174_Body
{
  HPAP174_Body <- CODEX_Final %>% filter(Image=="HPAP-174_CODEX_Body-of-pancreas_OCT.ome.tif")
  HPAP174_Body <- CreateSeuratObject(counts = t(HPAP174_Body[,-c(1:6)]), assay = "CODEX",meta.data = HPAP174_Body[,c(1:6)])
  HPAP174_Body <- NormalizeData(HPAP174_Body)
  HPAP174_Body <- FindVariableFeatures(HPAP174_Body)
  HPAP174_Body <- ScaleData(HPAP174_Body)
  HPAP174_Body <- RunPCA(HPAP174_Body)
}

#HPAP-175_Tail
{
  HPAP175_Tail <- CODEX_Final %>% filter(Image=="HPAP-175_CODEX_Tail-of-pancreas_OCT.ome.tif")
  HPAP175_Tail <- CreateSeuratObject(counts = t(HPAP175_Tail[,-c(1:6)]), assay = "CODEX",meta.data = HPAP175_Tail[,c(1:6)])
  HPAP175_Tail <- NormalizeData(HPAP175_Tail)
  HPAP175_Tail <- FindVariableFeatures(HPAP175_Tail)
  HPAP175_Tail <- ScaleData(HPAP175_Tail)
  HPAP175_Tail <- RunPCA(HPAP175_Tail)
}

# Find integration anchors and cluster----
anchors <- FindIntegrationAnchors(object.list = list(HPAP013_Body, 
                                                     HPAP016_Body,
                                                     HPAP017_Body,
                                                     HPAP027_Body,
                                                     HPAP040_Head,
                                                     HPAP045_Tail,
                                                     HPAP053_Head,
                                                     HPAP056_Head,
                                                     HPAP066_Tail,
                                                     HPAP093_Body,
                                                     HPAP094_Body,
                                                     HPAP095_Head,
                                                     HPAP112_Head,
                                                     HPAP155_Body,
                                                     HPAP169_Tail,
                                                     HPAP174_Body,
                                                     HPAP175_Tail), anchor.features = 4000, reduction = "rpca")
#Integrate data
Integrated_data <- IntegrateData(anchors = anchors, dims = 1:30)

#Seurat Pipeline - Normalize and Standardize
CODEX_seurat <- NormalizeData(Integrated_data, normalization.method = "CLR", margin = 2, verbose = FALSE)
CODEX_Seurat <- ScaleData(Integrated_data, features = rownames(Integrated_data), verbose = FALSE)

#PCA and Clustering
CODEX_Seurat <- RunPCA(CODEX_Seurat, features = rownames(CODEX_Seurat), npcs = 30, verbose = FALSE)
ElbowPlot(CODEX_Seurat, ndims = 30) #Decide how many dimensions to use
CODEX_Seurat <- FindNeighbors(CODEX_Seurat, dims = 1:20, verbose = TRUE) #15 minutes
CODEX_Seurat <- FindClusters(CODEX_Seurat, resolution = 0.6) #45 minutes
CODEX_Seurat <- RunUMAP(CODEX_Seurat, dims = 1:20) #40 minutes

#Check Dimplot for Clustering at 0.6 resolution
DimPlot(CODEX_Seurat,split.by = "Group", label = TRUE) + ggtitle("Resolution 0.6")

#Top_marker analysis
#Change Idents
Idents(CODEX_Seurat) <- "integrated_snn_res.0.6"
table(Idents(CODEX_Seurat))

CODEX_Markers <- FindAllMarkers(CODEX_Seurat,
                              only.pos = TRUE,
                              min.pct = 0.25,
                              logfc.threshold = 0.25)

#Export Top 10 markers per cluster
Markers_Seurat_top10 <- CODEX_Markers %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = 10, with_ties = FALSE) %>%
  arrange(cluster, desc(avg_log2FC))

#Level 1 marker identification (Broad markers for endothelial, immune, endocrine, epithelial/ductal, acinar, MSC, and SMA markers)
new.cluster.ids <- c("Endothelial",        #Cluster 0
                     "Endocrine",          #Cluster 1
                     "Endocrine",          #Cluster 10
                     "Endocrine",          #Cluster 11
                     "Immune",             #Cluster 12
                     "Immune",             #Cluster 13 
                     "Epithelial/Ductal",  #Cluster 14 
                     "MSC",                #Cluster 15 
                     "Immune",             #Cluster 16 
                     "Endocrine",          #Cluster 17 
                     "Endocrine",          #Cluster 18
                     "Endocrine",          #Cluster 19
                     "Acinar",             #Cluster 2
                     "Acinar",             #Cluster 20
                     "MSC",                #Cluster 21
                     "Immune",             #Cluster 22
                     "Endocrine",          #Cluster 23
                     "Endothelial",        #Cluster 24
                     "Immune",             #Cluster 3
                     "Endothelial",        #Cluster 4
                     "Acinar",             #Cluster 5
                     "Immune",             #Cluster 6
                     "SMA",                #Cluster 7
                     "SMA",                #Cluster 8
                     "Immune")             #Cluster 9

Idents(CODEX_Seurat) <- "integrated_snn_res.0.6"
table(Idents(CODEX_Seurat))
names(new.cluster.ids) <- levels(CODEX_Seurat)
CODEX_Seurat <- RenameIdents(CODEX_Seurat, new.cluster.ids)
DimPlot(CODEX_Seurat, 
        split.by = "Group",
        reduction = "umap", label = TRUE, repel = F)

#Supplementary Figure 6C: Overall UMAP for all cells
DimPlot(CODEX_Seurat, 
        reduction = "umap", 
        label = F, 
        repel = TRUE, 
        label.size = 4, 
        pt.size = 0.8, 
        cols = c(
          "#d62728","#2ca02c","#e377c2","#ff7f0e",
          "#1f77b4","#bcbd22","#17becf")
) +
  theme_void() +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    plot.title = element_blank(),
    plot.subtitle = element_text(size = 14, hjust = 0.5),
    axis.text = element_text(size=12),
    axis.ticks = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    legend.position ="none",
    plot.margin = margin(10, 10, 10, 10)
  ) +
  labs(
    title = "Pancreas (UMAP)")


##########################################################################################################  
######################################### Immune Seurat Markers #######################################
##########################################################################################################  
Immune_seurat <- subset(CODEX_Seurat, idents = "Immune")
#Normalization not necessary because preceding Seurat was already normalization prior to scaling
Immune_seurat <- ScaleData(Immune_seurat, features = rownames(Immune_seurat), verbose = FALSE)
Immune_seurat <- RunPCA(Immune_seurat, features = rownames(Immune_seurat), npcs = 30, verbose = FALSE)
ElbowPlot(Immune_seurat, ndims = 30) #Decide how many dimensions to use
Immune_seurat <- FindNeighbors(Immune_seurat, dims = 1:12, verbose = TRUE)
Immune_seurat <- FindClusters(Immune_seurat, resolution = c(0.6))
Immune_seurat <- RunUMAP(Immune_seurat, dims = 1:12)



#Find top markers at resolution =0.6
Idents(Immune_seurat) <- "integrated_snn_res.0.6"
table(Idents(Immune_seurat))
markers_Immune_seurat <- FindAllMarkers(Immune_seurat,only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25)

#Filter our the top 10 markers percluster at resolution=0.6
markers_Immune_seurat_top10 <- markers_Immune_seurat %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = 10, with_ties = FALSE) %>% arrange(cluster, desc(avg_log2FC))


new.cluster.ids <- c("cDC1",                   #Cluster 0
                     "CD8 T",                  #Cluster 1
                     "DC",                     #Cluster 2
                     "M1 Macrophage",          #Cluster 3
                     "CD8 T",                  #Cluster 4
                     "CD4 T",                  #Cluster 5
                     "M2 Macrophage",          #Cluster 6
                     "Granulocyte/Monocyte",   #Cluster 7
                     "B Cell",                 #Cluster 8
                     "Other",                  #Cluster 9
                     "Other",                  #Cluster 10
                     "CD8 T",                  #Cluster 11
                     "Granulocyte/Monocyte")   #Cluster 12

table(Idents(Immune_seurat))

#Rename Immune Seurat Clusters
names(new.cluster.ids) <- levels(Immune_seurat)
Immune_seurat <- RenameIdents(Immune_seurat, new.cluster.ids)
DimPlot(Immune_seurat, reduction = "umap",
        split.by = "Group",
        label = TRUE, repel = F)

#Filter out "Others"
Immune_seurat <- subset(Immune_seurat, idents = c("cDC1",
                                                  "CD8 T",
                                                  "DC",
                                                  "M1 Macrophage",
                                                  "CD4 T",
                                                  "M2 Macrophage",
                                                  "Granulocyte/Monocyte",
                                                  "B Cell"
                                                  ))
#Figure 4C
DimPlot(Immune_seurat, 
          reduction = "umap", 
          label = F, 
          repel = TRUE, 
          label.size = 4, 
          pt.size = 0.8, 
          cols = c(
            "#7f7f7f","#e377c2","#1f77b4","#ff7f0e","cyan","#d62728","#2ca02c","#bcbd22")
  ) +
    theme_void() +
    theme(
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      plot.title = element_blank(),
      plot.subtitle = element_text(size = 14, hjust = 0.5),
      axis.text = element_text(size=12),
      axis.ticks = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
      legend.position ="none",
      plot.margin = margin(10, 10, 10, 10)) +
    labs(
      title = "Immune Cells (UMAP)")

##########################################################################################################  
######################################### Endocrine Seurat Markers #######################################
##########################################################################################################  

Endocrine_seurat <- subset(Seurat, idents = "Endocrine")

#Normalization not necessary because preceding Seurat was already normalization prior to scaling
Endocrine_seurat <- ScaleData(Endocrine_seurat, features = rownames(Endocrine_seurat), verbose = FALSE)
Endocrine_seurat <- RunPCA(Endocrine_seurat, features = rownames(Endocrine_seurat), npcs = 30, verbose = FALSE)
ElbowPlot(Endocrine_seurat, ndims = 30) #Decide how many dimensions to use
Endocrine_seurat <- FindNeighbors(Endocrine_seurat, dims = 1:7, verbose = TRUE)
Endocrine_seurat <- FindClusters(Endocrine_seurat, resolution = c(0.2))
Endocrine_seurat <- RunUMAP(Endocrine_seurat, dims = 1:7)

DimPlot(Endocrine_seurat, split.by ="Group",label = TRUE) + ggtitle("Endocrine_res0.2")

#Find top markers at resolution =0.2 (resolution needs to be low to prevent over-clustering)
Idents(Endocrine_seurat) <- "CODEX_snn_res.0.2"
table(Idents(Endocrine_seurat))
markers_Endocrine_seurat <- FindAllMarkers(Endocrine_seurat,only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25)

#Filter our the top 10 markers
markers_Endocrine_seurat_top10 <- markers_Endocrine_seurat %>% 
  group_by(cluster) %>% 
  slice_max(order_by = avg_log2FC, n = 10, with_ties = FALSE) %>% 
  arrange(cluster, desc(avg_log2FC))

#Assign new cluster IDs
new.cluster.ids <- c("Alpha",        #Cluster 0
                     "Stromal",      #Cluster 1
                     "Delta",        #Cluster 2
                     "Beta",         #Cluster 3
                     "Beta",         #Cluster 4 
                     "Beta",         #Cluster 5 
                     "Mix",          #Cluster 6 
                     "Other",        #Cluster 7 
                     "PP",           #Cluster 8 
                     "Beta",         #Cluster 9 
                     "Other",        #Cluster 10 
                     "Ghrelin",      #Cluster 11 
                     "Other")        #Cluster 12


names(new.cluster.ids) <- levels(Endocrine_seurat)
Endocrine_seurat <- RenameIdents(Endocrine_seurat, new.cluster.ids)
DimPlot(Endocrine_seurat, reduction = "umap", label = TRUE, repel = F)

#Filter OUT all "Other" markers
Endocrine_seurat <- subset(Endocrine_seurat, idents = c("Alpha",
                                                        "Beta",
                                                        "Delta",
                                                        "Ghrelin",
                                                        "Mix",
                                                        "PP",
                                                        "Stromal"))
#Figure 4C
DimPlot(Endocrine_seurat, 
        reduction = "umap", 
        label = F, 
        repel = TRUE, 
        label.size = 4, 
        pt.size = 0.8, 
        cols = c(
          "#d62728","#1f77b4","#bcbd22","#2ca02c","#7f7f7f","#e377c2","cyan")
) +
  theme_void() +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    plot.title = element_blank(),
    plot.subtitle = element_text(size = 14, hjust = 0.5),
    axis.text = element_text(size=12),
    axis.ticks = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    legend.position ="none",
    plot.margin = margin(10, 10, 10, 10)
  ) +
  labs(
    title = "Endocrine Cells (UMAP)")



#Supplementary Figure S7A: Centroid Data for Endocrine Cells
#Convert Seurat Object to data.frame: Contains coordinate metadata and cluster names
df <- data.frame(Endocrine_seurat@meta.data,
                 Clusternames = Idents(Endocrine_seurat))

df <- df %>% filter(Clusternames %in% c("Alpha", "Beta","Delta", "PP", "Ghrelin"))
#Choose Young Adults
df <- df %>% filter(DonorID=="HPAP-013")
#Choose Old Adults
df <- df %>% filter(DonorID=="HPAP-174")


#Part1 - Whole Image
ggplot(df, aes(x = CentroidX, y = Centroidy, color = Clusternames)) +
  geom_point(size = 0.6, alpha = 0.85) +
  scale_color_manual(values=c( "#d62728","#bcbd22","#2ca02c","#e377c2","cyan")) +
  #scale_y_reverse() +                         # Flip if imaging coordinates are top-down
  coord_fixed() +
  scale_x_continuous(limits=c(1273,9299)) + #HPAP-013 Coordinates
  scale_y_continuous(limits = c(533,7798))+ #HPAP-013 Coordinates# Preserve spatial aspect ratio
  labs(
    title = NULL,
    x = "Spatial x (µm)",
    y = "Spatial y (µm)",
    color = "Cluster"
  ) +
  theme_classic(base_size = 16) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_text(size=10),
    legend.position = "none",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    #plot.margin = margin(10, 10, 10, 10),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 18)
  )



#Part 2: Inset
ggplot(df, aes(x = CentroidX, y = Centroidy, color = Clusternames)) +
  geom_point(size = 2.5, alpha = 0.85) +
  scale_color_manual(values=c( "#d62728","#bcbd22","#2ca02c","#e377c2","cyan")) +
  scale_y_reverse() +                         # Flip if imaging coordinates are top-down
  coord_fixed() +                             # Preserve spatial aspect ratio
  labs(
    title = NULL,
    x = NULL,
    y = NULL,
    color = "Cluster"
  ) +
  #scale_x_continuous(limits=c(5200,5400)) + #HPAP-013 Coordinates
  #scale_y_continuous(limits = c(3250,3450))+ #HPAP-013 Coordinates
  scale_x_continuous(limits=c(5400,5700)) + #HPAP-169 Coordinates
  scale_y_continuous(limits = c(2200,2500))+ #HPAP-169 Coordinates
  theme_classic(base_size = 16) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.position = "none",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    #plot.margin = margin(10, 10, 10, 10),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 18))



#Figure 4D: Centroid Data for Immune Cells
df <- rbind.data.frame(data.frame(Immune_seurat@meta.data,
                                  Clusternames = Idents(Immune_seurat)),
                       data.frame(Endocrine_seurat@meta.data,
                                  Clusternames = Idents(Endocrine_seurat)))
  
df <- df %>% filter(Clusternames %in% c("CD8 T", "CD4 T",
                                        "M2 Macrophage", "M1 Macrophage",
                                        "Granulocyte/Monocyte",
                                        "B Cell",
                                        "Alpha",
                                        "Beta",
                                        "Delta"))
#Choose Young Adults
df <- df %>% filter(DonorID=="HPAP-013")
#Choose Old Adults
df <- df %>% filter(DonorID=="HPAP-174")
  
#Part1 - Whole Image
ggplot(df, aes(x = CentroidX, y = Centroidy, color = Clusternames)) +
  geom_point(size = 0.3, alpha = 0.85) +
  scale_color_manual(values=c("#e377c2","#d62728","cyan","orange","#2ca02c","#bcbd22",
                              "gray50","gray50","gray50")) +
  coord_fixed() + 
  scale_x_continuous(limits=c(1273,9299)) + #MAX Coordinates
  scale_y_continuous(limits = c(533,7798))+ #MAX COORDINATES
  labs(
    title = NULL,
    x = "Spatial x (µm)",
    y = "Spatial y (µm)",
    color = "Cluster") +
  theme_classic(base_size = 16) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_text(size=10),
    legend.position = "none",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 18))

#Part2 - Create Inset  
  
ggplot(df, aes(x = CentroidX, y = Centroidy, color = Clusternames, alpha = Clusternames)) +
  geom_point(size = 2.5, alpha = 0.85) +
  scale_color_manual(values=c("#e377c2","#d62728","cyan","orange","#2ca02c","#bcbd22",
                              "gray90","gray90","gray90")) +
  scale_alpha_manual(values=c(1,1,1,1,1,1,0.01,0.01,0.01)) +
  scale_y_reverse() +                         # Flip if imaging coordinates are top-down
  coord_fixed() +                             # Preserve spatial aspect ratio
  labs(title = NULL,
       x = NULL,
       y = NULL,
       color = "Cluster") +
  #scale_x_continuous(limits=c(5300,5700)) + #HPAP-169 Coordinates
  #scale_y_continuous(limits = c(2200,2600))+ #HPAP-169 Coordinates
  scale_x_continuous(limits=c(3550,3900)) + #HPAP-174 Coordinates
  scale_y_continuous(limits = c(1150,1500))+ #HPAP-174 Coordinates
  theme_classic(base_size = 16) +
  theme(panel.background = element_rect(fill = "white", color = NA),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "none",
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        #plot.margin = margin(10, 10, 10, 10),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 18))



#Supplementary Figure S7B
#Endocrine_Group Statistics
Codex_population_Endocrine <- data.frame(table(Idents(Endocrine_seurat),
                                               Endocrine_seurat@meta.data$DonorID))
Codex_population_Endocrine$Group <- ifelse(Codex_population_Endocrine$Var2 %in% c("HPAP-053",
                                                                                  "HPAP-066",
                                                                                  "HPAP-093",
                                                                                  "HPAP-169",
                                                                                  "HPAP-174",
                                                                                  "HPAP-175"),
                                           "old","young")

names(Codex_population_Endocrine) <- c("Celltype","DonorID","value","Group")
Codex_population_Endocrine <- Codex_population_Endocrine %>% filter(Celltype %in% c("Alpha","Beta", "Delta","PP","Ghrelin"))
Codex_population_Endocrine$Group <- factor(Codex_population_Endocrine$Group, levels = c("young","old"))
Codex_population_Endocrine$Celltype <- factor(Codex_population_Endocrine$Celltype, levels = c("Alpha","Beta","Delta","PP","Ghrelin"))


Endocrine_proportions <- Codex_population_Endocrine %>%
  group_by(DonorID, Group) %>%
  mutate(TotalCells=sum(value))

Endocrine_proportions$Proportion <- (Endocrine_proportions$value/Endocrine_proportions$TotalCells) * 100

#Stacked Barplot (data not shown)
ggplot(Endocrine_proportions, aes(x = Group, y = Proportion, fill = Celltype)) +
  geom_bar(stat = "identity", position = "fill") +
 # facet_wrap(~ Group) +
  theme_classic(base_size = 12) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.25),              # 0, 0.25, 0.5, 0.75, 1
    labels = c(0, 25, 50, 75, 100)) +         # show as percent
  scale_fill_manual(values=c("#d62728","#2ca02c","#bcbd22","#e377c2","cyan")) + 
  labs(y = "Proportion", x = "Subject", fill = "Cell Type") +
  theme(legend.position = 'none',
        axis.title.x =element_blank())

#Bargraphs per cell type
Endocrine_Ind <- Endocrine_proportions %<% filter(Celltype)
Endocrine_summary <- Endocrine_proportions %>%
  group_by(Group, Celltype) %>%
  summarise(
    mean_prop = mean(Proportion, na.rm = TRUE),
    sem = sd(Proportion, na.rm = TRUE) / sqrt(n())
  ) %>%
  ungroup()

ggplot(Endocrine_summary, aes(x = Group, y = mean_prop, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.9,color='black') +
  geom_errorbar(aes(ymin = mean_prop - sem, ymax = mean_prop + sem),
                width = 0.2, position = position_dodge(0.7)) +
  geom_jitter(data=Endocrine_proportions, aes(x=Group,y=Proportion),shape=21,fill='white',width=0.3) +
  
  facet_wrap2(~ Celltype, ncol = 5, scales = "fixed",
              strip = strip_themed(
                background_x = elem_list_rect(
                  fill = c("#d62728","#2ca02c","#bcbd22","#e377c2","cyan"),
                  color = "black"))) +
  scale_fill_manual(values = c("#607C83","#E0DCB8")) +
  theme_classic(base_size = 16) +
  labs(
    y = "% endocrine cells",
    x = NULL,
    fill = "Group"
  ) +
  theme(
    strip.text = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),
    axis.title.y = element_text(face = "bold", size = 14),
    legend.position = "none",
    legend.title = element_blank())

#Statistics
Endocrine_ttest_results <- Endocrine_proportions %>%
  group_by(Celltype) %>%
  summarise(
    ttest = list(t.test(Proportion ~ Group, data = cur_data())),
    p.value = ttest[[1]]$p.value,
    t.statistic = ttest[[1]]$statistic,
    df = ttest[[1]]$parameter,
    mean_group1 = ttest[[1]]$estimate[1],
    mean_group2 = ttest[[1]]$estimate[2]) 


#Group Statistics - Immune

#Immune_Group
Codex_population_Immune$Group <- ifelse(Codex_population_Immune$Var2 %in% c("HPAP-040",
                                                                            "HPAP-017",
                                                                            "HPAP-053",
                                                                            "HPAP-066",
                                                                            "HPAP-155",
                                                                            "HPAP-174",
                                                                            "HPAP-175"),
                                 "old","young")  
  
Codex_population_Immune
names(Codex_population_Immune) <- c("Celltype","DonorID","value","Group")
Codex_population_Immune$Group <- factor(Codex_population_Immune$Group, levels = c("young","old"))
Codex_population_Immune$Celltype <- factor(Codex_population_Immune$Celltype, levels = c("M1 Macrophage","M2 Macrophage","Granulocyte/Monocyte","CD8 T","CD4 T",
                                                                                        "B Cell","DC","cDC1"))
  
Codex_population_Immune <- Codex_population_Immune %>%
  group_by(DonorID, Group) %>%
  mutate(TotalCells=sum(value))
  
Codex_population_Immune$Proportion <- (Codex_population_Immune$value/Codex_population_Immune$TotalCells) * 100
  
#Figure 4G - Stacked Barplot
ggplot(Codex_population_Immune, aes(x = Group, y = Proportion, fill = Celltype)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_classic(base_size = 12) +
  scale_fill_manual(values = c("#d62728","orange","#2ca02c","#e377c2","cyan","#bcbd22","#1f77b4","#7f7f7f")) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.25),              # 0, 0.25, 0.5, 0.75, 1
    labels = c(0, 25, 50, 75, 100)         # show as percent
  ) +    labs(y = "Percent", x = "Subject", fill = "Cell Type") +
  theme(legend.position = 'none',
        axis.title.x = element_blank())

#Bargraph per celltype
Immune_summary <- Codex_population_Immune %>%
  group_by(Group, Celltype) %>%
  summarise(mean_prop = mean(Proportion, na.rm = TRUE),
            sem = sd(Proportion, na.rm = TRUE) / sqrt(n())) %>%
  ungroup()
  
ggplot(Immune_summary, aes(x = Group, y = mean_prop, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.9,color='black') +
  geom_errorbar(aes(ymin = mean_prop - sem, ymax = mean_prop + sem),
                width = 0.2, position = position_dodge(0.7)) +
  geom_jitter(data=Codex_population_Immune, aes(x=Group,y=Proportion),shape=21,fill='white',width=0.3) +
  facet_wrap2(~ Celltype, ncol = 8, scales = "fixed",
              strip = strip_themed(
              background_x = elem_list_rect(
                fill = c("orange","#d62728","#2ca02c","#e377c2","cyan","#bcbd22","#1f77b4","#7f7f7f"),
                color = "black"))) +
  scale_fill_manual(values = c("#E0DCB8","#607C83")) +
  theme_classic(base_size = 14) +
  labs(y = "Proportion of Total Cells",
       x = NULL,
       fill = "Group") +
    theme(
      strip.text = element_text(face = "bold", size = 14),
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),
      axis.title.y = element_text(face = "bold", size = 14),
      legend.position = "none",
      legend.title = element_blank())

Immune_ttest_results <- Codex_population_Immune %>%
  group_by(Celltype) %>%
  summarise(
    ttest = list(t.test(Proportion ~ Group, data = cur_data())),
    p.value = ttest[[1]]$p.value,
    t.statistic = ttest[[1]]$statistic,
    df = ttest[[1]]$parameter,
    mean_group1 = ttest[[1]]$estimate[1],
    mean_group2 = ttest[[1]]$estimate[2]) 

#CD8 and CD4 T cell proportion
CD8_proportion <- Codex_population_Immune %>% 
  filter(Celltype %in% c("CD8 T", "CD4 T")) %>%
  group_by(DonorID, Group, Celltype) %>%
  summarise(Total = sum(value), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Celltype, values_from = Total) %>%
  mutate(CD3 = `CD8 T` + `CD4 T`,
         CD8_CD3 = (`CD8 T` / CD3) * 100,
         CD4_CD3 = (`CD4 T` / CD3) * 100)

CD8_proportion$Group <- ifelse(CD8_proportion$DonorID %in% c("HPAP-040",
                                                              "HPAP-017",
                                                              "HPAP-053",
                                                              "HPAP-066",
                                                              "HPAP-155",
                                                              "HPAP-174",
                                                              "HPAP-175"),
                               "old","young")


  
  CD8_summary <- CD8_proportion %>%
    group_by(Group) %>%
    summarise(
      mean_prop = mean(CD8_CD3, na.rm = TRUE),
      sem = sd(CD8_CD3, na.rm = TRUE) / sqrt(n())
    ) %>%
    ungroup()
  
  CD8_summary$Group <- factor(CD8_summary$Group, levels = c("young","old"))
  

#Figure 4G - % CD8 of CD4 
#Switch for CD4_CD3 when illustrating CD4 graph
ggplot(CD8_summary, aes(x = Group, y = mean_prop, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.9,color='black') +
  geom_errorbar(aes(ymin = mean_prop - sem, ymax = mean_prop + sem),
                width = 0.2, position = position_dodge(0.7)) +
  geom_jitter(data=CD8_proportion, aes(x=Group,y=CD8_CD3),shape=21,fill='white',width=0.1) +
  scale_fill_manual(values = c("#E0DCB8","#607C83")) +
  coord_cartesian(ylim = c(50,100)) +
  #scale_y_continuous(limits=c(0,50)) +
  theme_classic(base_size = 12) +
  labs(title="CD8+ proportion",
       y = bquote("%CD8"^"+" ~ "of CD3"^"+"),
       x = NULL,
       fill = "Group") +
  theme(plot.title = element_blank(),
        strip.text = element_text(face = "bold", size = 14),
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),
      axis.title.y = element_text(face = "bold", size = 14),
      legend.position = "none",
      legend.title = element_blank())

t.test(CD8_CD3 ~ Group, data = CD8_proportion) #%CD4 of CD3
t.test(CD4_CD3 ~ Group, data = CD8_proportion) #%CD8 of CD3

  

##### Supplementary Figure 6E-G ##### 
#Nebulosa plots for whole cell expression

#Figure 6E - All cells
#Cells used = Endocrine: CHGA..Cell..Mean, CPEP..Cell..Mean, SST..Cell..Mean,
#             Immune: CD45..Cell..Mean,
#             Lymphocyte: CD3..Cell..Mean,
#             Macrophage: CD68..Cell..Mean, CD163..Cell..Mean
#             Acinar: PNLIP..Cell..Mean, GP2..Cell..Mean, CD31..Cell..Mean
#             Smoothe Muscle Actin (SMA): ACTA2..Cell..Mean

Nebulosa::plot_density(CODEX_Seurat, features = "CHGA..Cell..Mean") +
  scale_color_viridis_c(option = "inferno") +  # Apply magma palette to color
  labs(title="Endpcrine",
       subtitle="Chromagranin-A") +
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5,size=14,face='bold'),
        plot.subtitle = element_text(hjust=0.5),
        legend.position='none')


# Supplementary Figure 6F: Endocrine Markers 
#Marker Identification = Endocrine: "CHGA..Cell..Mean",
#                        Beta cell: "CPEP..Cell..Mean", "NKX6.1..Nucleus..Mean", "PDX1..Nucleus..Mean"
#                        Alpha cell: "GCG..Cell..Mean",
#                        Delta cell: "SST..Cell..Mean",
#                        Epsilon cell: "GHRL..Cell..Mean"
#                        Pancreatic polypeptide cell: "PPY..Cell..Mean"


Nebulosa::plot_density(Endocrine_seurat, features = "CHGA..Cell..Mean") +
  scale_color_viridis_c(option = "inferno") +  # Apply magma palette to color
  labs(title="Endocrine",
       subtitle="Chromagranin-A") +
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5,size=14,face='bold'),
        plot.subtitle = element_text(hjust=0.5),
        legend.position='none')

# Supplementary Figure 6G: Immune Markers
#Marker Identification = Macrophage: "CD68..Cell..Mean", "CD163..Cell..Mean",
#                        Granulocyte/Monocyte: "CD66b..Cell..Mean",
#                        T cell: "CD3..Cell..Mean", "CD8..Cell..Mean", "CD4..Cell..Mean",
#                        Dendritic cell: "CD11c..Cell..Mean",
#                        B cell: ""CD16..Cell..Mean"

Nebulosa::plot_density(Immune_seurat, features = "GHRL..Cell..Mean") +
  scale_color_viridis_c(option = "inferno") +  # Apply magma palette to color
  labs(title="Ghrelin",
       subtitle="") +
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5,size=14,face='bold'),
        plot.subtitle = element_text(hjust=0.5),
        legend.position='none')





########### Supplementary Figure 6D, and Supplementary Figure 7C ###########
#Export all cells and present as scaled abundance
#Merge Acinar, Endothelial, Alpha, Beta, Delta, PP, Ghrelin, cDC1, CD8 T, DC, M1Mac, CD4 T, M2 Mac, Granulo/Monocyte, B Cell

#Setup CODEX DATA: Original CSV found in Zenodo (Zenodo (https://doi.org/10.5281/zenodo.20496994))
CODEX <- read.table(file="Figure 4_CODEX_Cell_measurements_merged",header = T,sep=",",check.names = FALSE) 

#Setup Raw Markers (repeat of lines 168-206 above)
raw_markers <- c("GCG..Cell..Mean",                #Alpha
                 "CD45..Cell..Mean",               #Immune
                 "CD19..Cell..Mean",               #B cells
                 "CD39L3..Cell..Mean",             #ENTPD
                 "GP2..Cell..Mean",                #Pancreatic Progenitor
                 "NKX2.2..Nucleus..Mean",          #Pancreas
                 "CHGA..Cell..Mean",               #Endocrine
                 "CD141..Cell..Mean",              #cD1Cs
                 "COL4A1..Cell..Mean",             #ECM
                 "COL6..Cell..Mean",               #ECM
                 "CD3..Cell..Mean",                #T cells
                 "CD8..Cell..Mean",                #CD8 T cells
                 "SST..Cell..Mean",                #Delta cells
                 "CD4..Cell..Mean",                #CD4T cells
                 "ACTA2..Cell..Mean",              #Smooth Muscle
                 "COL1A1..Cell..Mean",             #Connective Tissue 
                 "PAX6..Nucleus..Mean",            #Islet progenitors
                 "PPY..Cell..Mean",                #PP cells
                 "SOX9..Nucleus..Mean",            #Pancreas Development
                 "CD11c..Cell..Mean",              #Dendritic Cells
                 "NKX6.1..Nucleus..Mean",          #Beta Cells
                 "KRT..Cell..Mean",                #Epithelial
                 "EPCAM..Cell..Mean",              #Epithelial
                 "CD31..Cell..Mean",               #Endothelial
                 "PDX1..Nucleus..Mean",            #Pancreas
                 "GHRL..Cell..Mean",               #Ghrelin
                 "CD163..Cell..Mean",              #Macrophage
                 "HLA.DR..Cell..Mean",             #MHC-II
                 "CD68..Cell..Mean",               #Macrophage
                 "CD66b..Cell..Mean",              #Granulocytes
                 "PNLIP..Cell..Mean",              #Acinar cells
                 "CPEP..Cell..Mean",               #Beta Cells
                 "TUBB3..Cell..Mean",              #Neuron
                 "Ki67..Nucleus..Mean",            #Proliferation Status
                 "LYVE1..Cell..Mean")              #Lymphatic

CODEX <- CODEX[-1] %>% select(matches(raw_markers))
CODEX <- CODEX %>%
  tibble::rownames_to_column("CellID")  # ensure each row has a unique ID

#Rename Markers
names(CODEX)[2:36] <- c("GCG","CD45","CD19","CD39L3","GP2","NKX2.2","CHGA","CD141","COL4A1","COL6","CD3","CD8","SST","CD4","ACTA2","COL1A1","PAX6",
                        "PPY","SOX9","CD11c","NKX6.1","KRT","EPCAM","CD31","PDX1","GHRL","CD163","HLA.DR","CD68",
                        "CD66b","PNLIP","CPEP","TUBB3","Ki67","LYVE1")

#Bind OVERALL CODEX, Immune, and Endocrine Codex files
# NOTE.... General "Immune markers" from the Overall CODEX file is removed to limit overlap
scaled_df <- rbind.data.frame((data.frame(CODEX_Seurat@meta.data,
                             Clusternames=Idents(CODEX_Seurat))),
                 (data.frame(Immune_seurat@meta.data,
                             Clusternames=Idents(Immune_seurat))),
                 (data.frame(Endocrine_seurat@meta.data,
                             Clusternames=Idents(Endocrine_seurat))))

scaled_df <- scaled_df %>% filter(Clusternames %in% c("Acinar",
                                                      "Endothelial",
                                                      "Ghrelin",
                                                      "PP",
                                                      "Beta",
                                                      "Delta",
                                                      "Alpha",
                                                      "Granulocyte/Monocyte",
                                                      "M1 Macrophage",
                                                      "M2 Macrophage",
                                                      "DC",
                                                      "CD4 T",
                                                      "CD8 T"))

CODEX_Final <- right_join(CODEX, scaled_df, by="CellID")
CODEX_Final <- CODEX_Final %>% drop_na()


#Supplementary Figure 6D
Markers <- c("Alpha cell","Beta Cell","Cd4","Cd8","Delta Cell","EC","Epithelial Cell", "Ghrelin Cell","Macrophages","Other",
             "PP Cell")
Heatmap <- data.frame(CODEX_Final[,-c(1,37:47)])
Heatmap <- data.frame(t(scale(t(Heatmap))),
                      Celltype=CODEX_Final$Clusternames)
Heatmap <- melt(Heatmap,
                idvar = c("Clusternames"))

Heatmap <- Heatmap %>% 
  group_by(Celltype,variable) %>%
  summarise(expression = mean(value))
Heatmap <- Heatmap %>%
  pivot_wider(names_from = Celltype, values_from = expression)
Heatmap <- column_to_rownames(Heatmap, var = "variable")

pheatmap(Heatmap,
         scale = "row",  # Optional: z-score normalization across rows
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("#000004", "#BC3754", "#FDE725FF"))(50))


#Figure 6C - Aggregate endocrine counts across steps
#Endocrine Cells
DefaultAssay(Endocrine_seurmrat)
Endocrine_seurat$Clusternames <- Idents(Endocrine_seurat)
Endocrine_seurat_Split <- AggregateExpression(Endocrine_seurat,
                                           group.by = c("Clusternames","DonorID"),
                                           assays = 'CODEX',
                                           slot = "counts",
                                           return.seurat = FALSE)

Endocrine_seurat_Split <- Endocrine_seurat_Split$CODEX

# transpose
Endocrine_seurat_Split.t <- t(Endocrine_seurat_Split)

# convert to data.frame
Endocrine_seurat_Split.t <- as.data.frame(Endocrine_seurat_Split.t)

# get values where to split
splitRows <- gsub('_.*', '', rownames(Endocrine_seurat_Split.t))

# split data.frame
Endocrine_seurat_Split <- split.data.frame(Endocrine_seurat_Split.t,
                                        f = factor(splitRows))

# fix colnames and transpose

Endocrine_seurat_Split.modified <- lapply(Endocrine_seurat_Split, function(x){
  rownames(x) <- gsub('.*_(.*)', '\\1', rownames(x))
  t(x)
  
})


#Beta_Cells
Beta <- data.frame(Endocrine_seurat_Split$Beta)
Beta <- Beta %>% rownames_to_column("Identifier") %>%              # make rownames a column
  mutate(DonorID = str_sub(Identifier, -8))         # take last 8 characters
Beta$Group <- ifelse(Beta$DonorID %in% c("HPAP-053","HPAP-066","HPAP-093","HPAP-169","HPAP-174","HPAP-175"),
                     "old","young")
Beta$Group <- factor(Beta$Group, levels = c("young","old"))


#ttest
results <- Beta %>%
  select(-c("Identifier", "DonorID", "Group")) %>%
  map_df(~{
    t_res <- t.test(.x[Beta$Group == unique(Beta$Group)[1]],
                    .x[Beta$Group == unique(Beta$Group)[2]])
    tibble(
      p.value = t_res$p.value,
      statistic = t_res$statistic,
      mean_group1 = t_res$estimate[1],
      mean_group2 = t_res$estimate[2]
    )
  }, .id = "Variable")

print(results)

#Alpha_Cells
Alpha <- data.frame(Endocrine_seurat_Split$Alpha)
Alpha <- Alpha %>% rownames_to_column("Identifier") %>%              # make rownames a column
  mutate(DonorID = str_sub(Identifier, -8))         # take last 8 characters
Alpha$Group <- ifelse(Alpha$DonorID %in% c("HPAP-053","HPAP-066","HPAP-093","HPAP-169","HPAP-174","HPAP-175"),
                     "old","young")
#ttest
results <- Alpha %>%
  select(-c("Identifier", "DonorID", "Group")) %>%
  map_df(~{
    t_res <- t.test(.x[Alpha$Group == unique(Alpha$Group)[1]],
                    .x[Alpha$Group == unique(Alpha$Group)[2]])
    tibble(
      p.value = t_res$p.value,
      statistic = t_res$statistic,
      mean_group1 = t_res$estimate[1],
      mean_group2 = t_res$estimate[2]
    )
  }, .id = "Variable")



#Cell Summary
Beta_summary <- Beta[,-c(1,41)] %>%
  group_by(Group) %>%
  summarise(across(
    where(is.numeric),
    list(
      mean = ~mean(.x, na.rm = TRUE),
      sem  = ~sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x)))
    ),
    .names = "{.col}_{.fn}"
  ))

#Change for respective marker to quantify
ggplot(Beta_summary, aes(x = Group, y = PAX6..Nucleus..Mean_mean, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.9,color='black') +
  geom_errorbar(aes(ymin = PAX6..Nucleus..Mean_mean - PAX6..Nucleus..Mean_sem, ymax = PAX6..Nucleus..Mean_mean + PAX6..Nucleus..Mean_sem),
                 width = 0.2, position = position_dodge(0.7)) +
  geom_jitter(data=Beta, aes(x=Group,y=PAX6..Nucleus..Mean),shape=21,fill='white',width=0.1) +
  scale_fill_manual(values = c("#607C83", "#E0DCB8")) +
  #coord_cartesian(ylim = c(50,100)) +
  #scale_y_continuous(limits=c(0,50)) +
  theme_classic(base_size = 12) +
  labs(title="CODEX",
       y = "PAX6",
       x = NULL,
       fill = "Group"
  ) +
  theme(
    plot.title = element_text(hjust=0.5),
    strip.text = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),
    axis.title.y = element_text(face = "bold", size = 14),
    legend.position = "none",
    legend.title = element_blank())

#END