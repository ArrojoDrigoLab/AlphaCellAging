# =======================================================================================================================
# Alpha cell inflammation during human pancreas aging and type 2 diabetes and its reversal by calorie restriction in mice
# -----------------------------------------------------------------------------------------------------------------------
# Author:      Michael Schleh
# Repository:  https://github.com/ArrojoDrigoLab/AlphaCellAging
# Journal:     Science Advances

# Description:
#   Generates Figure 5A-G and Figure S8B
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



##################### Metadata #####################
#Figure 5A: Source Metadata_T2D and ND Cohort
PancDB_Metadata <- read.csv("T2D_PancDB_Metadata.csv",header=T,sep=",")
Metadata_Summary <- PancDB_Metadata %>%
  group_by(Group) %>%
  summarise(across(3:6,
                   list(
                     mean = ~round(mean(.x, na.rm = TRUE), 2),
                     sd   = ~round(sd(.x, na.rm = TRUE), 2),
                     sem  = ~round(sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x))), 2)
                   )))

#Figure 5A: i) Age
   ggplot() +
     geom_col(data=Metadata_Summary, aes(x=Group,y=age_years_mean,fill=Group),color='black') +
     geom_errorbar(data=Metadata_Summary,mapping = aes(x=Group, ymin = age_years_mean-age_years_sem, ymax = age_years_mean+age_years_sem,),width=0.5) +
     geom_jitter(data=PancDB_Metadata, aes(x=Group,y=age_years),width=0.2,alpha=0.5) +
    scale_fill_manual(values=c("#D0DCB8","#B4A1C1"))+
    labs(title="Age",
         y="years")+
    theme_classic() +
    theme(plot.title = element_text(hjust=0.5),
          legend.position = 'none',
          axis.title.x = element_blank(),
          axis.title.y = element_text(size=12))

#statistics   
t.test(PancDB_Metadata$age_years~PancDB_Metadata$Group)  
   

#Figure 5A: ii) BMI
ggplot() +
  geom_col(data=Metadata_Summary, aes(x=Group,y=bmi_mean,fill=Group),color='black') +
  geom_errorbar(data=Metadata_Summary,mapping = aes(x=Group, ymin = bmi_mean-bmi_sem, ymax = bmi_mean+bmi_sem,),width=0.5) +
  geom_jitter(data=PancDB_Metadata, aes(x=Group,y=bmi),width=0.2,alpha=0.5) +
  scale_fill_manual(values=c("#D0DCB8","#B4A1C1"))+
  labs(title="BMI",
       y=bquote("kg/m"^2))+
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5),
        legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=12))

#statistics   
t.test(PancDB_Metadata$bmi~PancDB_Metadata$Group)  
  

#Figure 5A: iii) C-Peptide
ggplot() +
  geom_col(data=Metadata_Summary, aes(x=Group,y=c_peptide_ng_ml_mean,fill=Group),color='black') +
  geom_errorbar(data=Metadata_Summary,mapping = aes(x=Group, ymin = c_peptide_ng_ml_mean-c_peptide_ng_ml_sem, ymax = c_peptide_ng_ml_mean+c_peptide_ng_ml_sem,),width=0.5) +
  geom_jitter(data=PancDB_Metadata, aes(x=Group,y=c_peptide_ng_ml),width=0.2,alpha=0.5) +
  scale_fill_manual(values=c("#D0DCB8","#B4A1C1"))+
  labs(title="C-Peptide",
       y="ng/mL")+
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5),
        legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=12))

#statistics   
t.test(PancDB_Metadata$c_peptide_ng_ml~PancDB_Metadata$Group)  

#Figure 5A: iv) HbA1c
ggplot() +
  geom_col(data=Metadata_Summary, aes(x=Group,y=hba1c_mean,fill=Group),color='black') +
  geom_errorbar(data=Metadata_Summary,mapping = aes(x=Group, ymin = hba1c_mean-hba1c_sem, ymax = hba1c_mean+hba1c_sem,),width=0.5) +
  geom_jitter(data=PancDB_Metadata, aes(x=Group,y=hba1c),width=0.2,alpha=0.5) +
  scale_fill_manual(values=c("#D0DCB8","#B4A1C1"))+
  labs(title="HbA1c",
       y="%")+
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5),
        legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=12))

#statistics   
t.test(PancDB_Metadata$hba1c~PancDB_Metadata$Group)  
  

#Data Processing - Combine Datasets
#See PancDB website (https://hpap.pmacs.upenn.edu/) and corresponding CODEX images for HPAP-### and download.
#Following image download, see QuPath code for cell segmentation and cell marker quantification.
#Exports are determined by marker-quantified csv exports from individual cell segmentation.

ND_Old_Full_Files <- c(
   "CODEX_Exports/HPAP-053_CODEX_Head-of-pancreas_OCT.ome.tif_measurements.csv",
   "CODEX_Exports/HPAP-066_CODEX_Tail-of-pancreas_OCT.ome.tif_measurements.csv",
   "CODEX_Exports/HPAP-093_CODEX_Body-of-pancreas_OCT.ome.tif_measurements.csv",
   "CODEX_Exports/HPAP-112_CODEX_Head-of-pancreas_OCT.ome.tif_measurements.csv",
   "CODEX_Exports/HPAP-169_CODEX_Tail-of-pancreas_OCT.ome.tif_measurements.csv",
   "CODEX_Exports/HPAP-174_CODEX_Body-of-pancreas_OCT.ome.tif_measurements.csv",
   "CODEX_Exports/HPAP-175_CODEX_Tail-of-pancreas_OCT.ome.tif_measurements.csv")

 ND_Old_Full <- ND_Old_Full_Files %>%
  lapply(function(f) read.csv(f, header = TRUE, sep = ",", check.names = FALSE)) %>%
  bind_rows() %>%
  filter(Classification != "Other") 

 ND_Young_Full_Files <- c(
  "CODEX_Exports/HPAP-016_CODEX_Body-of-pancreas_OCT.ome.tif_measurements.csv",
  "CODEX_Exports/HPAP-017_CODEX_Body-of-pancreas_OCT.ome.tif_measurements.csv",
  "CODEX_Exports/HPAP-018_CODEX_Head-of-pancreas_OCT.ome.tif_measurements.csv",
  "CODEX_Exports/HPAP-027_CODEX_Head-of-pancreas_OCT.ome.tif_measurements.csv",
  "CODEX_Exports/HPAP-037_CODEX_Head-of-pancreas_OCT.ome.tif_measurements.csv",
  "CODEX_Exports/HPAP-040_CODEX_Head-of-pancreas_OCT.ome.tif_measurements.csv",
  "CODEX_Exports/HPAP-045_CODEX_Tail-of-pancreas_OCT.ome.tif_measurements.csv",
  "CODEX_Exports/HPAP-056_CODEX_Head-of-pancreas_OCT.ome.tif_measurements.csv",
  "CODEX_Exports/HPAP-095_CODEX_Head-of-pancreas_OCT.ome.tif_measurements.csv",
  "CODEX_Exports/HPAP-155_CODEX_Body-of-pancreas_OCT.ome.tif_measurements.csv")
  
ND_Young_Full <- ND_Young_Full_Files %>%
  lapply(function(f) read.csv(f, header = TRUE, sep = ",", check.names = FALSE)) %>%
  bind_rows() %>%
  filter(Classification != "Other")    

T2D_Full_Files <- c(
  "CODEX_Exports/HPAP-013_CODEX_Body-of-pancreas_OCT.ome.tif_measurements.csv",
  "CODEX_Exports/HPAP-051_CODEX_Tail-of-Pancreas_FFPE.ome.tif_measurements.csv",
  "CODEX_Exports/HPAP-079_CODEX_Tail-of-Pancreas_FFPE.ome.tif_measurements.csv",
  "CODEX_Exports/HPAP-100_CODEX_Tail-of-Pancreas_FFPE.ome.tif_measurements.csv",
  "CODEX_Exports/HPAP-124_CODEX_Tail-of-Pancreas_FFPE.ome.tif_measurements.csv",
  "CODEX_Exports/HPAP-142_CODEX_Tail-of-Pancreas_FFPE.ome.tif_measurements.csv",
  "CODEX_Exports/HPAP-143_CODEX_Tail-of-Pancreas_FFPE.ome.tif_measurements.csv",
  "CODEX_Exports/HPAP-177_CODEX_Body-of-pancreas_OCT.ome.tif_measurements.csv")

T2D_Full <- T2D_Full_Files %>%
  lapply(function(f) read.csv(f, header = TRUE, sep = ",", check.names = FALSE)) %>%
  bind_rows() %>%
  filter(Classification != "Other")    

T2D_Full_Data <- rbind.data.frame(data.frame(ND_Young_Full,
                                                 Category="ND",
                                                 Group="Young",
                                                 Analysis="Full"),
                                      data.frame(ND_Old_Full,
                                                 Category="ND",
                                                 Group="Old",
                                                 Analysis="Full"),
                                      data.frame(T2D_Full,
                                                 Category="T2D",
                                                 Group="T2D",
                                                 Analysis="Full"))

#########################################################################################################
#Here is the node where CSV files are found for download in Zenodo (https://zenodo.org/records/20496994)#
#########################################################################################################
write.csv(T2D_Full_Data, file='Figure 5_CODEX_Cell_measurements_merged.csv')


#Metadata and pre-processing.  Metadata linked to individual cell barcodes
T2D_Metadata <- T2D_Full_Data %>%
  rownames_to_column("CellID") %>%     # ensure rownames are kept as CellID
  mutate(DonorID = substr(Image, 1, 8)) %>%     # extract donor ID from image name
  transmute(      # select/rename relevant columns cleanly
    CellID,
    Image,
    CentroidX = Centroid.X.µm,
    CentroidY = Centroid.Y.µm,
    DonorID,
    Category=Category,
    Group=Group,
    Analysis=Analysis)
      
#Full_Markers_ND  
Full_markers_ND <- c("GCG..Cell..Mean",         #Alpha
                 "CD45..Cell..Mean",            #Immune
                 "CD19..Cell..Mean",            #B cells
                 "CHGA..Cell..Mean",            #Endocrine
                 "PD.L1..Cell..Mean",           #Exhaustion
                 "ATP1A1..Cell..Mean",          #Na+/K= ATPase
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
                 "SOX9..Nucleus..Mean",         #Pancreas Development
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


#Full_Markers_T2D
Full_markers_T2D <- c("GCG..Cell..Mean",        #Alpha
                       "CD45..Cell..Mean",      #Immune
                       "CD19..Cell..Mean",      #B cells
                       "CD20..Cell..Mean",      #B cells
                       "CHGA..Cell..Mean",      #Endocrine
                       "PD.L1..Cell..Mean",     #Exhaustion
                       "CD141..Cell..Mean",     #cD1Cs
                       "COL4A1..Cell..Mean",    #ECM
                       "COL6..Cell..Mean",      #ECM
                       "CD3..Cell..Mean",       #T cells
                       "CD8..Cell..Mean",       #CD8 T cells
                       "SST..Cell..Mean",       #Delta cells
                       "CD4..Cell..Mean",       #CD4T cells
                       "ACTA2..Cell..Mean",     #Smooth Muscle
                       "PPY..Cell..Mean",       #PP cells
                       "SOX9..Nucleus..Mean",   #Pancreas Development
                       "CD11c..Cell..Mean",     #Dendritic Cells
                       "NKX6.1..Nucleus..Mean", #Beta Cells
                       "KRT..Cell..Mean",       #Epithelial
                       "CD31..Cell..Mean",      #Endothelial
                       "PDX1..Nucleus..Mean",   #Pancreas
                       "GHRL..Cell..Mean",      #Ghrelin
                       "CD163..Cell..Mean",     #Macrophage
                       "HLA.DR..Cell..Mean",    #MHC-II
                       "CD68..Cell..Mean",      #Macrophage
                       "CD66b..Cell..Mean",     #Granulocytes
                       "PNLIP..Cell..Mean",     #Acinar cells
                       "CPEP..Cell..Mean",      #Beta Cells
                       "TUBB3..Cell..Mean",     #Neuron
                       "Ki67..Nucleus..Mean")   #Proliferation Status



#T2D FULL
T2D_Full <- T2D_Full_Data %>% select(matches(Full_markers_T2D))
T2D_Full <- data.frame(T2D_Metadata,
                       T2D_Full)



#Integrate Dataframes, remove CD20 marker
#OLD_FULL
#HPAP053_Head - Level1
{
  HPAP053_Head <- T2D_Full %>% filter(Image=="HPAP-053_CODEX_Head-of-pancreas_OCT.ome.tif")
  HPAP053_Head <- HPAP053_Head %>% select(-c("CD20..Cell..Mean"))
  HPAP053_Head <- CreateSeuratObject(counts = t(HPAP053_Head[,-c(1:8)]), assay = "CODEX",meta.data = HPAP053_Head[,c(1:8)])
  HPAP053_Head <- NormalizeData(HPAP053_Head)
  HPAP053_Head <- FindVariableFeatures(HPAP053_Head)
  HPAP053_Head <- ScaleData(HPAP053_Head)
  HPAP053_Head <- RunPCA(HPAP053_Head,reduction.name = "pca")
}
#HPAP-066_Tail - Level1
{
  HPAP066_Tail <- T2D_Full %>% filter(Image=="HPAP-066_CODEX_Tail-of-pancreas_OCT.ome.tif")
  HPAP066_Tail <- HPAP066_Tail %>% select(-c("CD20..Cell..Mean"))
  HPAP066_Tail <- CreateSeuratObject(counts = t(HPAP066_Tail[,-c(1:8)]), assay = "CODEX",meta.data = HPAP066_Tail[,c(1:8)])
  HPAP066_Tail <- NormalizeData(HPAP066_Tail)
  HPAP066_Tail <- FindVariableFeatures(HPAP066_Tail)
  HPAP066_Tail <- ScaleData(HPAP066_Tail)
  HPAP066_Tail <- RunPCA(HPAP066_Tail,reduction.name = "pca")
}
#HPAP-093_Body - Level1
{
  HPAP093_Body <- T2D_Full %>% filter(Image=="HPAP-093_CODEX_Body-of-pancreas_OCT.ome.tif")
  HPAP093_Body <- HPAP093_Body %>% select(-c("CD20..Cell..Mean"))
  HPAP093_Body <- CreateSeuratObject(counts = t(HPAP093_Body[,-c(1:8)]), assay = "CODEX",meta.data = HPAP093_Body[,c(1:8)])
  HPAP093_Body <- NormalizeData(HPAP093_Body)
  HPAP093_Body <- FindVariableFeatures(HPAP093_Body)
  HPAP093_Body <- ScaleData(HPAP093_Body)
  HPAP093_Body <- RunPCA(HPAP093_Body)
}
#HPAP-112_Head - Level1
{
  HPAP112_Head <- T2D_Full %>% filter(Image=="HPAP-112_CODEX_Head-of-pancreas_OCT.ome.tif")
  HPAP112_Head <- HPAP112_Head %>% select(-c("CD20..Cell..Mean"))
  HPAP112_Head <- CreateSeuratObject(counts = t(HPAP112_Head[,-c(1:8)]), assay = "CODEX",meta.data = HPAP112_Head[,c(1:8)])
  HPAP112_Head <- NormalizeData(HPAP112_Head)
  HPAP112_Head <- FindVariableFeatures(HPAP112_Head)
  HPAP112_Head <- ScaleData(HPAP112_Head)
  HPAP112_Head <- RunPCA(HPAP112_Head,reduction.name = "pca")
}
#HPAP-169_Tail - Level1
{
  HPAP169_Tail <- T2D_Full %>% filter(Image=="HPAP-169_CODEX_Tail-of-pancreas_OCT.ome.tif")
  HPAP169_Tail <- HPAP169_Tail %>% select(-c("CD20..Cell..Mean"))
  HPAP169_Tail <- CreateSeuratObject(counts = t(HPAP169_Tail[,-c(1:8)]), assay = "CODEX",meta.data = HPAP169_Tail[,c(1:8)])
  HPAP169_Tail <- NormalizeData(HPAP169_Tail)
  HPAP169_Tail <- FindVariableFeatures(HPAP169_Tail)
  HPAP169_Tail <- ScaleData(HPAP169_Tail)
  HPAP169_Tail <- RunPCA(HPAP169_Tail,reduction.name = "pca")
}
#HPAP174_Body - Level1
{
  HPAP174_Body <- T2D_Full %>% filter(Image=="HPAP-174_CODEX_Body-of-pancreas_OCT.ome.tif")
  HPAP174_Body <- HPAP174_Body %>% select(-c("CD20..Cell..Mean"))
  HPAP174_Body <- CreateSeuratObject(counts = t(HPAP174_Body[,-c(1:8)]), assay = "CODEX",meta.data = HPAP174_Body[,c(1:8)])
  HPAP174_Body <- NormalizeData(HPAP174_Body)
  HPAP174_Body <- FindVariableFeatures(HPAP174_Body)
  HPAP174_Body <- ScaleData(HPAP174_Body)
  HPAP174_Body <- RunPCA(HPAP174_Body,reduction.name = "pca")
}
#HPAP175_Tail - Level1
{
  HPAP175_Tail <- T2D_Full %>% filter(Image=="HPAP-175_CODEX_Tail-of-pancreas_OCT.ome.tif")
  HPAP175_Tail <- HPAP175_Tail %>% select(-c("CD20..Cell..Mean"))
  HPAP175_Tail <- CreateSeuratObject(counts = t(HPAP175_Tail[,-c(1:8)]), assay = "CODEX",meta.data = HPAP175_Tail[,c(1:8)])
  HPAP175_Tail <- NormalizeData(HPAP175_Tail)
  HPAP175_Tail <- FindVariableFeatures(HPAP175_Tail)
  HPAP175_Tail <- ScaleData(HPAP175_Tail)
  HPAP175_Tail <- RunPCA(HPAP175_Tail,reduction.name = "pca")
}

#Young_Full
#HPAP016_Body - Level1
{
  HPAP016_Body <- T2D_Full %>% filter(Image=="HPAP-016_CODEX_Body-of-pancreas_OCT.ome.tif")
  HPAP016_Body <- HPAP016_Body %>% select(-c("CD20..Cell..Mean"))
  HPAP016_Body <- CreateSeuratObject(counts = t(HPAP016_Body[,-c(1:8)]), assay = "CODEX",meta.data = HPAP016_Body[,c(1:8)])
  HPAP016_Body <- NormalizeData(HPAP016_Body)
  HPAP016_Body <- FindVariableFeatures(HPAP016_Body)
  HPAP016_Body <- ScaleData(HPAP016_Body)
  HPAP016_Body <- RunPCA(HPAP016_Body,reduction.name = "pca")
}
#HPAP017_Body - Level1
{
  HPAP017_Body <- T2D_Full %>% filter(Image=="HPAP-017_CODEX_Body-of-pancreas_OCT.ome.tif")
  HPAP017_Body <- HPAP017_Body %>% select(-c("CD20..Cell..Mean"))
  HPAP017_Body <- CreateSeuratObject(counts = t(HPAP017_Body[,-c(1:8)]), assay = "CODEX",meta.data = HPAP017_Body[,c(1:8)])
  HPAP017_Body <- NormalizeData(HPAP017_Body)
  HPAP017_Body <- FindVariableFeatures(HPAP017_Body)
  HPAP017_Body <- ScaleData(HPAP017_Body)
  HPAP017_Body <- RunPCA(HPAP017_Body,reduction.name = "pca")
}
#HPAP018_Head - Level1
{
  HPAP018_Head <- T2D_Full %>% filter(Image=="HPAP-018_CODEX_Head-of-pancreas_OCT.ome.tif")
  HPAP018_Head <- HPAP018_Head %>% select(-c("CD20..Cell..Mean"))
  HPAP018_Head <- CreateSeuratObject(counts = t(HPAP018_Head[,-c(1:8)]), assay = "CODEX",meta.data = HPAP018_Head[,c(1:8)])
  HPAP018_Head <- NormalizeData(HPAP018_Head)
  HPAP018_Head <- FindVariableFeatures(HPAP018_Head)
  HPAP018_Head <- ScaleData(HPAP018_Head)
  HPAP018_Head <- RunPCA(HPAP018_Head,reduction.name = "pca")
}
#HPAP027_Head - Level1
{
  HPAP027_Head <- T2D_Full %>% filter(Image=="HPAP-027_CODEX_Head-of-pancreas_OCT.ome.tif")
  HPAP027_Head <- HPAP027_Head %>% select(-c("CD20..Cell..Mean"))
  HPAP027_Head <- CreateSeuratObject(counts = t(HPAP027_Head[,-c(1:8)]), assay = "CODEX",meta.data = HPAP027_Head[,c(1:8)])
  HPAP027_Head <- NormalizeData(HPAP027_Head)
  HPAP027_Head <- FindVariableFeatures(HPAP027_Head)
  HPAP027_Head <- ScaleData(HPAP027_Head)
  HPAP027_Head <- RunPCA(HPAP027_Head,reduction.name = "pca")
}
#HPAP037_Head - Level1
{
  HPAP037_Head <- T2D_Full %>% filter(Image=="HPAP-037_CODEX_Head-of-pancreas_OCT.ome.tif")
  HPAP037_Head <- HPAP037_Head %>% select(-c("CD20..Cell..Mean"))
  HPAP037_Head <- CreateSeuratObject(counts = t(HPAP037_Head[,-c(1:8)]), assay = "CODEX",meta.data = HPAP037_Head[,c(1:8)])
  HPAP037_Head <- NormalizeData(HPAP037_Head)
  HPAP037_Head <- FindVariableFeatures(HPAP037_Head)
  HPAP037_Head <- ScaleData(HPAP037_Head)
  HPAP037_Head <- RunPCA(HPAP037_Head,reduction.name = "pca")
}
#HPAP040_Head - Level1
{
  HPAP040_Head <- T2D_Full %>% filter(Image=="HPAP-040_CODEX_Head-of-pancreas_OCT.ome.tif")
  HPAP040_Head <- HPAP040_Head %>% select(-c("CD20..Cell..Mean"))
  HPAP040_Head <- CreateSeuratObject(counts = t(HPAP040_Head[,-c(1:8)]), assay = "CODEX",meta.data = HPAP040_Head[,c(1:8)])
  HPAP040_Head <- NormalizeData(HPAP040_Head)
  HPAP040_Head <- FindVariableFeatures(HPAP040_Head)
  HPAP040_Head <- ScaleData(HPAP040_Head)
  HPAP040_Head <- RunPCA(HPAP040_Head,reduction.name = "pca")
}
#HPAP045_Tail - Level1
{
  HPAP045_Tail <- T2D_Full %>% filter(Image=="HPAP-045_CODEX_Tail-of-pancreas_OCT.ome.tif")
  HPAP045_Tail <- HPAP045_Tail %>% select(-c("CD20..Cell..Mean"))
  HPAP045_Tail <- CreateSeuratObject(counts = t(HPAP045_Tail[,-c(1:8)]), assay = "CODEX",meta.data = HPAP045_Tail[,c(1:8)])
  HPAP045_Tail <- NormalizeData(HPAP045_Tail)
  HPAP045_Tail <- FindVariableFeatures(HPAP045_Tail)
  HPAP045_Tail <- ScaleData(HPAP045_Tail)
  HPAP045_Tail <- RunPCA(HPAP045_Tail,reduction.name = "pca")
}
#HPAP056_Head - Level1
{
  HPAP056_Head <- T2D_Full %>% filter(Image=="HPAP-056_CODEX_Head-of-pancreas_OCT.ome.tif")
  HPAP056_Head <- HPAP056_Head %>% select(-c("CD20..Cell..Mean"))
  HPAP056_Head <- CreateSeuratObject(counts = t(HPAP056_Head[,-c(1:8)]), assay = "CODEX",meta.data = HPAP056_Head[,c(1:8)])
  HPAP056_Head <- NormalizeData(HPAP056_Head)
  HPAP056_Head <- FindVariableFeatures(HPAP056_Head)
  HPAP056_Head <- ScaleData(HPAP056_Head)
  HPAP056_Head <- RunPCA(HPAP056_Head,reduction.name = "pca")
}
#HPAP095_Head - Level1
{
  HPAP095_Head <- T2D_Full %>% filter(Image=="HPAP-095_CODEX_Head-of-pancreas_OCT.ome.tif")
  HPAP095_Head <- HPAP095_Head %>% select(-c("CD20..Cell..Mean"))
  HPAP095_Head <- CreateSeuratObject(counts = t(HPAP095_Head[,-c(1:8)]), assay = "CODEX",meta.data = HPAP095_Head[,c(1:8)])
  HPAP095_Head <- NormalizeData(HPAP095_Head)
  HPAP095_Head <- FindVariableFeatures(HPAP095_Head)
  HPAP095_Head <- ScaleData(HPAP095_Head)
  HPAP095_Head <- RunPCA(HPAP095_Head,reduction.name = "pca")
}
#HPAP155_Body - Level1
{
  HPAP155_Body <- T2D_Full %>% filter(Image=="HPAP-155_CODEX_Body-of-pancreas_OCT.ome.tif")
  HPAP155_Body <- HPAP155_Body %>% select(-c("CD20..Cell..Mean"))
  HPAP155_Body <- CreateSeuratObject(counts = t(HPAP155_Body[,-c(1:8)]), assay = "CODEX",meta.data = HPAP155_Body[,c(1:8)])
  HPAP155_Body <- NormalizeData(HPAP155_Body)
  HPAP155_Body <- FindVariableFeatures(HPAP155_Body)
  HPAP155_Body <- ScaleData(HPAP155_Body)
  HPAP155_Body <- RunPCA(HPAP155_Body,reduction.name = "pca")
}

#T2D Cohort - LEVEL 1 AND 2
#HPAP013_BODY
{
  HPAP013_BODY <- T2D_Full %>% filter(Image=="HPAP-013_CODEX_Body-of-pancreas_OCT.ome.tif")
  HPAP013_BODY <- HPAP013_BODY %>% select(-c("CD20..Cell..Mean"))
  HPAP013_BODY <- CreateSeuratObject(counts = t(HPAP013_BODY[,-c(1:8)]), assay = "CODEX",meta.data = HPAP013_BODY[,c(1:8)])
  HPAP013_BODY <- NormalizeData(HPAP013_BODY)
  HPAP013_BODY <- FindVariableFeatures(HPAP013_BODY)
  HPAP013_BODY <- ScaleData(HPAP013_BODY)
  HPAP013_BODY <- RunPCA(HPAP013_BODY,reduction.name = "pca")
}
#HPAP051_TAIL
{
  HPAP051_TAIL <- T2D_Full %>% filter(Image=="HPAP-051_CODEX_Tail-of-Pancreas_FFPE.ome.tif")
  HPAP051_TAIL <- HPAP051_TAIL %>% select(-c("CD19..Cell..Mean"))
  names(HPAP051_TAIL)[11] <- "CD19..Cell..Mean"
  HPAP051_TAIL <- CreateSeuratObject(counts = t(HPAP051_TAIL[,-c(1:8)]), assay = "CODEX",meta.data = HPAP051_TAIL[,c(1:8)])
  HPAP051_TAIL <- NormalizeData(HPAP051_TAIL)
  HPAP051_TAIL <- FindVariableFeatures(HPAP051_TAIL)
  HPAP051_TAIL <- ScaleData(HPAP051_TAIL)
  HPAP051_TAIL <- RunPCA(HPAP051_TAIL,reduction.name = "pca")
}
#HPAP079_TAIL
{
  HPAP079_TAIL <- T2D_Full %>% filter(Image=="HPAP-079_CODEX_Tail-of-Pancreas_FFPE.ome.tif")
  HPAP079_TAIL <- HPAP079_TAIL %>% select(-c("CD19..Cell..Mean"))
  names(HPAP079_TAIL)[11] <- "CD19..Cell..Mean"
  HPAP079_TAIL <- CreateSeuratObject(counts = t(HPAP079_TAIL[,-c(1:8)]), assay = "CODEX",meta.data = HPAP079_TAIL[,c(1:8)])
  HPAP079_TAIL <- NormalizeData(HPAP079_TAIL)
  HPAP079_TAIL <- FindVariableFeatures(HPAP079_TAIL)
  HPAP079_TAIL <- ScaleData(HPAP079_TAIL)
  HPAP079_TAIL <- RunPCA(HPAP079_TAIL,reduction.name = "pca")
}
#HPAP100_TAIL
{
  HPAP100_TAIL <- T2D_Full %>% filter(Image=="HPAP-100_CODEX_Tail-of-Pancreas_FFPE.ome.tif")
  HPAP100_TAIL <- HPAP100_TAIL %>% select(-c("CD19..Cell..Mean"))
  names(HPAP100_TAIL)[11] <- "CD19..Cell..Mean"
  HPAP100_TAIL <- CreateSeuratObject(counts = t(HPAP100_TAIL[,-c(1:8)]), assay = "CODEX",meta.data = HPAP100_TAIL[,c(1:8)])
  HPAP100_TAIL <- NormalizeData(HPAP100_TAIL)
  HPAP100_TAIL <- FindVariableFeatures(HPAP100_TAIL)
  HPAP100_TAIL <- ScaleData(HPAP100_TAIL)
  HPAP100_TAIL <- RunPCA(HPAP100_TAIL,reduction.name = "pca")
 }
#HPAP124_TAIL
{
  HPAP124_TAIL <- T2D_Full %>% filter(Image=="HPAP-124_CODEX_Tail-of-Pancreas_FFPE.ome.tif")
  HPAP124_TAIL <- HPAP124_TAIL %>% select(-c("CD19..Cell..Mean"))
  names(HPAP124_TAIL)[11] <- "CD19..Cell..Mean"
  HPAP124_TAIL <- CreateSeuratObject(counts = t(HPAP124_TAIL[,-c(1:8)]), assay = "CODEX",meta.data = HPAP124_TAIL[,c(1:8)])
  HPAP124_TAIL <- NormalizeData(HPAP124_TAIL)
  HPAP124_TAIL <- FindVariableFeatures(HPAP124_TAIL)
  HPAP124_TAIL <- ScaleData(HPAP124_TAIL)
  HPAP124_TAIL <- RunPCA(HPAP124_TAIL,reduction.name = "pca")
}
#HPAP142_TAIL
{
  HPAP142_TAIL <- T2D_Full %>% filter(Image=="HPAP-142_CODEX_Tail-of-Pancreas_FFPE.ome.tif")
  HPAP142_TAIL <- HPAP142_TAIL %>% select(-c("CD19..Cell..Mean"))
  names(HPAP142_TAIL)[11] <- "CD19..Cell..Mean"
  HPAP142_TAIL <- CreateSeuratObject(counts = t(HPAP142_TAIL[,-c(1:8)]), assay = "CODEX",meta.data = HPAP142_TAIL[,c(1:8)])
  HPAP142_TAIL <- NormalizeData(HPAP142_TAIL)
  HPAP142_TAIL <- FindVariableFeatures(HPAP142_TAIL)
  HPAP142_TAIL <- ScaleData(HPAP142_TAIL)
  HPAP142_TAIL <- RunPCA(HPAP142_TAIL,reduction.name = "pca")
}
#HPAP143_TAIL
{
  HPAP143_TAIL <- T2D_Full %>% filter(Image=="HPAP-143_CODEX_Tail-of-Pancreas_FFPE.ome.tif")
  HPAP143_TAIL <- HPAP143_TAIL %>% select(-c("CD19..Cell..Mean"))
  names(HPAP143_TAIL)[11] <- "CD19..Cell..Mean"
  HPAP143_TAIL <- CreateSeuratObject(counts = t(HPAP143_TAIL[,-c(1:8)]), assay = "CODEX",meta.data = HPAP143_TAIL[,c(1:8)])
  HPAP143_TAIL <- NormalizeData(HPAP143_TAIL)
  HPAP143_TAIL <- FindVariableFeatures(HPAP143_TAIL)
  HPAP143_TAIL <- ScaleData(HPAP143_TAIL)
  HPAP143_TAIL <- RunPCA(HPAP143_TAIL,reduction.name = "pca")
}
#HPAP177_BODY
{
  HPAP177_BODY <- T2D_Full %>% filter(Image=="HPAP-143_CODEX_Tail-of-Pancreas_FFPE.ome.tif")
  HPAP177_BODY <- HPAP177_BODY %>% select(-c("CD19..Cell..Mean"))
  names(HPAP177_BODY)[11] <- "CD19..Cell..Mean"
  HPAP177_BODY <- CreateSeuratObject(counts = t(HPAP177_BODY[,-c(1:8)]), assay = "CODEX",meta.data = HPAP177_BODY[,c(1:8)])
  HPAP177_BODY <- NormalizeData(HPAP177_BODY)
  HPAP177_BODY <- FindVariableFeatures(HPAP177_BODY)
  HPAP177_BODY <- ScaleData(HPAP177_BODY)
  HPAP177_BODY <- RunPCA(HPAP177_BODY,reduction.name = "pca")
}



# Find integration anchors and cluster----
anchors <- FindIntegrationAnchors(
  object.list = list(HPAP013_BODY,HPAP016_Body,HPAP017_Body,
                     HPAP018_Head,HPAP027_Head,HPAP037_Head,
                     HPAP040_Head,HPAP045_Tail,HPAP051_TAIL,
                     HPAP053_Head,HPAP056_Head,HPAP066_Tail,
                     HPAP079_TAIL,HPAP093_Body,HPAP095_Head,
                     HPAP100_TAIL,HPAP112_Head,HPAP124_TAIL,
                     HPAP142_TAIL,HPAP143_TAIL,HPAP155_Body,
                     HPAP169_Tail,HPAP174_Body,HPAP175_Tail,
                     HPAP177_BODY),
  anchor.features = 30,
  reduction = "rpca",
  dims = 1:15
)

combined <- IntegrateData(
  anchorset = anchors,
  dims = 1:15
)




#Run standard Seurat Pipeline
CODEX_seurat <- NormalizeData(combined, normalization.method = "CLR", margin = 2, verbose = FALSE)
DefaultAssay(combined) <- "integrated"
combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
ElbowPlot(combined, ndims = 30) #Decide how many dimensions to use

combined <- FindNeighbors(combined, dims = 1:20, verbose = TRUE) #15minutes
combined <- FindClusters(combined, resolution = 0.6) #45minutes
combined <- RunUMAP(combined, dims = 1:20) #40minutes



# Check - DimPlot for each resolution
#resolution <- "CODEX_snn_res.0.6"

DimPlot(combined,
        split.by = "Category",
        label = TRUE) + ggtitle("Resolution 0.6_dim20")

#Top_marker analysis
#Change Idents
table(Idents(combined))

CODEX_Markers <- FindAllMarkers(combined,
                              only.pos = TRUE,
                              min.pct = 0.25,
                              logfc.threshold = 0.25
                              )

Markers_Seurat_top10 <- CODEX_Markers %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = 10, with_ties = FALSE) %>%
  arrange(cluster, desc(avg_log2FC))



#Validate by Featureplot
table(Idents(combined))

new.cluster.ids <- c(
  "Immune",          #Cluster 0
  "Endothelial",     #Cluster 1
  "Acinar",          #Cluster 2
  "Acinar",          #Cluster 3
  "Ductal",          #Cluster 4
  "Immune",          #Cluster 5
  "Endocrine",       #Cluster 6
  "Endocrine",       #Cluster 7
  "SMA",             #Cluster 8
  "Immune",          #Cluster 9
  "Endocrine",       #Cluster 10
  "Immune",          #Cluster 11
  "Immune",          #Cluster 12
  "Acinar",          #Cluster 13
  "Immune",          #Cluster 14
  "Endocrine",       #Cluster 15
  "MSC",             #Cluster 16
  "Other",           #Cluster 17
  "Endocrine",       #Cluster 18
  "Endocrine",       #Cluster 19
  "Acinar",          #Cluster 20
  "Immune",          #Cluster 21
  "Other")           #Cluster 22

table(Idents(combined))
names(new.cluster.ids) <- levels(combined)
CODEX_Seurat <- RenameIdents(combined, new.cluster.ids)

#Subset out "OTHER" cells
CODEX_Seurat <- subset(CODEX_Seurat, idents = c("Immune","Endothelial","Acinar","Ductal","Endocrine","SMA","MSC"))

DimPlot(CODEX_Seurat,
        split.by = "Group",
        reduction = "umap", label = TRUE, repel = F)

#Figure 5B: CODEX Dimplot from all cells
DimPlot(CODEX_Seurat,
        reduction = "umap",
        label = F,
        repel = TRUE,
        label.size = 4,
        pt.size = 0.8,
        cols = c(
          "#e377c2","#d62728","#1f77b4","#ff7f0e","#2ca02c","#17becf","#bcbd22")
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
    title = "Pancreas (UMAP)",
    #subtitle = "Dimensional reduction colored by Seurat cluster"
  )

# Subset Immune Markers
Immune_seurat <- subset(CODEX_Seurat, idents = "Immune")

#Re-scale and PCA data.  Normalization not necessary because this is subset from previously normalized data.
Immune_seurat <- ScaleData(Immune_seurat, features = rownames(Immune_seurat), verbose = FALSE)
Immune_seurat <- RunPCA(Immune_seurat, features = rownames(Immune_seurat), npcs = 30, verbose = FALSE)
ElbowPlot(Immune_seurat, ndims = 30) #Decide how many dimensions to use
Immune_seurat <- FindNeighbors(Immune_seurat, dims = 1:15, verbose = TRUE) #15 minutes
Immune_seurat <- FindClusters(Immune_seurat, resolution = c(0.6)) #45 minutes
Immune_seurat <- RunUMAP(Immune_seurat, dims = 1:15) #40 minutes

#Sanity Check
DimPlot(Immune_seurat,
        split.by = "Group",
        label = TRUE) + ggtitle("Immune_res_0.6")

#Find top markers at resolution =0.6
table(Idents(Immune_seurat))
markers_Immune_seurat <- FindAllMarkers(Immune_seurat,
                                     only.pos = TRUE,
                                     min.pct = 0.25,
                                     logfc.threshold = 0.25
                                     )
#Filter our the top 10 markers
markers_Immune_seurat_top10 <- markers_Immune_seurat %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = 10, with_ties = FALSE) %>% arrange(cluster, desc(avg_log2FC))

#Define Cell marker clusters for immune cells.
new.cluster.ids <- c("CD4 T",  #Cluster 0
  "M1 Macrophage",             #Cluster 1
  "Other",                     #Cluster 2
  "CD8 T",                     #Cluster 3
  "CD8 T",                     #Cluster 4
  "cDC",                       #Cluster 5
  "Granulocyte/Monocyte",      #Cluster 6
  "Other",                     #Cluster 7
  "M2 Macrophage",             #Cluster 8
  "Other",                     #Cluster 9
  "CD8 T",                     #Cluster 10
  "Other",                     #Cluster 11
  "CD8 T",                     #Cluster 12
  "Other",                     #Cluster 13
  "CD4 T",                     #Cluster 14
  "M1 Macrophage",             #Cluster 15
  "CD8 T",                     #Cluster 16
  "Other",                     #Cluster 17
  "Other")                     #Cluster 18


names(new.cluster.ids) <- levels(Immune_seurat)
Immune_seurat <- RenameIdents(Immune_seurat, new.cluster.ids)
DimPlot(Immune_seurat, reduction = "umap",
        split.by = "Category",
        label = TRUE, repel = F)

#re-run with new clusters.You will need to re-scale and PCA the UMAP at lines 707-712 before naming new clusters.
Immune_seurat <- subset(Immune_seurat, idents = c("CD4 T",
                                                  "M1 Macrophage",
                                                  "CD8 T",
                                                  "cDC",
                                                  "Granulocyte/Monocyte",
                                                  "M2 Macrophage"))
                        
#Figure 5B: Dimplot for immune filtered cells
DimPlot(Immune_seurat,
        reduction = "umap",
        label = F,
        repel = TRUE,
        label.size = 4,
        pt.size = 0.8,
        cols = c("cyan","#d62728","#e377c2","#7f7f7f","#2ca02c","#ff7f0e","#bcbd22")) +
    theme_void() +
    theme(legend.title = element_blank(),
          legend.text = element_text(size = 12),
          plot.title = element_blank(),
          plot.subtitle = element_text(size = 14, hjust = 0.5),
          axis.text = element_text(size=12),
          axis.ticks = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
          legend.position ="none",
          plot.margin = margin(10, 10, 10, 10)) +
    labs(title = "Immune Cells (UMAP)")




#Subset endocrine markers
Endocrine_seurat <- subset(CODEX_Seurat, idents = "Endocrine")

Endocrine_seurat <- ScaleData(Endocrine_seurat, features = rownames(Endocrine_seurat), verbose = FALSE)
Endocrine_seurat <- RunPCA(Endocrine_seurat, features = rownames(Endocrine_seurat), npcs = 30, verbose = FALSE)
ElbowPlot(Endocrine_seurat, ndims = 30) #Decide how many dimensions to use
Endocrine_seurat <- FindNeighbors(Endocrine_seurat, dims = 1:10, verbose = TRUE) #15 minutes
Endocrine_seurat <- FindClusters(Endocrine_seurat, resolution = c(0.6)) #45 minutes
Endocrine_seurat <- RunUMAP(Endocrine_seurat, dims = 1:10) #40 minutes

#Sanity check
DimPlot(Endocrine_seurat,split.by ="Group",label = TRUE) + ggtitle("Endocrine_res0.6")

#Identify cluster markers
table(Idents(Endocrine_seurat))
markers_Endocrine_seurat <- FindAllMarkers(Endocrine_seurat,
                                        only.pos = TRUE,
                                        min.pct = 0.25,
                                        logfc.threshold = 0.25)

#Filter our the top 10 markers
markers_Endocrine_seurat_top10 <- markers_Endocrine_seurat %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = 10, with_ties = FALSE) %>% arrange(cluster, desc(avg_log2FC))

#name new endocrine cluster IDs
new.cluster.ids <- c("Alpha",        #Cluster 0
                     "Beta",         #Cluster 1
                     "Delta",        #Cluster 2
                     "Alpha",        #Cluster 3
                     "Beta",         #Cluster 4
                     "Alpha",        #Cluster 5
                     "Beta",         #Cluster 6
                     "Other",        #Cluster 7
                     "Delta",        #Cluster 8
                     "Epsilon",      #Cluster 9
                     "Stromal",      #Cluster 10
                     "PP",           #Cluster 11
                     "Beta",         #Cluster 12
                     "Stromal",      #Cluster 13
                     "Other",        #Cluster 14
                     "Endothelial")  #Cluster 15

Idents(Endocrine_seurat) <- "integrated_snn_res.0.6"
table(Endocrine_seurat$seurat_clusters)


names(new.cluster.ids) <- levels(Endocrine_seurat)
Endocrine_seurat <- RenameIdents(Endocrine_seurat, new.cluster.ids)
#sanity check for new cluster IDs
DimPlot(Endocrine_seurat, reduction = "umap", label = TRUE, repel = F)

# subset out "Others" category
Endocrine_seurat <- subset(Endocrine_seurat, idents = c("Alpha",
                                                        "Beta",
                                                        "Delta",
                                                        "Epsilon",
                                                        "Stromal",
                                                        "PP",
                                                        "Endothelial"))

#Figure 5B: Dimplot for endocrine filtered cells
DimPlot(Endocrine_Seurat,
        reduction = "umap",
        label = F,
        repel = TRUE,
        label.size = 4,
        pt.size = 0.4,
        cols = c(
          "#d62728","#2ca02c","#bcbd22","cyan","#1f7a1f","#e377c2","#ff7f0e")) +
  theme_void() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 12),
        plot.title = element_blank(),
        plot.subtitle = element_text(size = 14, hjust = 0.5),
        axis.text = element_text(size=12),
        axis.ticks = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
        legend.position ="none",
        plot.margin = margin(10, 10, 10, 10)) +
  labs(title = "Endocrine Cells (UMAP)",)


####################################################################################
############################Group Statistics - Endocrine############################
####################################################################################

Codex_population_Endocrine <- data.frame(table(Idents(Endocrine_Seurat),
                                               Endocrine_Seurat@meta.data$DonorID))
Codex_population_Endocrine$Group <- ifelse(Codex_population_Endocrine$Var2 %in% c("HPAP-013","HPAP-051","HPAP-079","HPAP-100",
                                                                                  "HPAP-124","HPAP-142","HPAP-143","HPAP-177"),"T2D","ND")

names(Codex_population_Endocrine) <- c("Celltype","DonorID","value","Group")
Codex_population_Endocrine <- Codex_population_Endocrine %>% filter(Celltype %in% c("Alpha","Beta","Stromal","Delta","PP","Epsilon"))
Codex_population_Endocrine$Celltype <- factor(Codex_population_Endocrine$Celltype, levels = c("Alpha","Beta","Stromal","Delta","PP","Epsilon"))
## Note - the "Stromal" population was changed to 'fibrotic' beta cell due to Col6+Col4A1+ marker identification subsets. ##

#Endocrine Cell Proportions
Codex_population_Endocrine <- Codex_population_Endocrine %>%
  group_by(DonorID, Group) %>%
  mutate(TotalCells=sum(value))

Codex_population_Endocrine$Proportion <- (Codex_population_Endocrine$value/Codex_population_Endocrine$TotalCells) * 100

#Figure 5C - Stacked Barplot for endocrine cell populations
ggplot(Codex_population_Endocrine, aes(x = Group, y = Proportion, fill = Celltype)) +
  geom_bar(stat = "identity", position = "fill") +
 # facet_wrap(~ Group) +
  theme_classic(base_size = 10) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.25),              # 0, 0.25, 0.5, 0.75, 1
    labels = c(0, 25, 50, 75, 100)) +         # show as percent
  scale_fill_manual(values=c("#d62728","#2ca02c","#1f7a1f","#bcbd22","#e377c2","cyan")) +
  labs(y = "Proportion", x = "Subject", fill = "Cell Type") +
  theme(#legend.position = 'none',
        axis.title.x =element_blank())



#Bargraph per celltype
Endocrine_summary <- Codex_population_Endocrine %>%
  group_by(Group, Celltype) %>%
  summarise(
    mean_prop = mean(Proportion, na.rm = TRUE),
    sem = sd(Proportion, na.rm = TRUE) / sqrt(n())
  ) %>%
  ungroup()

#Summary statistics per population
Codex_population_Endocrine <- Codex_population_Endocrine %>% filter(Celltype!="NA")
Endocrine_summary <- Endocrine_summary %>% filter(Celltype!="NA")

#Figure 5C - Endocrine cell population by individual donor
ggplot(Endocrine_summary, aes(x = Group, y = mean_prop, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.9,color='black') +
  geom_errorbar(aes(ymin = mean_prop - sem, ymax = mean_prop + sem),width = 0.2, position = position_dodge(0.7)) +
  geom_jitter(data=Codex_population_Endocrine, aes(x=Group,y=Proportion),shape=21,fill='white',width=0.3) +
  facet_wrap2(~ Celltype, ncol = 6, scales = "fixed",
              strip = strip_themed(
                background_x = elem_list_rect(
                  fill = c("#d62728","#2ca02c","#1f7a1f","#bcbd22","#e377c2","cyan"),
                  color = "black"))) +
  scale_fill_manual(values = c("#D0DCB8","#B4A1C1")) +
  theme_classic(base_size = 16) +
  labs(y = "% Endocrine cells",
       x = NULL,
       fill = "Group") +
  theme(strip.text = element_text(face = "bold", size = 14),
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),
        axis.title.y = element_text(face = "bold", size = 14),
        legend.position = "none",
        legend.title = element_blank())

Endocrine_ttest_results <- Codex_population_Endocrine %>%
  group_by(Celltype) %>%
  summarise(
    ttest = list(t.test(Proportion ~ Group, data = cur_data())),
    p.value = ttest[[1]]$p.value,
    t.statistic = ttest[[1]]$statistic,
    df = ttest[[1]]$parameter,
    mean_group1 = ttest[[1]]$estimate[1],
    mean_group2 = ttest[[1]]$estimate[2]
  )

####################################################################################
############################Group Statistics - Endocrine############################
####################################################################################
Codex_population_Immune <- data.frame(table(Idents(Immune_seurat),
                                            Immune_seurat@meta.data$DonorID))
Codex_population_Immune$Group <- ifelse(Codex_population_Immune$Var2 %in% c("HPAP-013","HPAP-051","HPAP-079","HPAP-100",
                                                                            "HPAP-124","HPAP-142","HPAP-143","HPAP-177"),"T2D","ND")

#Filter out CDC1 and B Cell populations (too noisy)
Codex_population_Immune <- Codex_population_Immune %>% filter(!Var1 %in% c("cDC","B Cell"))
names(Codex_population_Immune) <- c("Celltype","DonorID","value","Group")
Codex_population_Immune$Celltype <- factor(Codex_population_Immune$Celltype, levels = c("M1 Macrophage","M2 Macrophage","Granulocyte/Monocyte","CD8 T","CD4 T"))

Codex_population_Immune <- Codex_population_Immune %>%
  group_by(DonorID, Group) %>%
  mutate(TotalCells=sum(value))

Codex_population_Immune$Proportion <- (Codex_population_Immune$value/Codex_population_Immune$TotalCells) * 100

#Figure 5F - Stacked Barplot for immune cell populations
ggplot(Codex_population_Immune, aes(x = Group, y = Proportion, fill = Celltype)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_classic(base_size = 12) +
  scale_fill_manual(values = c("#d62728","orange","#2ca02c","#e377c2","cyan","#bcbd22","#1f77b4","#7f7f7f")) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), labels = c(0, 25, 50, 75, 100)) +   # show as percent
  labs(y = "Percent", x = "Subject", fill = "Cell Type") +
  theme(legend.position = 'none',
        axis.title.x = element_blank())

#Bargraph per celltype
Immune_summary <- Codex_population_Immune %>%
  group_by(Group, Celltype) %>%
  summarise(mean_prop = mean(Proportion, na.rm = TRUE),
            sem = sd(Proportion, na.rm = TRUE) / sqrt(n())) %>%
  ungroup()

#Figure 5F - Immune cell population by individual donor
ggplot(Immune_summary, aes(x = Group, y = mean_prop, fill = Group)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.9,color='black') +
    geom_errorbar(aes(ymin = mean_prop - sem, ymax = mean_prop + sem), width = 0.2, position = position_dodge(0.7)) +
    geom_jitter(data=Codex_population_Immune, aes(x=Group,y=Proportion),shape=21,fill='white',width=0.3) +
    facet_wrap2(~ Celltype, ncol = 8, scales = "fixed",
                strip = strip_themed(
                background_x = elem_list_rect(fill = c("orange","#d62728","#2ca02c","#e377c2","cyan","#bcbd22","#1f77b4","#7f7f7f"),
                color = "black"))) +
  scale_fill_manual(values = c("#D0DCB8","#B4A1C1")) +
  theme_classic(base_size = 14) +
  labs(y = "% Immune Cells",
       x = NULL,
       fill = "Group") +
    theme(strip.text = element_text(face = "bold", size = 11),
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

#Immune Cell Density
#This spreadsheet will be used to present area per image, to present density for area.  Fill located in Zenodo (https://zenodo.org/records/20496994)
Immune_Density <- read.csv("Codex_population_Immune_Density.csv")
Immune_Density$Celltype <- factor(Immune_Density$Celltype, levels = c("M2 Macrophage","M1 Macrophage","Granulocyte/Monocyte",
                                                                        "CD8 T","CD4 T"))
    
Immune_Density <- Immune_Density %>%
  group_by(Group, Celltype) %>%
  summarise(mean_prop = mean(Cells.mm2, na.rm = TRUE),
            sem = sd(Cells.mm2, na.rm = TRUE) / sqrt(n())) %>%
  ungroup()
  
Immune_Density$Celltype <- factor(Immune_Density$Celltype, levels = c("M2 Macrophage","M1 Macrophage","Granulocyte/Monocyte",
                                                                      "CD8 T","CD4 T"))
#Figure 5G - Immune Density Plot
ggplot(Immune_Density, aes(x = Group, y = mean_prop, fill = Group)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.9,color='black') +
    geom_errorbar(aes(ymin = mean_prop - sem, ymax = mean_prop + sem), width = 0.2, position = position_dodge(0.7)) +
    geom_jitter(data=Immune_Density, aes(x=Group,y=Cells.mm2),shape=21,fill='white',width=0.3) +
    facet_wrap2(~ Celltype, ncol = 8, scales = "fixed",strip = strip_themed(background_x = elem_list_rect(
      fill = c("orange","#d62728","#2ca02c","#e377c2","cyan","#bcbd22","#1f77b4","#7f7f7f"),
      color = "black"))) +
  scale_fill_manual(values = c("#D0DCB8","#B4A1C1")) +
  theme_classic(base_size = 14) +
  labs(y = "Immune Cells per mm2",
       x = NULL,
       fill = "Group") +
  theme(strip.text = element_text(face = "bold", size = 11),
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),
        axis.title.y = element_text(face = "bold", size = 14),
        legend.position = "none",
        legend.title = element_blank())

Immune_Density_Stat <- Immune_Density %>%
  group_by(Celltype) %>%
  summarise(ttest = list(t.test(Cells.mm2 ~ Group, data = cur_data())),
            p.value = ttest[[1]]$p.value,
            t.statistic = ttest[[1]]$statistic,
            df = ttest[[1]]$parameter,
            mean_group1 = ttest[[1]]$estimate[1],
            mean_group2 = ttest[[1]]$estimate[2])




#Endocrine Cells
DefaultAssay(Endocrine_seurat)
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
Beta$Group <- ifelse(Beta$DonorID %in% c("HPAP-013",
                                         "HPAP-051",
                                         "HPAP-079",
                                         "HPAP-100",
                                         "HPAP-066",
                                         "HPAP-056",
                                         "HPAP-124",
                                         "HPAP-142",
                                         "HPAP-143"),"T2D","ND")
#ttest
results <- Beta %>%
  select(-c("Identifier", "DonorID", "Group")) %>%
  map_df(~{t_res <- t.test(.x[Beta$Group == unique(Beta$Group)[1]],
                    .x[Beta$Group == unique(Beta$Group)[2]])
  tibble(p.value = t_res$p.value,
         statistic = t_res$statistic,
         mean_group1 = t_res$estimate[1],
         mean_group2 = t_res$estimate[2])
  }, .id = "Variable")

print(results)

#Alpha_Cells - repeat from Endocrine seurat data, but the alpha cell fraction
Alpha <- data.frame(Endocrine_seurat_Split$Alpha)
Alpha <- Alpha %>% rownames_to_column("Identifier") %>%        # make rownames a column
  mutate(DonorID = str_sub(Identifier, -8))                    # take last 8 characters
Alpha$Group <- ifelse(Alpha$DonorID %in% c("HPAP-013",
                                           "HPAP-051",
                                           "HPAP-079",
                                           "HPAP-100",
                                           "HPAP-066",
                                           "HPAP-056",
                                           "HPAP-124",
                                           "HPAP-142",
                                           "HPAP-143"),"T2D","ND")
#ttest
results <- Alpha %>%
  select(-c("Identifier", "DonorID", "Group")) %>%
  map_df(~{t_res <- t.test(.x[Alpha$Group == unique(Alpha$Group)[1]],
                    .x[Alpha$Group == unique(Alpha$Group)[2]])
  tibble(p.value = t_res$p.value,
         statistic = t_res$statistic,
         mean_group1 = t_res$estimate[1],
         mean_group2 = t_res$estimate[2])
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


#Figure 5D - Beta cell scaled expression graph.  Below is example for NKX6-1.  Repeat with Chromagranin-A and C-peptide cell intensity markers where applicable   
ggplot(Beta_summary, aes(x = Group, y = NKX6.1..Nucleus..Mean, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.9,color='black') +
  geom_errorbar(aes(ymin = NKX6.1..Nucleus..Mean_mean - NKX6.1..Nucleus..Mean_sem, ymax = NKX6.1..Nucleus..Mean_mean + NKX6.1..Nucleus..Mean_sem),
                width = 0.2, position = position_dodge(0.7)) +
  geom_jitter(data=Beta, aes(x=Group,y=NKX6.1..Nucleus..Mean),shape=21,fill='white',width=0.1) +
  scale_fill_manual(values = c("#D0DCB8","#B4A1C1")) +
  theme_classic(base_size = 12) +
  labs(title="CODEX",y = "NKX6-1",x = NULL,fill = "Group") +
  theme(plot.title = element_text(hjust=0.5),
        strip.text = element_text(face = "bold", size = 14),
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),
        axis.title.y = element_text(face = "bold", size = 14),
        legend.position = "none",
        legend.title = element_blank())


#Figure 5E - Alpha cell SOX9 scaled nuclear expression graph.  
#Cell Summary
Alpha_summary <- Alpha[,-c(1,41)] %>%
  group_by(Group) %>%
  summarise(across(
    where(is.numeric),
    list(
      mean = ~mean(.x, na.rm = TRUE),
      sem  = ~sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x)))
    ),
    .names = "{.col}_{.fn}"))

ggplot(Alpha_summary, aes(x = Group, y = SOX9..Nucleus..Mean_mean, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.9,color='black') +
  geom_errorbar(aes(ymin = SOX9..Nucleus..Mean_mean - SOX9..Nucleus..Mean_sem, ymax = SOX9..Nucleus..Mean_mean + SOX9..Nucleus..Mean_sem),
                width = 0.2, position = position_dodge(0.7)) +
  geom_jitter(data=Beta, aes(x=Group,y=SOX9..Nucleus..Mean),shape=21,fill='white',width=0.1) +
  scale_fill_manual(values = c("#D0DCB8","#B4A1C1")) +
  theme_classic(base_size = 12) +
  labs(title="CODEX",
       y = "Ki67_Nucleus",
       x = NULL,
       fill = "Group") +
  theme(plot.title = element_text(hjust=0.5),
        strip.text = element_text(face = "bold", size = 14),
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),
        axis.title.y = element_text(face = "bold", size = 14),
        legend.position = "none",
        legend.title = element_blank())

#End



