# =======================================================================================================================
# Alpha cell inflammation during human pancreas aging and type 2 diabetes and its reversal by calorie restriction in mice
# -----------------------------------------------------------------------------------------------------------------------
# Author:      Michael Schleh
# Repository:  https://github.com/ArrojoDrigoLab/AlphaCellAging
# Journal:     Science Advances

# Description:
#   Generates Figure 4A-B,F,K-L 
#     - A) RDS setup, cleaning, and single cell clustering
#     - B) Second-level single-cell clustering
#     - C) Population analysis
#     - D) HLA-ABC and HLA-DR single cell expression

# Usage:
#   1. Download packages and set the working directory.
#   2. Ensure the required input files exist in `Figure 4/' (see README).
#   3. Run: RScript_Figure 4_ImagingMassCytometry (Young vs Old).R
#
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Setup: Working 
# ----------------------------------------------------------------------------
#1A. Working Directory 
setwd("")    #<--- Set to source file location

#1B. Download necessary packages
install.packages(c(
  "tidyverse", "stringr", "Seurat", "reshape2", "patchwork",
  "readxl", "Rtsne", "pheatmap", "RColorBrewer",
  "broom", "future", "scales"))

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("Nebulosa"))

#1B. Library packages
library(tidyverse)
library(stringr)
library(Seurat)
library(reshape2)
library(patchwork)
library(readxl)
library(Rtsne)
library(pheatmap)
library(RColorBrewer)
library(broom)
library(future)
library(scales)
library(Nebulosa)

#1C - Set future memory limit to 50 GB
options(future.globals.maxSize = 50 * 1024^3)
# Use 4 parallel workers
plan("multisession", workers = 4)


###################################################################
#2 - Upload Data and Assign MetaData (IMC = Imaging Mass Cytometry)
###################################################################
IMC <- read.csv(file="Figure 4_IMC_QuPath_Cell_measurements_merged_metadata_celltype.csv",header = T,sep=",")        #See zenodo for raw Qupath exported csv file
IMC <- IMC %>% mutate(group = if_else(str_starts(SampleAge, regex("^(2[2-9]yo|3[0-9]yo|40yo)")), "young","old"))     #Create Young/Old Identifiers, factor levels 
IMC$group <- factor(IMC$group, levels = c("young","old"))                                                            #Factor identifiers
IMC <- IMC %>% mutate(age_num = as.numeric(str_extract(SampleAge, "\\d+")))                                          #Cleanup


####### Supplementary MetaData ####### 
#2A Supplementary Figure 5B: MetaData_Summary
#Summary Statistics Setup
#Age
Metadata <- read.csv(file="IMC_Donors.csv",header = T,sep=",")   #Sample MetaData sourced in Github ("AlphaCellAging/Figure 4/Figure 4_IMC Donors.csv")
Metadata$Group <- factor(Metadata$Group, levels = c("Young","Old"))

Metadata_measure <- c("Age","BMI","HbA1c")
Metadata_summary <- Metadata %>%
  group_by(Group) %>%
  summarise(across(
    all_of(Metadata_measure),
    list(
      mean   = ~round(mean(.x, na.rm = TRUE), 2),
      sd     = ~round(sd(.x, na.rm = TRUE), 2),
      sem    = ~round(sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x))), 2))))

#fig S5B.i) Age
ggplot() +
  geom_col(data=Metadata_summary, aes(x=Group,y=Age_mean,fill=Group),color='black') +
  geom_errorbar(data=Metadata_summary,mapping = aes(x=Group, ymin = Age_mean-Age_sd, ymax = Age_mean+Age_sd),width=0.5) +
  geom_jitter(data=Metadata, aes(x=Group,y=Age),width=0.2,alpha=0.5) +
  scale_fill_manual(values=c("#E0DCB8","#607C83"))+
  labs(title="Age",
       y="years")+
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5),
        legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=12))
#fig S5B.ii) HbA1c
  ggplot() +
    geom_col(data=Metadata_summary, aes(x=Group,y=HbA1c_mean,fill=Group),color='black') +
    geom_errorbar(data=Metadata_summary,mapping = aes(x=Group, ymin = HbA1c_mean-HbA1c_sd, ymax = HbA1c_mean+HbA1c_sd),width=0.5) +
    geom_jitter(data=Metadata, aes(x=Group,y=HbA1c),width=0.2,alpha=0.5) +
    scale_fill_manual(values=c("#E0DCB8","#607C83"))+
    labs(title="HbA1c",
         y="%")+
    theme_classic() +
    theme(plot.title = element_text(hjust=0.5),
          legend.position = 'none',
          axis.title.x = element_blank(),
          axis.title.y = element_text(size=12))
#fig S5B.iii) BMI
ggplot() +
    geom_col(data=Metadata_summary, aes(x=Group,y=BMI_mean,fill=Group),color='black') +
    geom_errorbar(data=Metadata_summary,mapping = aes(x=Group, ymin = BMI_mean-BMI_sd, ymax = BMI_mean+BMI_sd),width=0.5) +
    geom_jitter(data=Metadata, aes(x=Group,y=BMI),width=0.2,alpha=0.5) +
    scale_fill_manual(values=c("#E0DCB8","#607C83"))+
    labs(title="BMI",
         y="Kg/m2")+
    theme_classic() +
    theme(plot.title = element_text(hjust=0.5),
          legend.position = 'none',
          axis.title.x = element_blank(),
          axis.title.y = element_text(size=12))


###################################################################
#2B - IMC Seurat Clusters ----
#2Bi) Setup MetaData for Seurat image convergence
#Setup Donor ID
IMC_Metadata <- IMC %>%
  mutate(ID = substr(Image, 1, 8)) %>%
  tibble::rownames_to_column("CellID")  # ensure each row has a unique ID

IMC_Metadata <- data.frame(CellID=rownames(IMC_Metadata),
                           Image=IMC_Metadata$Image,
                           CentroidX=IMC_Metadata$Centroid.X.px,
                           Centroidy=IMC_Metadata$Centroid.Y.px,
                           DonorID=IMC_Metadata$ID,
                           Age=IMC_Metadata$age_num,
                           Group=IMC_Metadata$group,
                           row.names = rownames(IMC_Metadata))

#Setup Markers to use for IMC Analysis
#raw markers
  {
raw_markers <- c("glucagon",              #Alpha Cell
                 "c.peptide", "nkx6.1",   #Beta Cells
                 "somatostatin",          #Delta Cells
                 "pp",                    #PP Cells
                 "ghrelin",               #Ghrelin
                 "ca2",                   #Acinar Cells
                 "cd31", "nestin",        #Endothelial Cells
                 "cd45",                  #Immune Cells
                 "cd3","cd8","cd4",       #T cells
                 "cd44",                  #Activated T cells
                 "collagentype1",         #ECM
                 "cd56",                  #NK
                 "cd14", "cd16",          #Monocytes
                 "cd11b",                 #myeloid cells
                 "cd68","cd163",          #Macrophages
                 "cd20",                  #B cells
                 "chga",                  #neuroendocrine cells
                 "cd57",                  #NK-T cells
                 "foxp3",                 #Tregs
                 "iapp",                  #Beta Cells
                 "pdx.1",                 #Beta cells and endocrine progenitors
                 "hla.dr","hla.abc",      #Antigens
                 "nfkb",                  #Inflammatory activation
                 "granzymeb",             #Cytotoxic lymphocytes
                 "cd45ro")                #Memory T cells
}
#Wrangle for for Seurat workflow  
IMC_selection <- data.frame(CellID=rownames(IMC),
                            IMC %>% select(c("Image","Centroid.X.px","Centroid.Y.px",
                                                    "DonorID","SampleSex","group","age_num")),
                            IMC %>% select(any_of(raw_markers)))
  
# Read the CSV file (contains Image column)
IMC_images <- read.csv("Figure 4_IMC_ImageNames.csv", stringsAsFactors = FALSE)            #Image names sourced in Github ("AlphaCellAging/Figure 4/Figure 4_IMC Donors.csv")

# Create a function to process each image
process_image <- function(image_name, data) {
  message("Processing: ", image_name)
    
# Filter for this image
df <- data %>% filter(Image == image_name)
    
# Create the Seurat object
seurat_obj <- CreateSeuratObject(counts = t(df[,-c(1:8)]), assay = "CODEX", meta.data = df[,c(1:8)])
    
    # Standard Seurat preprocessing
    seurat_obj <- NormalizeData(seurat_obj)
    seurat_obj <- FindVariableFeatures(seurat_obj)
    seurat_obj <- ScaleData(seurat_obj)
    seurat_obj <- RunPCA(seurat_obj)
    
    return(seurat_obj)
  }
  
  # Loop through image list
  image_names <- IMC_images$Image
  
  # Optional: clean names to valid R variable format
  clean_varnames <- image_names %>%
    str_replace_all("[^A-Za-z0-9]+", "_") %>%  # replace special chars with _
    str_replace("_ome_tiff_Image0$", "")       # optional cleanup
  
  # Process and assign each image as its own variable
  for (i in seq_along(image_names)) {
    img <- image_names[i]
    varname <- clean_varnames[i]
    
    assign(varname, process_image(img, IMC_selection))
  }
  
# Optional: save all Seurat objects as .rds
  for (nm in clean_varnames) {
    saveRDS(get(nm), file = paste0("SeuratObjects_Image/", nm, ".rds"))
  }

# Find integration anchors and cluster----
IMC_names <- ls(pattern = "^HPAP")   # list of character names
IMC_list <- mget(IMC_names)          # convert to actual Seurat objects

IMC_anchors <- FindIntegrationAnchors(object.list = IMC_list,anchor.features = 4000,reduction = "rpca")

rm(list = setdiff(ls(), "IMC_anchors"))


integrated_data <- IntegrateData(anchors = IMC_anchors,dims = 1:30)
integrated_data
  
#Run standard Seurat Pipeline
#IMC_seurat <- NormalizeData(integrated_data, normalization.method = "CLR", margin = 2, verbose = FALSE)
IMC_Seurat <- ScaleData(integrated_data, features = rownames(integrated_data), verbose = FALSE)
IMC_Seurat <- RunPCA(IMC_Seurat, features = rownames(IMC_Seurat), npcs = 30, verbose = FALSE)
ElbowPlot(IMC_Seurat, ndims = 30) #Decide how many dimensions to use
IMC_Seurat <- FindNeighbors(IMC_Seurat, dims = 1:12, verbose = TRUE) #15 minutes
IMC_Seurat <- FindClusters(IMC_Seurat, resolution = 0.6) #45 minutes
IMC_Seurat <- RunUMAP(IMC_Seurat, dims = 1:10) #40 minutes


  
# LEVEL 1 - Large Populations (Endocrine, Immune, Acinar, Stromal, Endothelial)
#Convert to features x cells (markers x cells)
  {
Seurat <- CreateSeuratObject(counts = IMC_selection, assay = "IMC", meta.data = IMC_Metadata)

#Run standard Seurat Pipeline
Seurat <- NormalizeData(Seurat, normalization.method = "CLR", margin = 2, verbose = FALSE)
Seurat <- ScaleData(Seurat, features = rownames(Seurat), verbose = FALSE)
Seurat <- RunPCA(Seurat, features = rownames(Seurat), npcs = 30, verbose = FALSE)
ElbowPlot(Seurat, ndims = 30) #Decide how many dimensions to use

Seurat <- FindNeighbors(Seurat, dims = 1:15, verbose = TRUE) #15 minutes
Seurat <- FindClusters(Seurat, resolution = c(0.4)) #45 minutes
Seurat <- RunUMAP(Seurat, dims = 1:15) #40 minutes

# Check - DimPlot for each resolution
resolution = "IMC_snn_res.0.4"
DimPlot(Seurat, group.by = resolution, label = TRUE) + ggtitle("Resolution 0.4")
}
#FindAllMarkers_All Cells
  {
Idents(Seurat) <- "IMC_snn_res.0.4"

markers_seurat <- FindAllMarkers(Seurat,
                              only.pos = TRUE,
                              min.pct = 0.25,
                              logfc.threshold = 0.25)
#Filter our the top 10 markers per cluster
markers_seu_top10 <- markers_seurat %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = 10, with_ties = FALSE) %>% arrange(cluster, desc(avg_log2FC))

  }

#Rename Cluster IDs
  {
new.cluster.ids <- c(
  "Acinar",                      #Cluster 0 
  "Immune",                      #Cluster 1 
  "Acinar",                      #Cluster 2
  "Endocrine",                   #Cluster 3
  "Immune",                      #Cluster 4
  "Immune",                      #Cluster 5
  "Endothelial",                 #Cluster 6
  "Endocrine",                   #Cluster 7
  "Immune",                      #Cluster 8
  "Endocrine",                   #Cluster 9
  "Immune",                      #Cluster 10
  "Endocrine",                   #Cluster 11
  "Immune",                      #Cluster 12
  "Immune",                      #Cluster 13 
  "Endocrine",                   #Cluster 14
  "Endocrine",                   #Cluster 15
  "Endocrine",                   #Cluster 16
  "Endocrine",                   #Cluster 17
  "Endocrine",                   #Cluster 18
  "Endocrine",                   #Cluster 19
  "Immune",                      #Cluster 20
  "Stromal",                     #Cluster 21
  "Endocrine",                   #Cluster 22
  "Immune",                      #Cluster 23
  "Immune",                      #Cluster 24
  "Immune",                      #Cluster 25
  "Immune")                      #Cluster 26


#Idents(Seurat) <- "IMC_snn_res.0.4"
table(Idents(Seurat))
names(new.cluster.ids) <- levels(Seurat)
Seurat <- RenameIdents(Seurat, new.cluster.ids)
table(Idents(Seurat))
DimPlot(Seurat, 
        split.by="Group",
        label = TRUE) + ggtitle("Resolution 0.4")
}

#Subset usable clusters (remove doublets/others)
  {
Seurat <- subset(Seurat, 
                 idents =  c("Acinar","Immune","Endocrine",
                           "Endothelial","Progenitor"))
}


#Supplementary Figiure S5C
#Final Dimplot
  DimPlot(Seurat, 
        reduction = "umap", 
        label = TRUE, 
        repel = TRUE, 
        label.size = 5, 
        pt.size = 0.8, 
        cols = c("#1f77b4","#d62728","#e377c2","#bcbd22","#2ca02c")) +
  #, , "#ff7f0e",,, , )
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  labs(
    title = "Whole Pancreas Clusters (UMAP)",
    #subtitle = "Dimensional reduction colored by Seurat cluster"
  )
  
#### Supplementary Figure S5D ####
#Marker Data_Heatmap----
Markers <- c("Alpha cell","Beta Cell","Cd4","Cd8","Delta Cell","EC","Epithelial Cell", "Ghrelin Cell","Macrophages","Other",
             "PP Cell")
Heatmap <- data.frame(IMC[,c(5:41,47)])
Heatmap <- Heatmap %>% select(-c("cd99","pan.keratin","beta.actin","cd11b"))
Heatmap <- data.frame(t(scale(t(Heatmap[,c(1:33)]))),
                      Celltype=Heatmap$CellType)
Heatmap <- melt(Heatmap,
                idvar = c("Celltype"))

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
           color = colorRampPalette(c("#440154FF", "#2A788EFF", "#FDE725FF"))(50))

#Supplementary Figure S5E
#Create dimplot for whole tissue fractions

#Markers used
#Acinar = "ca2"
#Endocrine = c("c.peptide", "chga", "somatostatin")
#Lymphocyte = c("cd8","cd4")
#Monocyte = c("cd14","cd16")
#Macrophage = c("cd68","cd163")
#Endothelial = c("cd31")
#Progenitor = c("nestin")

Nebulosa::plot_density(Seurat, features = "cd31") +
  labs(title="Endothelial",
       subtitle="CD31") +
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5,size=14,face='bold'),
        plot.subtitle = element_text(hjust=0.5),
        legend.position='none')


######################################################################
#LEVEL 2: subset Immune cells
#Subset Immune cells and run Seurat Workflow

Immune_Seurat <- subset(Seurat, idents = "Immune")
Immune_Seurat <- NormalizeData(Immune_Seurat, normalization.method = "CLR", margin = 2, verbose = FALSE)
Immune_Seurat <- ScaleData(Immune_Seurat, verbose = FALSE)
Immune_Seurat <- RunPCA(Immune_Seurat, features = rownames(Immune_Seurat), npcs = 30, verbose = FALSE)
ElbowPlot(Immune_Seurat, ndims = 30)

Immune_Seurat <- FindNeighbors(Immune_Seurat, dims = 1:15, verbose = TRUE) #15 minutes
Immune_Seurat <- FindClusters(Immune_Seurat, resolution = c(0.4)) #45 minutes
Immune_Seurat <- RunUMAP(Immune_Seurat, dims = 1:15) #40 minutes

resolution <- "IMC_snn_res.0.4"

#Find top markers at resolution =0.4
markers_Immune_Seurat <- FindAllMarkers(Immune_Seurat,
                                only.pos = TRUE,
                                min.pct = 0.25,
                                logfc.threshold = 0.25)
#Filter our the top 5 markers
markers_Immune_seu_top5 <- markers_Immune_Seurat %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = 5, with_ties = FALSE) %>% 
  arrange(cluster, desc(avg_log2FC))

new.cluster.ids <- c(
  "Monocyte",             #Cluster 0
  "Monocyte",             #Cluster 1
  "Macrophage",           #Cluster 2
  "CD8 T",                #Cluster 3
  #"Other",                Cluster 4
  "CD8 T",                #Cluster 5
  "NK",                   #Cluster 6
  "Macrophage",           #Cluster 7
  "CD8 T",                #Cluster 8
  "CD8 T",                #Cluster 9
  "CD8 T",                #Cluster 10
  #"Other",                Cluster 11
  "CD4 T",                #Cluster 12
  "CD8 T",                #Cluster 13
  "Macrophage",           #Cluster 14
  "B Cell",               #Cluster 15
  "CD4 T")                #Cluster 16
  
  
Idents(Immune_Seurat) <- "IMC_snn_res.0.4"
table(Idents(Immune_Seurat))
names(new.cluster.ids) <- levels(Immune_Seurat)
Immune_Seurat <- RenameIdents(Immune_Seurat, new.cluster.ids)


DimPlot(Immune_Seurat, reduction = "umap",    #<---- Use to check for appropriate clusters 
        split.by = "Group",
        label = TRUE, repel = F)


#Subset New Seurat and Save (repeat if necessary
Immune_Seurat <- subset(Immune_Seurat, idents = c("CD8 T",
                                                    "Monocyte",
                                                    "CD4 T",
                                                    "NK",
                                                    "B Cells",
                                                    "Macrophage"))
##### Figure 4A (bottom) #####
#Immune seurat Dimplot
  DimPlot(Immune_Seurat, 
          reduction = "umap", 
          label = F, 
          repel = TRUE, 
          label.size = 5, 
          pt.size = 0.8, 
          cols = c("#2ca02c","#d62728","#e377c2","#1f77b4","cyan","#bcbd22")) +
          #, , "#ff7f0e",,, , )
    theme(
      legend.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5),
      #axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
      plot.margin = margin(10, 10, 10, 10)
    ) +
    labs(
      title = "Immune Cell Clusters (UMAP)",
      #subtitle = "Dimensional reduction colored by Seurat cluster"
    )

#### Supplementary Figure S5G ####

#Markers used
#Monocytes = c("cd14","cd16")
#Macophage = "cd68"
#NK cell = "cd56"
#T cell = c("cd3","cd8,"cd4", "cd20)
ImmuneNebulosa::plot_density(Immune_Seurat, features = "cd3") +
  labs(title="T Cell",
       subtitle="CD3") +
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5,size=14,face='bold'),
        plot.subtitle = element_text(hjust=0.5),
        legend.position='none')
Immune

#Setup Metadata for cell populations
#Figure 4F (summary plot)
Immune_Metadata <- data.frame(Immune_Seurat@meta.data,
                                 Idents = Idents(Immune_Seurat))

Immune_Metadata_Table <- data.frame(table(Cluster = Immune_Metadata$Idents, 
                                          Image = Immune_Metadata$Image))

Immune_Metadata_Table <- data.frame(Immune_Metadata_Table %>% mutate(ID = substr(Image, 1, 8)))
Immune_Metadata_Table$Group <- ifelse(Immune_Metadata_Table$ID %in% c("HPAP-053","HPAP-066","HPAP-069","HPAP-093","HPAP-105","HPAP-118"),
                                      "old", "young")
Immune_Metadata_Summary <- Immune_Metadata_Table %>% group_by(Cluster, ID, Group) %>%
  summarise(mean_cells = mean(Freq))
Immune_Metadata_Summary$Group <- factor(Immune_Metadata_Summary$Group, levels = c("young","old"))


#Tcells
#Figure 4F (top right - CD8 T cells)
CD8T <- Immune_Metadata_Summary %>% filter(Cluster=="CD8 T")

#compute mean ± SEM
mean_sem <- function(x) {
  m <- mean(x, na.rm = TRUE)
  se <- sd(x, na.rm = TRUE) / sqrt(24)  # n = 24
  data.frame(y = m, ymin = m - se, ymax = m + se)
}

ggplot(CD8T, aes(x = Group, y = mean_cells, fill = Group)) +
  stat_summary(fun = mean, geom = "bar", color = "black", width = 0.6, alpha = 0.8) + #Base Layer
  stat_summary(fun.data = mean_sem, geom = "errorbar", width = 0.2, color = "black", linewidth = 0.6) +   #error bars (SEM) ---
  geom_jitter(color = "black", width = 0.15, size = 1, alpha = 0.7) + #--- jittered points ---
  scale_fill_manual(values = c("#E0DCB8", "#607C83")) +
  labs(title = "IMC",y = bquote("CD8+ T cells per mm"^2)) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text = element_text(color = "black"),
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16), 
    axis.title.y = element_text(face = "bold"))

#Statistics
t.test(log1p(CD8T$mean_cells)~CD8T$Group, unequal.var=T)

#CD4 T

CD4T <- Immune_Metadata_Summary %>% filter(Cluster=="CD4 T")
CD4T <- CD4T %>% filter(mean_cells<30)
  
#compute mean ± SEM
  mean_sem <- function(x) {
    m <- mean(x, na.rm = TRUE)
    se <- sd(x, na.rm = TRUE) / sqrt(24)  # n = 24
    data.frame(y = m, ymin = m - se, ymax = m + se)
  }
  
ggplot(CD4T, aes(x = Group, y = mean_cells, fill = Group)) +
  stat_summary(fun = mean, geom = "bar", color = "black", width = 0.6, alpha = 0.8) + #Base Layer
    stat_summary(fun.data = mean_sem, geom = "errorbar", width = 0.2, color = "black", linewidth = 0.6) +   #error bars (SEM) ---
    geom_jitter(color = "black", width = 0.15, size = 2, alpha = 0.7) + #--- jittered points ---
    scale_fill_manual(values = c("#E0DCB8", "#607C83")) +
    labs(title = "IMC",y = bquote("CD4+ T cells per mm"^2)) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.text = element_text(color = "black"),
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 16), 
      axis.title.y = element_text(face = "bold"))

#Statistics
t.test(log1p(CD4T$mean_cells)~CD4T$Group, unequal.var=T)


#Level2 - Subset Endocrine Cells from Original Seurat
Endocrine_Seurat <- subset(
  Seurat,
  idents = c("Endocrine"))

#Run Seurat Workflow
Endocrine_Seurat <- NormalizeData(Endocrine_Seurat, normalization.method = "CLR", margin = 2, verbose = FALSE)
Endocrine_Seurat <- ScaleData(Endocrine_Seurat, verbose = FALSE)
Endocrine_Seurat <- RunPCA(Endocrine_Seurat, features = rownames(Endocrine_Seurat), npcs = 30, verbose = FALSE)
#Decide how many dimensions to use
ElbowPlot(Endocrine_Seurat, ndims = 30)
Endocrine_Seurat <- FindNeighbors(Endocrine_Seurat, dims = 1:14, verbose = TRUE)
Endocrine_Seurat <- FindClusters(Endocrine_Seurat, resolution = c(0.6))
Endocrine_Seurat <- RunUMAP(Endocrine_Seurat, dims = 1:14)

Idents(Endocrine_Seurat) <- "IMC_snn_res.0.6"

DimPlot(Endocrine_Seurat,
        reduction = "umap", label = TRUE, repel = F)  #Check for appropriate clustering

#FindAllMarkers
Endocrine_markers <- FindAllMarkers(Endocrine_Seurat,
                                    only.pos = TRUE,
                                    min.pct = 0.25,
                                    logfc.threshold = 0.25)

Endocrine_markers_top5 <- Endocrine_markers %>% 
  group_by(cluster) %>% 
  slice_max(order_by = avg_log2FC, n = 5, with_ties = FALSE) %>% 
  arrange(cluster, desc(avg_log2FC))

#More Conservative Clustering
{
new.cluster.ids <- c(
  "Other",              #cluster 0
  "Other",              #cluster 1
  "Alpha",               #cluster 2
  "Mix",               #cluster 3
  "Other",              #cluster 4
  "Other",        #cluster 5
  "Other",              #cluster 6
  "Beta",              #cluster 7
  "Other",                 #cluster 8
  "Other",              #cluster 9
  "Other",              #cluster 10
  "Beta",              #cluster 11
  "Other",              #cluster 12
  "PP/Ghrelin",              #cluster 13
  "Other",              #cluster 14
  "Other",            #cluster 15
  "Delta",               #cluster 16
  "Other",              #cluster 17
  "Other",               #cluster 18
  "Other",               #cluster 19
  "Other",               #cluster 20
  "Other",              #cluster 21
  "Other",        #cluster 22
  "Other",              #cluster 23
  "Beta",              #cluster 24
  "Other",                 #cluster 25
  "Other",              #cluster 26
  "PP/Ghrelin",              #cluster 27
  "Other",              #cluster 28
  "Other",              #cluster 29
  "Other",              #cluster 30
  "Other"              #cluster 31
)
}


#Less Conservative Clustering
{
  new.cluster.ids <- c(
    "Beta",              #cluster 0
    "Beta",              #cluster 1
    "Alpha",               #cluster 2
    "Mix",               #cluster 3
    "SST",              #cluster 4
    "Beta",        #cluster 5
    "Alpha",              #cluster 6
    "Beta",              #cluster 7
    "PP/Ghrelin",                 #cluster 8
    "Other",              #cluster 9
    "Alpha",              #cluster 10
    "PP/Ghrelin",              #cluster 11
    "Other",              #cluster 12
    "PP/Ghrelin",              #cluster 13
    "Beta",              #cluster 14
    "Other",            #cluster 15
    "Delta",               #cluster 16
    "PP/Ghrelin",              #cluster 17
    "Alpha",               #cluster 18
    "Other",               #cluster 19
    "Beta",               #cluster 20
    "Other",              #cluster 21
    "Other",        #cluster 22
    "Other",              #cluster 23
    "Beta",              #cluster 24
    "Beta",                 #cluster 25
    "Other",              #cluster 26
    "PP/Ghrelin",              #cluster 27
    "Other",              #cluster 28
    "Other",              #cluster 29
    "Beta",              #cluster 30
    "Other"              #cluster 31
  )
}

#Final Round Clustering
{
new.cluster.ids <- c(
  "Other",              #cluster 0
  "Alpha",              #cluster 1
  "Delta",               #cluster 2
  "Beta",               #cluster 3
  "Beta",              #cluster 4
  "PP/Ghrelin",        #cluster 5
  "Beta",              #cluster 6
  "Beta",              #cluster 7
  "PP/Ghrelin",                 #cluster 8
  "Other",              #cluster 9
  "Alpha",              #cluster 10
  "Endothelial",              #cluster 11
  "PP/Ghrelin",              #cluster 12
  "PP/Ghrelin",              #cluster 13
  "Mix",              #cluster 14
  "Endothelial")            #cluster 15
}
Idents(Endocrine_Seurat) <- "IMC_snn_res.0.6"
table(Idents(Endocrine_Seurat))
names(new.cluster.ids) <- levels(Endocrine_Seurat)
Endocrine_Seurat <- RenameIdents(Endocrine_Seurat, new.cluster.ids)
  
DimPlot(Endocrine_Seurat, reduction = "umap",  #Check clustering.  There are several doublets/mixed cells that need to be filtered out.....
        label = TRUE, repel = F)
  
Endocrine_Seurat <- subset(Endocrine_Seurat,
                           idents=c("Alpha",
                                    "Mix",
                                    "Beta",
                                    "PP/Ghrelin",
                                    "Delta",
                                    "Endothelial"))
#Figure 4A (top)
DimPlot(Endocrine_Seurat, 
        reduction = "umap", 
        label = F, 
        repel = TRUE, 
        label.size = 5, 
        pt.size = 0.8, 
        cols = c("#d62728","#bcbd22","#2ca02c","#1f77b4","#ff7f0e","#9467bd")) +
    theme(
      legend.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5),
      #axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
      plot.margin = margin(10, 10, 10, 10)
    ) +
    labs(
      title = "Endocrine Cell Clusters (UMAP)",
      #subtitle = "Dimensional reduction colored by Seurat cluster"
    )

#Supplementary Figure S5F
##Markers used
#Beta = c("c.peptide", "nkx6.1")
#Alpha = "glucagon"
#Delta = "somatostatin"
#Epsilon = "ghrelin"
#Endocrine = "chga"
#Pancreatic polypeptide = "pp"
#Endothelial + "cd31"

Nebulosa::plot_density(Endocrine_Seurat, features = "nkx6.1") +
  labs(title="Beta",
       subtitle="NKX6-1") +
    #scale_x_continuous(limits = c(-10,2.5)) +
    theme_classic() +
    theme(plot.title = element_text(hjust=0.5,size=14,face='bold'),
          plot.subtitle = element_text(hjust=0.5),
          legend.position='none')

  
  #Save Endocrine MetaData
  Endocrine_Metadata <- data.frame(Endocrine_Seurat@meta.data,
                                   Idents = Idents(Endocrine_Seurat))
  head(Endocrine_Seurat)
  write.csv(Endocrine_Metadata, "D:/Mike/IMC/3_Endocrine_Cells/Endocrine_Markers_Metadata_LessConservative.csv")
  #Immune Clusters Per Donor
  write.csv(t(table(Endocrine_Metadata$Clusternames, Endocrine_Metadata$DonorID)),
            file="3_Endocrine_Cells/Endocrine_Cell_Proportions_xSubject_LessConservative.csv")
  
#Bind Dataframe
  
FinalDF <- rbind.data.frame(Immune_Metadata,
                            Endocrine_Metadata %>% select(-c("IMC_snn_res.0.2","IMC_snn_res.0.6")))

write.csv(FinalDF, file = "D:/Mike/IMC/IMC_Markers_Final_13Oct2025.csv")


################## Figure 4K ##################
#Compare Alpha and Beta Cell HLA-A Expression#

#Filter Alpha/Beta Cells from IMC_Selection
#Figure 4D - Alpha Cell HLA-ABC, HLA-DR

Alpha <- IMC %>% filter(CellType=="Alpha cell")
#Clean Outliers by 1.5-fold IQR methodd 
hla.abc_clean <- Alpha %>% filter(between(hla.abc,
                                            quantile(hla.abc, 0.25, na.rm=T) - 1.5 * IQR(hla.abc, na.rm=T),
                                            quantile(hla.abc, 0.75, na.rm=T) + 1.5 * IQR(hla.abc, na.rm=T)))
  
ggplot(hla.abc_clean, aes(x=group, y=hla.abc,fill=group)) +
  geom_violin(color='black') +
  # Add Median
  stat_summary(fun = median, geom = "crossbar", 
                 width = 0.3, color = "black", size = 0.5) +
    
    # Add Q1 and Q3 as horizontal lines
    stat_summary(fun.data = function(x) {
      q <- quantile(x, probs = c(0.25, 0.75))
      data.frame(y = q, ymin = q[1], ymax = q[2])
    }, geom = "errorbar", width = 0.2, color = "black", size = 0.5) +
    scale_fill_manual(values=c("#E0DCB8","#607C83")) +
    labs(title="HLA-ABC",
         y="Intensity") +
    theme_classic() +
    theme(legend.position = 'none',
          axis.title.x = element_blank(),
          plot.title = element_text(hjust=0.5))
  hla.abc_plot
  
  t.test(log1p(hla.abc_clean$hla.abc)~hla.abc_clean$group, var.equal = FALSE)

#HLA-DR
hla.dr_clean <- Alpha %>% filter(between(hla.dr,
                                           quantile(hla.dr, 0.25, na.rm=T) - 1.5 * IQR(hla.dr, na.rm=T),
                                           quantile(hla.dr, 0.75, na.rm=T) + 1.5 * IQR(hla.dr, na.rm=T)))
  
  hla.dr_plot <- ggplot(hla.dr_clean, aes(x=group, y=hla.dr,fill=group)) +
    geom_violin(color='black') +
    # Add Median
    stat_summary(fun = median, geom = "crossbar", 
                 width = 0.3, color = "black", size = 0.5) +
    
    # Add Q1 and Q3 as horizontal lines
    stat_summary(fun.data = function(x) {
      q <- quantile(x, probs = c(0.25, 0.75))
      data.frame(y = q, ymin = q[1], ymax = q[2])
    }, geom = "errorbar", width = 0.2, color = "black", size = 0.5) +
    scale_fill_manual(values=c("#E0DCB8","#607C83")) +
    labs(title="HLA-DR",
         y="expression") +
    theme_classic() +
    theme(legend.position = 'none',
          plot.title = element_text(hjust=0.5),
          axis.title = element_blank())
  hla.dr_plot
  
  t.test(hla.dr_clean$hla.dr~hla.dr_clean$group, var.equal=F)

#Figure 4L - Beta Cell HLA-ABC, HLA-DR
Beta <- IMC_Master %>% filter(CellType=="Beta Cell")
hla.abc_clean <- Beta %>% filter(between(hla.abc,
                                          quantile(hla.abc, 0.25, na.rm=T) - 1.5 * IQR(hla.abc, na.rm=T),
                                          quantile(hla.abc, 0.75, na.rm=T) + 1.5 * IQR(hla.abc, na.rm=T)))

ggplot(hla.abc_clean, aes(x=group, y=hla.abc,fill=group)) +
  geom_violin(color='black') +
  # Add Median
  stat_summary(fun = median, geom = "crossbar", 
               width = 0.3, color = "black", size = 0.5) +
      # Add Q1 and Q3 as horizontal lines
  stat_summary(fun.data = function(x) {
    q <- quantile(x, probs = c(0.25, 0.75))
    data.frame(y = q, ymin = q[1], ymax = q[2])
  }, geom = "errorbar", width = 0.2, color = "black", size = 0.5) +
  scale_fill_manual(values=c("#E0DCB8","#607C83")) +
  labs(title="HLA-ABC",
       y="Intensity") +
  scale_y_continuous(limits = c(0,50)) +
  theme_classic() +
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        plot.title = element_text(hjust=0.5))

t.test(log1p(hla.abc_clean$hla.abc)~hla.abc_clean$group, var.equal = FALSE)

#HLA-DR
hla.dr_clean <- Beta %>% filter(between(hla.dr,
                                         quantile(hla.dr, 0.25, na.rm=T) - 1.5 * IQR(hla.dr, na.rm=T),
                                         quantile(hla.dr, 0.75, na.rm=T) + 1.5 * IQR(hla.dr, na.rm=T)))

ggplot(hla.dr_clean, aes(x=group, y=hla.dr,fill=group)) +
  geom_violin(color='black') +
  # Add Median
  stat_summary(fun = median, geom = "crossbar", 
               width = 0.3, color = "black", size = 0.5) +
  scale_y_continuous(limits = c(0,6)) +
  # Add Q1 and Q3 as horizontal lines
  stat_summary(fun.data = function(x) {
    q <- quantile(x, probs = c(0.25, 0.75))
    data.frame(y = q, ymin = q[1], ymax = q[2])
  }, geom = "errorbar", width = 0.2, color = "black", size = 0.5) +
  scale_fill_manual(values=c("#E0DCB8","#607C83")) +
  labs(title="HLA-DR",
       y="expression") +
  theme_classic() +
  theme(legend.position = 'none',
        plot.title = element_text(hjust=0.5),
        axis.title = element_blank())
  hla.dr_plot

# END 
###############################################################################
######################### Data not shown in manuscript#########################
  
#Nebulso plot subset by age for endocrine marker expression
Young_seurat <- subset(Endocrine_Seurat, subset = Group == "young")
Old_seurat <- subset(Endocrine_Seurat, subset = Group == "old")

{
{
p_background <- Nebulosa::plot_density(Endocrine_Seurat, features = "hla.abc") +
  scale_color_gradient(low = "#FFF9C4", high = "#FFB300") +
  scale_fill_gradient(low = "#FFF9C4", high = "#FFB300") +
  scale_x_continuous(limits=c(-12,13))+
  scale_y_continuous(limits=c(-15,8))+
  theme_void() +
  labs(title="HLA-ABC",
       y="Young") +
  theme(plot.title = element_text(hjust=0.5),
        axis.line = element_line(linewidth = 1),
        axis.ticks = element_blank(),
        #axis.text = element_blank(),
        axis.title.x =  element_blank(),
        legend.position='none')
p_background
}
{

p1 <- Nebulosa::plot_density(Young_seurat, features = "hla.abc") +
    scale_fill_gradient(low = alpha("#FFF9C4", 0.3), high = alpha("red", 0.05)) +
    scale_color_gradient(low = alpha("#FFF9C4", 0.3), high = alpha("red", 0.05)) +
    scale_x_continuous(limits=c(-12,13))+
    scale_y_continuous(limits=c(-15,8))+
    labs(title="HLA-ABC",
         y="Young") +
    theme_void() +
    theme(plot.title = element_text(hjust=0.5),
          axis.line = element_line(linewidth = 1),
          axis.ticks = element_blank(),
          #axis.text = element_blank(),
          axis.title.x =  element_blank(),
          legend.position='none',
          panel.background = element_rect(fill = NA, color = NA) )
  p1
  
  p1 <- ggdraw() +
    draw_plot(p_background) +
    draw_plot(p1)  # defaults to same coordinates/size
  p1
  } 
}#p1
{
  {
    p_background <- Nebulosa::plot_density(Endocrine_seu, features = "hla.dr") +
      scale_color_gradient(low = "#FFF9C4", high = "#FFB300") +
      scale_fill_gradient(low = "#FFF9C4", high = "#FFB300") +
      scale_x_continuous(limits=c(-12,13))+
      scale_y_continuous(limits=c(-15,8))+
      theme_void() +
      labs(title="HLA-DR") +
      theme(plot.title = element_text(hjust=0.5),
            axis.line = element_line(linewidth = 1),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            #axis.text = element_blank(),
            axis.title.x =  element_blank(),
            legend.position='none')
    p_background
  }
  {
    p2 <- Nebulosa::plot_density(Young_seurat, features = "hla.dr") +
      scale_fill_gradient(low = alpha("#FFF9C4", 0.3), high = alpha("red", 0.05)) +
      scale_color_gradient(low = alpha("#FFF9C4", 0.3), high = alpha("red", 0.05)) +
      scale_x_continuous(limits=c(-12,13))+
      scale_y_continuous(limits=c(-15,8))+
      labs(title="HLA-DR") +
      theme_void() +
      theme(plot.title = element_text(hjust=0.5),
            axis.line = element_line(linewidth = 1),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            #axis.text = element_blank(),
            axis.title.x =  element_blank(),
            legend.position='none',
            panel.background = element_rect(fill = NA, color = NA) )
    p2
    
    p2 <- ggdraw() +
      draw_plot(p_background) +
      draw_plot(p2)  # defaults to same coordinates/size
    p2
  } 
}#p2
{
  {
    p_background <- Nebulosa::plot_density(Endocrine_seu, features = "hla.abc") +
      scale_color_gradient(low = "#FFF9C4", high = "#FDD835") +
      scale_fill_gradient(low = "#FFF9C4", high = "#FDD835") +
      scale_x_continuous(limits=c(-12,13))+
      scale_y_continuous(limits=c(-15,8))+
      theme_void() +
      labs(y="Old") +
      theme(plot.title = element_blank(),
            axis.line = element_line(linewidth = 1),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            #axis.text = element_blank(),
            axis.title.x =  element_blank(),
            legend.position='none')
    p_background
  }
  {
    p3 <- Nebulosa::plot_density(Old_seurat, features = "hla.abc") +
      scale_fill_gradient(low = alpha("#FFF9C4", 0.2), high = alpha("red", 0.5)) +
      scale_color_gradient(low = alpha("#FFF9C4", 0.2), high = alpha("red", 0.5)) +
      scale_x_continuous(limits=c(-12,13))+
      scale_y_continuous(limits=c(-15,8))+
      labs(y="Old") +
      theme_void() +
      theme(plot.title = element_blank(),
            axis.line = element_line(linewidth = 1),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            #axis.text = element_blank(),
            axis.title.x =  element_blank(),
            legend.position='none')
    p3
    
    p3.1 <- Nebulosa::plot_density(Young_seurat, features = "hla.abc") +
      scale_fill_gradient(low = alpha("#FFF9C4", 0.2), high = alpha("red", 0.01)) +
      scale_color_gradient(low = alpha("#FFF9C4", 0.2), high = alpha("red", 0.01)) +
      scale_x_continuous(limits=c(-12,13))+
      scale_y_continuous(limits=c(-15,8))+
      labs(y="Old") +
      theme_void() +
      theme(plot.title = element_blank(),
            axis.line = element_line(linewidth = 1),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            #axis.text = element_blank(),
            axis.title.x =  element_blank(),
            legend.position='none')
    p3.1
    
    p3 <- ggdraw() +
      draw_plot(p_background) +
      draw_plot(p3.1) +  # defaults to same coordinates/size
      draw_plot(p3)  # defaults to same coordinates/size
    p3
  } 
}#p3
{
  {
    p_background <- Nebulosa::plot_density(Endocrine_seu, features = "hla.dr") +
      scale_color_gradient(low = "#FFF9C4", high = "#FDD835") +
      scale_fill_gradient(low = "#FFF9C4", high = "#FDD835") +
      scale_x_continuous(limits=c(-12,13))+
      scale_y_continuous(limits=c(-15,8))+
      theme_void() +
      labs(y="Old") +
      theme(plot.title = element_blank(),
            axis.line = element_line(linewidth = 1),
            axis.title.y = element_blank(),
            axis.ticks = element_blank(),
            #axis.text = element_blank(),
            axis.title.x =  element_blank(),
            legend.position='none')
    p_background
  }
  {
    p4 <- Nebulosa::plot_density(Old_seurat, features = "hla.dr") +
      scale_fill_gradient(low = alpha("#FFF9C4", 0.2), high = alpha("red", 0.1)) +
      scale_color_gradient(low = alpha("#FFF9C4", 0.2), high = alpha("red", 0.1)) +
      scale_x_continuous(limits=c(-12,13))+
      scale_y_continuous(limits=c(-15,8))+
      labs(y="Old") +
      theme_void() +
      theme(plot.title = element_blank(),
            axis.title.y = element_blank(),
            axis.line = element_line(linewidth = 1),
            axis.ticks = element_blank(),
            #axis.text = element_blank(),
            axis.title.x =  element_blank(),
            legend.position='none',
            panel.background = element_rect(fill = NA, color = NA) )
    p4
    
    
    p4.1 <- Nebulosa::plot_density(Young_seurat, features = "hla.dr") +
      scale_fill_gradient(low = alpha("#FFF9C4", 0.3), high = alpha("red", 0.001)) +
      scale_color_gradient(low = alpha("#FFF9C4", 0.3), high = alpha("red", 0.001)) +
      scale_x_continuous(limits=c(-12,13))+
      scale_y_continuous(limits=c(-15,8))+
      labs(title="HLA-DR") +
      theme_void() +
      theme(plot.title = element_blank(),
            axis.title.y = element_blank(),
            axis.line = element_line(linewidth = 1),
            axis.ticks = element_blank(),
            #axis.text = element_blank(),
            axis.title.x =  element_blank(),
            legend.position='none',
            panel.background = element_rect(fill = NA, color = NA) )
    p4.1
    
    p4 <- ggdraw() +
      draw_plot(p_background) +
      draw_plot(p4) +
      draw_plot(p4.1)
      # defaults to same coordinates/size
    p4
  } 
}#p4
(p1+p2)/(p3+p4)







#### Beta Cell marker expression 
#### Not shown in manuscript
## iapp
iapp_clean <- Beta %>% filter(between(iapp,
                                      quantile(iapp, 0.25, na.rm=T) - 1.5 * IQR(iapp, na.rm=T),
                                      quantile(iapp, 0.75, na.rm=T) + 1.5 * IQR(iapp, na.rm=T)))

iapp_beta <- ggplot(iapp_clean, aes(x=group, y=log1p(iapp),fill=group)) +
  geom_violin(color='black', adjust =.5) +
  # Add Median
  stat_summary(fun = median, geom = "crossbar", 
               width = 0.3, color = "black", size = 0.5) +
  
  # Add Q1 and Q3 as horizontal lines
  stat_summary(fun.data = function(x) {
    q <- quantile(x, probs = c(0.25, 0.75))
    data.frame(y = q, ymin = q[1], ymax = q[2])
  }, geom = "errorbar", width = 0.2, color = "black", size = 0.5) +
  scale_fill_manual(values=c("#E0DCB8","#607C83")) +
  theme_classic() +
  labs(y="log(intensity) + 1")+
  theme(legend.position = 'none',
        axis.title.x = element_blank())
iapp_beta
t.test(log1p(iapp_clean$iapp)~iapp_clean$group, unequal.var=T)

#PDX-1
pdx.1_clean <- Beta %>% filter(between(pdx.1,
                                       quantile(pdx.1, 0.25, na.rm=T) - 1.5 * IQR(pdx.1, na.rm=T),
                                       quantile(pdx.1, 0.75, na.rm=T) + 1.5 * IQR(pdx.1, na.rm=T)))

pdx.1_plot <- ggplot(pdx.1_clean, aes(x=group, y=log1p(pdx.1),fill=group)) +
  geom_violin(color='black',adjust=.5,width=0.8) +
  
  # Add Median
  stat_summary(fun = median, geom = "crossbar", 
               width = 0.3, color = "black", size = 0.5) +
  
  # Add Q1 and Q3 as horizontal lines
  stat_summary(fun.data = function(x) {
    q <- quantile(x, probs = c(0.25, 0.75))
    data.frame(y = q, ymin = q[1], ymax = q[2])
  }, geom = "errorbar", width = 0.2, color = "black", size = 0.5) +
  
  
  scale_fill_manual(values=c("#E0DCB8","#607C83")) +
  labs(title="PDX-1",
       y="log(intensity) + 1") +
  theme_classic() +
  theme(legend.position = 'none',
        plot.title = element_text(hjust=0.5),
        axis.title.x = element_blank())
pdx.1_plot

t.test(log1p(pdx.1_clean$pdx.1)~pdx.1_clean$group)

#NKX6.1
nkx6.1_clean <- Beta %>% filter(between(nkx6.1,
                                        quantile(nkx6.1, 0.25, na.rm=T) - 1.5 * IQR(nkx6.1, na.rm=T),
                                        quantile(nkx6.1, 0.75, na.rm=T) + 1.5 * IQR(nkx6.1, na.rm=T)))

nkx6.1_plot <- ggplot(nkx6.1_clean, aes(x=group, y=log1p(nkx6.1),fill=group)) +
  geom_violin(color='black',adjust=.5,width=0.8) +
  # Add Median
  stat_summary(fun = median, geom = "crossbar", 
               width = 0.3, color = "black", size = 0.5) +
  
  # Add Q1 and Q3 as horizontal lines
  stat_summary(fun.data = function(x) {
    q <- quantile(x, probs = c(0.25, 0.75))
    data.frame(y = q, ymin = q[1], ymax = q[2])
  }, geom = "errorbar", width = 0.2, color = "black", size = 0.5) +
  scale_fill_manual(values=c("#E0DCB8","#607C83")) +
  labs(title="NKX6-1",
       y="log(intensity) + 1") +
  theme_classic() +
  theme(legend.position = 'none',
        plot.title = element_text(hjust=0.5),
        axis.title.x = element_blank())
nkx6.1_plot

t.test(log1p(nkx6.1_clean$nkx6.1)~nkx6.1_clean$group)

Indentity_Beta <- pdx.1_plot + nkx6.1_plot
Indentity_Beta

#### Endocrine Cell expression gating
#Filtering by Cell type (CD45->CD3->CD4/8)
CD45 <- IMC %>% filter(cd45>5)
#T Cells
CD3 <- CD45 %>% filter(cd3>2)
CD3 <- CD3 %>% filter(CellType %in% c("Cd4","Cd8"))



CD8 <- CD3 %>% filter(cd8>5)
CD4 <- CD3 %>% filter(cd4>3)

#Create CD8 Proportion relative to CD3+
#Workflow used to generate Figure 4F 
CD3.table <- data.frame(table(CD3$DonorID))
names(CD3.table)[2] <- "CD3" 

CD8.table <- data.frame(table(CD8$DonorID))
names(CD8.table)[2] <- "CD8"
CD8.table <- right_join(CD8.table, CD3.table, by="Var1")
CD8.table$proportion <- (CD8.table$CD8/CD8.table$CD3)*100
Old <- c("HPAP-118","HPAP-105","HPAP-093","HPAP-069","HPAP-066","HPAP-053")
CD8.table <- data.frame(CD8.table,
                        classification=ifelse(CD8.table$Var1 %in% Old,"Old","Young"))
CD8.table$classification <- factor(CD8.table$classification, levels = c("Young","Old"))

t.test(log(CD8.table$proportion)~CD8.table$classification)

#Create CD4 Proportion relative to CD3+
CD4.table <- data.frame(table(CD4$DonorID))
names(CD4.table)[2] <- "CD4"
CD4.table <- right_join(CD4.table, CD3.table, by="Var1")
CD4.table$proportion <- (CD4.table$CD4/CD4.table$CD3)*100
Old <- c("HPAP-118","HPAP-105","HPAP-093","HPAP-069","HPAP-066","HPAP-053")
CD4.table <- data.frame(CD4.table,
                        classification=ifelse(CD4.table$Var1 %in% Old,"Old","Young"))
CD4.table$classification <- factor(CD4.table$classification, levels = c("Young","Old"))
#Remove HPAP-027 Outlier
CD4.table <- CD4.table %>% filter(Var1!="HPAP-027")
  
t.test(log(CD4.table$proportion)~CD4.table$classification)

#Gating Strategy----
#### Not shown in manscuript
ggplot(IMC, aes(x=cd45)) +
  geom_density(fill = "darkgray", adjust=0.001) +
  scale_x_log10(limits=c(0.1,100)) +
  labs(title = "CD45",
       x = "CD45 intensity",
       y = "density of all cells") +
  #xlim(0.01,4)+
  geom_vline(xintercept =2,colour = 'red', linetype='dashed')+
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5))

ggplot(CD45, aes(x=cd3)) +
  geom_density(fill = "skyblue", adjust=0.001) +
  scale_x_log10() +
  labs(title = "CD3",
       x = "CD3 intensity",
       y = "density of CD45+ cells") +
  #xlim(0.01,4)+
  geom_vline(xintercept =2.5,colour = 'red', linetype='dashed')+
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5))

ggplot(CD3, aes(x=cd4)) +
  geom_density(fill = "skyblue", adjust=0.00001) +
  scale_x_log10(limits=c(0.1,100)) +
  labs(title = "CD4",
       x = "CD4 intensity",
       y = "density of CD3+ cells") +
  #xlim(0.01,4)+
  geom_vline(xintercept =3,colour = 'red', linetype='dashed')+
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5))

ggplot(CD3, aes(x=cd8)) +
  geom_density(fill = "skyblue", adjust=0.00001) +
  scale_x_log10(limits=c(0.1,100)) +
  labs(title = "CD8",
       x = "CD8 intensity",
       y = "density of CD3+ cells") +
  #xlim(0.01,4)+
  geom_vline(xintercept =3,colour = 'red', linetype='dashed')+
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5))

#Macrophage gating strategy
CD68 <- CD45 %>% filter(cd68>10)
CD163_hi <- CD68 %>% filter(cd163>10)
CD163_lo <- CD68 %>% filter(cd163>0.5 & cd163<5)

ggplot(CD45, aes(x=cd68)) +
  geom_density(fill = "orange", adjust=0.001) +
  scale_x_log10() +
  labs(title = "CD68",
       x = "CD68 intensity",
       y = "density of CD45+ cells") +
  #xlim(0.01,4)+
  geom_vline(xintercept =10,colour = 'red', linetype='dashed')+
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5))
CD68 <- CD68 %>% filter(cd163>0.5)

ggplot(CD68, aes(x=cd163)) +
  geom_density(fill = "orange", adjust=0.001) +
  scale_x_log10(limits=c(.5,300)) +
  labs(title = "CD163",
       x = "CD163 intensity",
       y = "density of CD68+ cells") +
  #xlim(0.01,4)+
  geom_vline(xintercept =c(5,10),colour = 'red', linetype='dashed')+
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5))

#Proportion CD163_lo cells----
CD68.table <- data.frame(table(CD68$DonorID))
names(CD68.table)[2] <- "CD68" 

CD163_lo.table <- data.frame(table(CD163_lo$DonorID))
names(CD163_lo.table)[2] <- "CD163"
CD163_lo.table <- right_join(CD163_lo.table, CD68.table, by="Var1")
CD163_lo.table$proportion <- (CD163_lo.table$CD163/CD163_lo.table$CD68)*100
Old <- c("HPAP-118","HPAP-105","HPAP-093","HPAP-069","HPAP-066","HPAP-053")
CD163_lo.table <- data.frame(CD163_lo.table,
                        classification=ifelse(CD163_lo.table$Var1 %in% Old,"Old","Young"))
CD163_lo.table$classification <- factor(CD163_lo.table$classification, levels = c("Young","Old"))
CD163_lo.table <- CD163_lo.table %>% filter(proportion>0)

t.test(log1p(CD163_lo.table$proportion)~CD163_lo.table$classification)

CD163lo_bar <- CD163_lo.table %>% 
  group_by(classification) %>% 
  summarise(mean=mean(proportion),
            SE=sd(proportion)/4.8)

ggplot() +
  geom_col(data=CD163lo_bar, aes(x=classification,y=mean,fill=classification),color='black') +
  geom_errorbar(data=CD163lo_bar,mapping = aes(x=classification, ymin = mean-SE, ymax = mean+SE),width=0.5) +
  geom_jitter(data=CD163_lo.table, aes(x=classification,y=proportion),width=0.2,alpha=0.5) +
  scale_fill_manual(values=c("#E0DCB8","#607C83"))+
  labs(y="%CD163Lo of CD68+ Cells")+
  theme_classic() +
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=12))

#CD163HI Cells
CD163_hi.table <- data.frame(table(CD163_hi$DonorID))
names(CD163_hi.table)[2] <- "CD163"
CD163_hi.table <- right_join(CD163_hi.table, CD68.table, by="Var1")
CD163_hi.table$proportion <- (CD163_hi.table$CD163/CD163_hi.table$CD68)*100
Old <- c("HPAP-118","HPAP-105","HPAP-093","HPAP-069","HPAP-066","HPAP-053")
CD163_hi.table <- data.frame(CD163_hi.table,
                             classification=ifelse(CD163_hi.table$Var1 %in% Old,"Old","Young"))
CD163_hi.table$classification <- factor(CD163_hi.table$classification, levels = c("Young","Old"))
CD163_hi.table <- CD163_hi.table %>% filter(proportion>0)


t.test(log1p(CD163_hi.table$proportion)~CD163_hi.table$classification)

CD163hi_bar <- CD163_hi.table %>% 
  group_by(classification) %>% 
  summarise(mean=mean(proportion),
            SE=sd(proportion)/4.8)

ggplot() +
  geom_col(data=CD163hi_bar, aes(x=classification,y=mean,fill=classification),color='black') +
  geom_errorbar(data=CD163hi_bar,mapping = aes(x=classification, ymin = mean-SE, ymax = mean+SE),width=0.5) +
  geom_jitter(data=CD163_hi.table, aes(x=classification,y=proportion),width=0.2,alpha=0.5) +
  scale_fill_manual(values=c("#E0DCB8","#607C83"))+
  labs(y="%CD163HI of CD68+ Cells")+
  theme_classic() +
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=12))

#Ratio
CD163ratio <- data.frame((CD163_hi.table %>% select(Var1,
                                                   CD163hi=CD163)),
                         CD163_lo.table %>% select(CD163lo=CD163))
CD163ratio$proportion <- (CD163ratio$CD163hi/CD163ratio$CD163lo)
CD163ratio <- data.frame(CD163ratio,
                             classification=ifelse(CD163ratio$Var1 %in% Old,"Old","Young"))

t.test(CD163ratio$proportion~CD163ratio$classification)


#Classify M1/M2 ratio
CD163ratio_bar <- CD163ratio %>% 
  group_by(classification) %>% 
  summarise(mean=mean(proportion),
            SE=sd(proportion)/4.8)

ggplot() +
  geom_col(data=CD163ratio_bar, aes(x=classification,y=mean,fill=classification),color='black') +
  geom_errorbar(data=CD163ratio_bar,mapping = aes(x=classification, ymin = mean-SE, ymax = mean+SE),width=0.5) +
  geom_jitter(data=CD163ratio, aes(x=classification,y=proportion),width=0.2,alpha=0.5) +
  scale_fill_manual(values=c("#E0DCB8","#607C83"))+
  labs(y="M1/M2 ratio")+
  theme_classic() +
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=12))




#Figure 1F -> Classify percent of CD8 % of CD3 cells
CD8_bar <- CD8.table %>% 
  group_by(classification) %>% 
  summarise(mean=mean(proportion),
            SE=sd(proportion)/4.8)

ggplot() +
  geom_col(data=CD8_bar, aes(x=classification,y=mean,fill=classification),color='black') +
  geom_errorbar(data=CD8_bar,mapping = aes(x=classification, ymin = mean-SE, ymax = mean+SE),width=0.5) +
  geom_jitter(data=CD8.table, aes(x=classification,y=proportion),width=0.2,alpha=0.5) +
  scale_fill_manual(values=c("#E0DCB8","#607C83"))+
  labs(y="%CD8 of CD3+ Cells")+
  theme_classic() +
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=12))


#Figure 1F Classify percent of CD4% of CD3 cells
CD4_bar <- CD4.table %>% 
  group_by(classification) %>% 
  summarise(mean=mean(proportion),
            SE=sd(proportion)/4.8)

ggplot() +
  geom_col(data=CD4_bar, aes(x=classification,y=mean,fill=classification),color='black') +
  geom_errorbar(data=CD4_bar,mapping = aes(x=classification, ymin = mean-SE, ymax = mean+SE),width=0.5) +
  geom_jitter(data=CD4.table, aes(x=classification,y=proportion),width=0.2,alpha=0.5) +
  scale_fill_manual(values=c("#E0DCB8","#607C83"))+
  labs(y="%CD4 of CD3+ Cells")+
  theme_classic() +
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=12))

#Classify CD8/CD4 ratio
#Create CD8 to CD4 Ratio
CD8_CD4_Ratio <- data.frame(right_join(CD8.table[1:2], CD4.table[1:2], by="Var1"))
CD8_CD4_Ratio$ratio <- (CD8_CD4_Ratio$CD8/CD8_CD4_Ratio$CD4)
CD8_CD4_Ratio$classification <- ifelse(CD8_CD4_Ratio$Var1 %in% Old,"Old","Young")
t.test(CD8_CD4_Ratio$ratio~CD8_CD4_Ratio$classification)

CD8_CD4_Ratio_bar <- CD8_CD4_Ratio %>% 
  group_by(classification) %>% 
  summarise(mean=mean(ratio),
            SE=sd(ratio)/4.8)

ggplot() +
  geom_col(data=CD8_CD4_Ratio_bar, aes(x=classification,y=mean,fill=classification),color='black') +
  geom_errorbar(data=CD8_CD4_Ratio_bar,mapping = aes(x=classification, ymin = mean-SE, ymax = mean+SE),width=0.5) +
  geom_jitter(data=CD8_CD4_Ratio, aes(x=classification,y=ratio),width=0.2,alpha=0.5) +
  scale_fill_manual(values=c("#E0DCB8","#607C83"))+
  labs(y="%CD4 of CD3+ Cells")+
  theme_classic() +
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=12))

#Classify Memory T Cells
ggplot(CD8, aes(x=cd45ro, group=group, fill=group)) +
  geom_density(adjust=0.1,alpha=0.1) +
  scale_x_log10() +
  labs(title = "CD8",
       x = "CD8 intensity",
       y = "density of CD45+ cells") +
  #xlim(0.01,4)+
  geom_vline(xintercept =5,colour = 'red', linetype='dashed')+
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5))

CD45RO <- CD8 %>% filter(cd45ro>3)


CD45RO <- data.frame(table(CD45RO$DonorID),
                     table(CD8$DonorID))
CD45RO <- CD45RO[-3]
CD45RO$proportion <- CD45RO$Freq/CD45RO$Freq.1*100
CD45RO <- data.frame(CD45RO,
                        classification=ifelse(CD45RO$Var1 %in% Old,"Old","Young"))
CD45RO$classification <- factor(CD45RO$classification, levels = c("Young","Old"))
#CD45RO_Porportion <- CD45RO_Porportion %>% filter(Var1!='HPAP-080')


CD45RO_bar <- CD45RO %>% 
  group_by(classification) %>% 
  summarise(mean=mean(proportion, na.rm=T),
            SE=sd(proportion,na.rm=T)/4.8)


ggplot() +
  geom_col(data=CD45RO_bar, aes(x=classification,y=mean,fill=classification),color='black') +
  geom_errorbar(data=CD45RO_bar,mapping = aes(x=classification, ymin = (mean-SE), ymax = (mean+SE)),width=0.5) +
  geom_jitter(data=CD45RO, aes(x=classification,y=proportion),width=0.2,alpha=0.5) +
  scale_fill_manual(values=c("#E0DCB8","#607C83"))+
  labs(y="% CD45RO of CD8+ T Cells")+
  theme_classic() +
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=12))
t.test(CD45RO$proportion~CD45RO$classification)



#Classify Activated CD44+ T cells
ggplot(CD8, aes(x=cd44, group=group, fill=group)) +
  geom_density(adjust=0.01,alpha=0.1) +
  scale_x_log10() +
  labs(title = "CD8",
       x = "CD8 intensity",
       y = "density of CD45+ cells") +
  #xlim(0.01,4)+
  geom_vline(xintercept =3,colour = 'red', linetype='dashed')+
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5))

#CD44-pos
CD44 <- CD8 %>% filter(cd44>3)
CD44_Low <- CD8 %>% filter(cd44<1)

CD44_Porportion <- data.frame(table(CD44$DonorID),
                                table(CD8$DonorID))
CD44_Porportion <- CD44_Porportion[-3]



CD44_Porportion <- data.frame(CD44_Porportion,
                              proportion=(CD44_Porportion$Freq/CD44_Porportion$Freq.1)*100)
  Old <- c("HPAP-118","HPAP-105","HPAP-093","HPAP-069","HPAP-066","HPAP-053")

CD44_Porportion <- data.frame(CD44_Porportion,
                                classification=ifelse(CD44_Porportion$Var1 %in% Old,"Old","Young"))
CD44_Porportion$classification <- factor(CD44_Porportion$classification, levels = c("Young","Old"))
#Remove Outlier
CD44_Porportion <- CD44_Porportion %>% filter(Var1!="HPAP-027")



CD44_bar <- CD44_Porportion %>% 
  group_by(classification) %>% 
  summarise(mean=mean(proportion, na.rm=T),
            SE=sd(proportion,na.rm=T)/4.8)


ggplot() +
  geom_col(data=CD44_bar, aes(x=classification,y=mean,fill=classification),color='black') +
  geom_errorbar(data=CD44_bar,mapping = aes(x=classification, ymin = (mean-SE), ymax = (mean+SE)),width=0.5) +
  geom_jitter(data=CD44_Porportion, aes(x=classification,y=proportion),width=0.2,alpha=0.5) +
  scale_fill_manual(values=c("#E0DCB8","#607C83"))+
  labs(y="% CD44 of CD8+ T Cells")+
  theme_classic() +
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=12))
t.test(CD44_Porportion$proportion~CD44_Porportion$classification)

#GranzymeB Expression on CD44+ Cells
CD44$group <- ifelse(CD44$DonorID %in% Old,"Young","Old")
CD44$group <- factor(CD44$group, levels = c("Young","Old"))

CD44_clean <- CD44 %>% filter(between(granzymeb,
                                         quantile(granzymeb, 0.25, na.rm=T) - 1.5 * IQR(granzymeb, na.rm=T),
                                         quantile(granzymeb, 0.75, na.rm=T) + 1.5 * IQR(granzymeb, na.rm=T)))

CD44_plot <- ggplot(CD44_clean, aes(x=group, y=log1p(granzymeb),fill=group)) +
  geom_violin(color='black') +
  # Add Median
  stat_summary(fun = median, geom = "crossbar", 
               width = 0.3, color = "black", size = 0.5) +
  
  # Add Q1 and Q3 as horizontal lines
  stat_summary(fun.data = function(x) {
    q <- quantile(x, probs = c(0.25, 0.75))
    data.frame(y = q, ymin = q[1], ymax = q[2])
  }, geom = "errorbar", width = 0.2, color = "black", size = 0.5) +
  scale_fill_manual(values=c("#E0DCB8","#607C83")) +
  labs(title="granzymeb",
       y="CD44Hi gzb expression") +
  theme_classic() +
  theme(legend.position = 'none',
        plot.title = element_text(hjust=0.5))
        #axis.title = element_blank())
CD44_plot

t.test(log1p(CD44_clean$granzymeb)~CD44_clean$group)

#Granzyme B on CD44LO Cells
CD44 <- CD44_Low
CD44$group <- ifelse(CD44$DonorID %in% Old,"Young","Old")
CD44$group <- factor(CD44$group, levels = c("Young","Old"))

CD44_clean <- CD44 %>% filter(between(granzymeb,
                                      quantile(granzymeb, 0.25, na.rm=T) - 1.5 * IQR(granzymeb, na.rm=T),
                                      quantile(granzymeb, 0.75, na.rm=T) + 1.5 * IQR(granzymeb, na.rm=T)))

CD44_plot <- ggplot(CD44_clean, aes(x=group, y=log1p(granzymeb),fill=group)) +
  geom_violin(color='black') +
  # Add Median
  stat_summary(fun = median, geom = "crossbar", 
               width = 0.3, color = "black", size = 0.5) +
  
  # Add Q1 and Q3 as horizontal lines
  stat_summary(fun.data = function(x) {
    q <- quantile(x, probs = c(0.25, 0.75))
    data.frame(y = q, ymin = q[1], ymax = q[2])
  }, geom = "errorbar", width = 0.2, color = "black", size = 0.5) +
  scale_fill_manual(values=c("#E0DCB8","#607C83")) +
  labs(title="granzymeb",
       y="CD44LO gzb expression") +
  theme_classic() +
  theme(legend.position = 'none',
        plot.title = element_text(hjust=0.5))
#axis.title = element_blank())
CD44_plot

#Classify percent of Marcrophages of CD45 cells
#Classify CD45 Table
CD45.table <- data.frame(table(CD45$DonorID))
names(CD45.table)[2] <- "CD45" 

CD68 <- Cytof %>% filter(CellType=="Macrophages")
CD68.table <- data.frame(table(CD68$DonorID))
names(CD68.table)[2] <- "CD68" 
CD68.table <- right_join(CD68.table, CD45.table, by="Var1")
CD68.table$proportion <- (CD68.table$CD68/CD68.table$CD45)*100

Old <- c("HPAP-118","HPAP-105","HPAP-093","HPAP-069","HPAP-066","HPAP-053")
CD68.table <- data.frame(CD68.table,
                        classification=ifelse(CD68.table$Var1 %in% Old,"Old","Young"))
CD68.table$classification <- factor(CD68.table$classification, levels = c("Young","Old"))

CD68_bar <- CD68.table %>% 
  group_by(classification) %>% 
  summarise(mean=mean(proportion, na.rm=T),
            SE=sd(proportion,na.rm=T)/4.8)
ggplot() +
  geom_col(data=CD68_bar, aes(x=classification,y=mean,fill=classification),color='black') +
  geom_errorbar(data=CD68_bar,mapping = aes(x=classification, ymin = mean-SE, ymax = mean+SE),width=0.5) +
  geom_jitter(data=CD68.table, aes(x=classification,y=proportion),width=0.2,alpha=0.5) +
  scale_fill_manual(values=c("#E0DCB8","#607C83"))+
  labs(y="% CD68 of CD45+ Cells")+
  theme_classic() +
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=12))
t.test(log1p(CD68.table$proportion)~CD68.table$classification)

#Classify percent of CD163 expressing macrophages
#Classify CD45 Table
CD45.table <- data.frame(table(CD45$DonorID))
names(CD45.table)[2] <- "CD45" 

Macrophages <- Cytof %>% filter(CellType=="Macrophages")
CD163.table <- Macrophages %>% 
  group_by(DonorID) %>%
  summarise(mean=mean(cd163))
names(CD163.table)[2] <- "CD163" 
Old <- c("HPAP-118","HPAP-105","HPAP-093","HPAP-069","HPAP-066","HPAP-053")
CD163.table$group <- ifelse(CD163.table$DonorID %in% Old,"old","young")
CD163.table$group <- factor(CD163.table$group, levels = c("young","old"))

CD163_bar <- Macrophages %>% 
  group_by(group) %>% 
  summarise(mean=mean(cd163, na.rm=T),
            SE=sd(cd163,na.rm=T)/4.8)


ggplot() +
  geom_col(data=CD163_bar, aes(x=group,y=mean,fill=group),color='black') +
  geom_errorbar(data=CD163_bar,mapping = aes(x=group, ymin = mean-SE, ymax = mean+SE),width=0.5) +
  geom_jitter(data=CD163.table, aes(x=group,y=CD163),width=0.2,alpha=0.5) +
  scale_fill_manual(values=c("#E0DCB8","#607C83"))+
  labs(y="CD163 expression on macrophages")+
  theme_classic() +
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=9.5))
t.test(log1p(CD163.table$CD163)~CD163.table$group)



#END


