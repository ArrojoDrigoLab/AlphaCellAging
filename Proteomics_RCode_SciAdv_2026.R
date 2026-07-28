# =======================================================================================================================
# Alpha cell inflammation during human pancreas aging and type 2 diabetes and its reversal by calorie restriction in mice
# -----------------------------------------------------------------------------------------------------------------------
# Author:      Michael Schleh
# Repository:  https://github.com/ArrojoDrigoLab/AlphaCellAging
# Journal:     Science Advances

# Description:
#   Generates Figure 1 panels for the islet proteomics analysis:
#     - A) Donor metadata plots (age, HbA1c (not shown))
#     - limma differential expression (Young vs. Old)
#     - Reactome pathway enrichment plot
#     - Individual marker violin plots (ER stress, autophagy/mTORC1, and inflammatory markers)
#
# Usage:
#   1. Download packages and set the working directory.
#   2. Ensure the required input files exist in `Figure 1/' (see README).
#   3. Run: Rscript_Figure1_Islet_Proteomics.R
#
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Setup: Working 
# ----------------------------------------------------------------------------
#1A. Working Directory 
setwd("") <---

#1B. Download necessary packages
install.packages(c("tidyverse", "RColorBrewer", "reshape2", "ggh4x", "uniprotREST"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install("limma")

#1C. Call into library
library(tidyverse)
library(limma)
library(RColorBrewer)
library(reshape2)
library(ggh4x)

# -----------------------------------------------------------------------------
# 2. Donor metadata plots (Figure 1A, 1B)
# -----------------------------------------------------------------------------



#Figure 1A_Metadata ----
MetaData <- read.csv("Proteomics_DonorMetadata.csv",header=T,sep=",")
MetaData$Category <- factor(MetaData$Category, levels = c("Young", "Old"))

#1A_AGE----
ggplot(MetaData, aes(x = Category,y=donorage, fill = Category)) + 
  geom_violin(aes(fill=Category)) + 
  geom_jitter(width=0.2,alpha=0.4) +
  labs(title = "Age", y = "years",x="") +
  scale_fill_brewer(palette = "Dark2")+
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(hjust=0.5),
        legend.position = "none",
        axis.text = element_text(color='black'))

#Hba1c (not shown)----
ggplot(MetaData, aes(x = Category,y=hba1c, fill = Category)) + 
  geom_violin(aes(fill=Category)) + 
  geom_jitter(width=0.2,alpha=0.4) +
  labs(title = "HbA1c", y = "%",x="") +
  scale_fill_brewer(palette = "Dark2")+
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(hjust=0.5),
        legend.position = "none",
        axis.text = element_text(color='black'))

# -----------------------------------------------------------------------------
# 3. Limma Differential Expression Analaysis
# -----------------------------------------------------------------------------

#LIMMA PROTEOMICS ----
#LIMMA
MetaData <- read.csv("GroupID_Metadata.csv",header=T,sep=',')
YoungID <- MetaData %>% filter(Group=="Young")
YoungID <- YoungID$ID

df<- read.csv("Proteomics_raw data.csv")
df <- na.omit(df)                                 #<--- Remove NA Proteins
df[-1] <- log2(df[-1])                            #<----Log2 Normalize before LIMMA

sample <- factor(c(rep("Young",30),rep("Old",19)))

design.mat <- model.matrix(~0+sample)
colnames(design.mat) <- levels(sample)
design.mat

contrast.mat <- makeContrasts(
  Diff = Young - Old,
  levels = design.mat)

fit <- lmFit(df, design.mat)
fit2 <- contrasts.fit(fit, contrast.mat)
fit3 <- eBayes(fit2)
deg<-topTable(fit3, coef = 'Diff',
             #p.value = 0.05,
             #adjust.method = 'fdr',
             #lfc =log2(1.5),
             number=nrow(df))


#Map proteins to Uniprot:
uniprot <- read.csv("idmapping_2026_02_18.csv")
DEProteins <- right_join(uniprot,deg, by="feature_id")
####################write.csv(df, "DifferentialAbundance_All.csv")

# -----------------------------------------------------------------------------
# 4. Pathway Analysis
# -----------------------------------------------------------------------------
#Increased and decreased DEGs input into DAVID selected Reactome Pathway
#https://davidbioinformatics.nih.gov/home.jsp

DAVID_UP <- read.delim("Pathways/UP_Old.txt")
DAVID_DOWN <- read.delim("Pathways/UP_YOUNG.txt")

DAVID_UP <- DAVID_UP %>% arrange(PValue)
DAVID_DOWN <- DAVID_DOWN %>% arrange(PValue)

#Filter redundant pathways 
DAVID <- rbind.data.frame(data.frame(DAVID_UP[c(2:3,5:7,9:10,12,14), ],Direction="Up Old"),
                          data.frame(DAVID_DOWN[c(1:3),],Direction="Down Old"))

DAVID <- data.frame(DAVID,
                    P.Value.New=ifelse(DAVID$Direction=="Down Old",(-log(DAVID$PValue)*-1),-log(DAVID$PValue)))

DAVID[c('GO', 'Term')] <- str_split_fixed(DAVID$Term, '~', 2)

#Reactome Pathway Plot (Figure 1B)----
ggplot(DAVID, aes(x=P.Value.New,y=reorder(Term,P.Value.New))) +
  geom_segment(aes(yend = reorder(Term,P.Value.New), xend = 0),color ="black",size=0.8,alpha=0.5) +
  geom_point(aes(size=Count,fill=Fold.Enrichment),shape=21,color="black") +
  labs(x=bquote(-log[10](FDR)),
       title="GO Biological process: M2 poliarzation") +
  scale_fill_viridis_c() +
  geom_vline(xintercept = 0,color="black",size=1) +
  theme(plot.title = element_blank(),
        axis.text = element_text(size=12,color="black"),
        axis.title = element_text(size=12,color="black"),
        axis.title.y = element_blank(),
        axis.line = element_line(size=1),
        legend.position = "none",
        strip.background = element_rect(fill = "white", color = "black"),
        strip.text = element_text(color = "black", face = "bold", size = 10),
        panel.background = element_blank()) +
  facet_grid(factor(Direction,levels=c("Up Old","Down Old")) ~ ., scales = "free_y", space = "free_y")


# -----------------------------------------------------------------------------
#Individual BARPLOTS----
# -----------------------------------------------------------------------------

#Figure 1C - Beta Cell Immaturity ----
#RBP4
RBP4 <- df %>% filter(feature_id=="P02753")
RBP4 <- melt(RBP4,by="feature_id")
RBP4$group <- ifelse(RBP4$variable %in% YoungID,"Young","Old")
RBP4$group <- factor(RBP4$group, levels = c("Young", "Old"))
RBP4$value <- scale(RBP4$value,center=T,scale=T)

RBP4.bar <- RBP4%>%
  group_by(group) %>%
  summarize(Mean = mean(value, na.rm=TRUE),
            SD=sd(value, na.rm=TRUE))
ggplot() +
  geom_violin(data=RBP4, aes(x=group,y=value,fill=group),color='black')+
  geom_jitter(data=RBP4, aes(x=group,y=value),width=0.2)+
  geom_errorbar(data=RBP4.bar, aes(x=group,ymin=Mean-SD, ymax=Mean+SD),width=0.2)+
  # geom_hline(yintercept = 6.5) +
  scale_fill_brewer(palette = "Dark2")+
  labs(title="RBP4",
       y="Log(Intensity)")+
  theme_classic(base_size = 12)+
  theme(plot.title = element_text(hjust=0.5),
        axis.text = element_text(color='black'),
        axis.title = element_blank(),
        legend.position = "none")

ggsave("RBP4.png",dpi=600,height = 2, width=1.5,units="in")

#1Figure 1D - Senescence ----
#CDKN1A
CDKN1A <- df %>% filter(feature_id=="P38936")
CDKN1A <- melt(CDKN1A,by="feature_id")
CDKN1A$group <- ifelse(CDKN1A$variable %in% YoungID,"Young","Old")
CDKN1A$group <- factor(CDKN1A$group, levels = c("Young", "Old"))
CDKN1A$value <- scale(CDKN1A$value,center=T,scale=T)


CDKN1A.bar <- CDKN1A%>%
  group_by(group) %>%
  summarize(Mean = mean(value, na.rm=TRUE),
            SD=sd(value, na.rm=TRUE))
ggplot() +
  geom_violin(data=CDKN1A, aes(x=group,y=value,fill=group),color='black')+
  geom_jitter(data=CDKN1A, aes(x=group,y=value),width=0.2)+
  geom_errorbar(data=CDKN1A.bar, aes(x=group,ymin=Mean-SD, ymax=Mean+SD),width=0.2)+
  # geom_hline(yintercept = 6.5) +
  scale_fill_brewer(palette = "Dark2")+
  labs(title="CDKN1A/p21",
         y="Log(Intensity)")+
  theme_classic(base_size = 12)+
  theme(plot.title = element_text(hjust=0.5),
        axis.text = element_text(color='black'),
        axis.title = element_blank(),
        legend.position = "none")
ggsave("CDKN1A.png",dpi=600,height = 2, width=1.5,units="in")

#1E - ER STRESS MARKERS----
#i) CHOP
CHOP <- df %>% filter(feature_id=="P35638")
CHOP <- melt(CHOP,by="feature_id")
CHOP$group <- ifelse(CHOP$variable %in% YoungID,"Young","Old")
CHOP$group <- factor(CHOP$group, levels = c("Young", "Old"))
CHOP$value <- scale(CHOP$value,center=T,scale=T)


CHOP.bar <- CHOP%>%
  group_by(group) %>%
  summarize(Mean = mean(value, na.rm=TRUE),
            SD=sd(value, na.rm=TRUE))
ggplot() +
  geom_violin(data=CHOP, aes(x=group,y=value,fill=group),color='black')+
  geom_jitter(data=CHOP, aes(x=group,y=value),width=0.2)+
  geom_errorbar(data=CHOP.bar, aes(x=group,ymin=Mean-SD, ymax=Mean+SD),width=0.2)+
  # geom_hline(yintercept = 6.5) +
  scale_fill_brewer(palette = "Dark2")+
  labs(title="DDIT3/CHOP",
       y="")+
  theme_classic(base_size = 12)+
  theme(plot.title = element_text(hjust=0.5),
        axis.text = element_text(color='black'),
        axis.title = element_blank(),
        legend.position = "none")
ggsave("CHOP.png",dpi=600,height = 2, width=1.5,units="in")

#ii) ATF6A
ATF6 <- df %>% filter(feature_id=="P18850")
ATF6 <- melt(ATF6,by="feature_id")
ATF6$group <- ifelse(ATF6$variable %in% YoungID,"Young","Old")
ATF6$group <- factor(ATF6$group, levels = c("Young", "Old"))
ATF6$value <- scale(ATF6$value,center=T,scale=T)


ATF6.bar <- ATF6%>%
  group_by(group) %>%
  summarize(Mean = mean(value, na.rm=TRUE),
            SD=sd(value, na.rm=TRUE))
ggplot() +
  geom_violin(data=ATF6, aes(x=group,y=value,fill=group),color='black')+
  geom_jitter(data=ATF6, aes(x=group,y=value),width=0.2)+
  geom_errorbar(data=ATF6.bar, aes(x=group,ymin=Mean-SD, ymax=Mean+SD),width=0.2)+
  # geom_hline(yintercept = 6.5) +
  scale_fill_brewer(palette = "Dark2")+
  labs(title="ATF6",
       y="scaled intensity")+
  theme_classic(base_size = 12)+
  theme(plot.title = element_text(hjust=0.5),
        axis.text = element_text(color='black'),
        axis.title = element_blank(),
        legend.position = "none")
ggsave("ATF6.png",dpi=600,height = 2, width=1.5,units="in")

#iii) PERK
PERK <- df %>% filter(feature_id=="Q9NZJ5")
PERK <- melt(PERK,by="feature_id")
PERK$group <- ifelse(PERK$variable %in% YoungID,"Young","Old")
PERK$group <- factor(PERK$group, levels = c("Young", "Old"))
PERK$value <- scale(PERK$value,center=T,scale=T)


PERK.bar <- PERK%>%
  group_by(group) %>%
  summarize(Mean = mean(value, na.rm=TRUE),
            SD=sd(value, na.rm=TRUE))
ggplot() +
  geom_violin(data=PERK, aes(x=group,y=value,fill=group),color='black')+
  geom_jitter(data=PERK, aes(x=group,y=value),width=0.2)+
  geom_errorbar(data=PERK.bar, aes(x=group,ymin=Mean-SD, ymax=Mean+SD),width=0.2)+
  # geom_hline(yintercept = 6.5) +
  scale_fill_brewer(palette = "Dark2")+
  labs(title="EIF2AK3/PERK",
       y="")+
  theme_classic(base_size = 12)+
  theme(plot.title = element_text(hjust=0.5,size=11),
        axis.text = element_text(color='black'),
        axis.title = element_blank(),
        legend.position = "none")
ggsave("PERK.png",dpi=600,height = 2, width=1.5,units="in")


#Figure 1F - mTORC1 activation----
#i) LAMP2
LAMP2 <- df %>% filter(feature_id=="P13473")
LAMP2 <- melt(LAMP2,by="feature_id")
LAMP2$group <- ifelse(LAMP2$variable %in% YoungID,"Young","Old")
LAMP2$group <- factor(LAMP2$group, levels = c("Young", "Old"))
LAMP2$value <- scale(LAMP2$value,center=T,scale=T)


LAMP2.bar <- LAMP2%>%
  group_by(group) %>%
  summarize(Mean = mean(value, na.rm=TRUE),
            SD=sd(value, na.rm=TRUE))
ggplot() +
  geom_violin(data=LAMP2, aes(x=group,y=value,fill=group),color='black')+
  geom_jitter(data=LAMP2, aes(x=group,y=value),width=0.2)+
  geom_errorbar(data=LAMP2.bar, aes(x=group,ymin=Mean-SD, ymax=Mean+SD),width=0.2)+
  # geom_hline(yintercept = 6.5) +
  scale_fill_brewer(palette = "Dark2")+
  labs(title="LAMP2",
       y="Log(Intensity)")+
  theme_classic(base_size = 12)+
  theme(plot.title = element_text(hjust=0.5,size=11),
        axis.text = element_text(color='black'),
        axis.title = element_blank(),
        legend.position = "none")
ggsave("LAMP2.png",dpi=600,height = 2, width=1.5,units="in")

#ii) LAMTOR
LTOR <- df %>% filter(feature_id=="Q6IAA8")
LTOR <- melt(LTOR,by="feature_id")
LTOR$group <- ifelse(LTOR$variable %in% YoungID,"Young","Old")
LTOR$group <- factor(LTOR$group, levels = c("Young", "Old"))
LTOR$value <- scale(LTOR$value,center=T,scale=T)


LTOR.bar <- LTOR%>%
  group_by(group) %>%
  summarize(Mean = mean(value, na.rm=TRUE),
            SD=sd(value, na.rm=TRUE))

ggplot() +
  geom_violin(data=LTOR, aes(x=group,y=value,fill=group),color='black')+
  geom_jitter(data=LTOR, aes(x=group,y=value),width=0.2)+
  geom_errorbar(data=LTOR.bar, aes(x=group,ymin=Mean-SD, ymax=Mean+SD),width=0.2)+
  # geom_hline(yintercept = 6.5) +
  scale_fill_brewer(palette = "Dark2")+
  labs(title="LAMTOR",
       y="")+
  theme_bw()+
  theme_classic(base_size = 12)+
  theme(plot.title = element_text(hjust=0.5,size=11),
        axis.text = element_text(color='black'),
        axis.title = element_blank(),
        legend.position = "none")
ggsave("LAMTOR.png",dpi=600,height = 2, width=1.5,units="in")


#iii) RAGA
RAGA <- df %>% filter(feature_id=="Q7L523")
RAGA <- melt(RAGA,by="feature_id")
RAGA$group <- ifelse(RAGA$variable %in% YoungID,"Young","Old")
RAGA$group <- factor(RAGA$group, levels = c("Young", "Old"))
RAGA$value <- scale(RAGA$value,center=T,scale=T)

RAGA.bar <- RAGA%>%
  group_by(group) %>%
  summarize(Mean = mean(value, na.rm=TRUE),
            SD=sd(value, na.rm=TRUE))

ggplot() +
  geom_violin(data=RAGA, aes(x=group,y=value,fill=group),color='black')+
  geom_jitter(data=RAGA, aes(x=group,y=value),width=0.2)+
  geom_errorbar(data=RAGA.bar, aes(x=group,ymin=Mean-SD, ymax=Mean+SD),width=0.2)+
  # geom_hline(yintercept = 6.5) +
  scale_fill_brewer(palette = "Dark2")+
  labs(title="RAGA",
       y="")+
  theme_classic(base_size = 12)+
  theme(plot.title = element_text(hjust=0.5,size=11),
        axis.text = element_text(color='black'),
        axis.title = element_blank(),
        legend.position = "none")
ggsave("RAGA.png",dpi=600,height = 2, width=1.5,units="in")

#iv) RAGC
RAGC <- df %>% filter(feature_id=="Q9HB90")
RAGC <- melt(RAGC,by="feature_id")
RAGC$group <- ifelse(RAGC$variable %in% YoungID,"Young","Old")
RAGC$group <- factor(RAGC$group, levels = c("Young", "Old"))
RAGC$value <- scale(RAGC$value,center=T,scale=T)


RAGC.bar <- RAGC%>%
  group_by(group) %>%
  summarize(Mean = mean(value, na.rm=TRUE),
            SD=sd(value, na.rm=TRUE))

ggplot() +
  geom_violin(data=RAGC, aes(x=group,y=value,fill=group),color='black')+
  geom_jitter(data=RAGC, aes(x=group,y=value),width=0.2)+
  geom_errorbar(data=RAGC.bar, aes(x=group,ymin=Mean-SD, ymax=Mean+SD),width=0.2)+
  # geom_hline(yintercept = 6.5) +
  scale_fill_brewer(palette = "Dark2")+
  labs(title="RAGC",
       y="")+
  theme_classic(base_size = 12)+
  theme(plot.title = element_text(hjust=0.5,size=11),
        axis.text = element_text(color='black'),
        axis.title = element_blank(),
        legend.position = "none")
ggsave("RAGC.png",dpi=600,height = 2, width=1.5,units="in")

#Figure 1G - Inflammatory Measures ----

#i) HLA-A
HLA_A <- df %>% filter(feature_id=="P04439")
HLA_A <- melt(HLA_A,by="feature_id")
HLA_A$group <- ifelse(HLA_A$variable %in% YoungID,"Young","Old")
HLA_A$group <- factor(HLA_A$group, levels = c("Young", "Old"))
HLA_A$value <- scale(HLA_A$value,center=T,scale=T)

HLA_A.bar <- HLA_A%>%
  group_by(group) %>%
  summarize(Mean = mean(value, na.rm=TRUE),
            SD=sd(value, na.rm=TRUE))
ggplot() +
  geom_violin(data=HLA_A, aes(x=group,y=value,fill=group),color='black')+
  geom_jitter(data=HLA_A, aes(x=group,y=value),width=0.2)+
  geom_errorbar(data=HLA_A.bar, aes(x=group,ymin=Mean-SD, ymax=Mean+SD),width=0.2)+
  # geom_hline(yintercept = 6.5) +
  scale_fill_brewer(palette = "Dark2")+
  labs(title="HLA-A",
       y="")+
  theme_classic(base_size = 12)+
  theme(plot.title = element_text(hjust=0.5,size=11),
        axis.text = element_text(color='black'),
        axis.title = element_blank(),
        legend.position = "none")
ggsave("HLA_A.png",dpi=600,height = 2, width=1.5,units="in")

#ii) TGFB1
TGFB1 <- df %>% filter(feature_id=="P01137")
TGFB1 <- melt(TGFB1,by="feature_id")
TGFB1$group <- ifelse(TGFB1$variable %in% YoungID,"Young","Old")
TGFB1$group <- factor(TGFB1$group, levels = c("Young", "Old"))
TGFB1$value <- scale(TGFB1$value,center=T,scale=T)


TGFB1.bar <- TGFB1%>%
  group_by(group) %>%
  summarize(Mean = mean(value, na.rm=TRUE),
            SD=sd(value, na.rm=TRUE))
ggplot() +
  geom_violin(data=TGFB1, aes(x=group,y=value,fill=group),color='black')+
  geom_jitter(data=TGFB1, aes(x=group,y=value),width=0.2)+
  geom_errorbar(data=TGFB1.bar, aes(x=group,ymin=Mean-SD, ymax=Mean+SD),width=0.2)+
  # geom_hline(yintercept = 6.5) +
  scale_fill_brewer(palette = "Dark2")+
  labs(title="TGFβ1",
       y="Log(Intensity)")+
  theme_classic(base_size = 12)+
  theme(plot.title = element_text(hjust=0.5,size=11),
        axis.text = element_text(color='black'),
        axis.title = element_blank(),
        legend.position = "none")
ggsave("TGFB1.png",dpi=600,height = 2, width=1.5,units="in")


#iii) TOM1
TOM1 <- df %>% filter(feature_id=="O60784")
TOM1 <- melt(TOM1,by="feature_id")
TOM1$group <- ifelse(TOM1$variable %in% YoungID,"Young","Old")
TOM1$group <- factor(TOM1$group, levels = c("Young", "Old"))
TOM1$value <- scale(TOM1$value,center=T,scale=T)

TOM1.bar <- TOM1%>%
  group_by(group) %>%
  summarize(Mean = mean(value, na.rm=TRUE),
            SD=sd(value, na.rm=TRUE))
ggplot() +
  geom_violin(data=TOM1, aes(x=group,y=value,fill=group),color='black')+
  geom_jitter(data=TOM1, aes(x=group,y=value),width=0.2)+
  geom_errorbar(data=TOM1.bar, aes(x=group,ymin=Mean-SD, ymax=Mean+SD),width=0.2)+
  # geom_hline(yintercept = 6.5) +
  scale_fill_brewer(palette = "Dark2")+
  labs(title="TOM1",
       y="")+
  theme_classic(base_size = 12)+
  theme(plot.title = element_text(hjust=0.5,size=11),
        axis.text = element_text(color='black'),
        axis.title = element_blank(),
        legend.position = "none")
ggsave("TOM1.png",dpi=600,height = 2, width=1.5,units="in")

#iv) CD63
CD63 <- df %>% filter(feature_id=="P08962")
CD63 <- melt(CD63,by="feature_id")
CD63$group <- ifelse(CD63$variable %in% YoungID,"Young","Old")
CD63$group <- factor(CD63$group, levels = c("Young", "Old"))
CD63$value <- scale(CD63$value,center=T,scale=T)


CD63.bar <- CD63%>%
  group_by(group) %>%
  summarize(Mean = mean(value, na.rm=TRUE),
            SD=sd(value, na.rm=TRUE))
ggplot() +
  geom_violin(data=CD63, aes(x=group,y=value,fill=group),color='black')+
  geom_jitter(data=CD63, aes(x=group,y=value),width=0.2)+
  geom_errorbar(data=CD63.bar, aes(x=group,ymin=Mean-SD, ymax=Mean+SD),width=0.2)+
  # geom_hline(yintercept = 6.5) +
  scale_fill_brewer(palette = "Dark2")+
  labs(title="CD63",
       y="Log(Intensity)")+
  theme_classic(base_size = 12)+
  theme(plot.title = element_text(hjust=0.5,size=11),
        axis.text = element_text(color='black'),
        axis.title = element_blank(),
        legend.position = "none")
ggsave("CD63.png",dpi=600,height = 2, width=1.5,units="in")


#v) CD44
CD44 <- df %>% filter(feature_id=="P16070")
CD44 <- melt(CD44,by="feature_id")
CD44$group <- ifelse(CD44$variable %in% YoungID,"Young","Old")
CD44$group <- factor(CD44$group, levels = c("Young", "Old"))
CD44$value <- scale(CD44$value,center=T,scale=T)


CD44.bar <- CD44%>%
  group_by(group) %>%
  summarize(Mean = mean(value, na.rm=TRUE),
            SD=sd(value, na.rm=TRUE))
ggplot() +
  geom_violin(data=CD44, aes(x=group,y=value,fill=group),color='black')+
  geom_jitter(data=CD44, aes(x=group,y=value),width=0.2)+
  geom_errorbar(data=CD44.bar, aes(x=group,ymin=Mean-SD, ymax=Mean+SD),width=0.2)+
  # geom_hline(yintercept = 6.5) +
  scale_fill_brewer(palette = "Dark2")+
  labs(title="CD44",
       y="Log(Intensity)")+
  theme_classic(base_size = 12)+
  theme(plot.title = element_text(hjust=0.5,size=11),
        axis.text = element_text(color='black'),
        axis.title = element_blank(),
        legend.position = "none")
ggsave("CD44.png",dpi=600,height = 2, width=1.5,units="in")


#vi) CD28
CD28 <- df %>% filter(feature_id=="P10747")
CD28 <- melt(CD28,by="feature_id")
CD28$group <- ifelse(CD28$variable %in% YoungID,"Young","Old")
CD28$group <- factor(CD28$group, levels = c("Young", "Old"))
CD28$value <- scale(CD28$value,center=T,scale=T)


CD28.bar <- CD28%>%
  group_by(group) %>%
  summarize(Mean = mean(value, na.rm=TRUE),
            SD=sd(value, na.rm=TRUE))
ggplot() +
  geom_violin(data=CD28, aes(x=group,y=value,fill=group),color='black')+
  geom_jitter(data=CD28, aes(x=group,y=value),width=0.2)+
  geom_errorbar(data=CD28.bar, aes(x=group,ymin=Mean-SD, ymax=Mean+SD),width=0.2)+
  # geom_hline(yintercept = 6.5) +
  scale_fill_brewer(palette = "Dark2")+
  labs(title="CD28",
       y="Log(Intensity)")+
  theme_classic(base_size = 12)+
  theme(plot.title = element_text(hjust=0.5,size=11),
        axis.text = element_text(color='black'),
        axis.title = element_blank(),
        legend.position = "none")
ggsave("CD28.png",dpi=600,height = 2, width=1.5,units="in")



#END




