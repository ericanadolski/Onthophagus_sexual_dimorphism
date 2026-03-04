##### load packages######
library(plyr)
library(tidyverse)
library(tximport)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(ggbreak)
library(gplots)
library(edgeR)
library(stats)
library(topGO)
library(RColorBrewer)
# I. RNAseq data analyses #############################################################################
###### Step 1: import RNAseq read abundance from salmon #####
### O taurus ######
# read in study design
Ot_RNA_info <- read.table("/Users/ericanadolski/GitHub/Onthophagus_sexual_dimorphism/Ot-RNA-sample-info.txt", header=TRUE)

# create file paths to abundance files using 'file.path' function
path <- file.path("/Volumes/T7_Drive_2/RNAseq/salmon-dsx", paste0(Ot_RNA_info$Sample, "_quant/quant.sf")) 
file.exists(path)

# import counts with tximport
Ot_salmon_tx <- tximport(path, 
                         type = "salmon",
                         txOut = TRUE,
                         countsFromAbundance = "lengthScaledTPM",
                         ignoreTxVersion = TRUE)
class(Ot_salmon_tx)
names(Ot_salmon_tx)

Ot_counts <- as.data.frame(Ot_salmon_tx$counts)
colnames(Ot_counts) <- Ot_RNA_info$Sample

### O sagittarius #####
Os_RNA_info <- read.table("/Users/ericanadolski/GitHub/Onthophagus_sexual_dimorphism/Os-RNA-sample-info.txt", header=TRUE)
path <- file.path("/Volumes/T7_Drive_2/RNAseq/salmon-dsx", paste0(Os_RNA_info$Sample, "_quant/quant.sf")) 
file.exists(path)

Os_salmon_tx <- tximport(path, 
                         type = "salmon",
                         txOut = TRUE,
                         countsFromAbundance = "lengthScaledTPM",
                         ignoreTxVersion = TRUE)

Os_counts <- as.data.frame(Os_salmon_tx$counts)
colnames(Os_counts) <- Os_RNA_info$Sample
### D gazella #####
Dg_RNA_info <- read.table("/Users/ericanadolski/GitHub/Onthophagus_sexual_dimorphism/Dg-RNA-sample-info.txt", header=TRUE)
path <- file.path("/Volumes/T7_Drive_2/RNAseq/salmon-dsx", paste0(Dg_RNA_info$Sample, "_quant/quant.sf")) 
file.exists(path)

Dg_salmon_tx <- tximport(path, 
                         type = "salmon",
                         txOut = TRUE,
                         countsFromAbundance = "lengthScaledTPM",
                         ignoreTxVersion = TRUE)

Dg_counts <- as.data.frame(Dg_salmon_tx$counts)
colnames(Dg_counts) <- Dg_RNA_info$Sample

###### Step 2: data filtering and visualization #####
##### O taurus ######
# generate DESeq dataset for exploratory analysis
Ot_dds_exp <- DESeqDataSetFromTximport(Ot_salmon_tx, Ot_RNA_info, design = ~ 1)
nrow(Ot_dds_exp) # 19500 

# filter out genes with very low read count across most samples
Ot_dds_exp <- Ot_dds_exp[rowSums(counts(Ot_dds_exp)) > 5, ]
nrow(Ot_dds_exp) # 15408

Ot_dds_rld <- rlog(Ot_dds_exp)
Ot_dds_rld <- as.data.frame(assay(Ot_dds_rld))
Ot_dds_rld <- as.data.frame(t(Ot_dds_rld))

Ot_sampleDists <- dist(Ot_dds_rld)
Ot_sampleDistMatrix <- as.matrix(Ot_sampleDists)

### heatmap all traits ####
rownames(Ot_sampleDistMatrix) <- paste(Ot_RNA_info$Sample)
colnames(Ot_sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
Ot_heatmap <- pheatmap(Ot_sampleDistMatrix,
         clustering_distance_rows=Ot_sampleDists,
         clustering_distance_cols=Ot_sampleDists,
         col=colors)

ggsave(Ot_heatmap, file = "/Users/ericanadolski/Desktop/beetle-sex-dimorph/figures/Ot_RNA_heatmap.pdf",
       width = 8, height = 8)

### PCA all traits ####
PCA_Ot_RNA <- prcomp(Ot_dds_rld, center = TRUE)
PCs_Ot_RNA <- as.data.frame(PCA_Ot_RNA$x)
PCs_Ot_RNA$Sample <- row.names(PCs_Ot_RNA)
PCs_Ot_RNA$Sample <- Ot_RNA_info$Sample
PCs_Ot_RNA <- inner_join(PCs_Ot_RNA, Ot_RNA_info, by='Sample') 
summary(PCA_Ot_RNA)

Ot_PCA <- ggplot(data = PCs_Ot_RNA, aes(x = PC1, y = PC2, color = Sex, shape=Trait)) + 
  geom_point(size = 6, alpha = 0.8) +
  scale_color_manual(values=c("#C1272D", "#2166AC")) + 
  #geom_text(aes(label = Replicate), nudge_y = 3) +
  labs(title="O. taurus RNA-seq",x = "PC1 (29.0%)", y = "PC2 (13.3%)") +
  theme_bw(base_size=16)
Ot_PCA

ggsave(Ot_PCA, file = "/Users/ericanadolski/Desktop/beetle-sex-dimorph/figures/Ot_RNA_PCA.pdf",
       width = 5, height = 4)

### PCA posterior head ####
Ot_PH_info <- filter(Ot_RNA_info, Trait == "PH")
Ot_PH_counts <- Ot_counts %>% dplyr::select(contains("PH")) 

dds_Ot_PH_dg <- DESeqDataSetFromMatrix(countData = round(Ot_PH_counts),
                                       colData = Ot_PH_info,
                                       design = ~ 1)
rlog_dds_Ot_PH <- rlog(dds_Ot_PH_dg) 
rlog_counts_Ot_PH <- as.data.frame(assay(rlog_dds_Ot_PH))
rlog_countsT_Ot_PH <- as.data.frame(t(rlog_counts_Ot_PH))

PCA_Ot_PH <- prcomp(rlog_countsT_Ot_PH, center = TRUE)
PCs_Ot_PH <- as.data.frame(PCA_Ot_PH$x)
PCs_Ot_PH$Sample <- row.names(PCs_Ot_PH)
PCs_Ot_PH <- inner_join(PCs_Ot_PH, Ot_PH_info, by='Sample') 
summary(PCA_Ot_PH)

PCA_Ot_PH <- ggplot(data = PCs_Ot_PH, 
                    aes(x = PC1, y = PC2, color=Sex)) + 
  geom_point(size = 6, shape=16) +
  scale_color_manual(values=c("#C1272D", "#2166AC")) +
  labs(title="O.taurus RNAseq PH",x = "PC1 (27.9%)", y = "PC2 (12.4%)") +
  theme_bw()

ggsave(PCA_Ot_PH, file = "/Users/ericanadolski/Desktop/beetle-sex-dimorph/figures/Ot_RNA_PCA_PH.pdf",
       width = 4.30, height = 3.77)

### PCA anterior head ####
Ot_AH_info <- filter(Ot_RNA_info, Trait == "AH")
Ot_AH_counts <- Ot_counts %>% dplyr::select(contains("AH")) 

dds_Ot_AH_dg <- DESeqDataSetFromMatrix(countData = round(Ot_AH_counts),
                                       colData = Ot_AH_info,
                                       design = ~ 1)
rlog_dds_Ot_AH <- rlog(dds_Ot_AH_dg) 
rlog_counts_Ot_AH <- as.data.frame(assay(rlog_dds_Ot_AH))
rlog_countsT_Ot_AH <- as.data.frame(t(rlog_counts_Ot_AH))

PCA_Ot_AH <- prcomp(rlog_countsT_Ot_AH, center = TRUE)
PCs_Ot_AH <- as.data.frame(PCA_Ot_AH$x)
PCs_Ot_AH$Sample <- row.names(PCs_Ot_AH)
PCs_Ot_AH <- inner_join(PCs_Ot_AH, Ot_AH_info, by='Sample') 
summary(PCA_Ot_AH)

PCA_Ot_AH <- ggplot(data = PCs_Ot_AH, 
                    aes(x = PC1, y = PC2, color=Sex)) + 
  geom_point(size = 6, shape=16) +
  scale_color_manual(values=c("#C1272D", "#2166AC")) +
  labs(title="O.taurus RNAseq AH",x = "PC1 (31.9%)", y = "PC2 (13.6%)") +
  theme_bw()

ggsave(PCA_Ot_AH, file = "/Users/ericanadolski/Desktop/beetle-sex-dimorph/figures/Ot_RNA_PCA_AH.pdf",
       width = 4.30, height = 3.77)

### PCA genitalia ####
Ot_G_info <- filter(Ot_RNA_info, Trait == "G")
Ot_G_counts <- Ot_counts %>% dplyr::select(contains("G")) 

dds_Ot_G_dg <- DESeqDataSetFromMatrix(countData = round(Ot_G_counts),
                                      colData = Ot_G_info,
                                      design = ~ 1)
rlog_dds_Ot_G <- rlog(dds_Ot_G_dg) 
rlog_counts_Ot_G <- as.data.frame(assay(rlog_dds_Ot_G))
rlog_countsT_Ot_G <- as.data.frame(t(rlog_counts_Ot_G))

PCA_Ot_G <- prcomp(rlog_countsT_Ot_G, center = TRUE)
PCs_Ot_G <- as.data.frame(PCA_Ot_G$x)
PCs_Ot_G$Sample <- row.names(PCs_Ot_G)
PCs_Ot_G <- inner_join(PCs_Ot_G, Ot_G_info, by='Sample') 
summary(PCA_Ot_G)

PCA_Ot_G <- ggplot(data = PCs_Ot_G, 
                   aes(x = PC1, y = PC2, color=Sex)) + 
  geom_point(size = 6, shape=16) +
  scale_color_manual(values=c("#C1272D", "#2166AC")) +
  labs(title="O.taurus RNAseq G",x = "PC1 (42.5%)", y = "PC2 (14.0%)") +
  theme_bw()

ggsave(PCA_Ot_G, file = "/Users/ericanadolski/Desktop/beetle-sex-dimorph/figures/Ot_RNA_PCA_G.pdf",
       width = 4.30, height = 3.77)

##### O sagittarius ######
# generate DESeq dataset for exploratory analysis
Os_dds_exp <- DESeqDataSetFromTximport(Os_salmon_tx, Os_RNA_info, design = ~ 1)
nrow(Os_dds_exp) # 20609

# filter out genes with very low read count across most samples
Os_dds_exp <- Os_dds_exp[rowSums(counts(Os_dds_exp)) > 5, ]
nrow(Os_dds_exp) # 15056

Os_dds_rld <- rlog(Os_dds_exp)
Os_dds_rld <- as.data.frame(assay(Os_dds_rld))
Os_dds_rld <- as.data.frame(t(Os_dds_rld))

Os_sampleDists <- dist(Os_dds_rld)
Os_sampleDistMatrix <- as.matrix(Os_sampleDists)

### heatmap all traits ####
rownames(Os_sampleDistMatrix) <- paste(Os_RNA_info$Sample)
colnames(Os_sampleDistMatrix) <- NULL
Os_heatmap <- pheatmap(Os_sampleDistMatrix,
                       clustering_distance_rows=Os_sampleDists,
                       clustering_distance_cols=Os_sampleDists,
                       col=colors)

ggsave(Os_heatmap, file = "/Users/ericanadolski/Desktop/beetle-sex-dimorph/figures/Os_RNA_heatmap.pdf",
       width = 8, height = 8)

### PCA all traits ####
PCA_Os_RNA <- prcomp(Os_dds_rld, center = TRUE)
PCs_Os_RNA <- as.data.frame(PCA_Os_RNA$x)
PCs_Os_RNA$Sample <- row.names(PCs_Os_RNA)
PCs_Os_RNA$Sample <- Os_RNA_info$Sample
PCs_Os_RNA <- inner_join(PCs_Os_RNA, Os_RNA_info, by='Sample') 
summary(PCA_Os_RNA)

Os_PCA <- ggplot(data = PCs_Os_RNA, aes(x = PC1, y = PC2, color = Sex, shape=Trait)) + 
  geom_point(size = 6, alpha = 0.8) +
  scale_color_manual(values=c("#C1272D", "#2166AC")) + 
  #geom_text(aes(label = Replicate), nudge_y = 3) +
  labs(title="O. sagittarius RNA-seq",x = "PC1 (35.6%)", y = "PC2 (14.6%)") +
  theme_bw(base_size=16)

ggsave(Os_PCA, file = "/Users/ericanadolski/Desktop/beetle-sex-dimorph/figures/Os_RNA_PCA.pdf",
       width = 5, height = 4)

### PCA posterior head ####
Os_PH_info <- filter(Os_RNA_info, Trait == "PH")
Os_PH_counts <- Os_counts %>% dplyr::select(contains("PH")) 

dds_Os_PH_dg <- DESeqDataSetFromMatrix(countData = round(Os_PH_counts),
                                       colData = Os_PH_info,
                                       design = ~ 1)
rlog_dds_Os_PH <- rlog(dds_Os_PH_dg) 
rlog_counts_Os_PH <- as.data.frame(assay(rlog_dds_Os_PH))
rlog_countsT_Os_PH <- as.data.frame(t(rlog_counts_Os_PH))

PCA_Os_PH <- prcomp(rlog_countsT_Os_PH, center = TRUE)
PCs_Os_PH <- as.data.frame(PCA_Os_PH$x)
PCs_Os_PH$Sample <- row.names(PCs_Os_PH)
PCs_Os_PH <- inner_join(PCs_Os_PH, Os_PH_info, by='Sample') 
summary(PCA_Os_PH)

PCA_Os_PH <- ggplot(data = PCs_Os_PH, 
                    aes(x = PC1, y = PC2, color=Sex)) + 
  geom_point(size = 6, shape=16) +
  scale_color_manual(values=c("#C1272D", "#2166AC")) +
  labs(title="O.sagittarius RNAseq PH",x = "PC1 (38.1%)", y = "PC2 (16.5%)") +
  theme_bw()

ggsave(PCA_Os_PH, file = "/Users/ericanadolski/Desktop/beetle-sex-dimorph/figures/Os_RNA_PCA_PH.pdf",
       width = 4.30, height = 3.77)

### PCA anterior head ####
Os_AH_info <- filter(Os_RNA_info, Trait == "AH")
Os_AH_counts <- Os_counts %>% dplyr::select(contains("AH")) 

dds_Os_AH_dg <- DESeqDataSetFromMatrix(countData = round(Os_AH_counts),
                                       colData = Os_AH_info,
                                       design = ~ 1)
rlog_dds_Os_AH <- rlog(dds_Os_AH_dg) 
rlog_counts_Os_AH <- as.data.frame(assay(rlog_dds_Os_AH))
rlog_countsT_Os_AH <- as.data.frame(t(rlog_counts_Os_AH))

PCA_Os_AH <- prcomp(rlog_countsT_Os_AH, center = TRUE)
PCs_Os_AH <- as.data.frame(PCA_Os_AH$x)
PCs_Os_AH$Sample <- row.names(PCs_Os_AH)
PCs_Os_AH <- inner_join(PCs_Os_AH, Os_AH_info, by='Sample') 
summary(PCA_Os_AH)

PCA_Os_AH <- ggplot(data = PCs_Os_AH, 
                    aes(x = PC1, y = PC2, color=Sex)) + 
  geom_point(size = 6, shape=16) +
  scale_color_manual(values=c("#C1272D", "#2166AC")) +
  labs(title="O.sagittarius RNAseq AH",x = "PC1 (48.1%)", y = "PC2 (9.3%)") +
  theme_bw()

ggsave(PCA_Os_AH, file = "/Users/ericanadolski/Desktop/beetle-sex-dimorph/figures/Os_RNA_PCA_AH.pdf",
       width = 4.30, height = 3.77)

### PCA genitalia ####
Os_G_info <- filter(Os_RNA_info, Trait == "G")
Os_G_counts <- Os_counts %>% dplyr::select(contains("G")) 

dds_Os_G_dg <- DESeqDataSetFromMatrix(countData = round(Os_G_counts),
                                      colData = Os_G_info,
                                      design = ~ 1)
rlog_dds_Os_G <- rlog(dds_Os_G_dg) 
rlog_counts_Os_G <- as.data.frame(assay(rlog_dds_Os_G))
rlog_countsT_Os_G <- as.data.frame(t(rlog_counts_Os_G))

PCA_Os_G <- prcomp(rlog_countsT_Os_G, center = TRUE)
PCs_Os_G <- as.data.frame(PCA_Os_G$x)
PCs_Os_G$Sample <- row.names(PCs_Os_G)
PCs_Os_G <- inner_join(PCs_Os_G, Os_G_info, by='Sample') 
summary(PCA_Os_G)

PCA_Os_G <- ggplot(data = PCs_Os_G, 
                   aes(x = PC1, y = PC2, color=Sex)) + 
  geom_point(size = 6, shape=16) +
  scale_color_manual(values=c("#C1272D", "#2166AC")) +
  labs(title="O.sagittarius RNAseq G",x = "PC1 (33.2%)", y = "PC2 (18.2%)") +
  theme_bw()

ggsave(PCA_Os_G, file = "/Users/ericanadolski/Desktop/beetle-sex-dimorph/figures/Os_RNA_PCA_G.pdf",
       width = 4.30, height = 3.77)
##### D gazella ######
# generate DESeq dataset for exploratory analysis
Dg_dds_exp <- DESeqDataSetFromTximport(Dg_salmon_tx, Dg_RNA_info, design = ~ 1)
nrow(Dg_dds_exp) # 31857 

# filter out genes with very low read count across most samples
Dg_dds_exp <- Dg_dds_exp[rowSums(counts(Dg_dds_exp)) > 5, ]
nrow(Dg_dds_exp) # 14953

Dg_dds_rld <- rlog(Dg_dds_exp) 
Dg_dds_rld <- as.data.frame(assay(Dg_dds_rld))
Dg_dds_rld <- as.data.frame(t(Dg_dds_rld))

Dg_sampleDists <- dist(Dg_dds_rld)
Dg_sampleDistMatrix <- as.matrix(Dg_sampleDists)

### heatmap all traits ####
rownames(Dg_sampleDistMatrix) <- paste(Dg_RNA_info$Sample)
colnames(Dg_sampleDistMatrix) <- NULL
Dg_heatmap <- pheatmap(Dg_sampleDistMatrix,
                       clustering_distance_rows=Dg_sampleDists,
                       clustering_distance_cols=Dg_sampleDists,
                       col=colors)

ggsave(Dg_heatmap, file = "/Users/ericanadolski/Desktop/beetle-sex-dimorph/figures/Dg_RNA_heatmap.pdf",
       width = 8, height = 8)

### PCA all traits ####
PCA_Dg_RNA <- prcomp(Dg_dds_rld, center = TRUE)
PCs_Dg_RNA <- as.data.frame(PCA_Dg_RNA$x)
PCs_Dg_RNA$Sample <- row.names(PCs_Dg_RNA)
PCs_Dg_RNA$Sample <- Dg_RNA_info$Sample
PCs_Dg_RNA <- inner_join(PCs_Dg_RNA, Dg_RNA_info, by='Sample') 
summary(PCA_Dg_RNA)

Dg_PCA <- ggplot(data = PCs_Dg_RNA, aes(x = PC1, y = PC2, color = Sex, shape=Trait)) + 
  geom_point(size = 6, alpha = 0.8) +
  scale_color_manual(values=c("#C1272D", "#2166AC")) + 
  #geom_text(aes(label = Replicate), nudge_y = 3) +
  labs(title="D. gazella RNA-seq",x = "PC1 (30.5%)", y = "PC2 (14.7%)") +
  theme_bw(base_size=16)

ggsave(Dg_PCA, file = "/Users/ericanadolski/Desktop/beetle-sex-dimorph/figures/Dg_RNA_PCA.pdf",
       width = 5, height = 4)

### PCA posterior head ####
Dg_PH_info <- filter(Dg_RNA_info, Trait == "PH")
Dg_PH_counts <- Dg_counts %>% dplyr::select(contains("PH")) 

dds_Dg_PH_dg <- DESeqDataSetFromMatrix(countData = round(Dg_PH_counts),
                                       colData = Dg_PH_info,
                                       design = ~ 1)
rlog_dds_Dg_PH <- rlog(dds_Dg_PH_dg) 
rlog_counts_Dg_PH <- as.data.frame(assay(rlog_dds_Dg_PH))
rlog_countsT_Dg_PH <- as.data.frame(t(rlog_counts_Dg_PH))

PCA_Dg_PH <- prcomp(rlog_countsT_Dg_PH, center = TRUE)
PCs_Dg_PH <- as.data.frame(PCA_Dg_PH$x)
PCs_Dg_PH$Sample <- row.names(PCs_Dg_PH)
PCs_Dg_PH <- inner_join(PCs_Dg_PH, Dg_PH_info, by='Sample') 
summary(PCA_Dg_PH)

PCA_Dg_PH <- ggplot(data = PCs_Dg_PH, 
                    aes(x = PC1, y = PC2, color=Sex)) + 
  geom_point(size = 6, shape=16) +
  scale_color_manual(values=c("#C1272D", "#2166AC")) +
  labs(title="D.gazella RNAseq PH",x = "PC1 (36..1%)", y = "PC2 (13.9%)") +
  theme_bw()

ggsave(PCA_Dg_PH, file = "/Users/ericanadolski/Desktop/beetle-sex-dimorph/figures/Dg_RNA_PCA_PH.pdf",
       width = 4.30, height = 3.77)

### PCA anterior head ####
Dg_AH_info <- filter(Dg_RNA_info, Trait == "AH")
Dg_AH_counts <- Dg_counts %>% dplyr::select(contains("AH")) 

dds_Dg_AH_dg <- DESeqDataSetFromMatrix(countData = round(Dg_AH_counts),
                                       colData = Dg_AH_info,
                                       design = ~ 1)
rlog_dds_Dg_AH <- rlog(dds_Dg_AH_dg) 
rlog_counts_Dg_AH <- as.data.frame(assay(rlog_dds_Dg_AH))
rlog_countsT_Dg_AH <- as.data.frame(t(rlog_counts_Dg_AH))

PCA_Dg_AH <- prcomp(rlog_countsT_Dg_AH, center = TRUE)
PCs_Dg_AH <- as.data.frame(PCA_Dg_AH$x)
PCs_Dg_AH$Sample <- row.names(PCs_Dg_AH)
PCs_Dg_AH <- inner_join(PCs_Dg_AH, Dg_AH_info, by='Sample') 
summary(PCA_Dg_AH)

PCA_Dg_AH <- ggplot(data = PCs_Dg_AH, 
                    aes(x = PC1, y = PC2, color=Sex)) + 
  geom_point(size = 6, shape=16) +
  scale_color_manual(values=c("#C1272D", "#2166AC")) +
  labs(title="D.gazella RNAseq AH",x = "PC1 (32.4%)", y = "PC2 (16.3%)") +
  theme_bw()

ggsave(PCA_Dg_AH, file = "/Users/ericanadolski/Desktop/beetle-sex-dimorph/figures/Dg_RNA_PCA_AH.pdf",
       width = 4.30, height = 3.77)

### PCA genitalia ####
Dg_G_info <- filter(Dg_RNA_info, Trait == "G")
Dg_G_counts <- Dg_counts %>% dplyr::select(contains("-G")) 

dds_Dg_G_dg <- DESeqDataSetFromMatrix(countData = round(Dg_G_counts),
                                      colData = Dg_G_info,
                                      design = ~ 1)
rlog_dds_Dg_G <- rlog(dds_Dg_G_dg) 
rlog_counts_Dg_G <- as.data.frame(assay(rlog_dds_Dg_G))
rlog_countsT_Dg_G <- as.data.frame(t(rlog_counts_Dg_G))

PCA_Dg_G <- prcomp(rlog_countsT_Dg_G, center = TRUE)
PCs_Dg_G <- as.data.frame(PCA_Dg_G$x)
PCs_Dg_G$Sample <- row.names(PCs_Dg_G)
PCs_Dg_G <- inner_join(PCs_Dg_G, Dg_G_info, by='Sample') 
summary(PCA_Dg_G)

PCA_Dg_G <- ggplot(data = PCs_Dg_G, 
                   aes(x = PC1, y = PC2, color=Sex)) + 
  geom_point(size = 6, shape=16) +
  scale_color_manual(values=c("#C1272D", "#2166AC")) +
  labs(title="D.gazella RNAseq G",x = "PC1 (42.2%)", y = "PC2 (10.2%)") +
  theme_bw()

ggsave(PCA_Dg_G, file = "/Users/ericanadolski/Desktop/beetle-sex-dimorph/figures/Dg_RNA_PCA_G.pdf",
       width = 4.30, height = 3.77)
###### Step 3: analyze differential gene expression #####
### O taurus ####
# generate DESeq dataset for DGE analysis with design - grouping variable
Ot_dds <- DESeqDataSetFromTximport(Ot_salmon_tx, Ot_RNA_info, design = ~ Sex + Trait)
Ot_dds$group <- factor(paste0(Ot_dds$Sex, Ot_dds$Trait))
design(Ot_dds) <- ~ group

nrow(Ot_dds) # 19500
Ot_dds <- Ot_dds[rowSums(counts(Ot_dds)) > 5, ] # filtering data
nrow(Ot_dds) # 15408

Ot_dds <- DESeq(Ot_dds)

# save normalized counts output
Ot_norm_counts <- counts(Ot_dds, normalized=TRUE)
colnames(Ot_norm_counts) <- Ot_RNA_info$Sample

### perform LFC shrink on results 
Ot_dds_shrink <- DESeq(Ot_dds, betaPrior=TRUE)
tally(as.data.frame(results(Ot_dds_shrink, contrast=c("group","MPH","FPH"))) %>% filter(padj <= 0.1)) 
# 1425
tally(as.data.frame(results(Ot_dds_shrink, contrast=c("group","MAH","FAH"))) %>% filter(padj <= 0.1)) 
# 1717
tally(as.data.frame(results(Ot_dds_shrink, contrast=c("group","MG","FG"))) %>% filter(padj <= 0.1)) 
# 1386

# DESeq2 independent filtering=true (default) removes genes with very low counts (with low power to be detected as significant in the first place because of high dispersion) to lower total number of tests which benefits multiple testing corrections

### save results in dataframe
# posterior head
Ot_PH_MvF_gene <- as.data.frame(results(Ot_dds_shrink, contrast=c("group","MPH","FPH")))
Ot_PH_MvF_gene <- rownames_to_column(Ot_PH_MvF_gene)
names(Ot_PH_MvF_gene)[1] <- "gene"
deg_Ot_PH_MvF <- Ot_PH_MvF_gene %>% filter(padj <= 0.1)

#AH
Ot_AH_MvF_gene <- as.data.frame(results(Ot_dds_shrink, contrast=c("group","MAH","FAH")))
Ot_AH_MvF_gene <- rownames_to_column(Ot_AH_MvF_gene)
names(Ot_AH_MvF_gene)[1] <- "gene"
deg_Ot_AH_MvF <- Ot_AH_MvF_gene %>% filter(padj <= 0.1)

#G
Ot_G_MvF_gene <- as.data.frame(results(Ot_dds_shrink, contrast=c("group","MG","FG")))
Ot_G_MvF_gene <- rownames_to_column(Ot_G_MvF_gene)
names(Ot_G_MvF_gene)[1] <- "gene"
deg_Ot_G_MvF <- Ot_G_MvF_gene %>% filter(padj <= 0.1)

### O sagittarius ####
# DESeq with group variable design
Os_dds <- DESeqDataSetFromTximport(Os_salmon_tx, Os_RNA_info, design = ~ Sex + Trait)
Os_dds$group <- factor(paste0(Os_dds$Sex, Os_dds$Trait))
design(Os_dds) <- ~ group

# filtering data
nrow(Os_dds) # 20609
Os_dds <- Os_dds[rowSums(counts(Os_dds)) > 5, ]
nrow(Os_dds) # 15056

# save normalized counts output
Os_dds <- DESeq(Os_dds)
Os_norm_counts <- counts(Os_dds, normalized=TRUE)
colnames(Os_norm_counts) <- Os_RNA_info$Sample

#### LFC shrink
Os_dds_shrink <- DESeq(Os_dds, betaPrior=TRUE)
tally(as.data.frame(results(Os_dds_shrink, contrast=c("group","MPH","FPH"))) %>% filter(padj <= 0.1)) # 481
tally(as.data.frame(results(Os_dds_shrink, contrast=c("group","MAH","FAH"))) %>% filter(padj <= 0.1)) # 215
tally(as.data.frame(results(Os_dds_shrink, contrast=c("group","MG","FG"))) %>% filter(padj <= 0.1)) # 1290

# save results in dataframe
# PH
Os_PH_MvF_gene <- as.data.frame(results(Os_dds_shrink, contrast=c("group","MPH","FPH")))
Os_PH_MvF_gene <- rownames_to_column(Os_PH_MvF_gene)
names(Os_PH_MvF_gene)[1] <- "gene"
deg_Os_PH_MvF <- Os_PH_MvF_gene %>% filter(padj <= 0.1)

#AH
Os_AH_MvF_gene <- as.data.frame(results(Os_dds_shrink, contrast=c("group","MAH","FAH")))
Os_AH_MvF_gene <- rownames_to_column(Os_AH_MvF_gene)
names(Os_AH_MvF_gene)[1] <- "gene"
deg_Os_AH_MvF <- Os_AH_MvF_gene  %>% filter(padj <= 0.1)

#G
Os_G_MvF_gene <- as.data.frame(results(Os_dds_shrink, contrast=c("group","MG","FG")))
Os_G_MvF_gene <- rownames_to_column(Os_G_MvF_gene)
names(Os_G_MvF_gene)[1] <- "gene"
deg_Os_G_MvF <- Os_G_MvF_gene %>% filter(padj <= 0.1)

### D gazella ####
# DESeq with group variable design
Dg_dds <- DESeqDataSetFromTximport(Dg_salmon_tx, Dg_RNA_info, design = ~ Sex + Trait)
Dg_dds$group <- factor(paste0(Dg_dds$Sex, Dg_dds$Trait))
design(Dg_dds) <- ~ group
# filtering data
nrow(Dg_dds) # 31857
Dg_dds <- Dg_dds[rowSums(counts(Dg_dds)) > 5, ]
nrow(Dg_dds) # 14953

Dg_dds <- DESeq(Dg_dds)

# save normalized counts output 
Dg_norm_counts <- counts(Dg_dds, normalized=TRUE)
colnames(Dg_norm_counts) <- Dg_RNA_info$Sample

#### LFC shrink 
Dg_dds_shrink <- DESeq(Dg_dds, betaPrior=TRUE)
tally(as.data.frame(results(Dg_dds_shrink, contrast=c("group","MPH","FPH"))) %>% filter(padj <= 0.1)) # 177
tally(as.data.frame(results(Dg_dds_shrink, contrast=c("group","MAH","FAH"))) %>% filter(padj <= 0.1)) # 34
tally(as.data.frame(results(Dg_dds_shrink, contrast=c("group","MG","FG"))) %>% filter(padj <= 0.1)) # 1188

# save results in dataframes 
# PH
Dg_PH_MvF_gene <- as.data.frame(results(Dg_dds_shrink, contrast=c("group","MPH","FPH")))
Dg_PH_MvF_gene <- rownames_to_column(Dg_PH_MvF_gene)
names(Dg_PH_MvF_gene)[1] <- "gene"
deg_Dg_PH_MvF <- Dg_PH_MvF_gene %>% filter(padj <= 0.1)

#AH
Dg_AH_MvF_gene <- as.data.frame(results(Dg_dds_shrink, contrast=c("group","MAH","FAH")))
Dg_AH_MvF_gene <- rownames_to_column(Dg_AH_MvF_gene)
names(Dg_AH_MvF_gene)[1] <- "gene"
deg_Dg_AH_MvF <- Dg_AH_MvF_gene %>% filter(padj <= 0.1)

#G
Dg_G_MvF_gene <- as.data.frame(results(Dg_dds_shrink, contrast=c("group","MG","FG")))
Dg_G_MvF_gene <- rownames_to_column(Dg_G_MvF_gene)
names(Dg_G_MvF_gene)[1] <- "gene"
deg_Dg_G_MvF <- Dg_G_MvF_gene %>% filter(padj <= 0.1)

###### Step 4: plotting differential gene expression results #####
### O taurus ####
### volcano plots ######
cols <- c("F" = "#C1272D", "M" = "#2166AC", "ns" = "grey") 
sizes <- c("F" = 2, "M" = 2, "ns" = 1) 
alphas <- c("F" = 0.7, "M" = 0.7, "ns" = 1)

### posterior head 
# annotate DEG groups 
Ot_PH_MvF_gene <- Ot_PH_MvF_gene %>% 
  mutate(upreg = ifelse(padj <= 0.1 & log2FoldChange < 0, "F",
                        ifelse(padj <= 0.1 & log2FoldChange > 0, "M",
                               "ns")))

Ot_PH_MvF_gene %>%
  dplyr::count(upreg)

Ot_PH_vol_plot <- Ot_PH_MvF_gene %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), fill = upreg, size = upreg, alpha = upreg)) +
  geom_point(shape = 21, # Specify shape and colour as fixed local parameters    
             colour = "black", stroke = 0.25) +
  scale_fill_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  ylim(0,12) +
  xlim(-5, 5) +
  theme_bw(base_size=14) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth= 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  labs(title="Otau posterior head")
Ot_PH_vol_plot

ggsave(Ot_PH_vol_plot, file = "/Users/ericanadolski/Desktop/beetle-sex-dimorph/figures/RNA_Ot_PH_volcano.pdf",
       width = 5, height = 4)

### anterior head 
Ot_AH_MvF_gene <- Ot_AH_MvF_gene %>% 
  mutate(upreg = ifelse(padj <= 0.1 & log2FoldChange < 0, "F",
                        ifelse(padj <= 0.1 & log2FoldChange > 0, "M",
                               "ns")))
Ot_AH_MvF_gene %>%
  dplyr::count(upreg)

Ot_AH_vol_plot <- Ot_AH_MvF_gene %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), fill = upreg, size = upreg, alpha = upreg)) +
  geom_point(shape = 22, 
             colour = "black", stroke = 0.25) +
  scale_fill_manual(values = cols) + 
  scale_size_manual(values = sizes) +
  scale_alpha_manual(values = alphas) + 
  ylim(0,12) +
  xlim(-5, 5) +
  theme_bw(base_size=14) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth= 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  labs(title="Otau anterior head")
Ot_AH_vol_plot

ggsave(Ot_AH_vol_plot, file = "/Users/ericanadolski/Desktop/beetle-sex-dimorph/figures/RNA_Ot_AH_volcano.pdf",
       width = 5, height = 4)

### genitalia 
Ot_G_MvF_gene <- Ot_G_MvF_gene %>% 
  mutate(upreg = ifelse(padj <= 0.1 & log2FoldChange < 0, "F",
                        ifelse(padj <= 0.1 & log2FoldChange > 0, "M",
                               "ns")))

Ot_G_MvF_gene %>%
  dplyr::count(upreg)

Ot_G_vol_plot <- Ot_G_MvF_gene %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), fill = upreg, size = upreg, alpha = upreg)) +
  geom_point(shape = 24, # Specify shape and colour as fixed local parameters    
             colour = "black", stroke = 0.25) +
  scale_fill_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  ylim(0,12) +
  xlim(-5, 5) +
  theme_bw(base_size=14) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth= 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  labs(title="Otau genitalia")
Ot_G_vol_plot

ggsave(Ot_G_vol_plot, file = "/Users/ericanadolski/Desktop/beetle-sex-dimorph/figures/RNA_Ot_G_volcano.pdf",
       width = 5, height = 4)


### clustered heatmap ####
myheatcolors <- rev(brewer.pal(name="RdBu", n=11))
### posterior head 
# subset count matrix by DEGs and trait 
Ot_PH_norm_counts <- dplyr::select(as.data.frame(Ot_norm_counts), contains("PH"))
Ot_PH_norm_counts <- rownames_to_column(Ot_PH_norm_counts)
names(Ot_PH_norm_counts)[1] <- "gene"
Ot_PH_deg <- inner_join(deg_Ot_PH_MvF, Ot_PH_norm_counts, by="gene")

clustRows <- hclust(as.dist(1-cor(t(Ot_PH_deg[,8:19]), method="pearson")), method="complete")
clustColumns <- hclust(as.dist(1-cor(Ot_PH_deg[,8:19], method="spearman")), method="complete")
module.assign <- cutree(clustRows, k=2)
module.color <- rainbow(length(unique(module.assign)), start=0.1, end=0.9) 
module.color <- module.color[as.vector(module.assign)] 
heatmap.2(as.matrix(Ot_PH_deg[,8:19]), 
          Rowv=as.dendrogram(clustRows), 
          Colv=as.dendrogram(clustColumns),
          RowSideColors=module.color,
          col=myheatcolors, scale='row', labRow=NA,
          density.info="none", trace="none",)

### anterior head 
Ot_AH_norm_counts <- dplyr::select(as.data.frame(Ot_norm_counts), contains("AH"))
Ot_AH_norm_counts <- rownames_to_column(Ot_AH_norm_counts)
names(Ot_AH_norm_counts)[1] <- "gene"
Ot_AH_deg <- inner_join(deg_Ot_AH_MvF, Ot_AH_norm_counts, by="gene")

clustRows <- hclust(as.dist(1-cor(t(Ot_AH_deg[,8:19]), method="pearson")), method="complete")
clustColumns <- hclust(as.dist(1-cor(Ot_AH_deg[,8:19], method="spearman")), method="complete")
module.assign <- cutree(clustRows, k=2)
module.color <- rainbow(length(unique(module.assign)), start=0.1, end=0.9) 
module.color <- module.color[as.vector(module.assign)] 
heatmap.2(as.matrix(Ot_AH_deg[,8:19]), 
          Rowv=as.dendrogram(clustRows), 
          Colv=as.dendrogram(clustColumns),
          RowSideColors=module.color,
          col=myheatcolors, scale='row', labRow=NA,
          density.info="none", trace="none",)

### genitalia
# subset count matrix by DEGs and trait 
Ot_G_norm_counts <- dplyr::select(as.data.frame(Ot_norm_counts), contains("G"))
Ot_G_norm_counts <- rownames_to_column(Ot_G_norm_counts)
names(Ot_G_norm_counts)[1] <- "gene"
Ot_G_deg <- inner_join(deg_Ot_G_MvF, Ot_G_norm_counts, by="gene")

clustRows <- hclust(as.dist(1-cor(t(Ot_G_deg[,8:19]), method="pearson")), method="complete")
clustColumns <- hclust(as.dist(1-cor(Ot_G_deg[,8:19], method="spearman")), method="complete")
module.assign <- cutree(clustRows, k=2)
module.color <- rainbow(length(unique(module.assign)), start=0.1, end=0.9) 
module.color <- module.color[as.vector(module.assign)] 
heatmap.2(as.matrix(Ot_G_deg[,8:19]), 
          Rowv=as.dendrogram(clustRows), 
          Colv=as.dendrogram(clustColumns),
          RowSideColors=module.color,
          col=myheatcolors, scale='row', labRow=NA,
          density.info="none", trace="none",)

### O sagittarius ####
### volcano plots ####
# annotate DEG groups 
Os_PH_MvF_gene <- Os_PH_MvF_gene %>% 
  mutate(upreg = ifelse(padj <= 0.1 & log2FoldChange < 0, "F",
                        ifelse(padj <= 0.1 & log2FoldChange > 0, "M",
                               "ns")))
Os_PH_MvF_gene %>%
  dplyr::count(upreg)

Os_PH_vol_plot <- Os_PH_MvF_gene %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), fill = upreg, size = upreg, alpha = upreg)) +
  geom_point(shape = 21, # Specify shape and colour as fixed local parameters    
             colour = "black", stroke = 0.25) +
  scale_fill_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  ylim(0,12) +
  xlim(-5, 5) +
  theme_bw(base_size=14) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth= 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  labs(title="Osag posterior head")
Os_PH_vol_plot

ggsave(Os_PH_vol_plot, file = "/Users/ericanadolski/Desktop/beetle-sex-dimorph/figures/RNA_Os_PH_volcano.pdf",
       width = 5, height = 4)

## AH
Os_AH_MvF_gene <- Os_AH_MvF_gene %>% 
  mutate(upreg = ifelse(padj <= 0.1 & log2FoldChange < 0, "F",
                        ifelse(padj <= 0.1 & log2FoldChange > 0, "M",
                               "ns")))
Os_AH_MvF_gene %>%
  dplyr::count(upreg)

Os_AH_vol_plot <- Os_AH_MvF_gene %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), fill = upreg, size = upreg, alpha = upreg)) +
  geom_point(shape = 22, # Specify shape and colour as fixed local parameters    
             colour = "black", stroke = 0.25) +
  scale_fill_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  ylim(0,12) +
  xlim(-5, 5) +
  theme_bw(base_size=14) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth= 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  labs(title="Osag anterior head")
Os_AH_vol_plot

ggsave(Os_AH_vol_plot, file = "/Users/ericanadolski/Desktop/beetle-sex-dimorph/figures/RNA_Os_AH_volcano.pdf",
       width = 5, height = 4)

## G
Os_G_MvF_gene <- Os_G_MvF_gene %>% 
  mutate(upreg = ifelse(padj <= 0.1 & log2FoldChange < 0, "F",
                        ifelse(padj <= 0.1 & log2FoldChange > 0, "M",
                               "ns")))
Os_G_MvF_gene %>%
  dplyr::count(upreg)

Os_G_vol_plot <- Os_G_MvF_gene %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), fill = upreg, size = upreg, alpha = upreg)) +
  geom_point(shape = 24, # Specify shape and colour as fixed local parameters    
             colour = "black", stroke = 0.25) +
  scale_fill_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  ylim(0,12) + 
  xlim(-5, 5) +
  theme_bw(base_size=14) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth= 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  labs(title="Osag genitalia")
Os_G_vol_plot

ggsave(Os_G_vol_plot, file = "/Users/ericanadolski/Desktop/beetle-sex-dimorph/figures/RNA_Os_G_volcano.pdf",
       width = 5, height = 4)

### clustered heatmap ####
### posterior head 
# subset count matrix by DEGs and trait 
Os_PH_norm_counts <- dplyr::select(as.data.frame(Os_norm_counts), contains("PH"))
Os_PH_norm_counts <- rownames_to_column(Os_PH_norm_counts)
names(Os_PH_norm_counts)[1] <- "gene"
Os_PH_deg <- inner_join(deg_Os_PH_MvF, Os_PH_norm_counts, by="gene")

clustRows <- hclust(as.dist(1-cor(t(Os_PH_deg[,8:19]), method="pearson")), method="complete")
clustColumns <- hclust(as.dist(1-cor(Os_PH_deg[,8:19], method="spearman")), method="complete")
module.assign <- cutree(clustRows, k=2)
module.color <- rainbow(length(unique(module.assign)), start=0.1, end=0.9) 
module.color <- module.color[as.vector(module.assign)] 
heatmap.2(as.matrix(Os_PH_deg[,8:19]), 
          Rowv=as.dendrogram(clustRows), 
          Colv=as.dendrogram(clustColumns),
          RowSideColors=module.color,
          col=myheatcolors, scale='row', labRow=NA,
          density.info="none", trace="none",)

### anterior head 
Os_AH_norm_counts <- dplyr::select(as.data.frame(Os_norm_counts), contains("AH"))
Os_AH_norm_counts <- rownames_to_column(Os_AH_norm_counts)
names(Os_AH_norm_counts)[1] <- "gene"
Os_AH_deg <- inner_join(deg_Os_AH_MvF, Os_AH_norm_counts, by="gene")

clustRows <- hclust(as.dist(1-cor(t(Os_AH_deg[,8:19]), method="pearson")), method="complete")
clustColumns <- hclust(as.dist(1-cor(Os_AH_deg[,8:19], method="spearman")), method="complete")
module.assign <- cutree(clustRows, k=2)
module.color <- rainbow(length(unique(module.assign)), start=0.1, end=0.9) 
module.color <- module.color[as.vector(module.assign)] 
heatmap.2(as.matrix(Os_AH_deg[,8:19]), 
          Rowv=as.dendrogram(clustRows), 
          Colv=as.dendrogram(clustColumns),
          RowSideColors=module.color,
          col=myheatcolors, scale='row', labRow=NA,
          density.info="none", trace="none",)

### genitalia 
# subset count matrix by DEGs and trait 
Os_G_norm_counts <- dplyr::select(as.data.frame(Os_norm_counts), contains("G"))
Os_G_norm_counts <- rownames_to_column(Os_G_norm_counts)
names(Os_G_norm_counts)[1] <- "gene"
Os_G_deg <- inner_join(deg_Os_G_MvF, Os_G_norm_counts, by="gene")

clustRows <- hclust(as.dist(1-cor(t(Os_G_deg[,8:19]), method="pearson")), method="complete")
clustColumns <- hclust(as.dist(1-cor(Os_G_deg[,8:19], method="spearman")), method="complete")
module.assign <- cutree(clustRows, k=2)
module.color <- rainbow(length(unique(module.assign)), start=0.1, end=0.9) 
module.color <- module.color[as.vector(module.assign)] 
heatmap.2(as.matrix(Os_G_deg[,8:19]), 
          Rowv=as.dendrogram(clustRows), 
          Colv=as.dendrogram(clustColumns),
          RowSideColors=module.color,
          col=myheatcolors, scale='row', labRow=NA,
          density.info="none", trace="none",)

### D gazella ####
### volcano plots ####
## PH
# annotate DEG groups 
Dg_PH_MvF_gene <- Dg_PH_MvF_gene %>% 
  mutate(upreg = ifelse(padj <= 0.1 & log2FoldChange < 0, "F",
                        ifelse(padj <= 0.1 & log2FoldChange > 0, "M",
                               "ns")))
Dg_PH_MvF_gene %>%
  dplyr::count(upreg)

Dg_PH_vol_plot <- Dg_PH_MvF_gene %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), fill = upreg, size = upreg, alpha = upreg)) +
  geom_point(shape = 21, 
             colour = "black", stroke = 0.25) +
  scale_fill_manual(values = cols) + 
  scale_size_manual(values = sizes) + 
  scale_alpha_manual(values = alphas) + 
  ylim(0,12) + 
  xlim(-5, 5) +
  theme_bw(base_size=14) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth= 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  labs(title="Dgaz posterior head")
Dg_PH_vol_plot

ggsave(Dg_PH_vol_plot, file = "/Users/ericanadolski/Desktop/beetle-sex-dimorph/figures/RNA_Dg_PH_volcano.pdf",
       width = 5, height = 4)

## AH
Dg_AH_MvF_gene <- Dg_AH_MvF_gene %>% 
  mutate(upreg = ifelse(padj <= 0.1 & log2FoldChange < 0, "F",
                        ifelse(padj <= 0.1 & log2FoldChange > 0, "M",
                               "ns")))
Dg_AH_MvF_gene %>%
  dplyr::count(upreg)

Dg_AH_vol_plot <- Dg_AH_MvF_gene %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), fill = upreg, size = upreg, alpha = upreg)) + 
  geom_point(shape = 22,  
             colour = "black", stroke = 0.25) +
  scale_fill_manual(values = cols) + 
  scale_size_manual(values = sizes) +
  scale_alpha_manual(values = alphas) + 
  ylim(0,12) + 
  xlim(-5, 5) +
  theme_bw(base_size=14) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth= 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  labs(title="Dgaz anterior head")
Dg_AH_vol_plot

ggsave(Dg_AH_vol_plot, file = "/Users/ericanadolski/Desktop/beetle-sex-dimorph/figures/RNA_Dg_AH_volcano.pdf",
       width = 5, height = 4)

## G
Dg_G_MvF_gene <- Dg_G_MvF_gene %>% 
  mutate(upreg = ifelse(padj <= 0.1 & log2FoldChange < 0, "F",
                        ifelse(padj <= 0.1 & log2FoldChange > 0, "M",
                               "ns")))
Dg_G_MvF_gene %>%
  dplyr::count(upreg)

Dg_G_vol_plot <- Dg_G_MvF_gene %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), fill = upreg, size = upreg, alpha = upreg)) +
  geom_point(shape = 24,  
             colour = "black", stroke = 0.25) +
  scale_fill_manual(values = cols) + 
  scale_size_manual(values = sizes) +
  scale_alpha_manual(values = alphas) + 
  ylim(0,12) + 
  xlim(-5, 5) +
  theme_bw(base_size=14) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth= 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  labs(title="Dgaz genitalia")
Dg_G_vol_plot

ggsave(Dg_G_vol_plot, file = "/Users/ericanadolski/Desktop/beetle-sex-dimorph/figures/RNA_Dg_G_volcano.pdf",
       width = 5, height = 4)

### clustered heatmap ####
### posterior head 
Dg_PH_norm_counts <- dplyr::select(as.data.frame(Dg_norm_counts), contains("PH"))
Dg_PH_norm_counts <- rownames_to_column(Dg_PH_norm_counts)
names(Dg_PH_norm_counts)[1] <- "gene"
Dg_PH_deg <- inner_join(deg_Dg_PH_MvF, Dg_PH_norm_counts, by="gene")

clustRows <- hclust(as.dist(1-cor(t(Dg_PH_deg[,8:19]), method="pearson")), method="complete")
clustColumns <- hclust(as.dist(1-cor(Dg_PH_deg[,8:19], method="spearman")), method="complete")
module.assign <- cutree(clustRows, k=2)
module.color <- rainbow(length(unique(module.assign)), start=0.1, end=0.9) 
module.color <- module.color[as.vector(module.assign)] 
heatmap.2(as.matrix(Dg_PH_deg[,8:19]), 
          Rowv=as.dendrogram(clustRows), 
          Colv=as.dendrogram(clustColumns),
          RowSideColors=module.color,
          col=myheatcolors, scale='row', labRow=NA,
          density.info="none", trace="none",)

### anterior head 
Dg_AH_norm_counts <- dplyr::select(as.data.frame(Dg_norm_counts), contains("AH"))
Dg_AH_norm_counts <- rownames_to_column(Dg_AH_norm_counts)
names(Dg_AH_norm_counts)[1] <- "gene"
Dg_AH_deg <- inner_join(deg_Dg_AH_MvF, Dg_AH_norm_counts, by="gene")

clustRows <- hclust(as.dist(1-cor(t(Dg_AH_deg[,8:19]), method="pearson")), method="complete")
clustColumns <- hclust(as.dist(1-cor(Dg_AH_deg[,8:19], method="spearman")), method="complete")
module.assign <- cutree(clustRows, k=2)
module.color <- rainbow(length(unique(module.assign)), start=0.1, end=0.9) 
module.color <- module.color[as.vector(module.assign)] 
heatmap.2(as.matrix(Dg_AH_deg[,8:19]), 
          Rowv=as.dendrogram(clustRows), 
          Colv=as.dendrogram(clustColumns),
          RowSideColors=module.color,
          col=myheatcolors, scale='row', labRow=NA,
          density.info="none", trace="none",)

### genitalia 
# subset count matrix by DEGs and trait 
Dg_G_norm_counts <- dplyr::select(as.data.frame(Dg_norm_counts), contains("-G"))
Dg_G_norm_counts <- rownames_to_column(Dg_G_norm_counts)
names(Dg_G_norm_counts)[1] <- "gene"
Dg_G_deg <- inner_join(deg_Dg_G_MvF, Dg_G_norm_counts, by="gene")

clustRows <- hclust(as.dist(1-cor(t(Dg_G_deg[,8:19]), method="pearson")), method="complete")
clustColumns <- hclust(as.dist(1-cor(Dg_G_deg[,8:19], method="spearman")), method="complete")
module.assign <- cutree(clustRows, k=2)
module.color <- rainbow(length(unique(module.assign)), start=0.1, end=0.9) 
module.color <- module.color[as.vector(module.assign)] 
heatmap.2(as.matrix(Dg_G_deg[,8:19]), 
          Rowv=as.dendrogram(clustRows), 
          Colv=as.dendrogram(clustColumns),
          RowSideColors=module.color,
          col=myheatcolors, scale='row', labRow=NA,
          density.info="none", trace="none",)

###### Step 5: plot bar chart of sex-responsive genes across species and traits ####
species <- c("Dgaz","Dgaz","Dgaz","Dgaz","Dgaz","Dgaz","Otau","Otau","Otau","Otau","Otau","Otau","Osag","Osag","Osag","Osag","Osag","Osag")
sex <- c('M','F', 'M','F', 'M','F', 'M','F', 'M','F', 'M','F', 'M','F', 'M','F','M','F') 
trait <- c('G','G','PH', 'PH','AH','AH','G','G','PH', 'PH','AH','AH','G','G','PH', 'PH','AH','AH')
upreg.genes <- c(568,620,100,77,23,11,448,938,648,777,771,946,431,859,300,181,179,36)

deg_numbers <- data.frame(species,sex,trait,upreg.genes)
head(deg_numbers)
deg_numbers$trait = factor(deg_numbers$trait, levels = c('G',"PH","AH"), ordered = TRUE)
deg_numbers$species = factor(deg_numbers$species, levels = c('Dgaz',"Otau","Osag"), ordered = TRUE)

bar_plot_deg <-ggplot(data=deg_numbers, aes(x=species, y=upreg.genes, fill=sex)) +
  geom_bar(stat="identity", position=position_dodge()) + 
  facet_grid(~trait) +
  ylim(0, 1000) +
  geom_text(aes(label=upreg.genes), vjust = -0.5, size=3.5, position = position_dodge(0.9))+ # outside bars
  scale_fill_manual(values=c("#C1272D","#2166AC"))+
  labs(y= "number of significantly upregulated genes")+
  theme_classic(base_size=14)
  #coord_flip()
bar_plot_deg

ggsave(bar_plot_deg, file = "/Users/ericanadolski/Desktop/beetle-sex-dimorph/figures/bar_plot_deg.pdf", width = 10, height = 5)

##### Step 6: GO term enrichment ###########
### generate gene to GO mapping by joining table of UniProt BLAST hits to UniProt GO table ######

# uniprot GO database
uniprot_go <- read.delim("/Users/ericanadolski/Documents/Sexual_dimorphism_project/GO_enrich/uniprot_reviewed_db.tsv")

# Otau
Ot_uniprot_hits <- read_delim("/Users/ericanadolski/Desktop/beetle-sex-dimorph/go_enrich/Ot-uniprot-blast-hits.txt", col_names = c("gene","Entry.Name"))
Ot_go_terms_full <- left_join(Ot_uniprot_hits, uniprot_go, by="Entry.Name") 
Ot_go_terms <- Ot_go_terms_full %>% dplyr::select(gene, Gene.Ontology.IDs)
head(Ot_go_terms)
write.table(Ot_go_terms, file = "/Users/ericanadolski/Desktop/beetle-sex-dimorph/go_enrich/Ot_go_terms.txt", sep = "\t", quote = FALSE, row.names = FALSE,)
# reformat to geneID2GO in BBEdit

# Osag
Os_uniprot_hits <- read_delim("/Users/ericanadolski/Desktop/beetle-sex-dimorph/go_enrich/Os-uniprot-blast-hits.txt", col_names = c("gene","Entry.Name"))
Os_go_terms_full <- left_join(Os_uniprot_hits, uniprot_go, by="Entry.Name") 
Os_go_terms <- Os_go_terms_full %>% dplyr::select(gene, Gene.Ontology.IDs)
head(Os_go_terms)
write.table(Os_go_terms, file = "/Users/ericanadolski/Desktop/beetle-sex-dimorph/go_enrich/Os_go_terms.txt", sep = "\t", quote = FALSE, row.names = FALSE,)

# Dgaz
Dg_uniprot_hits <- read_delim("/Users/ericanadolski/Desktop/beetle-sex-dimorph/go_enrich/Dg-uniprot-blast-hits.txt", col_names = c("gene","Entry.Name"))
Dg_go_terms_full <- left_join(Dg_uniprot_hits, uniprot_go, by="Entry.Name") 
Dg_go_terms <- Dg_go_terms_full %>% dplyr::select(gene, Gene.Ontology.IDs)
head(Dg_go_terms)
write.table(Dg_go_terms, file = "/Users/ericanadolski/Desktop/beetle-sex-dimorph/go_enrich/Dg_go_terms.txt", sep = "\t", quote = FALSE, row.names = FALSE,)

### read in custom GO mappings in geneID2GO format ####
OtgeneID2GO <- readMappings(file = "/Users/ericanadolski/Desktop/beetle-sex-dimorph/go_enrich/Ot_geneID2GO.txt")
str(head(OtgeneID2GO))

OsgeneID2GO <- readMappings(file = "/Users/ericanadolski/Desktop/beetle-sex-dimorph/go_enrich/Os_geneID2GO.txt")
str(head(OsgeneID2GO))

DggeneID2GO <- readMappings(file = "/Users/ericanadolski/Desktop/beetle-sex-dimorph/go_enrich/Dg_geneID2GO.txt")
str(head(DggeneID2GO))

##### run O taurus tests ####
### Ot PH female upreg #####
OtPHF_upreg <- Ot_PH_MvF_gene %>% 
  mutate(F.upreg = ifelse(padj <= 0.1 & log2FoldChange < 0, "1",
                          ifelse(padj <= 0.1 & log2FoldChange > 0, "0",
                                 "0")))
head(OtPHF_upreg)
OtPHF_geneList <- as.factor(OtPHF_upreg$F.upreg)
names(OtPHF_geneList) <- OtPHF_upreg$gene
length(OtPHF_geneList) # total number of genes
summary(OtPHF_geneList) # gives # of genes of interest

# build GOdata object
OtPHF_GOdata <- new("topGOdata",
                    description = "GO analysis of Otau PH Female upreg genes",
                    ontology = "BP", 
                    allGenes = OtPHF_geneList, 
                    annot = annFUN.gene2GO, 
                    gene2GO = OtgeneID2GO, 
                    nodeSize = 5)

# run test for significance using the weight01 algorithm (default) with fisher 
OtPHF_weight_fisher_result=runTest(OtPHF_GOdata, algorithm='weight01', statistic='fisher') 

# generate table of results
OtPHF_allGO=usedGO(OtPHF_GOdata)
OtPHF_all_res=GenTable(OtPHF_GOdata, weightFisher=OtPHF_weight_fisher_result, orderBy='weightFisher', topNodes=length(OtPHF_allGO))
OtPHF_all_res$weightFisher <- as.numeric(OtPHF_all_res$weightFisher)

#get list of significant GO terms
OtPHF_results.table.p= OtPHF_all_res[which(OtPHF_all_res$weightFisher<=0.001),]
dim(OtPHF_results.table.p) # 11

### Ot PH male upreg ####
OtPHM_upreg <- Ot_PH_MvF_gene %>% 
  mutate(M.upreg = ifelse(padj <= 0.1 & log2FoldChange < 0, "0", ifelse(padj <= 0.1 & log2FoldChange > 0, "1", "0")))

OtPHM_geneList <- as.factor(OtPHM_upreg$M.upreg)
names(OtPHM_geneList) <- OtPHM_upreg$gene
OtPHM_GOdata <- new("topGOdata", description = "GO analysis of Otau PH Male upreg genes",
                    ontology = "BP", allGenes = OtPHM_geneList, annot = annFUN.gene2GO, gene2GO = OtgeneID2GO, nodeSize = 5)

OtPHM_weight_fisher_result=runTest(OtPHM_GOdata, algorithm='weight01', statistic='fisher') # run test for significance

OtPHM_allGO=usedGO(OtPHM_GOdata) # generate table of results
OtPHM_all_res=GenTable(OtPHM_GOdata, weightFisher=OtPHM_weight_fisher_result, orderBy='weightFisher', topNodes=length(OtPHM_allGO))
OtPHM_all_res$weightFisher <- as.numeric(OtPHM_all_res$weightFisher)

OtPHM_results.table.p= OtPHM_all_res[which(OtPHM_all_res$weightFisher<=0.001),] #get list of significant GO 
dim(OtPHM_results.table.p) # 17

### Ot AH female upreg ####
OtAHF_upreg <- Ot_AH_MvF_gene %>% 
  mutate(F.upreg = ifelse(padj <= 0.1 & log2FoldChange < 0, "1", ifelse(padj <= 0.1 & log2FoldChange > 0, "0","0")))

OtAHF_geneList <- as.factor(OtAHF_upreg$F.upreg)
names(OtAHF_geneList) <- OtAHF_upreg$gene
summary(OtAHF_geneList)

# build GOdata object
OtAHF_GOdata <- new("topGOdata", description = "GO analysis of Otau AH Female upreg genes",
                    ontology = "BP", allGenes = OtAHF_geneList, annot = annFUN.gene2GO, gene2GO = OtgeneID2GO, nodeSize = 5)
OtAHF_weight_fisher_result=runTest(OtAHF_GOdata, algorithm='weight01', statistic='fisher') # run test for significance
OtAHF_allGO=usedGO(OtAHF_GOdata) # generate table of results
OtAHF_all_res=GenTable(OtAHF_GOdata, weightFisher=OtAHF_weight_fisher_result, orderBy='weightFisher', topNodes=length(OtAHF_allGO))
OtAHF_all_res$weightFisher <- as.numeric(OtAHF_all_res$weightFisher)
OtAHF_results.table.p= OtAHF_all_res[which(OtAHF_all_res$weightFisher<=0.001),] #get list of significant GO terms
dim(OtAHF_results.table.p) # 43

### Ot AH male upreg ####
OtAHM_upreg <- Ot_AH_MvF_gene %>% 
  mutate(M.upreg = ifelse(padj <= 0.1 & log2FoldChange < 0, "0", ifelse(padj <= 0.1 & log2FoldChange > 0, "1", "0")))
head(OtAHM_upreg)
OtAHM_geneList <- as.factor(OtAHM_upreg$M.upreg)
names(OtAHM_geneList) <- OtAHM_upreg$gene
OtAHM_GOdata <- new("topGOdata", description = "GO analysis of Otau AH Male upreg genes",
                    ontology = "BP", allGenes = OtAHM_geneList, annot = annFUN.gene2GO, gene2GO = OtgeneID2GO, nodeSize = 5)
OtAHM_weight_fisher_result=runTest(OtAHM_GOdata, algorithm='weight01', statistic='fisher') # run test for significance
OtAHM_allGO=usedGO(OtAHM_GOdata) # generate table of results
OtAHM_all_res=GenTable(OtAHM_GOdata, weightFisher=OtAHM_weight_fisher_result, orderBy='weightFisher', topNodes=length(OtAHM_allGO))
OtAHM_all_res$weightFisher <- as.numeric(OtAHM_all_res$weightFisher)
OtAHM_results.table.p= OtAHM_all_res[which(OtAHM_all_res$weightFisher<=0.001),] #get list of significant GO 
dim(OtAHM_results.table.p) # 11

### Ot G female upreg #####
OtGF_upreg <- Ot_G_MvF_gene %>% 
  mutate(F.upreg = ifelse(padj <= 0.1 & log2FoldChange < 0, "1", ifelse(padj <= 0.1 & log2FoldChange > 0, "0", "0")))

OtGF_geneList <- as.factor(OtGF_upreg$F.upreg)
names(OtGF_geneList) <- OtGF_upreg$gene
summary(OtGF_geneList)

OtGF_GOdata <- new("topGOdata", description = "GO analysis of Otau G Female upreg genes",
                   ontology = "BP", allGenes = OtGF_geneList, annot = annFUN.gene2GO, gene2GO = OtgeneID2GO, nodeSize = 5)
OtGF_weight_fisher_result=runTest(OtGF_GOdata, algorithm='weight01', statistic='fisher') # run test for significance 
OtGF_allGO=usedGO(OtGF_GOdata) # generate table of results
OtGF_all_res=GenTable(OtGF_GOdata, weightFisher=OtGF_weight_fisher_result, orderBy='weightFisher', topNodes=length(OtGF_allGO))
OtGF_all_res$weightFisher <- as.numeric(OtGF_all_res$weightFisher)
OtGF_results.table.p= OtGF_all_res[which(OtGF_all_res$weightFisher<=0.001),] #get list of significant GO terms
dim(OtGF_results.table.p) # 13

### Ot G male upreg ####
OtGM_upreg <- Ot_G_MvF_gene %>% 
  mutate(M.upreg = ifelse(padj <= 0.1 & log2FoldChange < 0, "0", ifelse(padj <= 0.1 & log2FoldChange > 0, "1", "0")))

OtGM_geneList <- as.factor(OtGM_upreg$M.upreg)
names(OtGM_geneList) <- OtGM_upreg$gene
summary(OtGM_geneList) 

OtGM_GOdata <- new("topGOdata", description = "GO analysis of Otau G Male upreg genes",
                   ontology = "BP", allGenes = OtGM_geneList, annot = annFUN.gene2GO, gene2GO = OtgeneID2GO, nodeSize = 5)
OtGM_weight_fisher_result=runTest(OtGM_GOdata, algorithm='weight01', statistic='fisher') # run test for significance
OtGM_allGO=usedGO(OtGM_GOdata) # generate table of results
OtGM_all_res=GenTable(OtGM_GOdata, weightFisher=OtGM_weight_fisher_result, orderBy='weightFisher', topNodes=length(OtGM_allGO))
OtGM_all_res$weightFisher <- as.numeric(OtGM_all_res$weightFisher)
OtGM_results.table.p= OtGM_all_res[which(OtGM_all_res$weightFisher<=0.001),] #get list of significant GO 
dim(OtGM_results.table.p) # 9

##### run O sagittarius tests ####
### Os PH female upreg #####
OsPHF_upreg <- Os_PH_MvF_gene %>% 
  mutate(F.upreg = ifelse(padj <= 0.1 & log2FoldChange < 0, "1",
                          ifelse(padj <= 0.1 & log2FoldChange > 0, "0",
                                 "0")))

OsPHF_geneList <- as.factor(OsPHF_upreg$F.upreg)
names(OsPHF_geneList) <- OsPHF_upreg$gene
summary(OsPHF_geneList) # gives # of genes of interest

# build GOdata object
OsPHF_GOdata <- new("topGOdata", description = "GO analysis of Osag PH Female upreg genes",
                    ontology = "BP", allGenes = OsPHF_geneList, annot = annFUN.gene2GO, gene2GO = OsgeneID2GO, nodeSize = 5)

# run test for significance using the weight01 algorithm (default) with fisher 
OsPHF_weight_fisher_result=runTest(OsPHF_GOdata, algorithm='weight01', statistic='fisher') 

# generate table of results
OsPHF_allGO=usedGO(OsPHF_GOdata)
OsPHF_all_res=GenTable(OsPHF_GOdata, weightFisher=OsPHF_weight_fisher_result, orderBy='weightFisher', topNodes=length(OsPHF_allGO))
OsPHF_all_res$weightFisher <- as.numeric(OsPHF_all_res$weightFisher)

#get list of significant GO terms
OsPHF_results.table.p= OsPHF_all_res[which(OsPHF_all_res$weightFisher<=0.001),]
dim(OsPHF_results.table.p) # 23

### Os PH male upreg ####
OsPHM_upreg <- Os_PH_MvF_gene %>% 
  mutate(M.upreg = ifelse(padj <= 0.1 & log2FoldChange < 0, "0", ifelse(padj <= 0.1 & log2FoldChange > 0, "1", "0")))

OsPHM_geneList <- as.factor(OsPHM_upreg$M.upreg)
names(OsPHM_geneList) <- OsPHM_upreg$gene
OsPHM_GOdata <- new("topGOdata", description = "GO analysis of Osag PH Male upreg genes",
                    ontology = "BP", allGenes = OsPHM_geneList, annot = annFUN.gene2GO, gene2GO = OsgeneID2GO, nodeSize = 5)

OsPHM_weight_fisher_result=runTest(OsPHM_GOdata, algorithm='weight01', statistic='fisher') # run test for significance
OsPHM_allGO=usedGO(OsPHM_GOdata) # generate table of results
OsPHM_all_res=GenTable(OsPHM_GOdata, weightFisher=OsPHM_weight_fisher_result, orderBy='weightFisher', topNodes=length(OsPHM_allGO))
OsPHM_all_res$weightFisher <- as.numeric(OsPHM_all_res$weightFisher)

OsPHM_results.table.p= OsPHM_all_res[which(OsPHM_all_res$weightFisher<=0.001),] #get list of significant GO 
dim(OsPHM_results.table.p) # 2

### Os AH female upreg ####
OsAHF_upreg <- Os_AH_MvF_gene %>% 
  mutate(F.upreg = ifelse(padj <= 0.1 & log2FoldChange < 0, "1", ifelse(padj <= 0.1 & log2FoldChange > 0, "0","0")))

OsAHF_geneList <- as.factor(OsAHF_upreg$F.upreg)
names(OsAHF_geneList) <- OsAHF_upreg$gene
summary(OsAHF_geneList)

# build GOdata object
OsAHF_GOdata <- new("topGOdata", description = "GO analysis of Osag AH Female upreg genes",
                    ontology = "BP", allGenes = OsAHF_geneList, annot = annFUN.gene2GO, gene2GO = OsgeneID2GO, nodeSize = 5)
OsAHF_weight_fisher_result=runTest(OsAHF_GOdata, algorithm='weight01', statistic='fisher') # run test for significance
OsAHF_allGO=usedGO(OsAHF_GOdata) # generate table of results
OsAHF_all_res=GenTable(OsAHF_GOdata, weightFisher=OsAHF_weight_fisher_result, orderBy='weightFisher', topNodes=length(OsAHF_allGO))
OsAHF_all_res$weightFisher <- as.numeric(OsAHF_all_res$weightFisher)
OsAHF_results.table.p= OsAHF_all_res[which(OsAHF_all_res$weightFisher<=0.001),] #get list of significant GO terms
dim(OsAHF_results.table.p) # 4

### Os AH male upreg ####
OsAHM_upreg <- Os_AH_MvF_gene %>% 
  mutate(M.upreg = ifelse(padj <= 0.1 & log2FoldChange < 0, "0", ifelse(padj <= 0.1 & log2FoldChange > 0, "1", "0")))
OsAHM_geneList <- as.factor(OsAHM_upreg$M.upreg)
names(OsAHM_geneList) <- OsAHM_upreg$gene

OsAHM_GOdata <- new("topGOdata", description = "GO analysis of Osag AH Male upreg genes",
                    ontology = "BP", allGenes = OsAHM_geneList, annot = annFUN.gene2GO, gene2GO = OsgeneID2GO, nodeSize = 5)
OsAHM_weight_fisher_result=runTest(OsAHM_GOdata, algorithm='weight01', statistic='fisher') # run test for significance
OsAHM_allGO=usedGO(OsAHM_GOdata) # generate table of results
OsAHM_all_res=GenTable(OsAHM_GOdata, weightFisher=OsAHM_weight_fisher_result, orderBy='weightFisher', topNodes=length(OsAHM_allGO))
OsAHM_all_res$weightFisher <- as.numeric(OsAHM_all_res$weightFisher)
OsAHM_results.table.p= OsAHM_all_res[which(OsAHM_all_res$weightFisher<=0.001),] #get list of significant GO 
dim(OsAHM_results.table.p) # 0

OsAHM_results.table.p0.05= OsAHM_all_res[which(OsAHM_all_res$weightFisher<=0.05),]
dim(OsAHM_results.table.p0.05) # 49

### Os G female upreg #####
OsGF_upreg <- Os_G_MvF_gene %>% 
  mutate(F.upreg = ifelse(padj <= 0.1 & log2FoldChange < 0, "1", ifelse(padj <= 0.1 & log2FoldChange > 0, "0", "0")))

OsGF_geneList <- as.factor(OsGF_upreg$F.upreg)
names(OsGF_geneList) <- OsGF_upreg$gene
summary(OsGF_geneList)

OsGF_GOdata <- new("topGOdata", description = "GO analysis of Osag G Female upreg genes",
                   ontology = "BP", allGenes = OsGF_geneList, annot = annFUN.gene2GO, gene2GO = OsgeneID2GO, nodeSize = 5)
OsGF_weight_fisher_result=runTest(OsGF_GOdata, algorithm='weight01', statistic='fisher') # run test for significance 
OsGF_allGO=usedGO(OsGF_GOdata) # generate table of results
OsGF_all_res=GenTable(OsGF_GOdata, weightFisher=OsGF_weight_fisher_result, orderBy='weightFisher', topNodes=length(OsGF_allGO))
OsGF_all_res$weightFisher <- as.numeric(OsGF_all_res$weightFisher)
OsGF_results.table.p= OsGF_all_res[which(OsGF_all_res$weightFisher<=0.001),] #get list of significant GO terms
dim(OsGF_results.table.p) # 8

### Os G male upreg ####
OsGM_upreg <- Os_G_MvF_gene %>% 
  mutate(M.upreg = ifelse(padj <= 0.1 & log2FoldChange < 0, "0", ifelse(padj <= 0.1 & log2FoldChange > 0, "1", "0")))

OsGM_geneList <- as.factor(OsGM_upreg$M.upreg)
names(OsGM_geneList) <- OsGM_upreg$gene
summary(OsGM_geneList) 

OsGM_GOdata <- new("topGOdata", description = "GO analysis of Osag G Male upreg genes",
                   ontology = "BP", allGenes = OsGM_geneList, annot = annFUN.gene2GO, gene2GO = OsgeneID2GO, nodeSize = 5)
OsGM_weight_fisher_result=runTest(OsGM_GOdata, algorithm='weight01', statistic='fisher') # run test for significance
OsGM_allGO=usedGO(OsGM_GOdata) # generate table of results
OsGM_all_res=GenTable(OsGM_GOdata, weightFisher=OsGM_weight_fisher_result, orderBy='weightFisher', topNodes=length(OsGM_allGO))
OsGM_all_res$weightFisher <- as.numeric(OsGM_all_res$weightFisher)
OsGM_results.table.p= OsGM_all_res[which(OsGM_all_res$weightFisher<=0.001),] #get list of significant GO 
dim(OsGM_results.table.p) # 23

##### run D gazella tests ####
### Dg PH female upreg #####
DgPHF_upreg <- Dg_PH_MvF_gene %>% 
  mutate(F.upreg = ifelse(padj <= 0.1 & log2FoldChange < 0, "1",
                          ifelse(padj <= 0.1 & log2FoldChange > 0, "0",
                                 "0")))
head(DgPHF_upreg)
DgPHF_geneList <- as.factor(DgPHF_upreg$F.upreg)
names(DgPHF_geneList) <- DgPHF_upreg$gene
length(DgPHF_geneList) # total number of genes
summary(DgPHF_geneList) # genes of interest

# build GOdata object
DgPHF_GOdata <- new("topGOdata", description = "GO analysis of Dgaz PH Female upreg genes", ontology = "BP", 
                    allGenes = DgPHF_geneList, annot = annFUN.gene2GO, gene2GO = DggeneID2GO, nodeSize = 5)

# run test for significance using the weight01 algorithm (default) with fisher 
DgPHF_weight_fisher_result=runTest(DgPHF_GOdata, algorithm='weight01', statistic='fisher') 

# generate table of results
DgPHF_allGO=usedGO(DgPHF_GOdata)
DgPHF_all_res=GenTable(DgPHF_GOdata, weightFisher=DgPHF_weight_fisher_result, orderBy='weightFisher', topNodes=length(DgPHF_allGO))
DgPHF_all_res$weightFisher <- as.numeric(DgPHF_all_res$weightFisher)

#get list of significant GO terms
DgPHF_results.table.p= DgPHF_all_res[which(DgPHF_all_res$weightFisher<=0.001),]
dim(DgPHF_results.table.p) # 9

### Dg PH male upreg ####
DgPHM_upreg <- Dg_PH_MvF_gene %>% 
  mutate(M.upreg = ifelse(padj <= 0.1 & log2FoldChange < 0, "0", ifelse(padj <= 0.1 & log2FoldChange > 0, "1", "0")))
DgPHM_geneList <- as.factor(DgPHM_upreg$M.upreg)
names(DgPHM_geneList) <- DgPHM_upreg$gene

DgPHM_GOdata <- new("topGOdata", description = "GO analysis of Dgaz PH Male upreg genes",
                    ontology = "BP", allGenes = DgPHM_geneList, annot = annFUN.gene2GO, gene2GO = DggeneID2GO, nodeSize = 5)
DgPHM_weight_fisher_result=runTest(DgPHM_GOdata, algorithm='weight01', statistic='fisher') # run test for significance
DgPHM_allGO=usedGO(DgPHM_GOdata) # generate table of results
DgPHM_all_res=GenTable(DgPHM_GOdata, weightFisher=DgPHM_weight_fisher_result, orderBy='weightFisher', topNodes=length(DgPHM_allGO))
DgPHM_all_res$weightFisher <- as.numeric(DgPHM_all_res$weightFisher)
DgPHM_results.table.p= DgPHM_all_res[which(DgPHM_all_res$weightFisher<=0.001),] #get list of significant GO 
dim(DgPHM_results.table.p) # 0

DgPHM_results.table.p0.05= DgPHM_all_res[which(DgPHM_all_res$weightFisher<=0.05),]
dim(DgPHM_results.table.p0.05) # 100

### Dg AH female upreg ####
DgAHF_upreg <- Dg_AH_MvF_gene %>% 
  mutate(F.upreg = ifelse(padj <= 0.1 & log2FoldChange < 0, "1", ifelse(padj <= 0.1 & log2FoldChange > 0, "0","0")))

DgAHF_geneList <- as.factor(DgAHF_upreg$F.upreg)
names(DgAHF_geneList) <- DgAHF_upreg$gene
summary(DgAHF_geneList)

DgAHF_GOdata <- new("topGOdata", description = "GO analysis of Dgaz AH Female upreg genes",
                    ontology = "BP", allGenes = DgAHF_geneList, annot = annFUN.gene2GO, gene2GO = DggeneID2GO, nodeSize = 5)
DgAHF_weight_fisher_result=runTest(DgAHF_GOdata, algorithm='weight01', statistic='fisher') # run test for significance
DgAHF_allGO=usedGO(DgAHF_GOdata) # generate table of results
DgAHF_all_res=GenTable(DgAHF_GOdata, weightFisher=DgAHF_weight_fisher_result, orderBy='weightFisher', topNodes=length(DgAHF_allGO))
DgAHF_all_res$weightFisher <- as.numeric(DgAHF_all_res$weightFisher)
DgAHF_results.table.p= DgAHF_all_res[which(DgAHF_all_res$weightFisher<=0.001),] #get list of significant GO terms
dim(DgAHF_results.table.p) # 7

### Dg AH male upreg ####
DgAHM_upreg <- Dg_AH_MvF_gene %>% 
  mutate(M.upreg = ifelse(padj <= 0.1 & log2FoldChange < 0, "0", ifelse(padj <= 0.1 & log2FoldChange > 0, "1", "0")))
DgAHM_geneList <- as.factor(DgAHM_upreg$M.upreg)
names(DgAHM_geneList) <- DgAHM_upreg$gene

DgAHM_GOdata <- new("topGOdata", description = "GO analysis of Dgaz AH Male upreg genes",
                    ontology = "BP", allGenes = DgAHM_geneList, annot = annFUN.gene2GO, gene2GO = DggeneID2GO, nodeSize = 5)
DgAHM_weight_fisher_result=runTest(DgAHM_GOdata, algorithm='weight01', statistic='fisher') # run test for significance
DgAHM_allGO=usedGO(DgAHM_GOdata) # generate table of results
DgAHM_all_res=GenTable(DgAHM_GOdata, weightFisher=DgAHM_weight_fisher_result, orderBy='weightFisher', topNodes=length(DgAHM_allGO))
DgAHM_all_res$weightFisher <- as.numeric(DgAHM_all_res$weightFisher)
DgAHM_results.table.p= DgAHM_all_res[which(DgAHM_all_res$weightFisher<=0.001),] #get list of significant GO 
dim(DgAHM_results.table.p) # 0

DgAHM_results.table.p0.05= DgAHM_all_res[which(DgAHM_all_res$weightFisher<=0.05),] 
dim(DgAHM_results.table.p0.05) # 38

### Dg G female upreg #####
DgGF_upreg <- Dg_G_MvF_gene %>% 
  mutate(F.upreg = ifelse(padj <= 0.1 & log2FoldChange < 0, "1", ifelse(padj <= 0.1 & log2FoldChange > 0, "0", "0")))

DgGF_geneList <- as.factor(DgGF_upreg$F.upreg)
names(DgGF_geneList) <- DgGF_upreg$gene
summary(DgGF_geneList)

DgGF_GOdata <- new("topGOdata", description = "GO analysis of Dgaz G Female upreg genes",
                   ontology = "BP", allGenes = DgGF_geneList, annot = annFUN.gene2GO, gene2GO = DggeneID2GO, nodeSize = 5)
DgGF_weight_fisher_result=runTest(DgGF_GOdata, algorithm='weight01', statistic='fisher') # run test for significance 
DgGF_allGO=usedGO(DgGF_GOdata) # generate table of results
DgGF_all_res=GenTable(DgGF_GOdata, weightFisher=DgGF_weight_fisher_result, orderBy='weightFisher', topNodes=length(DgGF_allGO))
DgGF_all_res$weightFisher <- as.numeric(DgGF_all_res$weightFisher)
DgGF_results.table.p= DgGF_all_res[which(DgGF_all_res$weightFisher<=0.001),] #get list of significant GO terms
dim(DgGF_results.table.p) # 22

### Dg G male upreg ####
DgGM_upreg <- Dg_G_MvF_gene %>% 
  mutate(M.upreg = ifelse(padj <= 0.1 & log2FoldChange < 0, "0", ifelse(padj <= 0.1 & log2FoldChange > 0, "1", "0")))

DgGM_geneList <- as.factor(DgGM_upreg$M.upreg)
names(DgGM_geneList) <- DgGM_upreg$gene
summary(DgGM_geneList) 

DgGM_GOdata <- new("topGOdata", description = "GO analysis of Dgaz G Male upreg genes",
                   ontology = "BP", allGenes = DgGM_geneList, annot = annFUN.gene2GO, gene2GO = DggeneID2GO, nodeSize = 5)
DgGM_weight_fisher_result=runTest(DgGM_GOdata, algorithm='weight01', statistic='fisher') # run test for significance
DgGM_allGO=usedGO(DgGM_GOdata) # generate table of results
DgGM_all_res=GenTable(DgGM_GOdata, weightFisher=DgGM_weight_fisher_result, orderBy='weightFisher', topNodes=length(DgGM_allGO))
DgGM_all_res$weightFisher <- as.numeric(DgGM_all_res$weightFisher)
DgGM_results.table.p= DgGM_all_res[which(DgGM_all_res$weightFisher<=0.001),] #get list of significant GO 
dim(DgGM_results.table.p) # 12


###### Step 6: get gene annotations via UniProt best hits ####
### read in protein best hits to Otau2.0 genome with annotations
Otau3_prot_hits <- read.delim("/Users/ericanadolski/Documents/Genomes/Otau3/Otau_OT2_top5_hits.txt")
dim(Otau3_prot_hits) # 242025 rows
head(Otau3_prot_hits)
Otau3_prot_hits_best <- Otau3_prot_hits %>% 
  group_by(Otau_ID) %>% 
  top_n(-1,OT2_evalue) %>% 
  dplyr::slice(which.max(OT2_pident))
head(Otau3_prot_hits_best)

Dgaz1_prot_hits <- read.delim("/Users/ericanadolski/Documents/Genomes/Dgaz1/Dgaz_OT2_top5_hits.txt")
Dgaz1_prot_hits_best <- Dgaz1_prot_hits %>% 
  group_by(Dgaz_ID) %>% 
  top_n(-1,OT2_evalue) %>% 
  dplyr::slice(which.max(OT2_pident))


Osag1_prot_hits <- read.delim("/Users/ericanadolski/Documents/Genomes/Osag1/Osag_OT2_top5_hits.txt")
Osag1_prot_hits_best <- Osag1_prot_hits %>% 
  group_by(Osag_ID) %>% 
  top_n(-1,OT2_evalue) %>% 
  dplyr::slice(which.max(OT2_pident))

dim(Otau3_prot_hits_best) # 17332 rows
dim(Dgaz1_prot_hits_best) # 19672 rows
dim(Osag1_prot_hits_best) # 16935 rows

### read in UniProt annotations ####
# go terms dataframe has Protein.names annotation
uniprot_anno <- go_terms_full %>% dplyr::select(gene, Protein.names, Organism)

### combine all annotations ####
deg_Ot_PH_MvF_anno <- join_all(list(deg_Ot_PH_MvF, uniprot_anno, Otau3_prot_anno), by='gene', type='left')
deg_Ot_AH_MvF_anno <- join_all(list(deg_Ot_AH_MvF, Otau3_prot_anno, uniprot_anno), by='gene', type='left')
deg_Ot_G_MvF_anno <- join_all(list(deg_Ot_G_MvF, Otau3_prot_anno, uniprot_anno), by='gene', type='left')

deg_Os_PH_MvF_anno <- left_join(deg_Os_PH_MvF, Osag1_prot_anno, by="gene")
deg_Os_AH_MvF_anno <- left_join(deg_Os_AH_MvF, Osag1_prot_anno, by="gene")
deg_Os_G_MvF_anno <- left_join(deg_Os_G_MvF, Osag1_prot_anno, by="gene")

deg_Dg_PH_MvF_anno <- left_join(deg_Dg_PH_MvF, Dgaz1_prot_anno, by="gene")
deg_Dg_AH_MvF_anno <- left_join(deg_Dg_AH_MvF, Dgaz1_prot_anno, by="gene")
deg_Dg_G_MvF_anno <- left_join(deg_Dg_G_MvF, Dgaz1_prot_anno, by="gene")

### Step 6.5: assign orthologs across genomes ######
### import Orthofinder analyses
orthogroup_count <- read.delim("/Users/ericanadolski/Documents/Sexual_dimorphism_project/Orthofinder/prot_results_May01/Orthogroups/Orthogroups.GeneCount.tsv")

head(orthogroup_count)
nrow(orthogroup_count)

orthogroup_counts <- orthogroup_count %>%
  filter(Dgaz <= 1,
         Osag <= 1,
         Otau <= 1)
nrow(orthogroup_counts)

orthogroup_IDs_full <- read.delim("/Users/ericanadolski/Desktop/beetle-sex-dimorph/Orthogroups.tsv")

orthogroup_IDs <- orthogroup_counts %>% 
  left_join(orthogroup_IDs_full, by = "Orthogroup") %>%
  dplyr::select(Orthogroup, Dgaz.y, Osag.y, Otau.y) %>%
  mutate_all(list(~na_if(.,"")))
nrow(orthogroup_IDs)

# clean in BBEdit
write.table(orthogroup_IDs, 
            file = "/Users/ericanadolski/Desktop/beetle-sex-dimorph/orthogroup_IDs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

orthogroup_IDs_final <- read.delim("/Users/ericanadolski/Desktop/beetle-sex-dimorph/orthogroup_IDs_final.txt")

# prep DGE columns 
Ot_PH_MvF_gene <- dplyr::rename(Ot_PH_MvF_gene, OtPH_upreg = upreg)
Ot_AH_MvF_gene <- dplyr::rename(Ot_AH_MvF_gene, OtAH_upreg = upreg)
Ot_G_MvF_gene <- dplyr::rename(Ot_G_MvF_gene, OtG_upreg = upreg)

Os_PH_MvF_gene <- dplyr::rename(Os_PH_MvF_gene, OsPH_upreg = upreg)
Os_AH_MvF_gene <- dplyr::rename(Os_AH_MvF_gene, OsAH_upreg = upreg)
Os_G_MvF_gene <- dplyr::rename(Os_G_MvF_gene, OsG_upreg = upreg)

Dg_PH_MvF_gene <- dplyr::rename(Dg_PH_MvF_gene, DgPH_upreg = upreg)
Dg_AH_MvF_gene <- dplyr::rename(Dg_AH_MvF_gene, DgAH_upreg = upreg)
Dg_G_MvF_gene <- dplyr::rename(Dg_G_MvF_gene, DgG_upreg = upreg)


orthogrp_sex_res_PH <- orthogroup_IDs_final %>%
  full_join(Ot_PH_MvF_gene, by = c("Otau" = "gene")) %>%
  full_join(Os_PH_MvF_gene, by = c("Osag" = "gene")) %>%
  full_join(Dg_PH_MvF_gene, by = c("Dgaz" = "gene")) %>%
  dplyr::select(Orthogroup,Otau,Dgaz,Osag,OtPH_upreg,OsPH_upreg,DgPH_upreg)

colSums(is.na(orthogrp_sex_res_PH))

orthogrp_sex_res_AH <- orthogroup_IDs_final %>%
  full_join(Ot_AH_MvF_gene, by = c("Otau" = "gene")) %>%
  full_join(Os_AH_MvF_gene, by = c("Osag" = "gene")) %>%
  full_join(Dg_AH_MvF_gene, by = c("Dgaz" = "gene")) %>%
  dplyr::select(Orthogroup,Otau,Dgaz,Osag,OtAH_upreg,OsAH_upreg,DgAH_upreg)

orthogrp_sex_res_G <- orthogroup_IDs_final %>%
  full_join(Ot_G_MvF_gene, by = c("Otau" = "gene")) %>%
  full_join(Os_G_MvF_gene, by = c("Osag" = "gene")) %>%
  full_join(Dg_G_MvF_gene, by = c("Dgaz" = "gene")) %>%
  dplyr::select(Orthogroup,Otau,Dgaz,Osag,OtG_upreg,OsG_upreg,DgG_upreg)


write.table(orthogrp_sex_res_PH, 
            file = "/Users/ericanadolski/Desktop/beetle-sex-dimorph/orthogrp_sex_res_PH.txt", sep = "\t", quote = FALSE, row.names = FALSE)

write.table(orthogrp_sex_res_AH, 
            file = "/Users/ericanadolski/Desktop/beetle-sex-dimorph/orthogrp_sex_res_AH.txt", sep = "\t", quote = FALSE, row.names = FALSE)

write.table(orthogrp_sex_res_G, 
            file = "/Users/ericanadolski/Desktop/beetle-sex-dimorph/orthogrp_sex_res_G.txt", sep = "\t", quote = FALSE, row.names = FALSE)

### parsing out the gene set totals #####
### POSTERIOR HEAD ####
### lineage-specific ####
nrow(orthogrp_sex_res_PH %>% filter(is.na(Orthogroup)) %>% filter(is.na(Osag)) %>% filter(is.na(Otau)) %>% filter(DgPH_upreg=="F")) 
# lineage-specific Dgaz F PH upreg genes = 23

nrow(orthogrp_sex_res_PH %>% filter(is.na(Orthogroup)) %>% filter(is.na(Osag)) %>% filter(is.na(Otau)) %>% filter(DgPH_upreg=="M")) 
# lineage-specific Dgaz M PH upreg genes = 34

nrow(orthogrp_sex_res_PH %>%  filter(is.na(Orthogroup)) %>% filter(is.na(Dgaz)) %>% filter(is.na(Otau)) %>% filter(OsPH_upreg=="F")) 
# lineage-specific Osag F PH upreg genes = 43

nrow(orthogrp_sex_res_PH %>%  filter(is.na(Orthogroup)) %>% filter(is.na(Dgaz)) %>% filter(is.na(Otau)) %>% filter(OsPH_upreg=="M")) 
# lineage-specific Osag M PH upreg genes = 103

nrow(orthogrp_sex_res_PH %>%  filter(is.na(Orthogroup)) %>% filter(is.na(Dgaz)) %>% filter(is.na(Osag)) %>% filter(OtPH_upreg=="F")) 
# lineage-specific Otau F PH upreg genes = 233

nrow(orthogrp_sex_res_PH %>%  filter(is.na(Orthogroup)) %>% filter(is.na(Dgaz)) %>% filter(is.na(Osag)) %>% filter(OtPH_upreg=="M")) 
# lineage-specific Otau M PH upreg genes = 205

### shared ####
nrow(orthogrp_sex_res_PH %>% filter(OtPH_upreg=="F") %>% filter(DgPH_upreg=="F")  %>% filter(OsPH_upreg=="F"))
# all 3 shared F PH upreg genes = 3

nrow(full_sex_res_PH %>% filter(OtPH_upreg=="F") %>% filter(OsPH_upreg=="F")) 
# Otau + Osag shared F PH upreg genes = 14

nrow(full_sex_res_PH %>% filter(DgPH_upreg=="F") %>% filter(OsPH_upreg=="F")) 
# Dgaz + Osag shared F PH upreg genes = 6

nrow(full_sex_res_PH %>% filter(OtPH_upreg=="F") %>% filter(DgPH_upreg=="F")) 
# Dgaz + Otau shared F PH upreg genes = 11

### posterior horn-promoting gene sets ###
nrow(full_sex_res_PH %>% filter(OtPH_upreg=="M") %>% filter(DgPH_upreg=="M")  %>% filter(OsPH_upreg=="F")) # 1

### GENITALIA ####
### lineage-specific ####
nrow(full_sex_res_G %>% filter(is.na(Orthogroup)) %>% filter(is.na(Osag)) %>% filter(is.na(Otau)) %>% filter(DgG_upreg=="F")) 
# lineage-specific Dgaz F G upreg genes = 178

nrow(full_sex_res_G %>% filter(is.na(Orthogroup)) %>% filter(is.na(Osag)) %>% filter(is.na(Otau)) %>% filter(DgG_upreg=="M")) 
# lineage-specific Dgaz M G upreg genes = 138

nrow(full_sex_res_G %>%  filter(is.na(Orthogroup)) %>% filter(is.na(Dgaz)) %>% filter(is.na(Otau)) %>% filter(OsG_upreg=="F")) 
# lineage-specific Osag F G upreg genes = 353

nrow(full_sex_res_G %>%  filter(is.na(Orthogroup)) %>% filter(is.na(Dgaz)) %>% filter(is.na(Otau)) %>% filter(OsG_upreg=="M")) 
# lineage-specific Osag M G upreg genes = 110

nrow(full_sex_res_G %>%  filter(is.na(Orthogroup)) %>% filter(is.na(Dgaz)) %>% filter(is.na(Osag)) %>% filter(OtG_upreg=="F")) 
# lineage-specific Otau F G upreg genes = 277

nrow(full_sex_res_G %>%  filter(is.na(Orthogroup)) %>% filter(is.na(Dgaz)) %>% filter(is.na(Osag)) %>% filter(OtG_upreg=="M")) 
# lineage-specific Otau M G upreg genes = 103

### shared ####
# female
nrow(full_sex_res_G %>% filter(OtG_upreg=="F") %>% filter(DgG_upreg=="F")  %>% filter(OsG_upreg=="F")) # 15

nrow(full_sex_res_G %>% filter(OtG_upreg=="F") %>% filter(OsG_upreg=="F")) # 70

nrow(full_sex_res_G %>% filter(DgG_upreg=="F") %>% filter(OsG_upreg=="F")) # 68

nrow(full_sex_res_G %>% filter(OtG_upreg=="F") %>% filter(DgG_upreg=="F")) # 85

# male
nrow(full_sex_res_G %>% filter(OtG_upreg=="M") %>% filter(DgG_upreg=="M")  %>% filter(OsG_upreg=="M"))# all 3 shared = 18

nrow(full_sex_res_G %>% filter(OtG_upreg=="M") %>% filter(OsG_upreg=="M")) # 55

nrow(full_sex_res_G %>% filter(DgG_upreg=="M") %>% filter(OsG_upreg=="M")) # 48

nrow(full_sex_res_G %>% filter(OtG_upreg=="M") %>% filter(DgG_upreg=="M")) # 66

### ANTERIOR HEAD ####
### lineage-specific ####
nrow(full_sex_res_AH %>% filter(is.na(Orthogroup)) %>% filter(is.na(Osag)) %>% filter(is.na(Otau)) %>% filter(DgAH_upreg=="F")) 
# lineage-specific Dgaz F AH upreg genes = 34

nrow(full_sex_res_AH %>% filter(is.na(Orthogroup)) %>% filter(is.na(Osag)) %>% filter(is.na(Otau)) %>% filter(DgAH_upreg=="M")) 
# lineage-specific Dgaz M AH upreg genes = 23

nrow(full_sex_res_AH %>%  filter(is.na(Orthogroup)) %>% filter(is.na(Dgaz)) %>% filter(is.na(Otau)) %>% filter(OsAH_upreg=="F")) 
# lineage-specific Osag F AH upreg genes = 103

nrow(full_sex_res_AH %>%  filter(is.na(Orthogroup)) %>% filter(is.na(Dgaz)) %>% filter(is.na(Otau)) %>% filter(OsAH_upreg=="M")) 
# lineage-specific Osag M AH upreg genes = 43

nrow(full_sex_res_AH %>%  filter(is.na(Orthogroup)) %>% filter(is.na(Dgaz)) %>% filter(is.na(Osag)) %>% filter(OtAH_upreg=="F")) 
# lineage-specific Otau F AH upreg genes = 234

nrow(full_sex_res_AH %>%  filter(is.na(Orthogroup)) %>% filter(is.na(Dgaz)) %>% filter(is.na(Osag)) %>% filter(OtAH_upreg=="M")) 
# lineage-specific Otau M AH upreg genes = 205

### shared ####
nrow(full_sex_res_AH %>% filter(OtAH_upreg=="F") %>% filter(DgAH_upreg=="F")  %>% filter(OsAH_upreg=="F"))
# all 3 shared F AH upreg genes = 0

nrow(full_sex_res_AH %>% filter(OtAH_upreg=="F") %>% filter(OsAH_upreg=="F")) 
# Otau + Osag shared F AH upreg genes = 10

nrow(full_sex_res_AH %>% filter(DgAH_upreg=="F") %>% filter(OsAH_upreg=="F")) 
# Dgaz + Osag shared F AH upreg genes = 0

nrow(full_sex_res_AH %>% filter(OtAH_upreg=="F") %>% filter(DgAH_upreg=="F")) 
# Dgaz + Otau shared F AH upreg genes = 0

###### Step 7: export sex-responsive genes ####
### Otau ####
write.table(deg_Ot_PH_MvF_anno, 
            file = "/Users/ericanadolski/Desktop/beetle-sex-dimorph/Ot_PH_MvF_deg.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(deg_Ot_AH_MvF_anno, 
            file = "/Users/ericanadolski/Desktop/beetle-sex-dimorph/Ot_AH_MvF_deg.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(deg_Ot_G_MvF_anno, 
            file = "/Users/ericanadolski/Desktop/beetle-sex-dimorph/Ot_G_MvF_deg.txt", sep = "\t", quote = FALSE, row.names = FALSE)

### Osag ####
write.table(deg_Os_PH_MvF_anno, 
            file = "/Users/ericanadolski/Desktop/beetle-sex-dimorph/Os_PH_MvF_deg.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(deg_Os_AH_MvF_anno, 
            file = "/Users/ericanadolski/Desktop/beetle-sex-dimorph/Os_AH_MvF_deg.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(deg_Os_G_MvF_anno, 
            file = "/Users/ericanadolski/Desktop/beetle-sex-dimorph/Os_G_MvF_deg.txt", sep = "\t", quote = FALSE, row.names = FALSE)
### Dgaz ####
write.table(deg_Dg_PH_MvF_anno, 
            file = "./GitHub/Beetle-sexual-dimorphism/Dg_PH_MvF_deg.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(deg_Dg_AH_MvF_anno, 
            file = "./GitHub/Beetle-sexual-dimorphism/Dg_AH_MvF_deg.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(deg_Dg_G_MvF_anno, 
            file = "./GitHub/Beetle-sexual-dimorphism/Dg_G_MvF_deg.txt",sep = "\t", quote = FALSE, row.names = FALSE)

### export upregulated gene sets for each sex, trait#####
### PH 
## Otau M upreg 
deg_Ot_PH_M_up <- deg_Ot_PH_MvF_anno %>% filter(log2FoldChange > 0)
nrow(deg_Ot_PH_M_up)
write.table(deg_Ot_PH_M_up, 
            file = "./GitHub/Beetle-sexual-dimorphism/Ot_PH_M_up.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

## Otau F upreg 
deg_Ot_PH_F_up <- deg_Ot_PH_MvF_anno %>% filter(log2FoldChange < 0)
nrow(deg_Ot_PH_F_up)
write.table(deg_Ot_PH_F_up, 
            file = "./GitHub/Beetle-sexual-dimorphism/Ot_PH_F_up.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

## Osag M upreg 
deg_Os_PH_M_up <- deg_Os_PH_MvF_anno %>% filter(log2FoldChange > 0)
nrow(deg_Os_PH_M_up)
write.table(deg_Os_PH_M_up, 
            file = "./GitHub/Beetle-sexual-dimorphism/Os_PH_M_up.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

## Osag F upreg 
deg_Os_PH_F_up <- deg_Os_PH_MvF_anno %>% filter(log2FoldChange < 0)
nrow(deg_Os_PH_F_up)
write.table(deg_Os_PH_F_up, 
            file = "./GitHub/Beetle-sexual-dimorphism/Os_PH_F_up.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

## Dgaz M upreg 
deg_Dg_PH_M_up <- deg_Dg_PH_MvF_anno %>% filter(log2FoldChange > 0)
nrow(deg_Dg_PH_M_up)
write.table(deg_Dg_PH_M_up, 
            file = "./GitHub/Beetle-sexual-dimorphism/Dg_PH_M_up.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

## Otau F upreg 
deg_Dg_PH_F_up <- deg_Dg_PH_MvF_anno %>% filter(log2FoldChange < 0)
nrow(deg_Dg_PH_F_up)
write.table(deg_Dg_PH_F_up, 
            file = "./GitHub/Beetle-sexual-dimorphism/Dg_PH_F_up.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)


## AH 
## Osag M  AH upreg 
deg_Os_AH_M_up <- deg_Os_AH_MvF_anno %>% filter(log2FoldChange > 0)
nrow(deg_Os_AH_M_up)
write.table(deg_Os_AH_M_up, 
            file = "./GitHub/Beetle-sexual-dimorphism/Os_AH_M_up.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

## Osag F AH upreg 
deg_Os_AH_F_up <- deg_Os_AH_MvF_anno %>% filter(log2FoldChange < 0)
nrow(deg_Os_AH_F_up)
write.table(deg_Os_AH_F_up, 
            file = "./GitHub/Beetle-sexual-dimorphism/Os_AH_F_up.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

## Otau M AH upreg 
deg_Ot_AH_M_up <- deg_Ot_AH_MvF_anno %>% filter(log2FoldChange > 0)
nrow(deg_Ot_AH_M_up)
write.table(deg_Ot_AH_M_up, 
            file = "./GitHub/Beetle-sexual-dimorphism/Ot_AH_M_up.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

## Otau F AH upreg 
deg_Ot_AH_F_up <- deg_Ot_AH_MvF_anno %>% filter(log2FoldChange < 0)
nrow(deg_Ot_AH_F_up)
write.table(deg_Ot_AH_F_up, 
            file = "./GitHub/Beetle-sexual-dimorphism/Ot_AH_F_up.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

## Dgaz M AH upreg 
deg_Dg_AH_M_up <- deg_Dg_AH_MvF_anno %>% filter(log2FoldChange > 0)
nrow(deg_Dg_AH_M_up)
write.table(deg_Dg_AH_M_up, 
            file = "./GitHub/Beetle-sexual-dimorphism/Dg_AH_M_up.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

## Dgaz F AH upreg 
deg_Dg_AH_F_up <- deg_Dg_AH_MvF_anno %>% filter(log2FoldChange < 0)
nrow(deg_Dg_AH_F_up)
write.table(deg_Dg_AH_F_up, 
            file = "./GitHub/Beetle-sexual-dimorphism/Dg_AH_F_up.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

### Genitalia
## Dgaz M G upreg 
deg_Dg_G_M_up <- deg_Dg_G_MvF_anno %>% filter(log2FoldChange > 0)
nrow(deg_Dg_G_M_up)
write.table(deg_Dg_G_M_up, 
            file = "./GitHub/Beetle-sexual-dimorphism/Dg_G_M_up.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

## Dgaz F G upreg 
deg_Dg_G_F_up <- deg_Dg_G_MvF_anno %>% filter(log2FoldChange < 0)
nrow(deg_Dg_G_F_up)
write.table(deg_Dg_G_F_up, 
            file = "./GitHub/Beetle-sexual-dimorphism/Dg_G_F_up.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

## Otau M G upreg 
deg_Ot_G_M_up <- deg_Ot_G_MvF_anno %>% filter(log2FoldChange > 0)
nrow(deg_Ot_G_M_up)
write.table(deg_Ot_G_M_up, 
            file = "./GitHub/Beetle-sexual-dimorphism/Ot_G_M_up.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

## Otau F G upreg 
deg_Ot_G_F_up <- deg_Ot_G_MvF_anno %>% filter(log2FoldChange < 0)
nrow(deg_Ot_G_F_up)
write.table(deg_Ot_G_F_up, 
            file = "./GitHub/Beetle-sexual-dimorphism/Ot_G_F_up.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

## Os M G upreg 
deg_Os_G_M_up <- deg_Os_G_MvF_anno %>% filter(log2FoldChange > 0)
nrow(deg_Os_G_M_up)
write.table(deg_Os_G_M_up, 
            file = "./GitHub/Beetle-sexual-dimorphism/Os_G_M_up.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

## Os F G upreg 
deg_Os_G_F_up <- deg_Os_G_MvF_anno %>% filter(log2FoldChange < 0)
nrow(deg_Os_G_F_up)
write.table(deg_Os_G_F_up, 
            file = "./GitHub/Beetle-sexual-dimorphism/Os_G_F_up.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

###### Step 9: Plotting doublesex expression levels #####
#### used updated GTF/fasta salmon pseudomapping run
### O taurus ####
##### use the 'cpm' function from EdgeR to get counts per million
# make the DGEList
Ot_eds_r <- DGEList(Ot_counts)
# calculate TMM normalization factors
Ot_eds_r <- calcNormFactors(Ot_eds_r)
#get the normalized counts
Ot_cpm <- cpm(Ot_eds_r, log=FALSE) # matrix output

# make new df 
Ot_dsx <- as.data.frame(t(Ot_cpm[c("jg10020.t1", "jg10020.t2", "jg10020.t3", "jg10020.t4", "jg10020.t5", "jg10020.t6", "jg10020.t7"),]))
Ot_dsx$Sex_Trait <- Ot_sample_table$Sex_Trait
Ot_dsx <- rename(Ot_dsx, "F1"="jg10020.t1", "F2"="jg10020.t2", "F3"="jg10020.t3", "F4"="jg10020.t4", "F5"="jg10020.t5", "M"="jg10020.t6", "B"="jg10020.t7")

ggplot(Ot_dsx, aes(x=Sex_Trait, y=jg10020.t1)) +
  geom_boxplot() +
  scale_x_discrete(limits=c("F_AH","M_AH","F_PH","M_PH","F_L","M_L","F_G","M_G","F_E","M_E")) +
  labs(title="O.tau dsxF1 read counts",x="Sample Type", y = "CPM") + theme_classic()

Ot_dsx_long <- Ot_dsx %>% rownames_to_column("sample") %>% pivot_longer(cols=c("F1", "F2", "F3", "F4", "F5", "M", "B"),
                    names_to='isoform',
                    values_to='counts')

Ot_plot <- ggplot(Ot_dsx_long) +
  geom_boxplot(aes(x=Sex_Trait, y=counts, color=isoform)) +
  scale_x_discrete(limits=c("F_AH","M_AH","F_PH","M_PH","F_L","M_L","F_G","M_G","F_E","M_E")) +
  labs(title="O.tau dsx isoform expression across sample types",x="Sample Type", y = "CPM") + theme_classic() + scale_color_manual(values=c("#999999", "#F564E3","#F8766D", "#B79F00", "#00BA38","#00BFC4","#619CFF"))

### O sagittarius ####
# make the DGEList
Os_eds_r <- DGEList(Os_counts)
# calculate TMM normalization factors
Os_eds_r <- calcNormFactors(Os_eds_r)
#get the normalized counts
Os_cpm <- cpm(Os_eds_r, log=FALSE) # matrix output
Os_dsx <- as.data.frame(t(Os_counts[c("jg15025.t1", "jg15025.t2", "jg15025.t3", "jg15025.t4", "jg15025.t5", "jg15025.t6"),]))
Os_dsx$Sex_Trait <- Os_sample_table$Sex_Trait
Os_dsx <- rename(Os_dsx, "F1"="jg15025.t1", "F2"="jg15025.t2", "F3"="jg15025.t3", "F4"="jg15025.t4", "F5"="jg15025.t5", "M"="jg15025.t6")

Os_dsx_long <- Os_dsx %>% rownames_to_column("sample") %>% pivot_longer(cols=c("F1", "F2", "F3", "F4", "F5", "M"), names_to='isoform', values_to='counts')
Os_plot <- ggplot(Os_dsx_long) +
  geom_boxplot(aes(x=Sex_Trait, y=counts, color=isoform)) +
  scale_x_discrete(limits=c("F_AH","M_AH","F_PH","M_PH","F_L","M_L","F_G","M_G","F_E","M_E")) +
  labs(title="Osag dsx isoform expression",x="Sample Type", y = "CPM") + theme_classic() +
  scale_color_manual(values=c("#F564E3","#F8766D", "#B79F00", "#00BA38","#00BFC4","#619CFF"))

ggplot_build(Os_plot)$dat

### D gazella #####
##### use the 'cpm' function from EdgeR to get counts per million
# make the DGEList
Dg_eds_r <- DGEList(Dg_counts)
# calculate TMM normalization factors
Dg_eds_r <- calcNormFactors(Dg_eds_r)
#get the normalized counts
Dg_cpm <- as.data.frame(cpm(Dg_eds_r, log=FALSE)) # matrix output

# make new df 
Dg_dsx <- as.data.frame(t(Dg_cpm[c("jg19569.t1", "jg19569.t2"),]))
Dg_dsx$Sex_Trait <- Dg_sample_table$Sex_Trait
Dg_dsx <- rename(Dg_dsx, "dsxF"="jg19569.t1", "dsxM"="jg19569.t2")

Dg_dsx_long <- Dg_dsx %>% rownames_to_column("sample") %>% pivot_longer(cols=c("dsxF", "dsxM"),
                                                                        names_to='isoform',
                                                                        values_to='counts')

ggplot(Dg_dsx_long) +
  geom_boxplot(aes(x=Sex_Trait, y=counts, color=isoform)) +
  scale_x_discrete(limits=c("F_AH","M_AH","F_PH","M_PH","F_L","M_L","F_G","M_G","F_E","M_E")) +
  labs(title="Dgaz dsx isoform expression across sample types",x="Sample Type", y = "CPM") + theme_classic() + scale_color_manual(values=c("#F564E3","#619CFF"))

# II. ATACseq data analyses ############################################################################
###### Step 1: import ATAC-seq read counts from bedtools multicov #####
# Otau
Ot_OCR_counts <- read.delim("/Users/ericanadolski/GitHub/Onthophagus_sexual_dimorphism/Ot_peak_counts.txt", header=TRUE)
Ot_info_ATAC <- read.delim("/Users/ericanadolski/GitHub/Onthophagus_sexual_dimorphism/Ot_sample_info_table.txt", row.names = NULL)

# Osag
Os_OCR_counts <- read.delim("/Users/ericanadolski/GitHub/Onthophagus_sexual_dimorphism/Os_peak_counts.txt", header=TRUE)
Os_info_ATAC <- read.delim("/Users/ericanadolski/GitHub/Onthophagus_sexual_dimorphism/Os_sample_info_table.txt")

# Dgaz
Dg_OCR_counts <- read.delim("/Users/ericanadolski/GitHub/Onthophagus_sexual_dimorphism/Dg_peak_counts.txt", header=TRUE)
Dg_info_ATAC <- read.delim("/Users/ericanadolski/GitHub/Onthophagus_sexual_dimorphism/Dg_sample_info_table.txt")

###### Step 2: normalize X chromosome counts and filter out low read counts #####
# write function to multiply by 2 ####
double <- function(observed) {
  result <- (observed * 2)
} 

### Otau ######
# subset out X chromosome rows
Ot_OCR_counts_X <- Ot_OCR_counts %>% filter(chr == "Scaffold8")

# subset female columns of X
Ot_OCR_counts_X_F <- Ot_OCR_counts_X[,1:19]

# subset male columns and multiply read counts by 2 to normalize
Ot_OCR_counts_X_M_2 <- Ot_OCR_counts_X[,20:34] %>% mutate_all(list(double))

# append normalized male to female columns 
Ot_OCR_counts_X_norm <- cbind(Ot_OCR_counts_X_M_2, Ot_OCR_counts_X_F)

# append X chr rows back
Ot_OCR_counts_norm <- Ot_OCR_counts %>% filter(!chr == "Scaffold8")
Ot_OCR_counts_norm <- rbind(Ot_OCR_counts_norm, Ot_OCR_counts_X_norm)

# filter out OCRs with extremely low read counts 
cpm_Ot_norm <- cpm(Ot_OCR_counts_norm[c(5:34)]) # 153897 OCRs
countcheck_Ot_norm <- cpm_Ot_norm > 3
keep_Ot_norm <- which(rowSums(countcheck_Ot_norm) >= 5)
Ot_filtered_OCR_counts_norm <- Ot_OCR_counts_norm[keep_Ot_norm,] # 78826 OCRs

### Osag ####
## subset out X chromosome rows
Os_OCR_counts_X <- Os_OCR_counts %>% filter(chr == "Scaffold_7")

# subset female columns of X
Os_OCR_counts_X_F <- Os_OCR_counts_X[,1:19]

# subset male columns and multiply read counts by 2 to normalize
Os_OCR_counts_X_M_2 <- Os_OCR_counts_X[,20:34] %>% mutate_all(list(double))

# append normalized male to female columns 
Os_OCR_counts_X_norm <- cbind(Os_OCR_counts_X_M_2, Os_OCR_counts_X_F)

# append X chr rows back
Os_OCR_counts_norm <- Os_OCR_counts %>% filter(!chr == "Scaffold_7")
Os_OCR_counts_norm <- rbind(Os_OCR_counts_norm, Os_OCR_counts_X_norm)

# filter out OCRs with extremely low read counts 
cpm_Os_norm <- cpm(Os_OCR_counts_norm[c(5:34)]) # 182967
countcheck_Os_norm <- cpm_Os_norm > 3
keep_Os_norm <- which(rowSums(countcheck_Os_norm) >= 5)
Os_filtered_OCR_counts_norm <- Os_OCR_counts_norm[keep_Os_norm,] # 97643

### Dgaz ######
# subset out X chromosome rows
Dg_OCR_counts_X <- Dg_OCR_counts %>% filter(chr == "ScIV947_12;HRSCAF=14")

# subset female columns of X
Dg_OCR_counts_X_F <- Dg_OCR_counts_X[,1:19]

# subset male columns and multiply read counts by 2 to normalize
Dg_OCR_counts_X_M_2 <- Dg_OCR_counts_X[,20:34] %>% mutate_all(list(double))

# append normalized male to female columns 
Dg_OCR_counts_X_norm <- cbind(Dg_OCR_counts_X_M_2, Dg_OCR_counts_X_F)

# append X chr rows back
Dg_OCR_counts_norm <- Dg_OCR_counts %>% filter(!chr == "ScIV947_12;HRSCAF=14")
Dg_OCR_counts_norm <- rbind(Dg_OCR_counts_norm, Dg_OCR_counts_X_norm)

# filter out OCRs with extremely low read counts 
cpm_Dg_norm <- cpm(Dg_OCR_counts_norm[c(5:34)]) # 199990
countcheck_Dg_norm <- cpm_Dg_norm > 3
keep_Dg_norm <- which(rowSums(countcheck_Dg_norm) >= 5)
Dg_filtered_OCR_counts_norm <- Dg_OCR_counts_norm[keep_Dg_norm,] # 90282 OCRs

###### Step 3: data visualization ######
############### heatmaps all traits  #####
# Otau 
# generate DESeq dataset for exploratory analysis
dds_Ot_ATAC <- DESeqDataSetFromMatrix(countData = Ot_filtered_OCR_counts_norm[5:34],
                                      colData = Ot_info_ATAC, design = ~ 1)
Ot_ATAC_dds_rld <- rlog(dds_Ot_ATAC)
Ot_ATAC_dds_rld <- as.data.frame(assay(Ot_ATAC_dds_rld))
Ot_ATAC_dds_rld <- as.data.frame(t(Ot_ATAC_dds_rld))
# heatmap
Ot_ATAC_sampleDists <- dist(Ot_ATAC_dds_rld)
Ot_ATAC_sampleDistMatrix <- as.matrix(Ot_ATAC_sampleDists)
rownames(Ot_ATAC_sampleDistMatrix) <- paste(Ot_info_ATAC$Sample) 
colnames(Ot_ATAC_sampleDistMatrix) <- NULL 
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) 
pheatmap(Ot_ATAC_sampleDistMatrix,
         clustering_distance_rows=Ot_ATAC_sampleDists,
         clustering_distance_cols=Ot_ATAC_sampleDists,
         col=colors)
# Osag 
dds_Os_ATAC <- DESeqDataSetFromMatrix(countData = Os_filtered_OCR_counts_norm[5:34],
                                      colData = Os_info_ATAC, design = ~ 1)
Os_ATAC_dds_rld <- rlog(dds_Os_ATAC)
Os_ATAC_dds_rld <- as.data.frame(assay(Os_ATAC_dds_rld))
Os_ATAC_dds_rld <- as.data.frame(t(Os_ATAC_dds_rld))
# heatmap
Os_ATAC_sampleDists <- dist(Os_ATAC_dds_rld)
Os_ATAC_sampleDistMatrix <- as.matrix(Os_ATAC_sampleDists)
rownames(Os_ATAC_sampleDistMatrix) <- paste(Os_info_ATAC$Sample) 
colnames(Os_ATAC_sampleDistMatrix) <- NULL 
pheatmap(Os_ATAC_sampleDistMatrix,
         clustering_distance_rows=Os_ATAC_sampleDists,
         clustering_distance_cols=Os_ATAC_sampleDists,
         col=colors)

# Dgaz
dds_Dg_ATAC <- DESeqDataSetFromMatrix(countData = Dg_filtered_OCR_counts_norm[5:34],
                                      colData = Dg_info_ATAC, design = ~ 1)
Dg_ATAC_dds_rld <- rlog(dds_Dg_ATAC)
Dg_ATAC_dds_rld <- as.data.frame(assay(Dg_ATAC_dds_rld))
Dg_ATAC_dds_rld <- as.data.frame(t(Dg_ATAC_dds_rld))
# heatmap
Dg_ATAC_sampleDists <- dist(Dg_ATAC_dds_rld)
Dg_ATAC_sampleDistMatrix <- as.matrix(Dg_ATAC_sampleDists)
rownames(Dg_ATAC_sampleDistMatrix) <- paste(Dg_info_ATAC$Sample) 
colnames(Dg_ATAC_sampleDistMatrix) <- NULL 
pheatmap(Dg_ATAC_sampleDistMatrix,
         clustering_distance_rows=Dg_ATAC_sampleDists,
         clustering_distance_cols=Dg_ATAC_sampleDists,
         col=colors)
############### PCA all traits ####
#Otau
PCA_Ot_ATAC <- prcomp(Ot_ATAC_dds_rld , center = TRUE)
PCs_Ot_ATAC <- as.data.frame(PCA_Ot_ATAC$x)
PCs_Ot_ATAC$Sample <- row.names(PCs_Ot_ATAC)
PCs_Ot_ATAC <- inner_join(PCs_Ot_ATAC, Ot_info_ATAC, by='Sample') 
summary(PCA_Ot_ATAC)

Ot_PCA_ATAC <- ggplot(data = PCs_Ot_ATAC, 
                      aes(x = PC1, y = PC2, color = Sex, shape=Trait)) + 
  geom_point(size = 6, alpha=0.8) +
  scale_color_manual(values=c("#C1272D", "#2166AC")) + 
  #geom_text(aes(label = Replicate), nudge_y = 5) +
  labs(title="O. taurus ATAC-seq",x = "PC1 (31.1%)", y = "PC2 (13.9%)") +
  theme_bw(base_size=16)
Ot_PCA_ATAC

ggsave(Ot_PCA_ATAC, file = "/Users/ericanadolski/Desktop/beetle-sex-dimorph/figures/ATAC_Ot_PCA.pdf",
       width = 5, height = 4)

#Osag
PCA_Os_ATAC <- prcomp(Os_ATAC_dds_rld , center = TRUE)
PCs_Os_ATAC <- as.data.frame(PCA_Os_ATAC$x)
PCs_Os_ATAC$Sample <- row.names(PCs_Os_ATAC)
PCs_Os_ATAC <- inner_join(PCs_Os_ATAC, Os_info_ATAC, by='Sample') 
summary(PCA_Os_ATAC)

Os_PCA_ATAC <- ggplot(data = PCs_Os_ATAC, 
                      aes(x = PC1, y = PC2, color = Sex, shape=Trait)) + 
  geom_point(size = 6, alpha=0.8) +
  scale_color_manual(values=c("#C1272D", "#2166AC")) + 
  #geom_text(aes(label = Replicate), nudge_y = 5) +
  labs(title="O. sagittarius ATAC-seq",x = "PC1 (32.0%)", y = "PC2 (12.2%)") +
  theme_bw(base_size=16)
Os_PCA_ATAC

ggsave(Os_PCA_ATAC, file = "/Users/ericanadolski/Desktop/beetle-sex-dimorph/figures/ATAC_Os_PCA.pdf",
       width = 5, height = 4)

#Dgaz
PCA_Dg_ATAC <- prcomp(Dg_ATAC_dds_rld , center = TRUE)
PCs_Dg_ATAC <- as.data.frame(PCA_Dg_ATAC$x)
PCs_Dg_ATAC$Sample <- row.names(PCs_Dg_ATAC)
PCs_Dg_ATAC <- inner_join(PCs_Dg_ATAC, Dg_info_ATAC, by='Sample') 
summary(PCA_Dg_ATAC)

Dg_PCA_ATAC <- ggplot(data = PCs_Dg_ATAC, 
                      aes(x = PC1, y = PC2, color = Sex, shape=Trait)) + 
  geom_point(size = 6, alpha=0.8) +
  scale_color_manual(values=c("#C1272D", "#2166AC")) + 
  #geom_text(aes(label = Replicate), nudge_y = 5) +
  labs(title="D. gazella ATAC-seq",x = "PC1 (23.2%)", y = "PC2 (11.1%)") +
  theme_bw(base_size=16)
Dg_PCA_ATAC

ggsave(Dg_PCA_ATAC, file = "/Users/ericanadolski/Desktop/beetle-sex-dimorph/figures/ATAC_Dg_PCA.pdf",
       width = 5, height = 4)
############### PCA subset by trait ####
########### D gazella ####
### PCA posterior head ####
Dg_PH_info_ATAC <- filter(Dg_info_ATAC, Trait == "PH")
Dg_PH_filtered_counts_ATAC <- Dg_filtered_OCR_counts_norm %>% dplyr::select(contains("PH")) 

dds_Dg_PH_ATAC <- DESeqDataSetFromMatrix(countData = Dg_PH_filtered_counts_ATAC,
                                         colData = Dg_PH_info_ATAC,
                                         design = ~ 1)
rlog_dds_Dg_PH_ATAC <- rlog(dds_Dg_PH_ATAC) 
rlog_counts_Dg_PH_ATAC <- as.data.frame(assay(rlog_dds_Dg_PH_ATAC))
rlog_countsT_Dg_PH_ATAC <- as.data.frame(t(rlog_counts_Dg_PH_ATAC))

PCA_Dg_PH_ATAC <- prcomp(rlog_countsT_Dg_PH_ATAC, center = TRUE)
PCs_Dg_PH_ATAC <- as.data.frame(PCA_Dg_PH_ATAC$x)
PCs_Dg_PH_ATAC$Sample <- row.names(PCs_Dg_PH_ATAC)
PCs_Dg_PH_ATAC <- inner_join(PCs_Dg_PH_ATAC, Dg_PH_info_ATAC, by='Sample') 
summary(PCA_Dg_PH_ATAC)

PCA_Dg_PH_ATAC <- ggplot(data = PCs_Dg_PH_ATAC, 
                         aes(x = PC1, y = PC2, color=Sex)) + 
  geom_point(size = 6, shape=16) +
  scale_color_manual(values = cols) +
  labs(title="D.gazella ATAC-seq Posterior Head",x = "PC1 (30.0%)", y = "PC2 (18.3%)") +
  theme_bw()
PCA_Dg_PH_ATAC

ggsave(PCA_Dg_PH_ATAC, file = "/Users/ericanadolski/Desktop/beetle-sex-dimorph/figures/ATAC_Dg_PH_PCA.pdf",
       width = 4.30, height = 3.77)

### PCA genitalia ####
Dg_G_info_ATAC <- filter(Dg_info_ATAC, Trait == "G")
Dg_G_filtered_counts_ATAC <- Dg_filtered_OCR_counts_norm %>% dplyr::select(contains("_G")) 

dds_Dg_G_ATAC <- DESeqDataSetFromMatrix(countData = Dg_G_filtered_counts_ATAC,
                                        colData = Dg_G_info_ATAC,
                                        design = ~ 1)
rlog_dds_Dg_G_ATAC <- rlog(dds_Dg_G_ATAC) 
rlog_counts_Dg_G_ATAC <- as.data.frame(assay(rlog_dds_Dg_G_ATAC))
rlog_countsT_Dg_G_ATAC <- as.data.frame(t(rlog_counts_Dg_G_ATAC))

PCA_Dg_G_ATAC <- prcomp(rlog_countsT_Dg_G_ATAC, center = TRUE)
PCs_Dg_G_ATAC <- as.data.frame(PCA_Dg_G_ATAC$x)
PCs_Dg_G_ATAC$Sample <- row.names(PCs_Dg_G_ATAC)
PCs_Dg_G_ATAC <- inner_join(PCs_Dg_G_ATAC, Dg_G_info_ATAC, by='Sample') 
summary(PCA_Dg_G_ATAC)

PCA_Dg_G_ATAC <- ggplot(data = PCs_Dg_G_ATAC, 
                        aes(x = PC1, y = PC2, color=Sex)) + 
  geom_point(size = 6, shape=16) +
  scale_color_manual(values = cols) +
  labs(title="D.gazella ATAC-seq Genitalia",x = "PC1 (34.2%)", y = "PC2 (13.5%)") +
  theme_bw()
PCA_Dg_G_ATAC

ggsave(PCA_Dg_G_ATAC, file = "/Users/ericanadolski/Desktop/beetle-sex-dimorph/figures/ATAC_Dg_G_PCA.pdf",
       width = 4.30, height = 3.77)

### PCA Anterior Head ####
Dg_AH_info_ATAC <- filter(Dg_info_ATAC, Trait == "AH")
Dg_AH_filtered_counts_ATAC <- Dg_filtered_OCR_counts_norm %>% dplyr::select(contains("AH")) 

dds_Dg_AH_ATAC <- DESeqDataSetFromMatrix(countData = Dg_AH_filtered_counts_ATAC,
                                         colData = Dg_AH_info_ATAC,
                                         design = ~ 1)
rlog_dds_Dg_AH_ATAC <- rlog(dds_Dg_AH_ATAC) 
rlog_counts_Dg_AH_ATAC <- as.data.frame(assay(rlog_dds_Dg_AH_ATAC))
rlog_countsT_Dg_AH_ATAC <- as.data.frame(t(rlog_counts_Dg_AH_ATAC))

PCA_Dg_AH_ATAC <- prcomp(rlog_countsT_Dg_AH_ATAC, center = TRUE)
PCs_Dg_AH_ATAC <- as.data.frame(PCA_Dg_AH_ATAC$x)
PCs_Dg_AH_ATAC$Sample <- row.names(PCs_Dg_AH_ATAC)
PCs_Dg_AH_ATAC <- inner_join(PCs_Dg_AH_ATAC, Dg_AH_info_ATAC, by='Sample') 
summary(PCA_Dg_AH_ATAC)

PCA_Dg_AH_ATAC <- ggplot(data = PCs_Dg_AH_ATAC, 
                         aes(x = PC1, y = PC2, color=Sex)) + 
  geom_point(size = 6, shape=16) +
  scale_color_manual(values = cols) +
  labs(title="D.gazella ATAC-seq Anterior Head",x = "PC1 (32.0%)", y = "PC2 (14.2%)") +
  theme_bw()
PCA_Dg_AH_ATAC

ggsave(PCA_Dg_AH_ATAC, file = "/Users/ericanadolski/Desktop/beetle-sex-dimorph/figures/ATAC_Dg_AH_PCA.pdf",
       width = 4.30, height = 3.77)
########### O taurus  #########
### PCA posterior head ####
Ot_PH_info_ATAC <- filter(Ot_info_ATAC, Trait == "PH")
Ot_PH_filtered_counts_ATAC <- Ot_filtered_OCR_counts_norm %>% dplyr::select(contains("PH")) 

dds_Ot_PH_ATAC <- DESeqDataSetFromMatrix(countData = Ot_PH_filtered_counts_ATAC,
                                         colData = Ot_PH_info_ATAC,
                                         design = ~ 1)
rlog_dds_Ot_PH_ATAC <- rlog(dds_Ot_PH_ATAC) 
rlog_counts_Ot_PH_ATAC <- as.data.frame(assay(rlog_dds_Ot_PH_ATAC))
rlog_countsT_Ot_PH_ATAC <- as.data.frame(t(rlog_counts_Ot_PH_ATAC))

PCA_Ot_PH_ATAC <- prcomp(rlog_countsT_Ot_PH_ATAC, center = TRUE)
PCs_Ot_PH_ATAC <- as.data.frame(PCA_Ot_PH_ATAC$x)
PCs_Ot_PH_ATAC$Sample <- row.names(PCs_Ot_PH_ATAC)
PCs_Ot_PH_ATAC <- inner_join(PCs_Ot_PH_ATAC, Ot_PH_info_ATAC, by='Sample') 
summary(PCA_Ot_PH_ATAC)

PCA_Ot_PH_ATAC <- ggplot(data = PCs_Ot_PH_ATAC, 
                         aes(x = PC1, y = PC2, color=Sex)) + 
  geom_point(size = 6, shape=16) +
  scale_color_manual(values = cols) +
  labs(title="O.taurus ATAC-seq Posterior Head",x = "PC1 (39.7%)", y = "PC2 (20.0%)") +
  theme_bw()
PCA_Ot_PH_ATAC

ggsave(PCA_Ot_PH_ATAC, file = "/Users/ericanadolski/Desktop/beetle-sex-dimorph/figures/ATAC_Ot_PH_PCA.pdf",
       width = 4.30, height = 3.77)

### PCA genitalia ####
Ot_G_info_ATAC <- filter(Ot_info_ATAC, Trait == "G")
Ot_G_filtered_counts_ATAC <- Ot_filtered_OCR_counts_norm %>% dplyr::select(contains("G")) 

dds_Ot_G_ATAC <- DESeqDataSetFromMatrix(countData = Ot_G_filtered_counts_ATAC,
                                        colData = Ot_G_info_ATAC,
                                        design = ~ 1)
rlog_dds_Ot_G_ATAC <- rlog(dds_Ot_G_ATAC) 
rlog_counts_Ot_G_ATAC <- as.data.frame(assay(rlog_dds_Ot_G_ATAC))
rlog_countsT_Ot_G_ATAC <- as.data.frame(t(rlog_counts_Ot_G_ATAC))

PCA_Ot_G_ATAC <- prcomp(rlog_countsT_Ot_G_ATAC, center = TRUE)
PCs_Ot_G_ATAC <- as.data.frame(PCA_Ot_G_ATAC$x)
PCs_Ot_G_ATAC$Sample <- row.names(PCs_Ot_G_ATAC)
PCs_Ot_G_ATAC <- inner_join(PCs_Ot_G_ATAC, Ot_G_info_ATAC, by='Sample') 
summary(PCA_Ot_G_ATAC)

PCA_Ot_G_ATAC <- ggplot(data = PCs_Ot_G_ATAC, 
                        aes(x = PC1, y = PC2, color=Sex)) + 
  geom_point(size = 6, shape=16) +
  scale_color_manual(values = cols) +
  labs(title="O.taurus ATAC-seq Genitalia",x = "PC1 (39.4%)", y = "PC2 (14.4%)") +
  theme_bw()
PCA_Ot_G_ATAC

ggsave(PCA_Ot_G_ATAC, file = "/Users/ericanadolski/Desktop/beetle-sex-dimorph/figures/ATAC_Ot_G_PCA.pdf",
       width = 4.30, height = 3.77)

### PCA Anterior Head ####
Ot_AH_info_ATAC <- filter(Ot_info_ATAC, Trait == "AH")
Ot_AH_filtered_counts_ATAC <- Ot_filtered_OCR_counts_norm %>% dplyr::select(contains("AH")) 

dds_Ot_AH_ATAC <- DESeqDataSetFromMatrix(countData = Ot_AH_filtered_counts_ATAC,
                                         colData = Ot_AH_info_ATAC,
                                         design = ~ 1)
rlog_dds_Ot_AH_ATAC <- rlog(dds_Ot_AH_ATAC) 
rlog_counts_Ot_AH_ATAC <- as.data.frame(assay(rlog_dds_Ot_AH_ATAC))
rlog_countsT_Ot_AH_ATAC <- as.data.frame(t(rlog_counts_Ot_AH_ATAC))

PCA_Ot_AH_ATAC <- prcomp(rlog_countsT_Ot_AH_ATAC, center = TRUE)
PCs_Ot_AH_ATAC <- as.data.frame(PCA_Ot_AH_ATAC$x)
PCs_Ot_AH_ATAC$Sample <- row.names(PCs_Ot_AH_ATAC)
PCs_Ot_AH_ATAC <- inner_join(PCs_Ot_AH_ATAC, Ot_AH_info_ATAC, by='Sample') 
summary(PCA_Ot_AH_ATAC)

PCA_Ot_AH_ATAC <- ggplot(data = PCs_Ot_AH_ATAC, 
                         aes(x = PC1, y = PC2, color=Sex)) + 
  geom_point(size = 6, shape=16) +
  scale_color_manual(values = cols) +
  labs(title="O.taurus ATAC-seq Anterior Head",x = "PC1 (47.8%)", y = "PC2 (16.0%)") +
  theme_bw()
PCA_Ot_AH_ATAC

ggsave(PCA_Ot_AH_ATAC, file = "/Users/ericanadolski/Desktop/beetle-sex-dimorph/figures/ATAC_Ot_AH_PCA.pdf",
       width = 4.30, height = 3.77)
########### O sagittarius  ############
### PCA posterior head ####
Os_PH_info_ATAC <- filter(Os_info_ATAC, Trait == "PH")
Os_PH_filtered_counts_ATAC <- Os_filtered_OCR_counts_norm %>% dplyr::select(contains("PH")) 

dds_Os_PH_ATAC <- DESeqDataSetFromMatrix(countData = Os_PH_filtered_counts_ATAC,
                                         colData = Os_PH_info_ATAC,
                                         design = ~ 1)
rlog_dds_Os_PH_ATAC <- rlog(dds_Os_PH_ATAC) 
rlog_counts_Os_PH_ATAC <- as.data.frame(assay(rlog_dds_Os_PH_ATAC))
rlog_countsT_Os_PH_ATAC <- as.data.frame(t(rlog_counts_Os_PH_ATAC))

PCA_Os_PH_ATAC <- prcomp(rlog_countsT_Os_PH_ATAC, center = TRUE)
PCs_Os_PH_ATAC <- as.data.frame(PCA_Os_PH_ATAC$x)
PCs_Os_PH_ATAC$Sample <- row.names(PCs_Os_PH_ATAC)
PCs_Os_PH_ATAC <- inner_join(PCs_Os_PH_ATAC, Os_PH_info_ATAC, by='Sample') 
summary(PCA_Os_PH_ATAC)

PCA_Os_PH_ATAC <- ggplot(data = PCs_Os_PH_ATAC, 
                         aes(x = PC1, y = PC2, color=Sex)) + 
  geom_point(size = 6, shape=16) +
  scale_color_manual(values = cols) +
  labs(title="O.sagittarius ATAC-seq Posterior Head",x = "PC1 (37.6%)", y = "PC2 (14.7%)") +
  theme_bw()
PCA_Os_PH_ATAC

ggsave(PCA_Os_PH_ATAC, file = "/Users/ericanadolski/Desktop/beetle-sex-dimorph/figures/ATAC_Os_PH_PCA.pdf",
       width = 4.30, height = 3.77)

### PCA genitalia ####
Os_G_info_ATAC <- filter(Os_info_ATAC, Trait == "G")
Os_G_filtered_counts_ATAC <- Os_filtered_OCR_counts_norm %>% dplyr::select(contains("G")) 

dds_Os_G_ATAC <- DESeqDataSetFromMatrix(countData = Os_G_filtered_counts_ATAC,
                                        colData = Os_G_info_ATAC,
                                        design = ~ 1)
rlog_dds_Os_G_ATAC <- rlog(dds_Os_G_ATAC) 
rlog_counts_Os_G_ATAC <- as.data.frame(assay(rlog_dds_Os_G_ATAC))
rlog_countsT_Os_G_ATAC <- as.data.frame(t(rlog_counts_Os_G_ATAC))

PCA_Os_G_ATAC <- prcomp(rlog_countsT_Os_G_ATAC, center = TRUE)
PCs_Os_G_ATAC <- as.data.frame(PCA_Os_G_ATAC$x)
PCs_Os_G_ATAC$Sample <- row.names(PCs_Os_G_ATAC)
PCs_Os_G_ATAC <- inner_join(PCs_Os_G_ATAC, Os_G_info_ATAC, by='Sample') 
summary(PCA_Os_G_ATAC)

PCA_Os_G_ATAC <- ggplot(data = PCs_Os_G_ATAC, 
                        aes(x = PC1, y = PC2, color=Sex)) + 
  geom_point(size = 6, shape=16) +
  scale_color_manual(values = cols) +
  labs(title="O.sagittarius ATAC-seq Genitalia",x = "PC1 (37.0%)", y = "PC2 (14.5%)") +
  theme_bw()
PCA_Os_G_ATAC

ggsave(PCA_Os_G_ATAC, file = "/Users/ericanadolski/Desktop/beetle-sex-dimorph/figures/ATAC_Os_G_PCA.pdf",
       width = 4.30, height = 3.77)

### PCA Anterior Head ####
Os_AH_info_ATAC <- filter(Os_info_ATAC, Trait == "AH")
Os_AH_filtered_counts_ATAC <- Os_filtered_OCR_counts_norm %>% dplyr::select(contains("AH")) 

dds_Os_AH_ATAC <- DESeqDataSetFromMatrix(countData = Os_AH_filtered_counts_ATAC,
                                         colData = Os_AH_info_ATAC,
                                         design = ~ 1)
rlog_dds_Os_AH_ATAC <- rlog(dds_Os_AH_ATAC) 
rlog_counts_Os_AH_ATAC <- as.data.frame(assay(rlog_dds_Os_AH_ATAC))
rlog_countsT_Os_AH_ATAC <- as.data.frame(t(rlog_counts_Os_AH_ATAC))

PCA_Os_AH_ATAC <- prcomp(rlog_countsT_Os_AH_ATAC, center = TRUE)
PCs_Os_AH_ATAC <- as.data.frame(PCA_Os_AH_ATAC$x)
PCs_Os_AH_ATAC$Sample <- row.names(PCs_Os_AH_ATAC)
PCs_Os_AH_ATAC <- inner_join(PCs_Os_AH_ATAC, Os_AH_info_ATAC, by='Sample') 
summary(PCA_Os_AH_ATAC)

PCA_Os_AH_ATAC <- ggplot(data = PCs_Os_AH_ATAC, 
                         aes(x = PC1, y = PC2, color=Sex)) + 
  geom_point(size = 6, shape=16) +
  scale_color_manual(values = cols) +
  labs(title="O.sagittarius ATAC-seq Anterior Head",x = "PC1 (39.9%)", y = "PC2 (13.3%)") +
  theme_bw()
PCA_Os_AH_ATAC

ggsave(PCA_Os_AH_ATAC, file = "/Users/ericanadolski/Desktop/beetle-sex-dimorph/figures/ATAC_Os_AH_PCA.pdf",
       width = 4.30, height = 3.77)
###### Step 4: differential chromatin accessibility analysis ####
### Dgaz #### 
### posterior head ####
Dg_PH_info_ATAC <- filter(Dg_info_ATAC, Trait == "PH")
Dg_PH_filtered_counts_ATAC <- Dg_filtered_OCR_counts_norm %>% dplyr::select(contains("PH")) 

dds_Dg_PH_ATAC <- DESeqDataSetFromMatrix(countData = Dg_PH_filtered_counts_ATAC,
                                         colData = Dg_PH_info_ATAC,
                                         design = ~ Sex)
dds_Dg_PH_ATAC <- DESeq(dds_Dg_PH_ATAC)
dds_Dg_PH_ATAC <- lfcShrink(dds_Dg_PH_ATAC, coef="Sex_M_vs_F", type="normal")
tally(as.data.frame(dds_Dg_PH_ATAC, independentFiltering=FALSE) %>% filter(padj <= 0.1)) # 33

### anterior head ####
Dg_AH_info_ATAC <- filter(Dg_info_ATAC, Trait == "AH")
Dg_AH_filtered_counts_ATAC <- Dg_filtered_OCR_counts_norm %>% dplyr::select(contains("AH")) 

dds_Dg_AH_ATAC <- DESeqDataSetFromMatrix(countData = Dg_AH_filtered_counts_ATAC,
                                         colData = Dg_AH_info_ATAC,
                                         design = ~ Sex)

dds_Dg_AH_ATAC <- DESeq(dds_Dg_AH_ATAC)
dds_Dg_AH_ATAC <- lfcShrink(dds_Dg_AH_ATAC, coef="Sex_M_vs_F", type="normal")
tally(as.data.frame(dds_Dg_AH_ATAC, independentFiltering=FALSE) %>% filter(padj <= 0.1)) # 34

### genitalia ####
Dg_G_info_ATAC <- filter(Dg_info_ATAC, Trait == "G")
Dg_G_filtered_counts_ATAC <- Dg_filtered_OCR_counts_norm %>% dplyr::select(contains("_G")) 
head(Dg_G_filtered_counts_ATAC)

dds_Dg_G_ATAC <- DESeqDataSetFromMatrix(countData = Dg_G_filtered_counts_ATAC,
                                        colData = Dg_G_info_ATAC,
                                        design = ~ Sex)
dds_Dg_G_ATAC <- DESeq(dds_Dg_G_ATAC)
dds_Dg_G_ATAC <- lfcShrink(dds_Dg_G_ATAC, coef="Sex_M_vs_F", type="normal")
tally(as.data.frame(dds_Dg_G_ATAC, independentFiltering=FALSE) %>% filter(padj <= 0.1)) # 13

### Otau #### 
### posterior head ####
Ot_PH_info_ATAC <- filter(Ot_info_ATAC, Trait == "PH")
Ot_PH_filtered_counts_ATAC <- Ot_filtered_OCR_counts_norm %>% dplyr::select(contains("PH")) 

dds_Ot_PH_ATAC <- DESeqDataSetFromMatrix(countData = Ot_PH_filtered_counts_ATAC,
                                         colData = Ot_PH_info_ATAC,
                                         design = ~ Sex)
nrow(dds_Ot_PH_ATAC) # 79738 (low read count OCRs already filtered out)

dds_Ot_PH_ATAC <- DESeq(dds_Ot_PH_ATAC)
dds_Ot_PH_ATAC <- lfcShrink(dds_Ot_PH_ATAC, coef="Sex_M_vs_F", type="normal")
tally(as.data.frame(dds_Ot_PH_ATAC, independentFiltering=FALSE) %>% filter(padj <= 0.1)) # 2560

### anterior head ####
Ot_AH_info_ATAC <- filter(Ot_info_ATAC, Trait == "AH")
Ot_AH_filtered_counts_ATAC <- Ot_filtered_OCR_counts_norm %>% dplyr::select(contains("AH")) 

dds_Ot_AH_ATAC <- DESeqDataSetFromMatrix(countData = Ot_AH_filtered_counts_ATAC,
                                         colData = Ot_AH_info_ATAC,
                                         design = ~ Sex)
dds_Ot_AH_ATAC <- DESeq(dds_Ot_AH_ATAC)
dds_Ot_AH_ATAC <- lfcShrink(dds_Ot_AH_ATAC, coef="Sex_M_vs_F", type="normal")
tally(as.data.frame(dds_Ot_AH_ATAC, independentFiltering=FALSE) %>% filter(padj <= 0.1)) # 20

### genitalia ####
Ot_G_info_ATAC <- filter(Ot_info_ATAC, Trait == "G")
Ot_G_filtered_counts_ATAC <- Ot_filtered_OCR_counts_norm %>% dplyr::select(contains("_G")) 

dds_Ot_G_ATAC <- DESeqDataSetFromMatrix(countData = Ot_G_filtered_counts_ATAC,
                                        colData = Ot_G_info_ATAC,
                                        design = ~ Sex)
dds_Ot_G_ATAC <- DESeq(dds_Ot_G_ATAC)
dds_Ot_G_ATAC <- lfcShrink(dds_Ot_G_ATAC, coef="Sex_M_vs_F", type="normal")
tally(as.data.frame(dds_Ot_G_ATAC, independentFiltering=FALSE) %>% filter(padj <= 0.1)) # 1104

### Osag #### 
### posterior head ####
Os_PH_info_ATAC <- filter(Os_info_ATAC, Trait == "PH")
Os_PH_filtered_counts_ATAC <- Os_filtered_OCR_counts_norm %>% dplyr::select(contains("PH")) 

dds_Os_PH_ATAC <- DESeqDataSetFromMatrix(countData = Os_PH_filtered_counts_ATAC,
                                         colData = Os_PH_info_ATAC,
                                         design = ~ Sex)
dds_Os_PH_ATAC <- DESeq(dds_Os_PH_ATAC)
dds_Os_PH_ATAC <- lfcShrink(dds_Os_PH_ATAC, coef="Sex_M_vs_F", type="normal")
tally(as.data.frame(dds_Os_PH_ATAC, independentFiltering=FALSE) %>% filter(padj <= 0.1)) # 162

### anterior head ####
Os_AH_info_ATAC <- filter(Os_info_ATAC, Trait == "AH")
Os_AH_filtered_counts_ATAC <- Os_filtered_OCR_counts_norm %>% dplyr::select(contains("AH")) 

dds_Os_AH_ATAC <- DESeqDataSetFromMatrix(countData = Os_AH_filtered_counts_ATAC,
                                         colData = Os_AH_info_ATAC,
                                         design = ~ Sex)

dds_Os_AH_ATAC <- DESeq(dds_Os_AH_ATAC)
dds_Os_AH_ATAC <- lfcShrink(dds_Os_AH_ATAC, coef="Sex_M_vs_F", type="normal")
tally(as.data.frame(dds_Os_AH_ATAC, independentFiltering=FALSE) %>% filter(padj <= 0.1)) # 213

### genitalia ####
Os_G_info_ATAC <- filter(Os_info_ATAC, Trait == "G")
Os_G_filtered_counts_ATAC <- Os_filtered_OCR_counts_norm %>% dplyr::select(contains("_G")) 

dds_Os_G_ATAC <- DESeqDataSetFromMatrix(countData = Os_G_filtered_counts_ATAC,
                                        colData = Os_G_info_ATAC,
                                        design = ~ Sex)
dds_Os_G_ATAC <- DESeq(dds_Os_G_ATAC)
dds_Os_G_ATAC <- lfcShrink(dds_Os_G_ATAC, coef="Sex_M_vs_F", type="normal")
tally(as.data.frame(dds_Os_G_ATAC, independentFiltering=FALSE) %>% filter(padj <= 0.1)) # 1122

### save results in dataframes ####
### Otau ####
#posterior head
Ot_PH_MvF_OCR <- as.data.frame(dds_Ot_PH_ATAC, independentFiltering=FALSE)
Ot_PH_MvF_OCR <- cbind(Ot_filtered_OCR_counts_norm[1:4], Ot_PH_MvF_OCR)
Ot_PH_MvF_OCR_sig <- Ot_PH_MvF_OCR %>% filter(padj <= 0.1)

#anterior head
Ot_AH_MvF_OCR <- as.data.frame(dds_Ot_AH_ATAC, independentFiltering=FALSE)
Ot_AH_MvF_OCR <- cbind(Ot_filtered_OCR_counts_norm[1:4], Ot_AH_MvF_OCR)
Ot_AH_MvF_OCR_sig <- Ot_AH_MvF_OCR %>% filter(padj <= 0.1)

#genitalia
Ot_G_MvF_OCR <- as.data.frame(dds_Ot_G_ATAC, independentFiltering=FALSE)
Ot_G_MvF_OCR <- cbind(Ot_filtered_OCR_counts_norm[1:4], Ot_G_MvF_OCR)
Ot_G_MvF_OCR_sig <- Ot_G_MvF_OCR %>% filter(padj <= 0.1)

### Osag ####
#posterior head
Os_PH_MvF_OCR <- as.data.frame(dds_Os_PH_ATAC, independentFiltering=FALSE)
Os_PH_MvF_OCR <- cbind(Os_filtered_OCR_counts_norm[1:4], Os_PH_MvF_OCR)
Os_PH_MvF_OCR_sig <- Os_PH_MvF_OCR %>% filter(padj <= 0.1)

#anterior head
Os_AH_MvF_OCR <- as.data.frame(dds_Os_AH_ATAC, independentFiltering=FALSE)
Os_AH_MvF_OCR <- cbind(Os_filtered_OCR_counts_norm[1:4], Os_AH_MvF_OCR)
Os_AH_MvF_OCR_sig <- Os_AH_MvF_OCR %>% filter(padj <= 0.1)

#genitalia
Os_G_MvF_OCR <- as.data.frame(dds_Os_G_ATAC, independentFiltering=FALSE)
Os_G_MvF_OCR <- cbind(Os_filtered_OCR_counts_norm[1:4], Os_G_MvF_OCR)
Os_G_MvF_OCR_sig <- Os_G_MvF_OCR %>% filter(padj <= 0.1)

### Dgaz ####
#posterior head
Dg_PH_MvF_OCR <- as.data.frame(dds_Dg_PH_ATAC, independentFiltering=FALSE)
Dg_PH_MvF_OCR <- cbind(Dg_filtered_OCR_counts_norm[1:4], Dg_PH_MvF_OCR)
Dg_PH_MvF_OCR_sig <- Dg_PH_MvF_OCR %>% filter(padj <= 0.1)

#anterior head
Dg_AH_MvF_OCR <- as.data.frame(dds_Dg_AH_ATAC, independentFiltering=FALSE)
Dg_AH_MvF_OCR <- cbind(Dg_filtered_OCR_counts_norm[1:4], Dg_AH_MvF_OCR)
Dg_AH_MvF_OCR_sig <- Dg_AH_MvF_OCR %>% filter(padj <= 0.1)

#genitalia
Dg_G_MvF_OCR <- as.data.frame(dds_Dg_G_ATAC, independentFiltering=FALSE)
Dg_G_MvF_OCR <- cbind(Dg_filtered_OCR_counts_norm[1:4], Dg_G_MvF_OCR)
Dg_G_MvF_OCR_sig <- Dg_G_MvF_OCR %>% filter(padj <= 0.1)

### volcano plots of sex-responsive peaks ####
cols <- c("F" = "#C1272D", "M" = "#2166AC", "ns" = "grey") 
sizes <- c("F" = 2, "M" = 2, "ns" = 1) 
alphas <- c("F" = 0.5, "M" = 0.5, "ns" = 1)

### Dgaz ####
### genitalia 
# annotate DA peaks as more open in females or males
Dg_G_MvF_OCR <- Dg_G_MvF_OCR %>% 
  mutate(DA = ifelse(padj <= 0.1 & log2FoldChange < 0, "F",
                     ifelse(padj <= 0.1 & log2FoldChange > 0, "M",
                            "ns")))
Dg_G_MvF_OCR %>% dplyr::count(DA)

Dg_G_DA_vol_plot <- Dg_G_MvF_OCR %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), fill = DA, size = DA, alpha = DA)) + 
  geom_point(shape = 24, colour = "black", stroke = 0.25) +
  scale_fill_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  ylim(0,10) + 
  xlim(-2, 2) +
  theme_bw(base_size=14) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth= 0.5),    
        panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  labs(title="Dgaz genitalia OCRs")
Dg_G_DA_vol_plot

ggsave(Dg_G_DA_vol_plot, file = "/Users/ericanadolski/Desktop/beetle-sex-dimorph/figures/ATAC-Dg-G-volcano.pdf",
       width = 5, height = 4)

### posterior head
Dg_PH_MvF_OCR <- Dg_PH_MvF_OCR %>% 
  mutate(DA = ifelse(padj <= 0.1 & log2FoldChange < 0, "F",
                     ifelse(padj <= 0.1 & log2FoldChange > 0, "M",
                            "ns")))
Dg_PH_MvF_OCR %>% dplyr::count(DA)

Dg_PH_DA_vol_plot <- Dg_PH_MvF_OCR %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), fill = DA, size = DA, alpha = DA)) + 
  geom_point(shape = 21, colour = "black", stroke = 0.25) +
  scale_fill_manual(values = cols) + 
  scale_size_manual(values = sizes) +
  scale_alpha_manual(values = alphas) +
  ylim(0,10) + 
  xlim(-2, 2) +
  theme_bw(base_size=14) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth= 0.5),    
        panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  labs(title="Dgaz posterior head OCRs")
Dg_PH_DA_vol_plot

ggsave(Dg_PH_DA_vol_plot, file = "/Users/ericanadolski/Desktop/beetle-sex-dimorph/figures/ATAC-Dg-PH-volcano.pdf",
       width = 5, height = 4)

### anterior head 
Dg_AH_MvF_OCR <- Dg_AH_MvF_OCR %>% 
  mutate(DA = ifelse(padj <= 0.1 & log2FoldChange < 0, "F",
                     ifelse(padj <= 0.1 & log2FoldChange > 0, "M",
                            "ns")))
Dg_AH_MvF_OCR %>% dplyr::count(DA)

Dg_AH_DA_vol_plot <- Dg_AH_MvF_OCR %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), fill = DA, size = DA, alpha = DA)) + 
  geom_point(shape = 22, colour = "black", stroke = 0.25) +
  scale_fill_manual(values = cols) + 
  scale_size_manual(values = sizes) +
  scale_alpha_manual(values = alphas) +
  ylim(0,10) + 
  xlim(-2, 2) +
  theme_bw(base_size=14) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth= 0.5),    
        panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  labs(title="Dgaz anterior head OCRs")
Dg_AH_DA_vol_plot

ggsave(Dg_AH_DA_vol_plot, file = "/Users/ericanadolski/Desktop/beetle-sex-dimorph/figures/ATAC-Dg-AH-volcano.pdf",
       width = 5, height = 4)

### Otau ####
### genitalia 
# annotate DA peaks as more open in females or males
Ot_G_MvF_OCR <- Ot_G_MvF_OCR %>% 
  mutate(DA = ifelse(padj <= 0.1 & log2FoldChange < 0, "F",
                     ifelse(padj <= 0.1 & log2FoldChange > 0, "M",
                            "ns")))
Ot_G_MvF_OCR %>% dplyr::count(DA)

Ot_G_DA_vol_plot <- Ot_G_MvF_OCR %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), fill = DA, size = DA, alpha = DA)) + 
  geom_point(shape = 24, colour = "black", stroke = 0.25) +
  scale_fill_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  ylim(0,10) + 
  xlim(-2, 2) +
  theme_bw(base_size=14) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth= 0.5),    
        panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  labs(title="Otau genitalia OCRs")
Ot_G_DA_vol_plot

ggsave(Ot_G_DA_vol_plot, file = "/Users/ericanadolski/Desktop/beetle-sex-dimorph/figures/ATAC-Ot-G-volcano.pdf",
       width = 5, height = 4)

### posterior head
Ot_PH_MvF_OCR <- Ot_PH_MvF_OCR %>% 
  mutate(DA = ifelse(padj <= 0.1 & log2FoldChange < 0, "F",
                     ifelse(padj <= 0.1 & log2FoldChange > 0, "M",
                            "ns")))
Ot_PH_MvF_OCR %>% dplyr::count(DA)

Ot_PH_DA_vol_plot <- Ot_PH_MvF_OCR %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), fill = DA, size = DA, alpha = DA)) + 
  geom_point(shape = 21, colour = "black", stroke = 0.25) +
  scale_fill_manual(values = cols) + 
  scale_size_manual(values = sizes) +
  scale_alpha_manual(values = alphas) +
  ylim(0,10) + 
  xlim(-2, 2) +
  theme_bw(base_size=14) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth= 0.5),    
        panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  labs(title="Otau posterior head OCRs")
Ot_PH_DA_vol_plot

ggsave(Ot_PH_DA_vol_plot, file = "/Users/ericanadolski/Desktop/beetle-sex-dimorph/figures/ATAC-Ot-PH-volcano.pdf",
       width = 5, height = 4)

### anterior head 
Ot_AH_MvF_OCR <- Ot_AH_MvF_OCR %>% 
  mutate(DA = ifelse(padj <= 0.1 & log2FoldChange < 0, "F",
                     ifelse(padj <= 0.1 & log2FoldChange > 0, "M",
                            "ns")))
Ot_AH_MvF_OCR %>% dplyr::count(DA)

Ot_AH_DA_vol_plot <- Ot_AH_MvF_OCR %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), fill = DA, size = DA, alpha = DA)) + 
  geom_point(shape = 22, colour = "black", stroke = 0.25) +
  scale_fill_manual(values = cols) + 
  scale_size_manual(values = sizes) +
  scale_alpha_manual(values = alphas) +
  ylim(0,10) + 
  xlim(-2, 2) +
  theme_bw(base_size=14) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth= 0.5),    
        panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  labs(title="Otau anterior head OCRs")
Ot_AH_DA_vol_plot

ggsave(Ot_AH_DA_vol_plot, file = "/Users/ericanadolski/Desktop/beetle-sex-dimorph/figures/ATAC-Ot-AH-volcano.pdf",
       width = 5, height = 4)

### Osag ####
### genitalia 
Os_G_MvF_OCR <- Os_G_MvF_OCR %>% 
  mutate(DA = ifelse(padj <= 0.1 & log2FoldChange < 0, "F",
                     ifelse(padj <= 0.1 & log2FoldChange > 0, "M",
                            "ns")))
Os_G_MvF_OCR %>% dplyr::count(DA)

Os_G_DA_vol_plot <- Os_G_MvF_OCR %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), fill = DA, size = DA, alpha = DA)) + 
  geom_point(shape = 24, colour = "black", stroke = 0.25) +
  scale_fill_manual(values = cols) +
  scale_size_manual(values = sizes) +
  scale_alpha_manual(values = alphas) + 
  ylim(0,10) + 
  xlim(-2, 2) +
  theme_bw(base_size=14) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth= 0.5),    
        panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  labs(title="Osag genitalia OCRs")
Os_G_DA_vol_plot

ggsave(Os_G_DA_vol_plot, file = "/Users/ericanadolski/Desktop/beetle-sex-dimorph/figures/ATAC-Os-G-volcano.pdf",
       width = 5, height = 4)

### posterior head
Os_PH_MvF_OCR <- Os_PH_MvF_OCR %>% 
  mutate(DA = ifelse(padj <= 0.1 & log2FoldChange < 0, "F",
                     ifelse(padj <= 0.1 & log2FoldChange > 0, "M",
                            "ns")))
Os_PH_MvF_OCR %>% dplyr::count(DA)

Os_PH_DA_vol_plot <- Os_PH_MvF_OCR %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), fill = DA, size = DA, alpha = DA)) + 
  geom_point(shape = 21, colour = "black", stroke = 0.25) +
  scale_fill_manual(values = cols) + 
  scale_size_manual(values = sizes) +
  scale_alpha_manual(values = alphas) +
  ylim(0,10) + 
  xlim(-2, 2) +
  theme_bw(base_size=14) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth= 0.5),    
        panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  labs(title="Osag posterior head OCRs")
Os_PH_DA_vol_plot

ggsave(Os_PH_DA_vol_plot, file = "/Users/ericanadolski/Desktop/beetle-sex-dimorph/figures/ATAC-Os-PH-volcano.pdf",
       width = 5, height = 4)

### anterior head 
Os_AH_MvF_OCR <- Os_AH_MvF_OCR %>% 
  mutate(DA = ifelse(padj <= 0.1 & log2FoldChange < 0, "F",
                     ifelse(padj <= 0.1 & log2FoldChange > 0, "M",
                            "ns")))
Os_AH_MvF_OCR %>% dplyr::count(DA)

Os_AH_DA_vol_plot <- Os_AH_MvF_OCR %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), fill = DA, size = DA, alpha = DA)) + 
  geom_point(shape = 22, colour = "black", stroke = 0.25) +
  scale_fill_manual(values = cols) + 
  scale_size_manual(values = sizes) +
  scale_alpha_manual(values = alphas) +
  ylim(0,10) + 
  xlim(-2, 2) +
  theme_bw(base_size=14) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth= 0.5),    
        panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  labs(title="Osag anterior head OCRs")
Os_AH_DA_vol_plot

ggsave(Os_AH_DA_vol_plot, file = "/Users/ericanadolski/Desktop/beetle-sex-dimorph/figures/ATAC-Os-AH-volcano.pdf",
       width = 5, height = 4)
### clustered heatmaps of sex-responsive peaks ####
myheatcolors <- rev(brewer.pal(name="RdBu", n=11))
### Dgaz ####
# genitalia 
# subset count matrix by DEGs and trait for heatmap 
Dg_G_filtered_counts_ATAC <- cbind(Dg_filtered_OCR_counts_norm[1:4], Dg_G_filtered_counts_ATAC)
Dg_G_DA_ocr <- inner_join(Dg_G_MvF_OCR_sig, Dg_G_filtered_counts_ATAC, by="peak")

clustRows <- hclust(as.dist(1-cor(t(Dg_G_DA_ocr[,14:23]), method="pearson")), method="complete")
clustColumns <- hclust(as.dist(1-cor(Dg_G_DA_ocr[,14:23], method="spearman")), method="complete")
module.assign <- cutree(clustRows, k=2)
module.color <- rainbow(length(unique(module.assign)), start=0.1, end=0.9) 
module.color <- module.color[as.vector(module.assign)] 
Dg_G_DA_heatmap <- heatmap.2(as.matrix(Dg_G_DA_ocr[,14:23]), 
                             Rowv=as.dendrogram(clustRows), 
                             Colv=as.dendrogram(clustColumns),
                             RowSideColors=module.color,
                             col=myheatcolors, scale='row', labRow=NA,
                             density.info="none", trace="none",)

### posterior head 
Dg_PH_filtered_counts_ATAC <- cbind(Dg_filtered_OCR_counts_norm[1:4], Dg_PH_filtered_counts_ATAC)
Dg_PH_DA_ocr <- inner_join(Dg_PH_MvF_OCR_sig, Dg_PH_filtered_counts_ATAC, by="peak")

clustRows <- hclust(as.dist(1-cor(t(Dg_PH_DA_ocr[,14:23]), method="pearson")), method="complete")
clustColumns <- hclust(as.dist(1-cor(Dg_PH_DA_ocr[,14:23], method="spearman")), method="complete")
module.assign <- cutree(clustRows, k=2)
module.color <- rainbow(length(unique(module.assign)), start=0.1, end=0.9) 
module.color <- module.color[as.vector(module.assign)] 
Dg_PH_DA_heatmap <- heatmap.2(as.matrix(Dg_PH_DA_ocr[,14:23]), 
                              Rowv=as.dendrogram(clustRows), 
                              Colv=as.dendrogram(clustColumns),
                              RowSideColors=module.color,
                              col=myheatcolors, scale='row', labRow=NA,
                              density.info="none", trace="none",)

### anterior head
Dg_AH_filtered_counts_ATAC <- cbind(Dg_filtered_OCR_counts_norm[1:4], Dg_AH_filtered_counts_ATAC)
Dg_AH_DA_ocr <- inner_join(Dg_AH_MvF_OCR_sig, Dg_AH_filtered_counts_ATAC, by="peak")

clustRows <- hclust(as.dist(1-cor(t(Dg_AH_DA_ocr[,14:23]), method="pearson")), method="complete")
clustColumns <- hclust(as.dist(1-cor(Dg_AH_DA_ocr[,14:23], method="spearman")), method="complete")
module.assign <- cutree(clustRows, k=2)
module.color <- rainbow(length(unique(module.assign)), start=0.1, end=0.9) 
module.color <- module.color[as.vector(module.assign)] 
Dg_AH_DA_heatmap <- heatmap.2(as.matrix(Dg_AH_DA_ocr[,14:23]), 
                              Rowv=as.dendrogram(clustRows), 
                              Colv=as.dendrogram(clustColumns),
                              RowSideColors=module.color,
                              col=myheatcolors, scale='row', labRow=NA,
                              density.info="none", trace="none",)
### Otau ####
# genitalia 
# subset count matrix by DEGs and trait for heatmap 
Ot_G_filtered_counts_ATAC <- cbind(Ot_filtered_OCR_counts_norm[1:4], Ot_G_filtered_counts_ATAC)
Ot_G_DA_ocr <- inner_join(Ot_G_MvF_OCR_sig, Ot_G_filtered_counts_ATAC, by="peak")

clustRows <- hclust(as.dist(1-cor(t(Ot_G_DA_ocr[,14:23]), method="pearson")), method="complete")
clustColumns <- hclust(as.dist(1-cor(Ot_G_DA_ocr[,14:23], method="spearman")), method="complete")
module.assign <- cutree(clustRows, k=2)
module.color <- rainbow(length(unique(module.assign)), start=0.1, end=0.9) 
module.color <- module.color[as.vector(module.assign)] 
Ot_G_DA_heatmap <- heatmap.2(as.matrix(Ot_G_DA_ocr[,14:23]), 
                             Rowv=as.dendrogram(clustRows), 
                             Colv=as.dendrogram(clustColumns),
                             RowSideColors=module.color,
                             col=myheatcolors, scale='row', labRow=NA,
                             density.info="none", trace="none",)

### posterior head 
Ot_PH_filtered_counts_ATAC <- cbind(Ot_filtered_OCR_counts_norm[1:4], Ot_PH_filtered_counts_ATAC)
Ot_PH_DA_ocr <- inner_join(Ot_PH_MvF_OCR_sig, Ot_PH_filtered_counts_ATAC, by="peak")

clustRows <- hclust(as.dist(1-cor(t(Ot_PH_DA_ocr[,14:23]), method="pearson")), method="complete")
clustColumns <- hclust(as.dist(1-cor(Ot_PH_DA_ocr[,14:23], method="spearman")), method="complete")
module.assign <- cutree(clustRows, k=2)
module.color <- rainbow(length(unique(module.assign)), start=0.1, end=0.9) 
module.color <- module.color[as.vector(module.assign)] 
Ot_PH_DA_heatmap <- heatmap.2(as.matrix(Ot_PH_DA_ocr[,14:23]), 
                              Rowv=as.dendrogram(clustRows), 
                              Colv=as.dendrogram(clustColumns),
                              RowSideColors=module.color,
                              col=myheatcolors, scale='row', labRow=NA,
                              density.info="none", trace="none",)

### anterior head
Ot_AH_filtered_counts_ATAC <- cbind(Ot_filtered_OCR_counts_norm[1:4], Ot_AH_filtered_counts_ATAC)
Ot_AH_DA_ocr <- inner_join(Ot_AH_MvF_OCR_sig, Ot_AH_filtered_counts_ATAC, by="peak")

clustRows <- hclust(as.dist(1-cor(t(Ot_AH_DA_ocr[,14:23]), method="pearson")), method="complete")
clustColumns <- hclust(as.dist(1-cor(Ot_AH_DA_ocr[,14:23], method="spearman")), method="complete")
module.assign <- cutree(clustRows, k=2)
module.color <- rainbow(length(unique(module.assign)), start=0.1, end=0.9) 
module.color <- module.color[as.vector(module.assign)] 
Ot_AH_DA_heatmap <- heatmap.2(as.matrix(Ot_AH_DA_ocr[,14:23]), 
                              Rowv=as.dendrogram(clustRows), 
                              Colv=as.dendrogram(clustColumns),
                              RowSideColors=module.color,
                              col=myheatcolors, scale='row', labRow=NA,
                              density.info="none", trace="none",)
### Osag ####
# genitalia 
# subset count matrix by DEGs and trait for heatmap 
Os_G_filtered_counts_ATAC <- cbind(Os_filtered_OCR_counts_norm[1:4], Os_G_filtered_counts_ATAC)
Os_G_DA_ocr <- inner_join(Os_G_MvF_OCR_sig, Os_G_filtered_counts_ATAC, by="peak")

clustRows <- hclust(as.dist(1-cor(t(Os_G_DA_ocr[,14:23]), method="pearson")), method="complete")
clustColumns <- hclust(as.dist(1-cor(Os_G_DA_ocr[,14:23], method="spearman")), method="complete")
module.assign <- cutree(clustRows, k=2)
module.color <- rainbow(length(unique(module.assign)), start=0.1, end=0.9) 
module.color <- module.color[as.vector(module.assign)] 
Os_G_DA_heatmap <- heatmap.2(as.matrix(Os_G_DA_ocr[,14:23]), 
                             Rowv=as.dendrogram(clustRows), 
                             Colv=as.dendrogram(clustColumns),
                             RowSideColors=module.color,
                             col=myheatcolors, scale='row', labRow=NA,
                             density.info="none", trace="none",)

### posterior head 
Os_PH_filtered_counts_ATAC <- cbind(Os_filtered_OCR_counts_norm[1:4], Os_PH_filtered_counts_ATAC)
Os_PH_DA_ocr <- inner_join(Os_PH_MvF_OCR_sig, Os_PH_filtered_counts_ATAC, by="peak")

clustRows <- hclust(as.dist(1-cor(t(Os_PH_DA_ocr[,14:23]), method="pearson")), method="complete")
clustColumns <- hclust(as.dist(1-cor(Os_PH_DA_ocr[,14:23], method="spearman")), method="complete")
module.assign <- cutree(clustRows, k=2)
module.color <- rainbow(length(unique(module.assign)), start=0.1, end=0.9) 
module.color <- module.color[as.vector(module.assign)] 
Os_PH_DA_heatmap <- heatmap.2(as.matrix(Os_PH_DA_ocr[,14:23]), 
                              Rowv=as.dendrogram(clustRows), 
                              Colv=as.dendrogram(clustColumns),
                              RowSideColors=module.color,
                              col=myheatcolors, scale='row', labRow=NA,
                              density.info="none", trace="none",)

### anterior head
Os_AH_filtered_counts_ATAC <- cbind(Os_filtered_OCR_counts_norm[1:4], Os_AH_filtered_counts_ATAC)
Os_AH_DA_ocr <- inner_join(Os_AH_MvF_OCR_sig, Os_AH_filtered_counts_ATAC, by="peak")

clustRows <- hclust(as.dist(1-cor(t(Os_AH_DA_ocr[,14:23]), method="pearson")), method="complete")
clustColumns <- hclust(as.dist(1-cor(Os_AH_DA_ocr[,14:23], method="spearman")), method="complete")
module.assign <- cutree(clustRows, k=2)
module.color <- rainbow(length(unique(module.assign)), start=0.1, end=0.9) 
module.color <- module.color[as.vector(module.assign)] 
Os_AH_DA_heatmap <- heatmap.2(as.matrix(Os_AH_DA_ocr[,14:23]), 
                              Rowv=as.dendrogram(clustRows), 
                              Colv=as.dendrogram(clustColumns),
                              RowSideColors=module.color,
                              col=myheatcolors, scale='row', labRow=NA,
                              density.info="none", trace="none",)
###### Step 5: plot bar chart of sex-responsive OCRs across traits ####
species <- c("Dgaz","Dgaz","Dgaz","Dgaz","Dgaz","Dgaz","Otau","Otau","Otau","Otau","Otau","Otau","Osag","Osag","Osag","Osag","Osag","Osag")
sex <- c('F', 'M','F', 'M','F', 'M','F', 'M','F', 'M','F', 'M','F', 'M','F', 'M','F', 'M') 
trait <- c('G','G','PH', 'PH','AH','AH','G','G','PH', 'PH','AH','AH','G','G','PH', 'PH','AH','AH')
DA_peaks <- c(8,5,8,25,5,29,212,892,981,1579,4,16,245,877,13,149,43,170)

DA_peak_numbers <- data.frame(species,sex,trait,DA_peaks)
head(DA_peak_numbers)
DA_peak_numbers$trait = factor(DA_peak_numbers$trait, levels = c('G',"PH","AH"), ordered = TRUE)
DA_peak_numbers$species = factor(DA_peak_numbers$species, levels = c('Dgaz',"Otau","Osag"), ordered = TRUE)

bar_plot_OCR <-ggplot(data=DA_peak_numbers, aes(x=species, y=DA_peaks, fill=sex)) +
  geom_bar(stat="identity", position=position_dodge()) + 
  facet_grid(~trait) +
  scale_y_break(c(1000,1500)) +
  geom_text(aes(label=DA_peaks), vjust = -0.5, size=3.5, position = position_dodge(0.9))+ # outside bars
  scale_fill_manual(values=c("#C1272D","#2166AC"))+
  labs(y= "number of differentially accessible OCRs")+
  theme_classic(base_size=14)
#coord_flip()
bar_plot_OCR

ggsave(bar_plot_OCR, file = "/Users/ericanadolski/Desktop/beetle-sex-dimorph/figures/bar_plot_OCR.pdf", width = 10, height = 5)

###### Step 5.5: differential accessibility across traits - needed for paper? ########
### Osag both sexes, trait responsive  ######
dds_Os<- DESeqDataSetFromMatrix(countData = Os_filtered_counts[5:54], 
                                colData = Os_info, 
                                design = ~Trait)
dds_Os<- DESeq(dds_Os)

# PH vs AH
Os_AHPH <- as.data.frame(results(dds_Os, contrast=c("Trait","AH","PH"), independentFiltering=FALSE))
Os_AHPH <- cbind(Os_filtered_counts[1:4], Os_AHPH)
Os_AHPH_sig1 <- Os_AHPH %>% filter(padj <= 0.01)
tally(as.data.frame(results(dds_Os, contrast=c("Trait","AH", "PH"), independentFiltering=FALSE)) %>%
        filter(padj <= 0.01 )) # 705

# PH vs L
Os_LPH <- as.data.frame(results(dds_Os, contrast=c("Trait","L","PH"), independentFiltering=FALSE))
Os_LPH <- cbind(Os_filtered_counts[1:4], Os_LPH)
Os_LPH_sig1 <- Os_LPH %>% filter(padj <= 0.01)
tally(as.data.frame(results(dds_Os, contrast=c("Trait","L","PH"), independentFiltering=FALSE)) %>% 
        filter(padj <= 0.01 )) # 5517

# PH vs E
Os_EPH <- as.data.frame(results(dds_Os, contrast=c("Trait","E","PH"), independentFiltering=FALSE))
Os_EPH <- cbind(Os_filtered_counts[1:4], Os_EPH)
Os_EPH_sig1 <- Os_EPH %>% filter(padj <= 0.01)
tally(as.data.frame(results(dds_Os, contrast=c("Trait","E","PH"), independentFiltering=FALSE)) %>% 
        filter(padj <= 0.01 )) # 9781

# PH vs G
Os_GPH <- as.data.frame(results(dds_Os, contrast=c("Trait","G","PH"), independentFiltering=FALSE))
Os_GPH <- cbind(Os_filtered_counts[1:4], Os_GPH)
Os_GPH_sig1 <- Os_GPH %>% filter(padj <= 0.01)
tally(as.data.frame(results(dds_Os, contrast=c("Trait","G","PH"), independentFiltering=FALSE)) %>% 
        filter(padj <= 0.01 )) # 10516

# AH vs G
Os_GAH <- as.data.frame(results(dds_Os, contrast=c("Trait","G","AH"), independentFiltering=FALSE))
Os_GAH <- cbind(Os_filtered_counts[1:4], Os_GAH)
Os_GAH_sig1 <- Os_GAH %>% filter(padj <= 0.01)
tally(as.data.frame(results(dds_Os, contrast=c("Trait","G","AH"), independentFiltering=FALSE)) %>% 
        filter(padj <= 0.01 )) # 18401

# AH vs E
Os_EAH <- as.data.frame(results(dds_Os, contrast=c("Trait","E","AH"), independentFiltering=FALSE))
Os_EAH <- cbind(Os_filtered_counts[1:4],Os_EAH)
Os_EAH_sig1 <- Os_EAH %>% filter(padj <= 0.01)
tally(as.data.frame(results(dds_Os, contrast=c("Trait","E","AH"), independentFiltering=FALSE)) %>% 
        filter(padj <= 0.01 )) # 16177

# AH vs L
Os_LAH <- as.data.frame(results(dds_Os, contrast=c("Trait","L","AH"), independentFiltering=FALSE))
Os_LAH <- cbind(Os_filtered_counts[1:4], Os_LAH)
Os_LAH_sig1 <- Os_LAH %>% filter(padj <= 0.01)
tally(as.data.frame(results(dds_Os, contrast=c("Trait","L","AH"), independentFiltering=FALSE)) %>% 
        filter(padj <= 0.01 )) # 8528

# E vs L
Os_LE <- as.data.frame(results(dds_Os, contrast=c("Trait","L","E"), independentFiltering=FALSE))
Os_LE <- cbind(Os_filtered_counts[1:4], Os_LE)
Os_LE_sig1 <- Os_LE %>% filter(padj <= 0.01)
tally(as.data.frame(results(dds_Os, contrast=c("Trait","L","E"), independentFiltering=FALSE)) %>% 
        filter(padj <= 0.01 )) # 9821

# E vs G
Os_GE <- as.data.frame(results(dds_Os, contrast=c("Trait","G","E"), independentFiltering=FALSE))
Os_GE <- cbind(Os_filtered_counts[1:4], Os_GE)
Os_GE_sig1 <- Os_GE %>% filter(padj <= 0.01)
tally(as.data.frame(results(dds_Os, contrast=c("Trait","G","E"), independentFiltering=FALSE)) %>% 
        filter(padj <= 0.01 )) # 13617

# G vs L
Os_LG <- as.data.frame(results(dds_Os, contrast=c("Trait","L","G"), independentFiltering=FALSE))
Os_LG <- cbind(Os_filtered_counts[1:4], Os_LG)
Os_LG_sig1 <- Os_LG %>% filter(padj <= 0.01)
tally(as.data.frame(results(dds_Os, contrast=c("Trait","L","G"), independentFiltering=FALSE)) %>% 
        filter(padj <= 0.01 )) # 11487

### Dgaz both sexes, trait responsive  ######
dds_Dg<- DESeqDataSetFromMatrix(countData = Dg_filtered_counts[5:54], 
                                colData = Dg_info, 
                                design = ~Trait)
dds_Dg<- DESeq(dds_Dg)

# PH vs AH
Dg_AHPH <- as.data.frame(results(dds_Dg, contrast=c("Trait","AH","PH"), independentFiltering=FALSE))
Dg_AHPH <- cbind(Dg_filtered_counts[1:4], Dg_AHPH)
Dg_AHPH_sig1 <- Dg_AHPH %>% filter(padj <= 0.01)
tally(as.data.frame(results(dds_Dg, contrast=c("Trait","AH", "PH"), independentFiltering=FALSE)) %>%
        filter(padj <= 0.01 )) # 17

# PH vs L
Dg_LPH <- as.data.frame(results(dds_Dg, contrast=c("Trait","L","PH"), independentFiltering=FALSE))
Dg_LPH <- cbind(Dg_filtered_counts[1:4], Dg_LPH)
Dg_LPH_sig1 <- Dg_LPH %>% filter(padj <= 0.01)
tally(as.data.frame(results(dds_Dg, contrast=c("Trait","L","PH"), independentFiltering=FALSE)) %>% 
        filter(padj <= 0.01 )) # 10542

# PH vs E
Dg_EPH <- as.data.frame(results(dds_Dg, contrast=c("Trait","E","PH"), independentFiltering=FALSE))
Dg_EPH <- cbind(Dg_filtered_counts[1:4], Dg_EPH)
Dg_EPH_sig1 <- Dg_EPH %>% filter(padj <= 0.01)
tally(as.data.frame(results(dds_Dg, contrast=c("Trait","E","PH"), independentFiltering=FALSE)) %>% 
        filter(padj <= 0.01 )) # 4911

# PH vs G
Dg_GPH <- as.data.frame(results(dds_Dg, contrast=c("Trait","G","PH"), independentFiltering=FALSE))
Dg_GPH <- cbind(Dg_filtered_counts[1:4], Dg_GPH)
Dg_GPH_sig1 <-Dg_GPH %>% filter(padj <= 0.01)
tally(as.data.frame(results(dds_Dg, contrast=c("Trait","G","PH"), independentFiltering=FALSE)) %>% 
        filter(padj <= 0.01 )) # 380

# AH vs G
Dg_GAH <- as.data.frame(results(dds_Dg, contrast=c("Trait","G","AH"), independentFiltering=FALSE))
Dg_GAH <- cbind(Dg_filtered_counts[1:4], Dg_GAH)
Dg_GAH_sig1 <- Dg_GAH %>% filter(padj <= 0.01)
tally(as.data.frame(results(dds_Dg, contrast=c("Trait","G","AH"), independentFiltering=FALSE)) %>% 
        filter(padj <= 0.01 )) # 679

# AH vs E
Dg_EAH <- as.data.frame(results(dds_Dg, contrast=c("Trait","E","AH"), independentFiltering=FALSE))
Dg_EAH <- cbind(Dg_filtered_counts[1:4], Dg_EAH)
Dg_EAH_sig1 <- Dg_EAH %>% filter(padj <= 0.01)
tally(as.data.frame(results(dds_Dg, contrast=c("Trait","E","AH"), independentFiltering=FALSE)) %>% 
        filter(padj <= 0.01 )) # 9854

# AH vs L
Dg_LAH <- as.data.frame(results(dds_Dg, contrast=c("Trait","L","AH"), independentFiltering=FALSE))
Dg_LAH <- cbind(Dg_filtered_counts[1:4], Dg_LAH)
Dg_LAH_sig1 <- Dg_LAH %>% filter(padj <= 0.01)
tally(as.data.frame(results(dds_Dg, contrast=c("Trait","L","AH"), independentFiltering=FALSE)) %>% 
        filter(padj <= 0.01 )) # 18432

# E vs L
Dg_LE <- as.data.frame(results(dds_Dg, contrast=c("Trait","L","E"), independentFiltering=FALSE))
Dg_LE <- cbind(Dg_filtered_counts[1:4], Dg_LE)
Dg_LE_sig1 <- Dg_LE %>% filter(padj <= 0.01)
tally(as.data.frame(results(dds_Dg, contrast=c("Trait","L","E"), independentFiltering=FALSE)) %>% 
        filter(padj <= 0.01 )) # 52

# E vs G
Dg_GE <- as.data.frame(results(dds_Dg, contrast=c("Trait","G","E"), independentFiltering=FALSE))
Dg_GE <- cbind(Dg_filtered_counts[1:4], Dg_GE)
Dg_GE_sig1 <- Dg_GE %>% filter(padj <= 0.01)
tally(as.data.frame(results(dds_Dg, contrast=c("Trait","G","E"), independentFiltering=FALSE)) %>% 
        filter(padj <= 0.01 )) # 4003

# G vs L
Dg_LG <- as.data.frame(results(dds_Dg, contrast=c("Trait","L","G"), independentFiltering=FALSE))
Dg_LG <- cbind(Dg_filtered_counts[1:4], Dg_LG)
Dg_LG_sig1 <- Dg_LG %>% filter(padj <= 0.01)
tally(as.data.frame(results(dds_Dg, contrast=c("Trait","L","G"), independentFiltering=FALSE)) %>% 
        filter(padj <= 0.01 )) # 5801

### Otau both sexes, trait responsive  ######
dds_Ot<- DESeqDataSetFromMatrix(countData = Ot_filtered_counts[5:54], 
                                colData = Ot_info, design = ~Trait)
dds_Ot<- DESeq(dds_Ot)

# PH vs AH
Ot_AHPH <- as.data.frame(results(dds_Ot, contrast=c("Trait","AH","PH"), independentFiltering=FALSE))
Ot_AHPH <- cbind(Ot_filtered_counts[1:4], Ot_AHPH)
Ot_AHPH_sig1 <- Ot_AHPH %>% filter(padj <= 0.01)
tally(as.data.frame(results(dds_Ot, contrast=c("Trait","AH", "PH"), independentFiltering=FALSE)) %>% 
        filter(padj <= 0.01 )) # 1208

# PH vs L
Ot_LPH <- as.data.frame(results(dds_Ot, contrast=c("Trait","L","PH"), independentFiltering=FALSE))
Ot_LPH <- cbind(Ot_filtered_counts[1:4], Ot_LPH)
Ot_LPH_sig1 <- Ot_LPH %>% filter(padj <= 0.01)
tally(as.data.frame(results(dds_Ot, contrast=c("Trait","L","PH"), independentFiltering=FALSE)) %>% 
        filter(padj <= 0.01 )) # 2877

# PH vs E
Ot_EPH <- as.data.frame(results(dds_Ot, contrast=c("Trait","E","PH"), independentFiltering=FALSE))
Ot_EPH <- cbind(Ot_filtered_counts[1:4], Ot_EPH)
Ot_EPH_sig1 <- Ot_EPH %>% filter(padj <= 0.01)
tally(as.data.frame(results(dds_Ot, contrast=c("Trait","E","PH"), independentFiltering=FALSE)) %>% 
        filter(padj <= 0.01 )) # 3240

# PH vs G
Ot_GPH <- as.data.frame(results(dds_Ot, contrast=c("Trait","G","PH"), independentFiltering=FALSE))
Ot_GPH <- cbind(Ot_filtered_counts[1:4],Ot_GPH)
Ot_GPH_sig1 <- Ot_GPH %>% filter(padj <= 0.01)
tally(as.data.frame(results(dds_Ot, contrast=c("Trait","G","PH"), independentFiltering=FALSE)) %>% 
        filter(padj <= 0.01 )) # 8630

# AH vs G
Ot_GAH <- as.data.frame(results(dds_Ot, contrast=c("Trait","G","AH"), independentFiltering=FALSE))
Ot_GAH <- cbind(Ot_filtered_counts[1:4], Ot_GAH)
Ot_GAH_sig1 <- Ot_GAH %>% filter(padj <= 0.01)
tally(as.data.frame(results(dds_Ot, contrast=c("Trait","G","AH"), independentFiltering=FALSE)) %>% 
        filter(padj <= 0.01 )) # 9369

# AH vs E
Ot_EAH <- as.data.frame(results(dds_Ot, contrast=c("Trait","E","AH"), independentFiltering=FALSE))
Ot_EAH <- cbind(Ot_filtered_counts[1:4],Ot_EAH)
Ot_EAH_sig1 <- Ot_EAH %>% filter(padj <= 0.01)
tally(as.data.frame(results(dds_Ot, contrast=c("Trait","E","AH"), independentFiltering=FALSE)) %>% 
        filter(padj <= 0.01 )) # 3512

# AH vs L
Ot_LAH <- as.data.frame(results(dds_Ot, contrast=c("Trait","L","AH"), independentFiltering=FALSE))
Ot_LAH <- cbind(Ot_filtered_counts[1:4], Ot_LAH)
Ot_LAH_sig1 <- Ot_LAH %>% filter(padj <= 0.01)
tally(as.data.frame(results(dds_Ot, contrast=c("Trait","L","AH"), independentFiltering=FALSE)) %>% 
        filter(padj <= 0.01 )) # 3192

# E vs L
Ot_LE <- as.data.frame(results(dds_Ot, contrast=c("Trait","L","E"), independentFiltering=FALSE))
Ot_LE <- cbind(Ot_filtered_counts[1:4], Ot_LE)
Ot_LE_sig1 <- Ot_LE %>% filter(padj <= 0.01)
tally(as.data.frame(results(dds_Ot, contrast=c("Trait","L","E"), independentFiltering=FALSE)) %>% 
        filter(padj <= 0.01 )) # 3836

# E vs G
Ot_GE <- as.data.frame(results(dds_Ot, contrast=c("Trait","G","E"), independentFiltering=FALSE))
Ot_GE <- cbind(Ot_filtered_counts[1:4], Ot_GE)
Ot_GE_sig1 <- Ot_GE %>% filter(padj <= 0.01)
tally(as.data.frame(results(dds_Ot, contrast=c("Trait","G","E"), independentFiltering=FALSE)) %>% 
        filter(padj <= 0.01 )) # 6070

# G vs L
Ot_LG <- as.data.frame(results(dds_Ot, contrast=c("Trait","L","G"), independentFiltering=FALSE))
Ot_LG <- cbind(Ot_filtered_counts[1:4], Ot_LG)
Ot_LG_sig1 <- Ot_LG %>% filter(padj <= 0.01)
tally(as.data.frame(results(dds_Ot, contrast=c("Trait","L","G"), independentFiltering=FALSE)) %>% 
        filter(padj <= 0.01 )) # 6146

###### Step 6: annotate peaks with information on nearest genes ######
#### read in genomic transcript coordinates ###
Ot_transcript_coords <- read.delim("/Users/ericanadolski/Documents/Genomes/Otau3/Otau_transcript_coords.txt")
# correct start and end positions for the negative strand transcripts
Ot_transcript_starts <- Ot_transcript_coords %>% 
  mutate(start_correct = ifelse(orientation == "+", start, end),
         end_correct = ifelse(orientation == "+",end, start)) %>% 
  select(gene, chr, start_gene = start_correct, orientation)

### anterior head ####
# denote which condition showed higher accessibility
OtAH_sex_res_OCR <- OtAH_sex_res_OCR %>% 
  mutate(DA = ifelse(padj <= 0.1 & log2FoldChange < 0, "F",
                     ifelse(padj <= 0.1 & log2FoldChange > 0, "M",
                            "ns")))

OtAH_OCR_sig_location  <- OtAH_sex_res_OCR %>%
  filter(padj <= 0.1) %>% mutate(mid = (end-start)/2 + start) %>% select(peak, chr, DA, start, mid, end)

# combine peak list with nearby genes
OtAH_DA_peak_gene_map <- 
  left_join(OtAH_OCR_sig_location, Ot_transcript_starts, by = "chr", relationship = "many-to-many") %>% 
  mutate(peakGeneDist = mid-start_gene) %>% 
  filter(abs(peakGeneDist) <= 25000) %>% # peak should be within 25kb of gene
  group_by(peak) %>% top_n(-2, abs(peakGeneDist)) #selecting 2 closest genes

# combine peak-gene map with protein annotations
OtAH_DA_peak_gene_map_anno <- left_join(OtAH_DA_peak_gene_map, Otau3_prot_anno, by = "gene")
OtAH_F_peak_map_anno <- OtAH_DA_peak_gene_map_anno %>% filter(DA == "F")
OtAH_M_peak_map_anno <- OtAH_DA_peak_gene_map_anno %>% filter(DA == "M")
###### PRE-ANNOTATION STEPS #####
### read in genome transcript coordinates ######

### Dgaz1 ### 
Dg_transcript_coords <- read.delim("/Users/ericanadolski/Documents/Genomes/Dgaz1/Dgaz_transcript_coords.txt")
# correct start and end positions for the negative strand transcripts
Dg_transcript_starts <- Dg_transcript_coords %>% 
  mutate(start_correct = ifelse(orientation == "+", start, end),
         end_correct = ifelse(orientation == "+",end, start)) %>% 
  select(gene, chr, start_gene = start_correct, orientation)
# now the Dg_transcript_starts file is ready for all annotation uses

### Osag1 ### 
Os_transcript_coords <- read.delim("/Users/ericanadolski/Documents/Genomes/Osag1/Osag_transcript_coords.txt")
# correct start and end positions for the negative strand transcripts
Os_transcript_starts <- Os_transcript_coords %>% 
  mutate(start_correct = ifelse(orientation == "+", start, end),
         end_correct = ifelse(orientation == "+",end, start)) %>% 
  select(gene, chr, start_gene = start_correct, orientation)
# now file is ready for all annotation uses

### Otau3 ### 
Ot_transcript_coords <- read.delim("/Users/ericanadolski/Documents/Genomes/Otau3/Otau_transcript_coords.txt")
# correct start and end positions for the negative strand transcripts
Ot_transcript_starts <- Ot_transcript_coords %>% 
  mutate(start_correct = ifelse(orientation == "+", start, end),
         end_correct = ifelse(orientation == "+",end, start)) %>% 
  select(gene, chr, start_gene = start_correct, orientation)
# now file is ready for all annotation uses

### read in annotated protein best hits ####

# get annotations of Dgaz1 protein best hits
Dgaz1_prot_hits_clean <- read.delim("/Users/ericanadolski/Documents/Genomes/Dgaz1/Dgaz1_prot_hits_clean.txt")
Dgaz1_prot_hits_clean$eval <- as.numeric(Dgaz1_prot_hits_clean$eval)
Dgaz1_prot_hits_best <- Dgaz1_prot_hits_clean %>% 
  group_by(DG3_ID) %>% 
  top_n(-1,eval) %>% 
  dplyr::slice(which.max(pident))

# get annotations of Osag1 protein best hits
Osag1_prot_hits_clean <- read.delim("/Users/ericanadolski/Documents/Genomes/Osag1/Osag1_prot_hits_clean.txt")
Osag1_prot_hits_clean$eval <- as.numeric(Osag1_prot_hits_clean$eval)
Osag1_prot_hits_best <- Osag1_prot_hits_clean %>% 
  group_by(OS1_ID) %>% 
  top_n(-1,eval) %>% 
  dplyr::slice(which.max(pident))


# get annotations of Otau3 protein best hits
Otau3_prot_hits_clean <- read.delim("/Users/ericanadolski/Documents/Genomes/Otau3/Otau3_prot_hits_clean.txt")
Otau3_prot_hits_clean$eval <- as.numeric(Otau3_prot_hits_clean$eval)
Otau3_prot_hits_best <- Otau3_prot_hits_clean %>% 
  group_by(OT3_ID) %>% 
  top_n(-1,eval) %>% 
  dplyr::slice(which.max(pident))

# combine with annotations from Otau2 genome
Otau2_prot_names <- read.delim("/Users/ericanadolski/Documents/Genomes/Otau3/Otau2_prot_names.txt")

Dgaz1_prot_anno <- left_join(Dgaz1_prot_hits_best, Otau2_prot_names)
colnames(Dgaz1_prot_anno)<- c("gene", "OT2_ID", "pident", "eval", "OT2_description")

Osag1_prot_anno <- left_join(Osag1_prot_hits_best, Otau2_prot_names)
colnames(Osag1_prot_anno)<- c("gene", "OT2_ID", "pident", "evalue","eval", "OT2_description")

Otau3_prot_anno <- left_join(Otau3_prot_hits_best, Otau2_prot_names)
colnames(Otau3_prot_anno)<- c("OT3_ID", "OT2_ID", "pident", "eval", "OT2_description")

###### PEAK ANNOTATION ######
###### Dgaz ###### 
### Dg AH sex res ####

### denote which condition showed higher accessibility - ordered in DESeq2
DgAH_sex_anno_all <- DgAH_sex_res %>% 
  mutate(DA = ifelse(padj <= 0.05 & log2FoldChange < 0, "F",
                     ifelse(padj <= 0.05 & log2FoldChange > 0, "M",
                            "notDA")))

DgAH_sex_res_sig_anno <- DgAH_sex_res_sig %>% 
  mutate(DA = ifelse(padj <= 0.05 & log2FoldChange < 0, "F",
                     ifelse(padj <= 0.05 & log2FoldChange > 0, "M",
                            "notDA")))

DgAH_sex_peaks <- DgAH_sex_res_sig_anno %>% mutate(mid = (end-start)/2 + start) %>% 
  select(peak, chr, DA, start, mid, end)

# fix peaks df so chr names are matching
#### can use this code
#df$chr <- gsub(';HRSCAF=\\d{1,2}','', df$chr)
DgAH_sex_peaks$chr <- gsub('ScIV947_','chr', DgAH_sex_peaks$chr)
DgAH_sex_peaks$chr <- gsub(';HRSCAF=\\d{1,2}','', DgAH_sex_peaks$chr)

# combine gene info with peak info 
DgAH_sex_peak_gene_map <- 
  left_join(DgAH_sex_peaks, Dg_transcript_starts, by = "chr") %>% 
  mutate(peakGeneDist = mid-start_gene) %>% 
  filter(abs(peakGeneDist) <= 25000) %>% # peak should be within 25kb of gene
  group_by(peak) %>% 
  top_n(-2, abs(peakGeneDist)) #selecting 2 closest genes

# combine with protein best hits
DgAH_sex_peak_gene_map_anno <- left_join(DgAH_sex_peak_gene_map, Dgaz1_prot_anno, by = "gene")
DgAH_F_peak_map_anno <- DgAH_sex_peak_gene_map_anno %>% filter(DA == "F")
DgAH_M_peak_map_anno <- DgAH_sex_peak_gene_map_anno %>% filter(DA == "M")


### Dg G sex res ############################
### denote which condition showed higher accessibility - ordered in DESeq2
DgG_sex_anno_all <- DgG_sex_res %>% 
  mutate(DA = ifelse(padj <= 0.05 & log2FoldChange < 0, "F",
                     ifelse(padj <= 0.05 & log2FoldChange > 0, "M",
                            "notDA")))

DgG_sex_res_sig_anno <- DgG_sex_res_sig %>% 
  mutate(DA = ifelse(padj <= 0.05 & log2FoldChange < 0, "F",
                     ifelse(padj <= 0.05 & log2FoldChange > 0, "M",
                            "notDA")))

DgG_sex_peaks <- DgG_sex_res_sig_anno %>% mutate(mid = (end-start)/2 + start) %>% 
  select(peak, chr, DA, start, mid, end)

# fix peaks df so chr names are matching
DgG_sex_peaks$chr <- gsub('ScIV947_','chr', DgG_sex_peaks$chr)
DgG_sex_peaks$chr <- gsub(';HRSCAF=\\d{1,2}','', DgG_sex_peaks$chr)

# combine gene info with peak info 
DgG_sex_peak_gene_map <- 
  left_join(DgG_sex_peaks, Dg_transcript_starts, by = "chr") %>% 
  mutate(peakGeneDist = mid-start_gene) %>% 
  filter(abs(peakGeneDist) <= 25000) %>% # peak should be within 25kb of gene
  group_by(peak) %>% 
  top_n(-2, abs(peakGeneDist)) #selecting 2 closest genes

# combine with protein best hits
DgG_sex_peak_gene_map_anno <- left_join(DgG_sex_peak_gene_map, Dgaz1_prot_anno, by = "gene")
DgG_F_peak_map_anno <- DgG_sex_peak_gene_map_anno %>% filter(DA == "F")
DgG_M_peak_map_anno <- DgG_sex_peak_gene_map_anno %>% filter(DA == "M")

### Dg PH sex_res ####
### denote which condition showed higher accessibility - ordered in DESeq2
DgPH_sex_anno_all <- DgPH_sex_res %>% 
  mutate(DA = ifelse(padj <= 0.05 & log2FoldChange < 0, "F",
                     ifelse(padj <= 0.05 & log2FoldChange > 0, "M",
                            "notDA")))

DgPH_sex_res_sig_anno <- DgPH_sex_res_sig %>% 
  mutate(DA = ifelse(padj <= 0.05 & log2FoldChange < 0, "F",
                     ifelse(padj <= 0.05 & log2FoldChange > 0, "M",
                            "notDA")))

DgPH_sex_peaks <- DgPH_sex_res_sig_anno %>% mutate(mid = (end-start)/2 + start) %>% 
  select(peak, chr, DA, start, mid, end)

# fix peaks df so chr names are matching
DgPH_sex_peaks$chr <- gsub('ScIV947_','chr', DgPH_sex_peaks$chr)
DgPH_sex_peaks$chr <- gsub(';HRSCAF=\\d{1,2}','', DgPH_sex_peaks$chr)

# combine gene info with peak info 
DgPH_sex_peak_gene_map <- 
  left_join(DgPH_sex_peaks, Dg_transcript_starts, by = "chr") %>% 
  mutate(peakGeneDist = mid-start_gene) %>% 
  filter(abs(peakGeneDist) <= 25000) %>% # peak should be within 25kb of gene
  group_by(peak) %>% 
  top_n(-2, abs(peakGeneDist)) #selecting 2 closest genes

# combine with protein best hits
DgPH_sex_peak_gene_map_anno <- left_join(DgPH_sex_peak_gene_map, Dgaz1_prot_anno, by = "gene")
DgPH_F_peak_map_anno <- DgPH_sex_peak_gene_map_anno %>% filter(DA == "F")
DgPH_M_peak_map_anno <- DgPH_sex_peak_gene_map_anno %>% filter(DA == "M")
### DgF_AHPH_sig1 ####
# add midpoint of peak
DgF_AHPH_peaks <- DgF_AHPH_sig %>% mutate(mid = (end-start)/2 + start) %>% 
  select(peak, chr, start, mid, end)

# have to fix so chr names are matching
# regex search and replace ;HRSCAF=\d{1,2}
write.table(DgF_AHPH_peaks, 
            file = "./GitHub/Beetle-sexual-dimorphism/DgF_AHPH_peaks.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
DgF_AHPH_peaks <- read.delim("./GitHub/Beetle-sexual-dimorphism/DgF_AHPH_peaks.txt")

# combine gene info with peak info 
DgF_AHPH_peak_gene_map <- 
  left_join(DgF_AHPH_peaks, Dg_transcript_starts, by = "chr") %>% 
  mutate(peakGeneDist = mid-start_gene) %>% 
  filter(abs(peakGeneDist) <= 25000) %>% # peak should be within 25kb of gene
  group_by(peak) %>% 
  top_n(-2, abs(peakGeneDist)) #selecting 2 closest genes

# combine with Dgaz protein best hits
DgF_AHPH_peak_gene_map_anno <- left_join(DgF_AHPH_peak_gene_map, Dgaz1_prot_anno, by = "gene")

write.table(DgF_AHPH_res_peak_gene_map_anno, 
            file = "./GitHub/Beetle-sexual-dimorphism/DgF_AHPH_peak_gene_map_anno.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

####### denote which condition showed higher accessibility 
## if large fold change was negative and condition A had higher pileups, A is negative

DgF_AHPH_anno_all <- DgF_AHPH %>% 
  mutate(DA = ifelse(padj <= 0.05 & log2FoldChange < 0, "AH",
                     ifelse(padj <= 0.05 & log2FoldChange > 0, "PH",
                            "notDA")
  )
  ) # %>% inner_join(., DgAH_sex_res_peak_gene_map_anno, by = "peak")
DgF_AHPH_anno_sig <- DgF_AHPH_anno_all %>% 

### DgM_AHPH_sig1 #### 
# add midpoint of peak
DgM_AHPH_peaks <- DgM_AHPH_sig %>% mutate(mid = (end-start)/2 + start) %>% 
  select(peak, chr, start, mid, end)

# have to fix so chr names are matching
# regex search and replace ;HRSCAF=\d{1,2}
write.table(DgM_AHPH_peaks, 
            file = "./GitHub/Beetle-sexual-dimorphism/DgM_AHPH_peaks.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
DgM_AHPH_peaks <- read.delim("./GitHub/Beetle-sexual-dimorphism/DgM_AHPH_peaks.txt")

# combine gene info with peak info 
DgM_AHPH_peak_gene_map <- 
  left_join(DgM_AHPH_peaks, Dg_transcript_starts, by = "chr") %>% 
  mutate(peakGeneDist = mid-start_gene) %>% 
  filter(abs(peakGeneDist) <= 25000) %>% # peak should be within 25kb of gene
  group_by(peak) %>% 
  top_n(-2, abs(peakGeneDist)) #selecting 2 closest genes

# combine with Dgaz protein best hits
DgM_AHPH_peak_gene_map_anno <- left_join(DgM_AHPH_peak_gene_map, Dgaz1_prot_anno, by = "gene")

write.table(DgM_AHPH_res_peak_gene_map_anno, 
            file = "./GitHub/Beetle-sexual-dimorphism/DgM_AHPH_peak_gene_map_anno.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

###### Otau ######
### OtAH_sex_res_sig1 #### 
### denote which condition showed higher accessibility - ordered in DESeq2
OtAH_sex_anno_all <- OtAH_sex_res %>% 
  mutate(DA = ifelse(padj <= 0.05 & log2FoldChange < 0, "F",
                     ifelse(padj <= 0.05 & log2FoldChange > 0, "M",
                            "notDA")))

OtAH_sex_res_sig_anno <- OtAH_sex_res_sig %>% 
  mutate(DA = ifelse(padj <= 0.05 & log2FoldChange < 0, "F",
                     ifelse(padj <= 0.05 & log2FoldChange > 0, "M",
                            "notDA")))

OtAH_sex_peaks <- OtAH_sex_res_sig_anno %>% mutate(mid = (end-start)/2 + start) %>% 
  select(peak, chr, DA, start, mid, end)

# combine gene info with peak info 
OtAH_sex_peak_gene_map <- 
  left_join(OtAH_sex_peaks, Ot_transcript_starts, by = "chr") %>% 
  mutate(peakGeneDist = mid-start_gene) %>% 
  filter(abs(peakGeneDist) <= 25000) %>% # peak should be within 25kb of gene
  group_by(peak) %>% 
  rename("OT3_ID"="gene") %>% 
  top_n(-2, abs(peakGeneDist)) #selecting 2 closest genes

# combine with protein best hits
OtAH_sex_peak_gene_map_anno <- left_join(OtAH_sex_peak_gene_map, Otau3_prot_anno, by = "OT3_ID")
OtAH_F_peak_map_anno <- OtAH_sex_peak_gene_map_anno %>% filter(DA == "F")
OtAH_M_peak_map_anno <- OtAH_sex_peak_gene_map_anno %>% filter(DA == "M")
### OtG_sex_res ######
### denote which condition showed higher accessibility - ordered in DESeq2
OtG_sex_anno_all <- OtG_sex_res %>% 
  mutate(DA = ifelse(padj <= 0.05 & log2FoldChange < 0, "F",
                     ifelse(padj <= 0.05 & log2FoldChange > 0, "M",
                            "notDA")))

OtG_sex_res_sig_anno <- OtG_sex_res_sig %>% 
  mutate(DA = ifelse(padj <= 0.05 & log2FoldChange < 0, "F",
                     ifelse(padj <= 0.05 & log2FoldChange > 0, "M",
                            "notDA")))

OtG_sex_peaks <- OtG_sex_res_sig_anno %>% mutate(mid = (end-start)/2 + start) %>% 
  select(peak, chr, DA, start, mid, end)

# combine gene info with peak info 
OtG_sex_peak_gene_map <- 
  left_join(OtG_sex_peaks, Ot_transcript_starts, by = "chr") %>% 
  mutate(peakGeneDist = mid-start_gene) %>% 
  filter(abs(peakGeneDist) <= 25000) %>% # peak should be within 25kb of gene
  group_by(peak) %>% 
  rename("OT3_ID"="gene") %>% 
  top_n(-2, abs(peakGeneDist)) #selecting 2 closest genes

# combine with protein best hits
OtG_sex_peak_gene_map_anno <- left_join(OtG_sex_peak_gene_map, Otau3_prot_anno, by = "OT3_ID")
OtG_F_peak_map_anno <- OtG_sex_peak_gene_map_anno %>% filter(DA == "F")
OtG_M_peak_map_anno <- OtG_sex_peak_gene_map_anno %>% filter(DA == "M")

### OtPH_sex_res_sig #####
### denote which condition showed higher accessibility - ordered in DESeq2
OtPH_sex_anno_all <- OtPH_sex_res %>% 
  mutate(DA = ifelse(padj <= 0.05 & log2FoldChange < 0, "F",
                     ifelse(padj <= 0.05 & log2FoldChange > 0, "M",
                            "notDA")))

OtPH_sex_res_sig_anno <- OtPH_sex_res_sig %>% 
  mutate(DA = ifelse(padj <= 0.05 & log2FoldChange < 0, "F",
                     ifelse(padj <= 0.05 & log2FoldChange > 0, "M",
                            "notDA")))

OtPH_sex_peaks <- OtPH_sex_res_sig_anno %>% mutate(mid = (end-start)/2 + start) %>% 
  select(peak, chr, DA, start, mid, end)

# combine gene info with peak info 
OtPH_sex_peak_gene_map <- 
  left_join(OtPH_sex_peaks, Ot_transcript_starts, by = "chr") %>% 
  mutate(peakGeneDist = mid-start_gene) %>% 
  filter(abs(peakGeneDist) <= 25000) %>% # peak should be within 25kb of gene
  group_by(peak) %>% 
  rename("OT3_ID"="gene") %>% 
  top_n(-2, abs(peakGeneDist)) #selecting 2 closest genes

# combine with protein best hits
OtPH_sex_peak_gene_map_anno <- left_join(OtPH_sex_peak_gene_map, Otau3_prot_anno, by = "OT3_ID")
OtPH_F_peak_map_anno <- OtPH_sex_peak_gene_map_anno %>% filter(DA == "F")
OtPH_M_peak_map_anno <- OtPH_sex_peak_gene_map_anno %>% filter(DA == "M")

### OtF_AHPH_sig1 #### 
# add midpoint of peak
OtF_AHPH_peaks <- OtF_AHPH_sig %>% mutate(mid = (end-start)/2 + start) %>% 
  select(peak, chr, start, mid, end)

# combine gene info with peak info 
OtF_AHPH_peak_gene_map <- 
  left_join(OtF_AHPH_peaks, Ot_transcript_starts, by = "chr") %>% 
  mutate(peakGeneDist = mid-start_gene) %>% 
  filter(abs(peakGeneDist) <= 25000) %>% # peak should be within 25kb of gene
  group_by(peak) %>% 
  rename("OT3_ID"="gene") %>% 
  top_n(-2, abs(peakGeneDist)) #selecting 2 closest genes

# combine with protein best hits
OtF_AHPH_peak_gene_map_anno <- left_join(OtF_AHPH_peak_gene_map, Otau3_prot_anno, by = "gene")

write.table(OtF_AHPH_peak_gene_map_anno, 
            file = "./GitHub/Beetle-sexual-dimorphism/OtF_AHPH_peak_gene_map_anno.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

## denote which condition showed higher accessibility - ordered in DESeq2
OtF_AHPH_anno_all <- OtM_AHPH %>% 
  mutate(DA = ifelse(padj <= 0.05 & log2FoldChange < 0, "AH",
                     ifelse(padj <= 0.05 & log2FoldChange > 0, "PH",
                            "notDA")
  )
  ) # %>% inner_join(., DgAH_sex_res_peak_gene_map_anno, by = "peak")
OtF_AHPH_anno_sig <- OtF_AHPH_anno_all %>% 
  filter(padj <= 0.05)
### OtM_AHPH_sig1 #####
# add midpoint of peak
OtM_AHPH_peaks <- OtM_AHPH_sig %>% mutate(mid = (end-start)/2 + start) %>% 
  select(peak, chr, start, mid, end)

# combine gene info with peak info 
OtM_AHPH_peak_gene_map <- 
  left_join(OtM_AHPH_peaks, Ot_transcript_starts, by = "chr") %>% 
  mutate(peakGeneDist = mid-start_gene) %>% 
  filter(abs(peakGeneDist) <= 25000) %>% # peak should be within 25kb of gene
  group_by(peak) %>% 
  rename("OT3_ID"="gene") %>% 
  top_n(-2, abs(peakGeneDist)) #selecting 2 closest genes

# combine with protein best hits
OtM_AHPH_peak_gene_map_anno <- left_join(OtM_AHPH_peak_gene_map, Otau3_prot_anno, by = "gene")


write.table(OtM_AHPH_peak_gene_map_anno, 
            file = "./GitHub/Beetle-sexual-dimorphism/OtM_AHPH_peak_gene_map_anno.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)



### denote which condition showed higher accessibility - ordered in DESeq2
OtM_AHPH_anno_all <- OtM_AHPH %>% 
  mutate(DA = ifelse(padj <= 0.05 & log2FoldChange < 0, "AH",
                     ifelse(padj <= 0.05 & log2FoldChange > 0, "PH",
                            "notDA")
  )
  ) # %>% inner_join(., DgAH_sex_res_peak_gene_map_anno, by = "peak")
OtM_AHPH_anno_sig <- OtM_AHPH_anno_all %>% 
  filter(padj <= 0.05)

## get separate peak bed files 
OtM_AHPH_anno_sig_PH_up <- OtM_AHPH_anno_sig %>%
  filter(DA == "PH") %>% 
  mutate(mid = (end-start)/2 + start)

OtM_AHPH_anno_sig_AH_up <- OtM_AHPH_anno_sig %>%
  filter(DA == "AH") %>% 
  mutate(mid = (end-start)/2 + start)

write.table(OtM_AHPH_anno_sig_PH_up, 
            file = "./GitHub/Beetle-sexual-dimorphism/OtM_AHPH_peaks_PH_up.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(OtM_AHPH_anno_sig_AH_up, 
            file = "./GitHub/Beetle-sexual-dimorphism/OtM_AHPH_peaks_AH_up.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

###### Step 9: export results ######
### export separate male open and female open OCR list for motif enrichment analysis ########
### Dgaz 
Dg_PH_F_up <- Dg_PH_MvF_OCR %>% filter(DA == "F")
write.table(Dg_PH_F_up, file = "/Users/ericanadolski/Documents/Sexual_dimorphism_project/homer/onthophagus-paper/DgPH_peaks_F_up.txt", sep = "\t", quote = FALSE, row.names = FALSE)

Dg_PH_M_up <- Dg_PH_MvF_OCR %>% filter(DA == "M")
write.table(Dg_PH_M_up, file = "/Users/ericanadolski/Documents/Sexual_dimorphism_project/homer/onthophagus-paper/DgPH_peaks_M_up.txt", sep = "\t", quote = FALSE, row.names = FALSE)

Dg_AH_F_up <- Dg_AH_MvF_OCR %>% filter(DA == "F")
write.table(Dg_AH_F_up, file = "/Users/ericanadolski/Documents/Sexual_dimorphism_project/homer/onthophagus-paper/DgAH_peaks_F_up.txt", sep = "\t", quote = FALSE, row.names = FALSE)

Dg_AH_M_up <- Dg_AH_MvF_OCR %>% filter(DA == "M")
write.table(Dg_AH_M_up, file = "/Users/ericanadolski/Documents/Sexual_dimorphism_project/homer/onthophagus-paper/DgAH_peaks_M_up.txt", sep = "\t", quote = FALSE, row.names = FALSE)

Dg_G_F_up <- Dg_G_MvF_OCR %>% filter(DA == "F")
write.table(Dg_G_F_up, file = "/Users/ericanadolski/Documents/Sexual_dimorphism_project/homer/onthophagus-paper/DgG_peaks_F_up.txt", sep = "\t", quote = FALSE, row.names = FALSE)

Dg_G_M_up <- Dg_G_MvF_OCR %>% filter(DA == "M")
write.table(Dg_G_M_up, file = "/Users/ericanadolski/Documents/Sexual_dimorphism_project/homer/onthophagus-paper/DgG_peaks_M_up.txt", sep = "\t", quote = FALSE, row.names = FALSE)

### Otau
Ot_PH_F_up <- Ot_PH_MvF_OCR %>% filter(DA == "F")
write.table(Ot_PH_F_up, file = "/Users/ericanadolski/Documents/Sexual_dimorphism_project/homer/onthophagus-paper/OtPH_peaks_F_up.txt", sep = "\t", quote = FALSE, row.names = FALSE)

Ot_PH_M_up <- Ot_PH_MvF_OCR %>% filter(DA == "M")
write.table(Ot_PH_M_up, file = "/Users/ericanadolski/Documents/Sexual_dimorphism_project/homer/onthophagus-paper/OtPH_peaks_M_up.txt", sep = "\t", quote = FALSE, row.names = FALSE)

Ot_AH_F_up <- Ot_AH_MvF_OCR %>% filter(DA == "F")
write.table(Ot_AH_F_up, file = "/Users/ericanadolski/Documents/Sexual_dimorphism_project/homer/onthophagus-paper/OtAH_peaks_F_up.txt", sep = "\t", quote = FALSE, row.names = FALSE)

Ot_AH_M_up <- Ot_AH_MvF_OCR %>% filter(DA == "M")
write.table(Ot_AH_M_up, file = "/Users/ericanadolski/Documents/Sexual_dimorphism_project/homer/onthophagus-paper/OtAH_peaks_M_up.txt", sep = "\t", quote = FALSE, row.names = FALSE)

Ot_G_F_up <- Ot_G_MvF_OCR %>% filter(DA == "F")
write.table(Ot_G_F_up, file = "/Users/ericanadolski/Documents/Sexual_dimorphism_project/homer/onthophagus-paper/OtG_peaks_F_up.txt", sep = "\t", quote = FALSE, row.names = FALSE)

Ot_G_M_up <- Ot_G_MvF_OCR %>% filter(DA == "M")
write.table(Ot_G_M_up, file = "/Users/ericanadolski/Documents/Sexual_dimorphism_project/homer/onthophagus-paper/OtG_peaks_M_up.txt", sep = "\t", quote = FALSE, row.names = FALSE)

### Osag 
Os_PH_F_up <- Os_PH_MvF_OCR %>% filter(DA == "F")
write.table(Os_PH_F_up, file = "/Users/ericanadolski/Documents/Sexual_dimorphism_project/homer/onthophagus-paper/OsPH_peaks_F_up.txt", sep = "\t", quote = FALSE, row.names = FALSE)

Os_PH_M_up <- Os_PH_MvF_OCR %>% filter(DA == "M")
write.table(Os_PH_M_up, file = "/Users/ericanadolski/Documents/Sexual_dimorphism_project/homer/onthophagus-paper/OsPH_peaks_M_up.txt", sep = "\t", quote = FALSE, row.names = FALSE)

Os_AH_F_up <- Os_AH_MvF_OCR %>% filter(DA == "F")
write.table(Os_AH_F_up, file = "/Users/ericanadolski/Documents/Sexual_dimorphism_project/homer/onthophagus-paper/OsAH_peaks_F_up.txt", sep = "\t", quote = FALSE, row.names = FALSE)

Os_AH_M_up <- Os_AH_MvF_OCR %>% filter(DA == "M")
write.table(Os_AH_M_up, file = "/Users/ericanadolski/Documents/Sexual_dimorphism_project/homer/onthophagus-paper/OsAH_peaks_M_up.txt", sep = "\t", quote = FALSE, row.names = FALSE)

Os_G_F_up <- Os_G_MvF_OCR %>% filter(DA == "F")
write.table(Os_G_F_up, file = "/Users/ericanadolski/Documents/Sexual_dimorphism_project/homer/onthophagus-paper/OsG_peaks_F_up.txt", sep = "\t", quote = FALSE, row.names = FALSE)

Os_G_M_up <- Os_G_MvF_OCR %>% filter(DA == "M")
write.table(Os_G_M_up, file = "/Users/ericanadolski/Documents/Sexual_dimorphism_project/homer/onthophagus-paper/OsG_peaks_M_up.txt", sep = "\t", quote = FALSE, row.names = FALSE)

### export annotated significant peak tables ######
#Otau
write.table(OtAH_sex_peak_gene_map_anno, 
            file = "./GitHub/Beetle-sexual-dimorphism/DAP-tables-annotated/OtAH_sex_peak_gene_map_anno.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(OtAH_F_peak_map_anno, 
            file = "./GitHub/Beetle-sexual-dimorphism/DAP-tables-annotated/OtAH_F_peak_map_anno.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(OtAH_M_peak_map_anno, 
            file = "./GitHub/Beetle-sexual-dimorphism/DAP-tables-annotated/OtAH_M_peak_map_anno.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(OtG_sex_peak_gene_map_anno, 
            file = "./GitHub/Beetle-sexual-dimorphism/DAP-tables-annotated/OtG_sex_peak_gene_map_anno.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(OtG_F_peak_map_anno, 
            file = "./GitHub/Beetle-sexual-dimorphism/DAP-tables-annotated/OtG_F_peak_map_anno.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(OtG_M_peak_map_anno, 
            file = "./GitHub/Beetle-sexual-dimorphism/DAP-tables-annotated/OtG_M_peak_map_anno.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(OtPH_sex_peak_gene_map_anno, 
            file = "./GitHub/Beetle-sexual-dimorphism/DAP-tables-annotated/OtPH_sex_peak_gene_map_anno.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(OtPH_F_peak_map_anno, 
            file = "./GitHub/Beetle-sexual-dimorphism/DAP-tables-annotated/OtPH_F_peak_map_anno.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(OtPH_M_peak_map_anno, 
            file = "./GitHub/Beetle-sexual-dimorphism/DAP-tables-annotated/OtPH_M_peak_map_anno.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
#Osag
write.table(OsAH_sex_peak_gene_map_anno, 
            file = "./GitHub/Beetle-sexual-dimorphism/DAP-tables-annotated/OsAH_sex_peak_gene_map_anno.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(OsAH_F_peak_map_anno, 
            file = "./GitHub/Beetle-sexual-dimorphism/DAP-tables-annotated/OsAH_F_peak_map_anno.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(OsAH_M_peak_map_anno, 
            file = "./GitHub/Beetle-sexual-dimorphism/DAP-tables-annotated/OsAH_M_peak_map_anno.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(OsG_sex_peak_gene_map_anno, 
            file = "./GitHub/Beetle-sexual-dimorphism/DAP-tables-annotated/OsG_sex_peak_gene_map_anno.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(OsG_F_peak_map_anno, 
            file = "./GitHub/Beetle-sexual-dimorphism/DAP-tables-annotated/OsG_F_peak_map_anno.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(OsG_M_peak_map_anno, 
            file = "./GitHub/Beetle-sexual-dimorphism/DAP-tables-annotated/OsG_M_peak_map_anno.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
rite.table(OsPH_sex_peak_gene_map_anno, 
            file = "./GitHub/Beetle-sexual-dimorphism/DAP-tables-annotated/OsPH_sex_peak_gene_map_anno.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(OsPH_F_peak_map_anno, 
            file = "./GitHub/Beetle-sexual-dimorphism/DAP-tables-annotated/OsPH_F_peak_map_anno.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(OsPH_M_peak_map_anno, 
            file = "./GitHub/Beetle-sexual-dimorphism/DAP-tables-annotated/OsPH_M_peak_map_anno.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
#Dgaz
write.table(DgAH_sex_peak_gene_map_anno, 
            file = "./GitHub/Beetle-sexual-dimorphism/DAP-tables-annotated/DgAH_sex_peak_gene_map_anno.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(DgAH_F_peak_map_anno, 
            file = "./GitHub/Beetle-sexual-dimorphism/DAP-tables-annotated/DgAH_F_peak_map_anno.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(DgAH_M_peak_map_anno, 
            file = "./GitHub/Beetle-sexual-dimorphism/DAP-tables-annotated/DgAH_M_peak_map_anno.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(DgG_sex_peak_gene_map_anno, 
            file = "./GitHub/Beetle-sexual-dimorphism/DAP-tables-annotated/DgG_sex_peak_gene_map_anno.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(DgG_F_peak_map_anno, 
            file = "./GitHub/Beetle-sexual-dimorphism/DAP-tables-annotated/DgG_F_peak_map_anno.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(DgG_M_peak_map_anno, 
            file = "./GitHub/Beetle-sexual-dimorphism/DAP-tables-annotated/DgG_M_peak_map_anno.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(DgPH_sex_peak_gene_map_anno, 
            file = "./GitHub/Beetle-sexual-dimorphism/DAP-tables-annotated/DgPH_sex_peak_gene_map_anno.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(DgPH_F_peak_map_anno, 
            file = "./GitHub/Beetle-sexual-dimorphism/DAP-tables-annotated/DgPH_F_peak_map_anno.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(DgPH_M_peak_map_anno, 
            file = "./GitHub/Beetle-sexual-dimorphism/DAP-tables-annotated/DgPH_M_peak_map_anno.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
### tobias motifs ######
homer <- universalmotif::read_homer("/Users/ericanadolski/Documents/Sexual_dimorphism_project/homer/OtPH-all-motifs/out_test/homerResults/nonRedundant2.motif")

homer <- universalmotif::read_homer("/Users/ericanadolski/Documents/Sexual_dimorphism_project/homer/OtPH-all-motifs/homerMotifs.all.motifs")

universalmotif::write_jaspar(homer,"/Users/ericanadolski/Documents/Sexual_dimorphism_project/tobias/jasparOtPH2.motifs", overwrite = TRUE)