###############################################################
# Transcriptomic Stability During Long-Term Cell Culture
# Differential Expression Analysis of RNA-seq Data
#
# Dataset: GSE268437
# Author: Mahdis Hamidi
#
# This script reproduces the downstream differential
# expression analysis performed on the processed count matrix.
###############################################################

###############################
# 1. Load packages
###############################

library(readxl)
library(DESeq2)
library(EnhancedVolcano)
library(pheatmap)

###############################
# 2. Create output directories
###############################

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

###############################
# 3. Import processed count matrix
###############################

counts <- read_excel(
  "data/processed/GSE268437_Processed_Combated_countmatrix_240521 (1).xlsx"
)

counts <- as.data.frame(counts)

rownames(counts) <- counts$Sample

counts$Sample <- NULL

counts <- counts[rownames(counts) != "Identifier", ]

###############################
# 4. Select samples used in the analysis
###############################

time_counts <- counts[, c(
  
  "A673.ITP","A673.MTP","A673.FTP",
  
  "MHH.ES1.ITP","MHH.ES1.MTP","MHH.ES1.FTP",
  
  "SKES.ITP","SKES.MTP","SKES.FTP",
  
  "SKNMC.ITP","SKNMC.MTP","SKNMC.FTP",
  
  "TC71.ITP","TC71.MTP","TC71.FTP",
  
  "MCF7.0.ITP...15","MCF7.0.MTP...16","MCF7.0.FTP",
  
  "Hela.0.ITP...13","Hela.0.MTP...14","Hela.0.FTP"
  
)]

###############################
# 5. Sample metadata
###############################

sample_info <- data.frame(
  
  sample = colnames(time_counts),
  
  cellline = factor(c(
    
    rep("A673",3),
    
    rep("MHHES1",3),
    
    rep("SKES",3),
    
    rep("SKNMC",3),
    
    rep("TC71",3),
    
    rep("MCF7",3),
    
    rep("HeLa",3)
    
  )),
  
  time = factor(
    
    rep(c("m0","m6","m12"),7),
    
    levels = c("m0","m6","m12")
    
  )
  
)

rownames(sample_info) <- sample_info$sample

###############################
# 6. Differential expression
###############################

dds <- DESeqDataSetFromMatrix(
  
  countData = round(time_counts),
  
  colData = sample_info,
  
  design = ~ cellline + time
  
)

dds <- DESeq(dds)

###############################
# 7. Extract results
###############################

res_m6 <- results(dds, name = "time_m6_vs_m0")

res_m12 <- results(dds, name = "time_m12_vs_m0")

res_m12 <- res_m12[order(res_m12$padj), ]

write.csv(
  
  as.data.frame(res_m12),
  
  file = "results/DEGs_m12_vs_m0.csv"
  
)

###############################
# 8. Variance stabilizing transformation
###############################

vsd <- vst(dds, blind = FALSE)

###############################
# 9. PCA by cell line
###############################

png(
  
  "figures/pca_cellline.png",
  
  width = 2200,
  
  height = 1800,
  
  res = 300
  
)

plotPCA(
  
  vsd,
  
  intgroup = "cellline"
  
)

dev.off()

###############################
# 10. PCA by time point
###############################

png(
  
  "figures/pca_time.png",
  
  width = 2200,
  
  height = 1800,
  
  res = 300
  
)

plotPCA(
  
  vsd,
  
  intgroup = "time"
  
)

dev.off()

###############################
# 11. Volcano plot
###############################

png(
  
  "figures/volcano_m12_vs_m0.png",
  
  width = 2400,
  
  height = 1800,
  
  res = 300
  
)

EnhancedVolcano(
  
  res_m12,
  
  lab = rownames(res_m12),
  
  x = "log2FoldChange",
  
  y = "padj",
  
  title = "Differential Expression: m12 vs m0",
  
  pCutoff = 0.05,
  
  FCcutoff = 1
  
)

dev.off()

###############################
# 12. Heatmap of top DE genes
###############################

top_genes <- head(
  
  rownames(res_m12),
  
  50
  
)

heatmap_matrix <- assay(vsd)[top_genes, ]

heatmap_matrix <- heatmap_matrix -
  
  rowMeans(heatmap_matrix)

png(
  
  "figures/heatmap_top50_genes.png",
  
  width = 2200,
  
  height = 2200,
  
  res = 300
  
)

pheatmap(
  
  heatmap_matrix,
  
  scale = "none",
  
  cluster_rows = TRUE,
  
  cluster_cols = TRUE,
  
  show_rownames = TRUE,
  
  show_colnames = TRUE,
  
  fontsize_row = 8,
  
  fontsize_col = 10
  
)

dev.off()

###############################
# 13. Session information
###############################

writeLines(
  
  capture.output(sessionInfo()),
  
  "results/sessionInfo.txt"
  
)

###############################################################
# Analysis completed successfully
###############################################################