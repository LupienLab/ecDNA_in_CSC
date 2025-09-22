#******** step 1: gene activity of CSC ********#

rm(list = ls())

library(ArchR)
library(data.table)
library(dplyr)
library(tidyverse)
library(caret)
library(cowplot)
library(readxl)
library(ggplot2)
library(ggnewscale)
library(gridExtra)
library(tidyr)
library(BSgenome.Hsapiens.UCSC.hg38)
set.seed(2025)

addArchRGenome("hg38")
addArchRThreads(threads = floor(parallel::detectCores())/2, force = F)

obj <- readRDS(paste0("./outputs.gro/mergedObject.rds"))
obj <- addImputeWeights(obj) # use MAGIC to impute gene scores by smoothing signal across nearby cells

matrixs <- getMatrixFromProject(obj)
genescoreMatrix <- matrixs@assays@data$GeneScoreMatrix 
geneName <- rowData(matrixs)$name

#~~~ load CSC information
SVD <- readRDS("/Users/chupan/Documents/data_su2c/scATAC/CSC/integrate_SVD.rds") %>% dplyr::filter(type == "T")
CSC <- readRDS("/Users/chupan/Documents/data_su2c/scATAC/CSC/overlap_cells.rds") %>% dplyr::filter(type == "T") 

selected_genes <- c("SOX2", "NANOG", "POU5F1", "PROM1", "NES", "KLF4",     # CSC marker 
                    "MKI67", "RRM2", "PCNA", "TOP2A", "MCM2", "MCM6", "AURKA", "AURKB", "CCNB1", "CDK1", "BIRC5",
                    "GFAP", "S100B", "MBP", "MOG", "TUBB3")
gene_idx <- which(geneName %in% selected_genes)
sub_genescore <- genescoreMatrix[gene_idx, ]
rownames(sub_genescore) <- geneName[gene_idx]
colnames(sub_genescore) <- colnames(genescoreMatrix)
sub_df <- t(as.data.frame(as.matrix(sub_genescore))) %>% as.data.frame() %>% tibble::rownames_to_column(var = "barcode")

#~~~ malignant cells
SVD_df <- SVD %>% inner_join(sub_df, by = "barcode") 

#~~~ CSC-positive cells
CSC_positive <- SVD_df %>% inner_join(CSC, by = "barcode")
CSC_negative <- SVD_df %>% anti_join(CSC, by = "barcode")

#~~~~ two-sample permutation test
gene <- "MCM2"

x_pos <-CSC_positive %>% select(paste0(gene)) %>% pull()   # CSC-positive cells
length(x_pos)
x_neg <- CSC_negative %>% select(paste0(gene)) %>% pull()  # CSC-negative cells
length(x_neg)

obs_diff <- mean(x_pos) - mean(x_neg)
obs_diff

all_values <- c(x_pos, x_neg)
labels <- c(rep(1, length(x_pos)), rep(0, length(x_neg)))

R <- 10000
perm_diffs <- replicate(R, {
  perm_labels <- sample(labels)
  mean(all_values[perm_labels == 1]) - mean(all_values[perm_labels == 0])
})

pval <- mean(perm_diffs >= obs_diff)
pval

b <- sum(perm_diffs >= obs_diff)      # exceedances for right tail
pval_corrected <- (b + 1) / (R + 1)   # bias-corrected one-sided p

df_plot <- bind_rows(
  data.frame(value = perm_diffs, type = "perm_diff"),
  data.frame(value = obs_diff, type = "obs_diff")
)

# plot
dir.create("./outputs.CSC")
pdf(paste0("./outputs.CSC/", gene, "_perm_test.pdf"), width = 3, height = 3)

ggplot(df_plot, aes(x = "permutation", y = value)) +
  geom_boxplot(fill = "orange", outlier.shape = NA, width = 0.4) +
  geom_hline(yintercept = obs_diff, color = "blue", linetype = "dashed", linewidth = 0.25) +
  annotate("text", x = 0.5, y = obs_diff, 
           label = paste0("P = ", signif(pval_corrected, 10), 
                          "\nObserved = ", round(obs_diff, 3)),  
           hjust = 0, color = "blue", size = 4.5) +
  labs(title = paste0(gene),
       y = "Difference in means", x = NULL) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black"))

dev.off()


