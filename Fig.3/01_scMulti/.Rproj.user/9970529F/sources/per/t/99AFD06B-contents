rm(list = ls())

library(Seurat) #~~~ 5.2.1
library(Signac)
library(dplyr)
library(Matrix)
library(ggplot2)
library(patchwork)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

#~~~ step1: load all Seurat objects into a named list
rds_files <- list.files("./object/", pattern = "\\.rds$", recursive = FALSE, full.names = TRUE)

samples <- list()
for(file_path in rds_files){
  seurat_object <- readRDS(file_path)
  file_name <- tools::file_path_sans_ext(basename(file_path))
  samples[[file_name]] <- seurat_object
}

#~~~ step2: merge all Seurat objects using Reduce function
for (i in seq_along(samples)) {
  DefaultAssay(samples[[i]]) <- "RNA"
  sample_id <- names(samples)[i]
  samples[[i]] <- RenameCells(samples[[i]], new.names = paste0(sample_id, "|", Cells(samples[[i]])))
}

#~~~ step3: merge seurat objects with merge.data = FALSE
merged_seurat <- Reduce(function(x, y) merge(x, y, merge.data = FALSE), samples)

#~~~ step4: explicitly add counts layer (sum across objects)
counts_list <- lapply(samples, function(s) GetAssayData(s, slot = "counts"))
merged_counts <- do.call(cbind, counts_list)
merged_seurat[["RNA"]]["counts"] <- merged_counts

#~~~ confirm result
counts_matrix <- GetAssayData(merged_seurat, layer = "counts")

saveRDS(merged_seurat, file = "./object_merged/merged_seurat.rds")

#~~~ checks after merging
merged_seurat <- readRDS("./object_merged/merged_seurat.rds")
Assays(merged_seurat)
table(merged_seurat$orig.ident)

#~~~ extract PCA embeddings (coordinates of each cell in PC space)
DefaultAssay(merged_seurat) <- "RNA"
merged_seurat <- NormalizeData(merged_seurat)
merged_seurat <- FindVariableFeatures(merged_seurat)
merged_seurat <- ScaleData(merged_seurat)
merged_seurat <- RunPCA(merged_seurat) %>% RunUMAP(dims = 1:30, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
pca_df <- as.data.frame(Embeddings(merged_seurat, "pca"))
pca_df$sample <- merged_seurat$id

pca_var <- Stdev(merged_seurat, "pca")^2
pca_var_exp <- round(pca_var / sum(pca_var) * 100, 1)  # percentage variance explained

pdf(paste0("./object_merged/PCA_integrate_RNA.pdf"), width = 3.5, height = 3.5)

p2 <- ggplot(pca_df, aes(x = PC_1, y = PC_2, color = sample)) + theme_bw() +
  geom_point(size = 1, alpha = 1) +
  labs(x = paste0("PC1:", pca_var_exp[1], "% variance"),
       y = paste0("PC2:", pca_var_exp[2], "% variance"),
       title = NULL) +
  scale_color_manual(values = c("G523_L" = "#1f77b4", "G620_L" = "#2ca02c", "G828_L" = "#9467bd",
                                "G583_L" = "#ff7f0e", "G797_L" = "#d62728", "G837_L" = "#8c564b")) +
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

p2

dev.off()
