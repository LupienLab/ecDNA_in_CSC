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

#~~~ save data sets for velocity analysis
dir.create("outputs.velocity" , showWarnings = FALSE)
merged_seurat_obj <- readRDS("./object_merged/merged_seurat.rds")
tmp_obj <- merged_seurat_obj
tmp_obj$barcode <- colnames(tmp_obj)

SVD_df <- readRDS(file = paste0("/Users/chupan/Documents/data_su2c/scATAC/CSC/integrate_SVD.rds")) %>% dplyr::filter(type == "L") %>% mutate(barcode = gsub("#", "|", barcode))
SVD_df_ranked <- SVD_df[match(colnames(tmp_obj), SVD_df$barcode), ]

tmp_obj$SV_1 <- SVD_df_ranked$SV2
tmp_obj$SV_2 <- SVD_df_ranked$SV3

write.csv(tmp_obj@meta.data, file = "./outputs.velocity/seurat_obj_metadata.csv", quote = F, row.names = F)

DefaultAssay(tmp_obj) <- "RNA"
counts_matrix <- GetAssayData(tmp_obj, assay='RNA', layer = 'counts')
writeMM(counts_matrix, file =  "./outputs.velocity/seurat_obj_counts.mtx")

# write gene names
write.table(data.frame('gene' = rownames(counts_matrix)), file = "./outputs.velocity/seurat_obj_genenames.csv", quote = F, row.names = F, col.names = F)


