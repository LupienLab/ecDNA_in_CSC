rm(list = ls())

library(scatools)
library(tidyverse)
library(patchwork)
library(dittoSeq)
library(purrr)
library(ggplotify)
library(Seurat)
library(ArchR)
library(data.table)
library(zellkonverter)
library(BSgenome.Hsapiens.UCSC.hg38)
'%ni%' <- Negate('%in%')

id = "G958_T"

#~~~ input the data
fragment <- paste0("/Users/chupan/Documents/data_su2c/scATAC/fragment/", id, ".fragments.tsv.gz")
barcode <- fread(paste0("/Users/chupan/Documents/ArchR/outputs.", id, "/", id, ".filtercell.txt"), header = F, data.table = F) %>% {colnames(.) = "barcode";.} %>% pull(barcode)
metacluster <- fread(paste0("/Users/chupan/Documents/ArchR/outputs.", id, "/", id, ".metacluster.txt"), header = F, data.table = F) %>% {colnames(.) = c("barcode", "cluster") ;.}

#~~~ generate bins object
bins <- get_tiled_bins(
  bs_genome = BSgenome.Hsapiens.UCSC.hg38,
  tilewidth = 1e7,
  respect_chr_arms = TRUE
)

#~~~ get blacklist regions
blacklist <- get_blacklist(genome = "hg38")

#~~~ run scatools to identify CNA at single-cell resolution 
scCNA <- run_scatools(
  sample_id = id,
  fragment_file = fragment,
  cells = barcode,
  bins = bins,
  blacklist = blacklist,
  min_cell_counts = 0,
  min_cell_prop = 0,
  outdir = paste0("outputs.", id),
  overwrite = TRUE,
  verbose = TRUE,
  ncores = 16,
  save_h5ad = F
)

#~~~ rearrange the clusters based on the metacluster 
metacluster <- metacluster %>% dplyr::arrange(scCNA$barcode)
scCNA$clusters <- metacluster$cluster

#~~~ plot the results usign the dittoSeq package
col_clones <- function(clones) {
  c(dittoColors()[1:length(unique(clones[clones != "N"]))], "black")
}

logr_min <- -0.5
logr_max <- 0.5
logr_breaks <- seq(logr_min, logr_max, length.out = 3)
logr_colors <- colorRampPalette(c("blue", "white", "red"))(3)

logr <- assay(scCNA, "logr_modal")
logr[!is.finite(logr)] <- NA
logr[is.na(logr)] <- 0
assay(scCNA, "logr_modal") <- logr

ht1 <- scatools::cloneCnaHeatmap(scCNA, assay_name = "logr_modal", clone_name = "clusters",
                                 cluster_clones = TRUE,
                                 col_clones = col_clones(scCNA$clusters), legend_name = "logr", 
                                 col_fun = logr_col_fun(breaks = logr_breaks, colors = logr_colors)) %>%
  as.ggplot()

ht1

saveRDS(scCNA, file = paste0("./outputs.", id, "/", id, ".scCNA.rds"))

