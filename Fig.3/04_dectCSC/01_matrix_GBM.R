#~~~~ step1: generate peak-by-cell matrix from the GBM side ~~~#

rm(list = ls())

library(ggplot2)
library(dplyr)
library(Signac)
library(GenomicRanges)
library(rtracklayer)
library(data.table)
library(Seurat)

dir.create(paste0("outputs"))

#~~~ fragment
fragFileDir <- list.files(paste0("/Users/chupan/Documents/data_su2c/scATAC/fragment/"), full.names = TRUE)
selectedfragFileDir <- grep("fragments\\.tsv\\.gz$", fragFileDir, value = T)

fragment_list <- lapply(selectedfragFileDir, function(path) {
  CreateFragmentObject(path = path)
})

#~~~ union peak across the eight GBM tumors and six GSC cultures
unifiedPeak <- fread(paste0("/Users/chupan/Documents/data_su2c/scMulti/unifiedPeak/malign.unifiedPeakSet.bed"), header = T, data.table = F, stringsAsFactors = F) %>% dplyr::select(seqnames, start, end) %>% GRanges()

#~~~ cell barcode 
cellFileDir <- list.files(paste0("/Users/chupan/Documents/data_su2c/scATAC/maligncell/"), full.names = TRUE)
selectedcellFileDir <- grep("maligncell\\.txt", cellFileDir, value = T)

cell_ids <- sub(".maligncell\\.txt", "", basename(selectedcellFileDir))
cell_list <- lapply(selectedcellFileDir, function(path) {
  fread(path, header = F, stringsAsFactors = F) %>% as.data.frame() %>% {colnames(.) = "barcode";.}
})
names(cell_list) <- cell_ids

#~~~ generate the peak_cell_matrix
peak_cell_list <- lapply(seq_along(fragment_list), function(i) {
  
  frag_id <- sub("\\.fragments\\.tsv\\.gz$", "", basename(selectedfragFileDir[i]))
  cell_id <- sub("\\.maligncell\\.txt$", "", basename(selectedcellFileDir[i]))
 
  if(frag_id == cell_id){
    featureMatrix <- FeatureMatrix(
      fragments = fragment_list[[i]],                        # the i-th fragment file
      features = unifiedPeak,                                # common feature set
      cells = cell_list[[i]]$barcode)                        # the i-th cell list
    colnames(featureMatrix) <- paste0(cell_id, "#", colnames(featureMatrix))}
  return(featureMatrix)
})

peak_cell_matrix <- do.call(cbind, peak_cell_list)

dir.create("outputs")
saveRDS(peak_cell_matrix, file = paste0("outputs/peak_by_cell_GBM.rds"))


