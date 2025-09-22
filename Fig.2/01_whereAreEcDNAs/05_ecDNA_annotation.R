rm(list = ls())

library(dplyr)
library(tibble)
library(data.table)

sample <- "A7TK_T"
nth <- 1

obj <- readRDS(file = paste0(sample, ".outputs", "/", sample, "-", nth, ".obj.rds"))

cell_df <- attr(obj, "cellColData")
cell_df <- as.data.frame(cell_df) %>% mutate(species = paste0(sample, "-", nth)) %>% rownames_to_column(var = "barcode") %>%
  dplyr::select(barcode, species, ecDNAcontent, celltype)

malign_cell <- read.table(paste0("/Users/chupan/Documents/data_su2c/scATAC/maligncell/", sample, ".maligncell.txt"), header = F, stringsAsFactors = F) %>% 
  {colnames(.) = "barcode";.} %>% mutate(barcode = paste0(sample, "#", barcode))

cell_df <- cell_df %>% mutate(malignancy = if_else(barcode %in% malign_cell$barcode, "malignant", "normal"))

write.table(cell_df, file = paste0(sample, "-", nth, "_cell_info.txt"), row.names = F, col.names = T, sep = "\t", quote = F)

