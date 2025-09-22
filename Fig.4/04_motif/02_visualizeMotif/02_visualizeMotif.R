rm(list = ls())

library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(data.table)
library(purrr)

#~~~ define TF family
tf_family <- list(
  bHLH  = c("TWIST1", "TWIST2", "OLIG2", "ASCL1", "ASCL2", "NEUROD1", "NEUROG1", "TCF4", "TCF7", "TCF12", "TCF3", "TCF4", "TCF21"),
  SOX   = c("SOX1", "SOX2", "SOX3", "SOX4", "SOX5", "SOX6", "SOX7", "SOX8", "SOX9", "SOX10", "SOX11", "SOX12", "SOX14", "SOX15", "SOX17", "SOX21", "SOX30"),
  POU   = c("POU3F2", "BRN2", "POU5F1", "OCT4", "POU5F1::SOX2", "POU3F1", "OCT6", "POU3F3", "BRN1", "POU2F2", "OCT2", "POU2F1", "OCT1", "POU6F2", "POU4F3", "POU6F1"),
  KLF   = c("KLF1", "KLF4", "KLF5", "KLF7", "KLF9", "KLF10", "KLF12", "KLF14", "KLF16"),
  RUNX  = c("RUNX1", "RUNX2"),
  SMAD  = c("SMAD2", "SMAD3", "SMAD4"),
  Tbox  = c("TBX1", "TBX4", "TBX5", "TBX15", "TBX20", "TBOX:SMAD", "EOMES"),
  AP1   = c("AP-1", "FOSL1", "FOSL2"),
  AP2   = c("TFAP2A", "TFAP2B", "TFAP2C", "TFAP2D", "TFAP2E", "TCFAP2B", "TCFAP2E")
)

#~~~ load and process files
files <- list.files(path = "./outputs/", pattern = "*.RData", full.names = TRUE)

loadFiles <- function(file) {
  load(file)
  sample_name <- gsub("\\.motif\\.RData$", "", basename(file))
  motif_df %>% dplyr::select(motif, neglogPvalue) %>% mutate(sample = sample_name)
}

motif_data <- bind_rows(lapply(files, loadFiles)) %>% filter(motif != "") %>% as.data.frame()

#~~~ collapse motifs into TF families
motif_to_family <- map2_dfr(names(tf_family), tf_family, ~ data.frame(motif = .y, family = .x, stringsAsFactors = FALSE))

# Left join and replace motif names with family name where applicable
motif_data <- motif_data %>% left_join(motif_to_family, by = "motif") %>% mutate(motif = if_else(!is.na(family), family, motif)) %>% select(-family) %>% filter(motif != "")
unique(motif_data$motif)

filtered_df <- motif_data %>% group_by(motif, sample) %>% slice_max(order_by = neglogPvalue, n = 1, with_ties = FALSE) %>% ungroup() %>% as.data.table()

data <- data.table::dcast(filtered_df, sample ~ motif, value.var = "neglogPvalue") %>% magrittr::set_rownames(.$sample) %>% select(-sample)

heatmap_data <- motif_data %>% group_by(sample, motif) %>%
  summarise(neglogPvalue = mean(neglogPvalue, na.rm = TRUE)) %>%  # handle duplicates
  pivot_wider(names_from = motif, values_from = neglogPvalue, values_fill = 0) %>%  # fill NAs with 0
  column_to_rownames("sample")  %>% select(sort(names(.))) # set rownames as sample names

#~~~ convert to matrix for heatmap
heatmap_matrix <- as.matrix(heatmap_data)

pdf(paste0("./plts/GSC_motif_family.pdf"), width = 6, height = 1.5)

p <- pheatmap::pheatmap(heatmap_matrix, cluster_rows = F, cluster_cols = F, border_color = "black", lwd = 0.5,
                          color = colorRampPalette(c("white", "#0066FF33", "#3300FF33"))(200))
                          
print(p)

dev.off()

