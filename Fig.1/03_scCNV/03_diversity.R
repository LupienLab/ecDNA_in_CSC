rm(list = ls())

library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(data.table)
library(DescTools)
library(Informeasure)
library(SingleCellExperiment)
'%ni%' <- Negate('%in%')

#~~~ define region of interest
# EGFR: chr7_50000001_60000000
# AGAP2 and CDK4: chr12_55500001_65500000
# MDM2:chr12_65500001_75500000
# chr19_1_10000000, chr6_89800001_99800000, chr3_150900001_160900000, chr3_20000001_30000000, chr5_30000001_40000000

#id = "A7TK_T"
#scCNA <- readRDS(file = paste0("./outputs.", id, "/", id, ".scCNA.rds"))
#regions <- ROWNAMES(scCNA)

region = "chr7_50000001_60000000"

sample_ids <- c("A7TK_T", "G933_T", "G946_T", "G958_T", "AA9S_T", "G797_T", "G837_T", "G900_T")  # replace with your actual sample names
  
data_list <- lapply(sample_ids, function(id) {
  cat("id: ", id, "\n")
  rds_path <- file.path(paste0("outputs.", id), paste0(id, ".scCNA.rds"))
  scCNA_obj <- readRDS(rds_path)
    
  # extract malignant cells
  malign_cell <- fread(paste0("/Users/chupan/Documents/ArchR/outputs.", id, "/", id, ".maligncell.txt"), header = F, data.table = F) %>% {colnames(.) = "barcode";.} %>% pull(barcode)
  malign_cells <- scCNA_obj$Barcode[scCNA_obj$Barcode %in% malign_cell]
    
  cat("id: ", id, " length(malign_cells),", length(malign_cells), "\n")
  # extract copy numbers for the given region
  cna_mat <- assay(scCNA_obj, "counts_gc_modal_smoothed_ratios")
    
  # check that the region and cells exist in the matrix
  if (!(region %in% rownames(cna_mat)) || length(malign_cells) == 0) {
    return(NULL)  # skip this sample if data is missing
  }
    
  values <- as.numeric(cna_mat[region, malign_cells]) #%>% sample(., size = 750, replace = F)
    
  # check again that values are non-empty
  if (length(values) == 0) {
    return(NULL)
  }
    
  # return a data.frame with cna and sample label
  data.frame(cna = values, sample = id) %>% drop_na()
})
  
# combine non-null results
data_df <- bind_rows(data_list) 
data_df$sample <- factor(data_df$sample, levels = sample_ids)
  
pal = c("G958_T" = "purple", 
        "G946_T" = "orchid",
        "G933_T" = "orange",
        "A7TK_T" = "salmon", 
        "G900_T" = "yellowgreen",
        "G837_T" = "olivedrab",
        "G797_T" = "navy",
        "AA9S_T" = "darkgreen")
  
diversity <- lapply(data_list, function(df) {
  counts <- df$cna
  IQR(counts)})
  
diversity <- do.call(rbind, diversity) %>% as.data.frame() %>% {colnames(.) = c("IRQ");.} %>% mutate(id = sample_ids)
diversity$id <- factor(diversity$id, levels = sample_ids)
  
# visualize the copy number variation result
dir.create("outputs.diversity")
pdf(file = paste0("outputs.diversity", "/", region, ".pdf"), width = 10, height = 4)
  
p <- ggplot(data_df, aes(x = sample, y = cna)) + 
  theme_bw() +
  geom_jitter(aes(colour = sample), width = 0.3, size = 0.05) +
  geom_violin(aes(fill = sample), trim = F) +
  geom_text(data = diversity, aes(x = id, y = 6, label = IRQ), inherit.aes = FALSE, size = 5) +
  labs(title = paste0(region), x = NULL, y = "single-cell CNA") +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.25) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
    axis.text.y = element_text(colour = "black"),
    legend.position = "none",
    panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(p)
  
dev.off()

