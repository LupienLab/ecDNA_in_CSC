rm(list = ls())

library(scatools)
library(tidyverse)
library(purrr)
library(ggpubr)
library(data.table)
library(BSgenome.Hsapiens.UCSC.hg38)
'%ni%' <- Negate('%in%')

#~~~ define region of interest
# EGFR: chr7_50000001_60000000 |G933_T-2, G946_T-1, G958_T-1
# AGAP2 and CDK4: chr12_55500001_65500000  |A7TK_T-1, G933_T-1, G958_T-1
# MDM2:chr12_65500001_75500000 |A7TK_T-2, G933_T-3, G958_T-1
# chr17_45100001_55100000 , chr5_78800001_88800000 chr1_30000001_40000000

id = "G958_T"
n = 1
region <- "chr12_65500001_75500000"

scCNA <- readRDS(file = paste0("./outputs.", id, "/", id, ".scCNA.rds"))

positive_cell <- fread(
  paste0("/Users/chupan/Documents/data_su2c/scATAC/celltype/", id, "-", n, ".ecDNA_positive_cell.txt"), 
  header = F, data.table = F) %>%
  {colnames(.) = "barcode";.} %>% pull(barcode)

negative_cell <- fread(
  paste0("/Users/chupan/Documents/data_su2c/scATAC/celltype/", id, "-", n, ".ecDNA_negative_cell.txt"), 
  header = F, data.table = F) %>%
  {colnames(.) = "barcode";.} %>% pull(barcode)

cna_df <- scCNA@assays@data$counts_gc_modal_smoothed_ratios %>% as.data.frame()

cna_values <- data.frame(
 value = c(as.numeric(cna_df[region, positive_cell]), as.numeric(cna_df[region, negative_cell])),
  group = c(rep("positive_cell", length(positive_cell)), rep("negative_cell", length(negative_cell)))
)

if(id == "A7TK_T") pal = c("positive_cell" = "salmon", "negative_cell" = "#FFA07A")
if(id == "G933_T") pal = c("positive_cell" = "orange", "negative_cell" = "#FFDAB9")
if(id == "G946_T") pal = c("positive_cell" = "orchid", "negative_cell" = "#F2D4E6")
if(id == "G958_T") pal = c("positive_cell" = "purple", "negative_cell" = "#E6E6FA")

pval <- wilcox.test(value ~ group, data = cna_values)$p.value
pval_text <- paste0("p = ", signif(pval, 3))

dir.create(paste0("outputs.", id))
pdf(paste0("./outputs.", id, "/", id, "-", n, "-", region, ".pdf"), width = 2, height = 3)

ggplot(cna_values, aes(x = group, y = value)) + 
  theme_bw() +
  geom_jitter(aes(colour = group), width = 0.3, size = 0.05) +
  geom_violin(aes(fill = group), trim = F) +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.25) +
  labs(title = region, x = NULL, y = "single-cell CNA") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
    axis.text.y = element_text(colour = "black"),
    legend.position = "none",
    panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  stat_compare_means(method = "wilcox.test", label = "p.format")

dev.off()

