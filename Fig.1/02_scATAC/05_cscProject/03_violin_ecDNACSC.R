#******** step 3: create a umap plot per sample ********#

rm(list = ls())

library(ArchR)
library(Seurat)
library(dplyr)
library(tidydr)
library(data.table)
library(caret)
library(knitr)
library(patchwork)
library(BSgenome.Hsapiens.UCSC.hg38)
set.seed(2025)

addArchRGenome("hg38")
addArchRThreads(threads = floor(parallel::detectCores()/2), force = F)

id = "G958_T"
n = 1

obj <- readRDS(paste0("./outputs.", id, "/", id, ".obj.rds"))
cell_cluster_number <- table(obj$Clusters) %>% as.data.frame() %>% {colnames(.) = c("cluster", "count");.}

uscore <- fread(paste0("/Users/chupan/Documents/whereAreEcDNAs/package/", id, ".outputs/", id, "-", n, ".ecDNA_content.txt"), header = F, stringsAsFactors = F) %>% 
  as.data.frame() %>% {colnames(.) = c("barcode", "uscore");.} %>% mutate(barcode = paste0(id, "#", barcode))

umap_data <- getEmbedding(ArchRProj = obj, embedding = "UMAP") %>% 
  as.data.frame() %>% {colnames(.) = c("umap_1", "umap_2");.} %>% mutate(cluster = obj$Clusters) %>% tibble::rownames_to_column(var = "barcode")

#~~~ cell + cluster + uscore
cell_cluster_uscore <- umap_data %>% dplyr::inner_join(uscore, by = "barcode") %>% dplyr::select(barcode, cluster, uscore) %>% 
  mutate(cluster = factor(cluster, levels = mixedsort(unique(cluster))))

#~~~ CSCs
overlap_df <- readRDS(file = paste0("/Users/chupan/Documents/data_su2c/scATAC/CSC/overlap_cells.rds"))
stem_cell <- overlap_df %>% dplyr::filter(type == "T") 
stem_cell_cluster_uscore <- stem_cell %>% inner_join(cell_cluster_uscore, by = "barcode")

#~~~ ecDNA-positive cells
ecDNA_cell <- fread(paste0("/Users/chupan/Documents/whereAreEcDNAs/package/", id, ".outputs/", id, "-", n, ".ecDNA_positive_cell.txt"), header = F, stringsAsFactors = F) %>% 
  as.data.frame() %>% {colnames(.) = "barcode";.} %>% mutate(barcode = paste0(id, "#", barcode))
cancer_cell <- fread(paste0("/Users/chupan/Documents/whereAreEcDNAs/package/", id, ".outputs/", id, "-", n, ".ecDNA_negative_cell.txt"), header = F, stringsAsFactors = F) %>% 
  as.data.frame() %>% {colnames(.) = "barcode";.} %>% mutate(barcode = paste0(id, "#", barcode))
ecDNA_cell_thr <- ecDNA_cell %>% inner_join(uscore, by = "barcode") %>% select(uscore) %>% pull(uscore) %>% min(.)

#~~~~ ecDNA, malignant + cluster
ecDNA_cell_cluster <- ecDNA_cell %>% inner_join(umap_data, by = "barcode")
ecDNA_cell_cluster_number <- table(ecDNA_cell_cluster$cluster) %>% as.data.frame() %>% {colnames(.) = c("cluster", "number");.} 
ecDNA_cell_cluster_number <- full_join(ecDNA_cell_cluster_number, cell_cluster_number, by = "cluster") %>% 
  tidyr::replace_na(list(number = 0)) %>% mutate(cluster = factor(cluster, levels = mixedsort(unique(cluster)))) 

#~~~ CSC +/- ecDNA-positive cells
ecDNA_stem_cell <- stem_cell %>% inner_join(ecDNA_cell, by = "barcode")

#~~~ CSC +/- ecDNA-positive cell + cluster + number
ecDNA_stem_cell_cluster <- ecDNA_stem_cell %>% inner_join(umap_data, by = "barcode")
ecDNA_stem_cell_cluster_number <- table(ecDNA_stem_cell_cluster$cluster) %>% as.data.frame() %>% {colnames(.) = c("cluster", "number");.}
ecDNA_stem_cell_cluster_number <- full_join(ecDNA_stem_cell_cluster_number, cell_cluster_number, by = "cluster") %>% 
  tidyr::replace_na(list(number = 0)) %>% mutate(cluster = factor(cluster, levels = mixedsort(unique(cluster)))) %>%
  mutate(percentage = as.integer(number)/as.integer(count))

#~~~ CSC +/- ecDNA cell + cluster + uscore
cell_cluster_uscore_df <- cell_cluster_uscore %>% mutate(celltype = ifelse(barcode %in% ecDNA_stem_cell$barcode, "ecDNA_CSC", "other")) %>%
  mutate(cluster = factor(cluster, levels = mixedsort(unique(cluster))))

width <- if(id == "A7TK_T") 4.00 else 4
width <- if(id == "G933_T") 4.50 else 4
width <- if(id == "G946_T") 5.75 else 4
width <- if(id == "G958_T") 5.50 else 4

ecDNA_stem_cell_cluster_percentage_df <- ecDNA_stem_cell_cluster_number %>% dplyr::select(cluster, percentage) 
ecDNA_stem_cell_cluster_percentage_df$cluster <- factor(
  ecDNA_stem_cell_cluster_percentage_df$cluster,
  levels = paste0("C", sort(as.numeric(gsub("C", "", levels(ecDNA_stem_cell_cluster_percentage_df$cluster)))))
)

#~~~ ecDNA-positive CSC per cluster
dir.create("outputs.CSC")
pdf(paste0("./outputs.CSC/", id, "-", n, "_CSC_cluster.pdf"), width = width, height = 2)

violin_plot <- ggplot(cell_cluster_uscore_df, aes(x = cluster, y = uscore)) +
  theme_bw() +
  geom_violin(trim = T, width = 2.5) +  
  geom_text(data = ecDNA_stem_cell_cluster_number, aes(x = cluster, y = 0.6, label = number), inherit.aes = FALSE, angle = 45, size = 2, color = "black") +
  geom_jitter(data = cell_cluster_uscore_df %>% filter(celltype == "other"), aes(color = "grey"), width = 0.1, size = 0.0000015) + 
  geom_jitter(data = cell_cluster_uscore_df %>% filter(celltype == "ecDNA_CSC"), aes(color = celltype), width = 0.15, size = 0.000005) + 
  geom_hline(yintercept = ecDNA_cell_thr, linetype = "dashed", color = "red", size = 0.15) + 
  scale_color_manual(values = c("ecDNA_CSC" = "blue", "other" = "black")) +  
  labs( x = NULL, y = "ecDNA content score", title = paste0(id, "-", n)) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"), 
        axis.text.y = element_text(colour = "black"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
violin_plot

dev.off()
