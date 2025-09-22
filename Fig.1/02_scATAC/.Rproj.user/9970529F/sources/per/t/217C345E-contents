#******** step 2: create a umap plot for stem cells ********#

rm(list = ls())

library(ArchR)
library(Seurat)
library(tidydr)
library(dplyr)
library(caret)
library(knitr)
library(data.table)
library(mascarade)

library(BSgenome.Hsapiens.UCSC.hg38)
set.seed(2024)

addArchRGenome("hg38")
addArchRThreads(threads = floor(parallel::detectCores()), force = F)

id = "A7TK_T"
n = 1

obj <- readRDS(paste0("./outputs.", id, "/", id, ".obj.rds"))
metacluster <- data.frame(barcode = gsub(paste0(id, "#"), "", obj$cellNames), cluster = obj$Clusters)

ecDNA_content <- fread(paste0("/Users/chupan/Documents/whereAreEcDNAs/package/", id, ".outputs/", id, "-", n, ".ecDNA_content.txt"), header = F, stringsAsFactors = F) %>% 
  as.data.frame() %>% {colnames(.) = c("barcode", "uscore");.} %>% mutate(barcode = paste0(id, "#", barcode))

CSC <- readRDS("/Users/chupan/Documents/data_su2c/scATAC/CSC/overlap_cells.rds") %>% filter(type == "T" & sample == id) %>% 
  dplyr::select(barcode)
CSC_ecDNA_content <- CSC %>% inner_join(ecDNA_content, by = "barcode")

dir.create("outputs.CSC")
pdf(paste0("./outputs.CSC/", id, "-", n, ".umap.CSC.cluster.pdf"), width = 6, height = 5) #width = 4, height = 5

umap_data <- getEmbedding(ArchRProj = obj, embedding = "UMAP") %>% 
  as.data.frame() %>% {colnames(.) = c("umap_1", "umap_2");.} %>% mutate(cluster = obj$Clusters) %>% tibble::rownames_to_column(var = "barcode") 
merged_data <- CSC_ecDNA_content %>% left_join(umap_data , by = "barcode")
centroids <- umap_data %>% group_by(cluster) %>% summarize(umap_1 = mean(umap_1), umap_2 = mean(umap_2))
maskTable <- generateMask(dims = umap_data[, c("umap_1", "umap_2")], clusters = umap_data$cluster, minDensity = 1.5, smoothSigma = 0.05)

solarExtra = c('#14B3FF', '#FDB31A', '#E42A2A', '#A31D1D') 

p0 <- ggplot(merged_data, aes(x = umap_1, y = umap_2)) + 
  theme_bw() +
  geom_point(aes(color = uscore), size = 1) + 
  geom_path(data = maskTable, aes(group = group), linewidth = 0.25, linetype = 2) +
  geom_text(data = centroids, aes(label = cluster, x = umap_1, y = umap_2), size = 3) + 
  scale_color_gradientn(colours = solarExtra) + 
  theme(panel.grid = element_blank()) + labs(x = "UMAP1", y = "UMAP2", title = c()) +
  coord_fixed() 

p0

dev.off()


