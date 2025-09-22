#******** step 2: create an UMAP plot for each tumor ********#

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

id = "AA9S_T"

obj <- readRDS(paste0("./outputs.", id, "/", id, ".obj.rds"))

metacluster <- data.frame(barcode = gsub(paste0(id, "#"), "", obj$cellNames), cluster = obj$Clusters)
write.table(metacluster, file = paste0("outputs.", id, "/", id, ".metacluster.txt"), row.names = F, col.names = F, sep = "\t", quote = F)

pdf(paste0("./outputs.", id, "/", id, ".umap.cluster.pdf"), width = 4, height = 6)

umap_data <- getEmbedding(ArchRProj = obj, embedding = "UMAP") %>% as.data.frame() %>% {colnames(.) = c("umap_1", "umap_2");.} %>% mutate(cluster = obj$Clusters)
centroids <- umap_data %>% group_by(cluster) %>% summarize(umap_1 = mean(umap_1), umap_2 = mean(umap_2))
maskTable <- generateMask(dims = umap_data[, c("umap_1", "umap_2")], clusters = umap_data$cluster, minDensity = 1.5, smoothSigma = 0.05)

p1 <- ggplot(umap_data, aes(x = umap_1, y = umap_2)) + 
  geom_point(aes(color = cluster), size = 0.05) + 
  geom_path(data = maskTable, aes(group = group), linewidth = 0.25, linetype = 2) +
  geom_text(data = centroids, aes(label = cluster, x = umap_1, y = umap_2), size = 2) + 
  coord_fixed()+ 
  theme_classic()

p1 <- p1 + theme_classic() + guides(color = guide_legend(title = "cell type", override.aes = list(size = 1.5))) +
  theme_dr(xlength = 0.25, ylength = 0.25, arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) + 
  theme(panel.grid = element_blank()) + labs(x = "UMAP1", y = "UMAP2", title = c()) 
p1

dev.off()



