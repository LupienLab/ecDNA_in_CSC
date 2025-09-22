#******** Chapter 6: create a umap plot per sample ********#

rm(list = ls())

library(ArchR)
library(Seurat)
library(tidydr)
library(data.table)
library(caret)
library(knitr)
library(dplyr)
library(spatstat.core)
library(mascarade)
library(spatstat.core)
library(BSgenome.Hsapiens.UCSC.hg38)
set.seed(2024)

addArchRGenome("hg38")
addArchRThreads(threads = floor(parallel::detectCores()/2), force = F)

proj <- readRDS(paste0("./outputs/", "object.rds"))

pdf(paste0("./outputs/umap.ecDNA.sample.pdf"), width = 5, height = 5)

celltype <- ifelse(proj$Clusters %in% c("C1", "C2"), "Myeloid", 
                   ifelse(proj$Clusters %in% c("C9"), "Oligodendrocyte",
                          ifelse(proj$Clusters %in% c("C10"), "Vascular",
                                 ifelse(proj$Clusters %in% c("C11", "C12", "C13"), "Neuron",
                                        ifelse(proj$Clusters %in% c("C16"), "OPC",
                                               ifelse(proj$Clusters %in% c("C15", "C17", "C18", "C19", "C20", "C21"), "ecDNA(-)", "ecDNA(+)"))))))
proj@cellColData$celltype <- celltype

umap_data <- getEmbedding(ArchRProj = proj, embedding = "UMAP") %>% as.data.frame() %>% {colnames(.) = c("umap_1", "umap_2");.} %>% 
mutate(cluster = proj$Clusters) %>% mutate(supcluster = proj@cellColData$celltype) %>% tibble::rownames_to_column(var = "barcode") %>%
mutate(sample = word(barcode, 1, sep = "#"))

centroids <- umap_data %>% group_by(supcluster) %>% summarize(umap_1 = mean(umap_1), umap_2 = mean(umap_2))
maskTable <- generateMask(dims = umap_data[, c("umap_1", "umap_2")], clusters = umap_data$supcluster, minDensity = 1.25, smoothSigma = 0.05)

custom_colors <- c("G958" = "purple", 
                   "G946" = "orchid",
                   "G933" = "orange",
                   "A7TK" = "salmon", 
                   "G900" = "yellowgreen",
                   "G837" = "olivedrab",
                   "G797" = "navy",
                   "AA9S" = "darkgreen")

p0 <- ggplot(umap_data, aes(x = umap_1, y = umap_2)) + 
  theme_classic() +
  geom_point(aes(color = sample), size = 0.05) + 
  scale_color_manual(values = custom_colors) +
  geom_path(data = maskTable, aes(group = group), linewidth = 0.25, linetype = 2) +
  geom_text(data = centroids, aes(label = supcluster, x = umap_1, y = umap_2), size = 2, fontface = "bold") + 
  coord_fixed()

p0 <- p0 + theme_classic() + 
  theme_dr(xlength = 0.25, ylength = 0.25, arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) + 
  theme(legend.key.size = unit(0.4, "inches"), panel.grid = element_blank()) + labs(x = "UMAP1", y = "UMAP2", title = c()) 
p0

dev.off()

#~~~
#~~~
#~~~
#~~~
#~~~
#~~~ only keep malignant cancer cells in the UMAP
pdf(paste0("./outputs/umap.malign.sample.pdf"), width = 5, height = 5)

celltype <- ifelse(proj$Clusters %in% c("C1", "C2"), "Myeloid", 
                   ifelse(proj$Clusters %in% c("C9"), "Oligodendrocyte",
                          ifelse(proj$Clusters %in% c("C10"), "Vascular",
                                 ifelse(proj$Clusters %in% c("C11", "C12", "C13"), "Neuron",
                                        ifelse(proj$Clusters %in% c("C16"), "OPC",
                                               ifelse(proj$Clusters %in% c("C14"), "A7TK-malign",
                                                      ifelse(proj$Clusters %in% c("C5", "C6", "C7", "C8", "C22"), "G958-malign",
                                                             ifelse(proj$Clusters %in% c("C3", "C4"), "G933-malign",
                                                                    ifelse(proj$Clusters %in% c("C23", "C24", "C25"), "G946-malign",
                                                                           ifelse(proj$Clusters %in% c("C19"), "AA9S-malign",
                                                                                  ifelse(proj$Clusters %in% c("C15", "C18"), "G900-malign",
                                                                                         ifelse(proj$Clusters %in% c("C17", "C21"), "G797-malign", "G837-malign"))))))))))))
proj@cellColData$celltype <- celltype

umap_data <- getEmbedding(ArchRProj = proj, embedding = "UMAP") %>% as.data.frame() %>% {colnames(.) = c("umap_1", "umap_2");.} %>% 
mutate(cluster = proj$Clusters) %>% mutate(supcluster = proj@cellColData$celltype) %>% tibble::rownames_to_column(var = "barcode") %>%
mutate(sample = word(barcode, 1, sep = "#"))

centroids <- umap_data %>% group_by(supcluster) %>% summarize(umap_1 = mean(umap_1), umap_2 = mean(umap_2))
maskTable <- generateMask(dims = umap_data[, c("umap_1", "umap_2")], clusters = umap_data$supcluster, minDensity = 1.25, smoothSigma = 0.05)

custom_colors <- c("G958-malign" = "purple", 
                   "G946-malign" = "orchid",
                   "G933-malign" = "orange",
                   "A7TK-malign" = "salmon", 
                   "G900-malign" = "yellowgreen",
                   "G837-malign" = "olivedrab",
                   "G797-malign" = "navy",
                   "AA9S-malign" = "darkgreen", 
                   "Myeloid" = "white", "Oligodendrocyte" = "white", "Vascular" = "white", "Neuron" = "white", "OPC" = "white")

p1 <- ggplot(umap_data, aes(x = umap_1, y = umap_2)) + 
  theme_classic() +
  geom_point(aes(color = celltype), size = 0.05) + 
  scale_color_manual(values = custom_colors) +
  geom_path(data = maskTable, aes(group = group), linewidth = 0.25, linetype = 2) +
  geom_text(data = centroids, aes(label = supcluster, x = umap_1, y = umap_2), size = 2, fontface = "bold") + 
  coord_fixed()

p1 <- p1 + theme_classic() + 
  theme_dr(xlength = 0.25, ylength = 0.25, arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) + 
  theme(legend.key.size = unit(0.4, "inches"), panel.grid = element_blank()) + labs(x = "UMAP1", y = "UMAP2", title = c()) 
p1

dev.off()










