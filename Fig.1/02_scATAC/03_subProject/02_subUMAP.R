#******** step 2: create a UMAP plot for sub objects ********#

rm(list = ls())

library(ArchR) #package.version("ArchR")
library(Seurat) # package.version("Seurat") 
library(tidydr)
library(data.table)
library(caret)
library(knitr)
library(dplyr)

set.seed(2024)

addArchRGenome("hg38")
addArchRThreads(threads = floor(parallel::detectCores()), force = F)

obj <- readRDS(paste0("./outputs.sub/subObject.rds"))

pdf(paste0("./outputs.sub/umap.cluster.subObject.pdf"), width = 6, height = 6)
p0 <- plotEmbedding(obj, colorBy = "cellColData", continuousSet = NULL, name = "Clusters", embedding = "UMAP")
p0 <- p0 + theme_classic() + guides(color = guide_legend(title = "cell type", override.aes = list(size = 1.5))) +
  theme_dr(xlength = 0.25, ylength = 0.25, arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) + 
  theme(panel.grid = element_blank()) + labs(x = "UMAP1", y = "UMAP2", title = c()) 
p0
dev.off()

pdf(paste0("./outputs.sub/umap.sample.subObject.pdf"), width = 6, height = 6)
p1 <- plotEmbedding(ArchRProj = obj, colorBy = "cellColData", 
                    pal = c("G900_T" = "yellowgreen",
                            "G837_T" = "olivedrab",
                            "G797_T" = "navy",
                            "AA9S_T" = "darkgreen"), 
                    continuousSet = NULL, name = "Sample", embedding = "UMAP", labelMeans = F) #, labelMeans = FALSE
p1 <- p1 + theme_classic() + guides(color = guide_legend(title = "cell type", override.aes = list(size = 1.5))) +
  theme_dr(xlength = 0.25, ylength = 0.25, arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) + 
  theme(panel.grid = element_blank()) + labs(x = "UMAP1", y = "UMAP2", title = c()) 
p1
dev.off()


