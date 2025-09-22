#******** step 3: gene activity score for marker genes ********#

rm(list = ls())

library(ArchR)
library(data.table)
library(dplyr)
library(tidyverse)
library(caret)
library(cowplot)
library(readxl)
library(ggplot2)
library(ggnewscale)
library(gridExtra)
library(tidyr)
library(BSgenome.Hsapiens.UCSC.hg38)

addArchRGenome("hg38")

addArchRThreads(threads = floor(parallel::detectCores()), force = F)

obj <- readRDS(paste0("./outputs.gro/mergedObject.rds"))
obj <- addImputeWeights(obj) # use MAGIC to impute gene scores by smoothing signal across nearby cells

markerGenes  <- c(
  "EGFR",   # malignant
  "MDM2",   # malignant
  "CDK4",   # malignant
  "AGAP2",  # malignant
  "CSF1R",  # myeloid
  "PECAM1", # myeloid
  "CD14",   # myeloid
  "ESAM",   # vascular
  "CSPG4",  # OPC
  "ASCL1",  # OPC
  "RBFOX3", # neuron
  "DLGAP1", # neuron
  "CLDN11", # oligodendrocyte
  "SOX10"   # oligodendrocyte
)
   
p <- plotEmbedding(
  ArchRProj = obj, 
  colorBy = "GeneScoreMatrix", 
  name = markerGenes, 
  embedding = "UMAP",
  imputeWeights = getImputeWeights(obj)
) 

p

