#******** step 1: create an object for each tumor ********#

rm(list = ls())

id = "G958_T"

library(ArchR) #~~~ packageVersion("ArchR") 1.0.2
library(Seurat) # packageVersion("Seurat") 5.2.1
library(data.table)
library(ggplot2)
library(Matrix)
library(dplyr)

addArchRGenome("hg38")
addArchRThreads(threads = floor(parallel::detectCores())/2, force = F)
inputFiles <- getInputFiles(paths = "/Users/chupan/Documents/data_su2c/scATAC/fragment/")
inputFiles <- inputFiles[which(names(inputFiles) == id)]

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 5,
  minFrags = 3000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  threads = getArchRThreads(),
  force = FALSE
)

doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10,
  knnMethod = "UMAP",
  LSIMethod = 1
)

obj <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = paste0("outputs.", id),
  copyArrows = FALSE,
  geneAnnotation = getGeneAnnotation(),
  genomeAnnotation = getGenomeAnnotation(),
  showLogo = TRUE,
  threads = getArchRThreads()
)

obj <- filterDoublets(obj) %>% addIterativeLSI(., useMatrix = "TileMatrix", name = "IterativeLSI", force = T) %>% addClusters(., reducedDims = "IterativeLSI", force = T) %>% addUMAP(., reducedDims = "IterativeLSI", force = T)

filtercell <- str_replace(obj$cellNames, paste0(id, "#"), "")
write.table(filtercell, file = paste0("./outputs.", id, "/", id, ".filtercell.txt"), row.names = F, col.names = F, quote = F)
saveRDS(obj, file = paste0("./outputs.", id, "/", id, ".obj.rds"))

#>>> TSS = 5, Frags = 3000 <<<#
#--- A7TK:  1905
#--- AA9S:  1518
#--- G958:  7042
#--- G946:  9096 
#--- G933:  4749
#--- G900:  2824 
#--- G837:  4477
#--- G797:  7829
