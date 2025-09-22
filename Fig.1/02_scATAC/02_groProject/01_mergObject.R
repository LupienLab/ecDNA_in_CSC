#******** step 1: merge multiple objects into a single merged object ********#

rm(list = ls())

library(ArchR) # 1.0.2
library(Seurat) # 5.2.1
library(data.table)
library(Matrix)
library(parallel)

addArchRGenome("hg38")
addArchRThreads(threads = floor(parallel::detectCores())/2, force = F)

tumorNames <- c("AA9S_T", "G797_T", "G837_T", "G900_T", "A7TK_T", "G933_T", "G946_T", "G958_T")

#~~~ arrow file paths 
arrow_dir <- "/Users/chupan/Documents/ArchR"
arrowFiles <- file.path(arrow_dir, paste0(tumorNames, ".arrow"))
names(arrowFiles) <- tumorNames

#~~~ generate a single merged ArchR project
obj <- ArchRProject(
  ArrowFiles = arrowFiles,
  outputDirectory = "outputs.gro",
  copyArrows = F,
  geneAnnotation = getGeneAnnotation(),
  genomeAnnotation = getGenomeAnnotation(),
  threads = getArchRThreads()
)

#~~~ filtercell file paths
filtercell_dir <- file.path(paste0("outputs.", tumorNames), paste0(tumorNames, ".filtercell.txt"))
names(filtercell_dir) <- tumorNames

filtercell_list <- lapply(names(filtercell_dir), function(id) {
  fp <- filtercell_dir[[id]]
  barcode <- paste0(id, "#", fread(fp, header = F, data.table = F)[, 1])
})
filtercell <- unlist(filtercell_list)

obj <- subsetArchRProject(
  ArchRProj       = obj,
  cells           = filtercell,
  outputDirectory = "outputs.gro",
  dropCells       = TRUE,
  threads = getArchRThreads(),
  force           = TRUE
)

obj <- addIterativeLSI(obj, useMatrix = "TileMatrix", name = "IterativeLSI", force = T) %>% 
  addHarmony(., reducedDims = "IterativeLSI", name = "Harmony", groupBy = "Sample") %>%
  addClusters(., reducedDims = "IterativeLSI", force = T) %>% addUMAP(., reducedDims = "IterativeLSI", force = T)

saveRDS(obj, file = paste0("./outputs.gro/mergedObject.rds"))
