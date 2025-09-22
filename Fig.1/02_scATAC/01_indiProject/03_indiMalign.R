#******** step 3: extract malignant cells for each tumor ********#

rm(list = ls())

library(ArchR)
library(Seurat)
library(tidydr)
library(dplyr)
library(caret)
library(knitr)
library(data.table)
library(mascarade)
'%ni%' <- Negate('%in%')

library(BSgenome.Hsapiens.UCSC.hg38)
set.seed(2024)

addArchRGenome("hg38")

id = "A7TK_T"
obj <- readRDS(paste0("./outputs.", id, "/", id, ".obj.rds"))

if(id == "A7TK_T")
  malign_cell <- obj$cellNames[which(obj$Clusters %ni% c("C1", "C2", "C3"))] 
if(id == "AA9S_T")
  malign_cell <- obj$cellNames[which(obj$Clusters %ni% c("C1", "C2", "C3", "C4"))] 
if(id == "G797_T")
  malign_cell <- obj$cellNames[which(obj$Clusters %ni% c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9"))] 
if(id == "G837_T")
  malign_cell <- obj$cellNames[which(obj$Clusters %ni% c("C1", "C2", "C3", "C5", "C12", "C13", "C14", "C15"))] 
if(id == "G900_T")
  malign_cell <- obj$cellNames[which(obj$Clusters %ni% c("C1", "C2", "C3", "C4", "C5", "C11", "C12", "C13"))] 
if(id == "G933_T")
  malign_cell <- obj$cellNames[which(obj$Clusters %ni% c("C1", "C2", "C3"))] 
if(id == "G946_T")
  malign_cell <- obj$cellNames[which(obj$Clusters %ni% c("C1", "C2", "C11", "C12", "C13"))] 
if(id == "G958_T")
  malign_cell <- obj$cellNames[which(obj$Clusters %ni% c("C1", "C2"))] 

write.table(gsub(paste0(id,"#"), "", malign_cell), file = paste0("./outputs.", id, "/", id, ".maligncell.txt"), row.names = F, col.names = F, quote = F)


