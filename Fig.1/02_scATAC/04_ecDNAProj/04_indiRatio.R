rm(list = ls())

library(ArchR)
library(Seurat)
library(tidydr)
library(data.table)
library(caret)
library(knitr)
library(dplyr)
library(spatstat.core)
set.seed(2024)

addArchRGenome("hg38")
addArchRThreads(threads = floor(parallel::detectCores()/2), force = F)

id = "G933_T-1"

#~~~ input file
fileDir <- list.files("/Users/cpan/OneDrive/rawresult/celltype/")
fileDir <- fileDir[grep(paste0(id), fileDir)]
fileDir.positive <- fileDir[grep("ecDNAcell", fileDir)]
fileDir.negative <- fileDir[grep("normalLikecell", fileDir)]

#~~~ input cell types
ecDNAPoscell <- fread(paste0("/Users/cpan/OneDrive/rawresult/celltype/", fileDir.positive), header = F, stringsAsFactors = F) %>% {colnames(.) = "barcode";.}
ecDNANegcell <- fread(paste0("/Users/cpan/OneDrive/rawresult/celltype/", fileDir.negative), header = F, stringsAsFactors = F) %>% {colnames(.) = "barcode";.}

#~~~ input the project
proj <- readRDS(paste0("./outputs/", "object.rds"))
table(proj$Sample, proj$Clusters)

if(grepl("G958_T", id)){
  maligncell <- proj$cellNames[proj$Clusters %in% c("C5", "C6", "C7", "C8", "C22")]
  maligncell <- sub(".*#", "", maligncell)
}

if(grepl("G946_T", id)){
  maligncell <- proj$cellNames[proj$Clusters %in% c("C23", "C24", "C25")]
  maligncell <- sub(".*#", "", maligncell)
}

if(grepl("G933_T", id)){
  maligncell <- proj$cellNames[proj$Clusters %in% c("C3", "C4")]
  maligncell <- sub(".*#", "", maligncell)
}

if(grepl("A7TK_T", id)){
  maligncell <- proj$cellNames[proj$Clusters %in% c("C14")]
  maligncell <- sub(".*#", "", maligncell)
}

ecDNAPoscells <- intersect(maligncell, ecDNAPoscell$barcode)
ecDNANegcells <- intersect(maligncell, ecDNANegcell$barcode)

write.table(ecDNAPoscells, file = paste0("./outputs/", id, ".ecDNAcell.txt"), row.names = F, col.names = F, quote = F)
write.table(ecDNANegcells, file = paste0("./outputs/", id, ".cancercell.txt"), row.names = F, col.names = F, quote = F)


data <- data.frame(
  sample = rep(id, each = 1),
  category = rep(c("1", "2"), times = 1),
  number = c(length(ecDNAPoscells), length(ecDNANegcells))
)

data_percentage <- data %>% group_by(sample) %>% mutate(percentage = number / sum(number) * 100)
data_percentage <- data_percentage %>% filter(sample == id)

if(id == "A7TK_T-1" | id == "A7TK_T-2") colorPalette = c("1" = "salmon", "2" = "#FFA07A")
if(id == "G933_T-1" | id == "G933_T-2" | id == "G933_T-3") colorPalette = c("1" = "orange", "2" = "#FFDAB9")
if(id == "G946_T-1") colorPalette = c("1" = "orchid", "2" = "#F2D4E6")
if(id == "G958_T-1") colorPalette = c("1" = "purple", "2" = "#E6E6FA")
celltype = c("ecDNA(+)", "ecDNA(-)", "normal-like") 

if( id == "G958_T-1"){
  ggplot(data_percentage, aes(x = sample, y = percentage, fill = category)) + theme_bw() +
    geom_bar(stat = "identity") + coord_flip() +
    scale_fill_manual(values = colorPalette, labels = celltype) +
    geom_hline(yintercept = data_percentage$percentage[2], linetype = "dashed", color = "black", linewidth = 0.25) + 
    geom_text(aes(label = paste0(round(percentage, 1), "%")), position = position_stack(vjust = 0.5), size = 3) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, color = "black"),
          axis.text.y = element_text(color = "black"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    labs(title = NULL, x = NULL, y = "Percentage (%)", fill = "cell type") 
}
if(id == "A7TK_T-1" | id == "A7TK_T-2" |id == "G933_T-1" | id == "G933_T-2" | id == "G933_T-3" | id == "G946_T-1" ){
  ggplot(data_percentage, aes(x = sample, y = percentage, fill = category)) + theme_bw() +
    geom_bar(stat = "identity") + coord_flip() +
    scale_fill_manual(values = colorPalette, labels = celltype) +
    geom_hline(yintercept = data_percentage$percentage[2], linetype = "dashed", color = "black", linewidth = 0.25) + 
    geom_text(aes(label = paste0(round(percentage, 1), "%")), position = position_stack(vjust = 0.5), size = 3) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, color = "black"),
          axis.text.y = element_text(color = "black"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    labs(title = NULL, x = NULL, y = "Percentage (%)", fill = "cell type") 
}

