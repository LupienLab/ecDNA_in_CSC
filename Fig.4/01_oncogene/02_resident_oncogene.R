rm(list = ls())

library(data.table)
library(magrittr)
library(tidyr)
library(dplyr)
library(pheatmap)

#~~~ input the ecDNA list data
ecDNAList <- fread("./outputs/ecDNAList.bed") %>% as.data.frame() 

#--- input the gene list data
fileList <- list.files("inputs/oncogeneList")
splitelements <- strsplit(fileList, "_")
sampleName <- sapply(splitelements, function(x) paste0(x[1], "_", x[2]))

#~~~ the first file 
geneList <- fread(paste0("./inputs/oncogeneList/", fileList[1])) %>% as.data.frame() %>% filter(grepl("ecDNA", feature) & is_canonical_oncogene == TRUE) %>% dplyr::select(sample_name, gene)
if(dim(geneList)[1] == 0){
  geneList <- data.frame(sample_name =  paste0(unlist(strsplit(fileList[1], "_"))[1], "_", unlist(strsplit(fileList[1], "_"))[2]), gene = "")
} 

for(i in 2:length(sampleName)){
  cat("sampleName: ", sampleName[i], "\n")
  geneListNew <- fread(paste0("inputs/oncogeneList/", fileList[i])) %>% as.data.frame() %>% filter(grepl("ecDNA", feature) & is_canonical_oncogene == TRUE) %>% dplyr::select(sample_name, gene)
  if(dim(geneListNew)[1] == 0){
    geneListNew <- data.frame(sample_name = paste0(unlist(strsplit(fileList[i], "_"))[1], "_", unlist(strsplit(fileList[i], "_"))[2]), gene = "") 
  } 
  geneList <- rbind(geneList, geneListNew)
}

geneTable <- table(geneList)
names(attr(geneTable, "dimnames")) <- NULL
class(geneTable) <- "matrix"
geneTable <- t(geneTable) 

geneTable <- geneTable %>% as.data.frame() %>% dplyr::select(ecDNAList$sample_name) %>% as.matrix()

annotation_col = data.frame(
  ecDNA = as.character( ecDNAList$count )
)
rownames(annotation_col) = colnames(geneTable)
head(annotation_col)

pdf(file = paste0("./outputs/ecDNA_oncogene.pdf"), height = 6, width = 5)

p <- pheatmap(geneTable, cellwidth = 10, cellheight = 15, cluster_row = TRUE, cluster_col = FALSE,
              border = T, border_color = "black", color = c("white", "red"), legend_breaks = c(0, 1), legend_labels = c("0", "1"), 
              annotation_col = annotation_col, annotation_colors = list(ecDNA = c("1" = "#e5f5e0", "2" = "#a1d99b", "3" = "#41ab5d")), 
              gaps_col = c(2,4,6,8), angle_col = "45")

p

dev.off()


