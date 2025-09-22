rm(list = ls())

library(data.table)
library(magrittr)
library(tidyr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(ComplexHeatmap)

#--- input the gene list data
fileList <- list.files("inputs/ecDNAList")
splitelements <- strsplit(fileList, "_")
sampleName <- sapply(splitelements, function(x) paste0(x[1], "_", x[2]))

#--- the first file 
ecDNAList <- fread(paste0("inputs/ecDNAList/", fileList[1])) %>% as.data.frame() %>% magrittr::set_colnames(c("sample_name", "count"))
if(dim(ecDNAList)[1] == 0){
  ecDNAList <- data.frame(sample_name =  paste0(unlist(strsplit(fileList[1], "_"))[1], "_", unlist(strsplit(fileList[1], "_"))[2]), count = 0)
} 

for(i in 2:length(sampleName)){
  cat("sampleName: ", sampleName[i], "\n")
  ecDNAListNew <- fread(paste0("inputs/ecDNAList/", fileList[i])) %>% as.data.frame() %>% magrittr::set_colnames(c("sample_name", "count"))
  if(dim(ecDNAListNew)[1] == 0){
    ecDNAListNew <- data.frame(sample_name = paste0(unlist(strsplit(fileList[i], "_"))[1], "_", unlist(strsplit(fileList[i], "_"))[2]), count = 0) 
  } 
  ecDNAList <- rbind(ecDNAList, ecDNAListNew)
}

ecDNAList.df <- ecDNAList %>% group_by(substr(sample_name, 1, nchar(sample_name) - 2)) %>% filter(any(count != 0)) %>% as.data.frame() %>% dplyr::select(sample_name, count)

dir.create("outputs")
write.table(ecDNAList.df, file = paste0("./outputs/", "ecDNAList.bed"), sep = "\t", row.names = F, col.names = T)

