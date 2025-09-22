rm(list = ls())

library(data.table)
library(tidydr)
library(tidyverse)
library(dplyr)
library(reshape2)
library(stringr)
library(ggrepel)
library(magrittr)
library(MASS)
library(ggplot2)

#~~~ oncogenes of interest
gene <- "CDK4"

#~~~~ oncogene list
oncogene <- sort(c("EGFR", "AGAP2", "CDK4", "MDM2", "MYCN"))

#~~~ input the CNA data
fileList <- list.files("./inputs/CNA")
sample <- unlist(lapply(strsplit(fileList, split = "_"), function(x) x[1]))

#~~~ read the first file 
geneCNA <- fread(paste0("./inputs/CNA/", fileList[1])) %>% as.data.frame() 
geneCNA.split <- geneCNA %>% separate_rows(gene) %>% group_by(log2) %>% mutate(gene = trimws(gene)) 
geneCNA <- geneCNA.split %>% filter(gene %in% oncogene) %>% group_by(gene) %>% filter(log2 == max(log2)) %>% distinct() %>% mutate(log2 = 2^log2) %>% dplyr::select(gene, log2) %>% magrittr::set_colnames(c("gene", sample[1]))

for(i in 2:length(sample)){
  cat("i: ", i, ", sample: ", sample[i],  "\n")
  geneCNANew <- fread(paste0("./inputs/CNA/", fileList[i])) %>% as.data.frame() 
  geneCNA.split <- geneCNANew %>% separate_rows(gene) %>% group_by(log2) %>% mutate(gene = trimws(gene)) 
  geneCNANew <- geneCNA.split %>% filter(gene %in% oncogene) %>% group_by(gene) %>% filter(log2 == max(log2)) %>% distinct() %>% mutate(log2 = 2^log2) %>% dplyr::select(gene, log2) %>% magrittr::set_colnames(c("gene", sample[i]))
  geneCNA <- full_join(geneCNA, geneCNANew, by = "gene")
}
geneCNA <- geneCNA %>% mutate_at(vars(-gene), ~ifelse(is.na(.), 0, .))


#~~~ input the FPKM data
fileList <- list.files("./inputs/FPKM/")
id <- unlist(lapply(strsplit(fileList, split = "_"), function(x) x[2]))

geneFPKM <- fread(paste0("./inputs/FPKM/", fileList[1])) %>% as.data.frame() %>% dplyr::select(c("gene_id", "FPKM")) %>% mutate(FPKM = as.numeric(FPKM)) %>% set_colnames(c("gene", id[1]))
for(i in 2:length(id)){
  cat("i: ", i, ", id: ", id[i],  "\n")
  geneFPKMNew <- fread(paste0("./inputs/FPKM/", fileList[i])) %>% as.data.frame() %>% dplyr::select(c("gene_id", "FPKM")) %>% mutate(FPKM = as.numeric(FPKM)) %>% set_colnames(c("gene", id[i]))
  geneFPKM <- merge(geneFPKM, geneFPKMNew, by = "gene")
}
geneFPKM <- geneFPKM %>% filter(gene %in% oncogene)
dim(geneFPKM)

#~~~ find the overlapping between CNA and FPKM 
geneCNA.df <- melt(geneCNA) %>% set_colnames(c("oncogene", "sample", "CNA")) 
geneFPKM.df <- melt(geneFPKM) %>% set_colnames(c("oncogene", "sample", "FPKM"))

cnvfpkm.df <- geneCNA.df %>% full_join(geneFPKM.df, by = c("oncogene", "sample")) %>% mutate(sample = as.character(sample)) %>% mutate(FPKM = log2(FPKM)) %>% mutate(CNA = log2(CNA))

#~~~ ggplot plot ""EGFR", "AGAP2", "CDK4", "MDM2", "MYCN")

tmp.cnvfpkm.df <- cnvfpkm.df %>% filter(oncogene == "MYCN")
tmp.cnvfpkm.df <- tmp.cnvfpkm.df %>% mutate(labels = ifelse(sample %in% c("G523", "G549", "G583", "G620", "G828"), sample, "")) 
if(tmp.cnvfpkm.df$oncogene[1] == "EGFR")
  sample = c("G523", "G549","G620", "G828")
if(tmp.cnvfpkm.df$oncogene[1] == "AGAP2" | tmp.cnvfpkm.df$oncogene[1] == "CDK4")
  sample = c("G523", "G583")
if(tmp.cnvfpkm.df$oncogene[1] == "MDM2")
  sample = c("G523", "G583", "G828")
if(tmp.cnvfpkm.df$oncogene[1] == "MYCN")
  sample = c("G583")

dir.create("outputs")
pdf(file = paste0("outputs", "/", gene, "_FPKM_CNA.pdf"), height = 3.5, width = 2.5)

p <- ggplot(tmp.cnvfpkm.df, aes(x = CNA, y = FPKM, fill = oncogene)) + theme_bw() +
  geom_text_repel(aes(label = labels)) + geom_smooth(method = 'lm', linewidth = 0.5, fill = "lightblue") + 
  geom_point(data = tmp.cnvfpkm.df[-which(tmp.cnvfpkm.df$labels %in% sample), ], aes(x = CNA, y = FPKM), colour="grey", size = 1) +
  geom_point(data = tmp.cnvfpkm.df[which(tmp.cnvfpkm.df$labels %in% sample), ],  aes(x = CNA, y = FPKM), colour="#fc4e2a", size = 1) + 
  theme(legend.position = "none",panel.grid = element_blank()) +
  xlab("Gene copy number\nlog2(CNA)") + ylab("Gene expression\nlog2(FPKM)") 

p <- p + facet_grid(. ~ oncogene) 

p
 
dev.off()


