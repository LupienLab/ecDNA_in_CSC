#******** step 1: comparison of ecDNA-resident oncogene activity between ecDNA(+) and ecDNA(-) cells ********#

rm(list = ls())

library(ArchR)
library(data.table)
library(ggplot2)
library(magrittr)
library(dplyr)
library(ggpubr)
library(reshape2)
library(rstatix)
Sys.info()[['sysname']]

# setting default genome to hg38
addArchRGenome("hg38")
# setting default number of parallel threads to ** 
addArchRThreads(threads = floor(parallel::detectCores())/2, force = F)

id = "G946_T"
n = 1

if(id == "A7TK_T" & n == 1) genelist <- c("AGAP2", "CDK4")
if(id == "A7TK_T" & n == 2) genelist <- c("MDM2")
if(id == "G933_T" & n == 1) genelist <- c("AGAP2", "CDK4")
if(id == "G933_T" & n == 2) genelist <- c("EGFR")
if(id == "G933_T" & n == 3) genelist <- c("MDM2")
if(id == "G946_T" & n == 1) genelist <- c("EGFR")
if(id == "G958_T" & n == 1) genelist <- c("AGAP2", "CDK4", "EGFR", "MDM2")

obj <- readRDS(paste0("./outputs.", id, "/", id, ".obj.rds"))
obj <- addImputeWeights(obj)

#~~~ get gene activity score matrix
matrixs <- getMatrixFromProject(obj)
genescoreMatrix <- matrixs@assays@data$GeneScoreMatrix 
if(length(genelist) == 1)
  genescoreMatrix <- matrixs@assays@data$GeneScoreMatrix[which(matrixs@elementMetadata$name %in% genelist), ] %>% as.matrix() %>% t() %>% set_rownames(c(rowData(matrixs)$name[which(matrixs@elementMetadata$name %in% genelist)]))
if(length(genelist) > 1)
  genescoreMatrix <- matrixs@assays@data$GeneScoreMatrix[which(matrixs@elementMetadata$name %in% genelist), ] %>% as.matrix() %>% set_rownames(c(rowData(matrixs)$name[which(matrixs@elementMetadata$name %in% genelist)]))
  
#~~~ load the ecDNA-positive/negative cells
positive_cell <- fread(paste0("/Users/chupan/Documents/data_su2c/scATAC/celltype/", id, "-", n, ".ecDNA_positive_cell.txt"), header = F, data.table = F) %>% set_colnames("barcode") %>% mutate(barcode = paste0(id, "#", barcode))
negative_cell <- fread(paste0("/Users/chupan/Documents/data_su2c/scATAC/celltype/", id, "-", n, ".ecDNA_negative_cell.txt"), header = F, data.table = F) %>% set_colnames("barcode") %>% mutate(barcode = paste0(id, "#", barcode))

#~~~ data frame for plot
if(length(genelist) == 1){
  positive_cell.numeric <- genescoreMatrix[, positive_cell$barcode] 
  positive_cell.df <- data.frame(barcode = names(positive_cell.numeric), gene = genelist, value = as.numeric(positive_cell.numeric)) %>% mutate(type = "a")
  negative_cell.numeric <- genescoreMatrix[, negative_cell$barcode] 
  negative_cell.df <- data.frame(barcode = names(negative_cell.numeric), gene = genelist, value = as.numeric(negative_cell.numeric)) %>% mutate(type = "b")
}

if(length(genelist) > 1){
  negative_cell.df <- genescoreMatrix[, negative_cell$barcode] %>% t() %>% reshape2::melt() %>% set_colnames(c("barcode", "gene", "value")) %>% mutate(type = rep("b", dim(.)[1]))
  positive_cell.df <- genescoreMatrix[, positive_cell$barcode] %>% t() %>% reshape2::melt() %>% set_colnames(c("barcode", "gene", "value")) %>% mutate(type = rep("a", dim(.)[1]))
}

cell.df <- rbind(negative_cell.df, positive_cell.df) %>% mutate(value = log(value + 1)) %>% mutate(gene = factor(gene)) 
cell.df %>% sample_n_by(type, size = 3)

stat.test <- cell.df %>% group_by(gene) %>% wilcox_test(value ~ type, alternative = "two.sided") # t_test(value ~ type, ref.group = "a") 
stat.test <- stat.test %>% add_xy_position(x = "gene", dodge = 0.8)

if(id == "A7TK_T") pal = c("a" = "salmon", "b" = "#FFA07A")
if(id == "G933_T") pal = c("a" = "orange", "b" = "#FFDAB9")
if(id == "G946_T") pal = c("a" = "orchid", "b" = "#F2D4E6")
if(id == "G958_T") pal = c("a" = "purple", "b" = "#E6E6FA")

#--- plot a box plot for gene activity score (i.e., log(copy number + 1)) #### 2.2, 2.85, 4
if((id == "A7TK_T" & n == "2") | (id == "G933_T" & n == "2" ) | (id == "G933_T" & n == "3") |(id == "G946_T" & n == "1")) {width = 2.2; height = 2.25} # 1 oncogene
if((id == "A7TK_T" & n == "1") | (id == "G933_T" & n == "1" )) {width = 2.95; height = 2.25}  # 2 oncogene
if(id == "G958_T" & n == "1"){width = 3.75; height = 2.25} # 4 oncogenes

dir.create("outputs.ecDNA")
pdf(paste0("./outputs.ecDNA/", id, "-", n, ".oncogene_activity.pdf"), width = width, height = height)

ggplot(cell.df, aes(x = gene, y = value, color = type, facet.by = "gene")) + 
  geom_boxplot(outlier.shape = NA, alpha = 1, width = 0.8) + 
  theme_bw() + 
  theme_test() + 
  scale_color_manual(values = pal, name = "Cell type", labels = c("ecDNA(+)", "ecDNA(-)")) +
  scale_x_discrete(guide = guide_axis(angle = 45)) + 
  labs(title = paste0(id, "-", n)) +
  theme(text = element_text(size = 6, colour = "black"),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        panel.grid.major = element_blank(),
        legend.key.size = unit(0.35, 'cm'),
        legend.text.align = 0) +
  xlab(" ") + ylab("Gene activity score") +
  stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, size = 1.5)

dev.off() 




