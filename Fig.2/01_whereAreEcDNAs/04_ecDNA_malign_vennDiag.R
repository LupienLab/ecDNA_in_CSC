rm(list = ls())

library(dplyr)
library(ggplot2)
library(data.table)
library(Vennerable)
library(VennDiagram)
"%ni%" <- Negate("%in%")

#~~~ inputs 
id = "A7TK_T"
nth = 1
ecDNA_malign_cell_1 <- read.table(paste0(id, ".outputs/", id, "-", nth, "_ecDNA_positive_malign_cell.txt"), header = F, stringsAsFactors = F) %>% {colnames(.) = "barcode";.} %>% pull("barcode")

nth = 2
ecDNA_malign_cell_2 <- read.table(paste0(id, ".outputs/", id, "-", nth, "_ecDNA_positive_malign_cell.txt"), header = F, stringsAsFactors = F) %>% {colnames(.) = "barcode";.} %>% pull("barcode")

hit12 <- intersect(ecDNA_malign_cell_1, ecDNA_malign_cell_2)
length(hit12)

venn_plot <- venn.diagram(
  x = list(A7TK_T_1 = ecDNA_malign_cell_1, A7TK_T_2 = ecDNA_malign_cell_2),
  filename = NULL,  
  fill = "white",
  alpha = 1,
  cex = 2,
  cat.cex = 1.5
)

grid::grid.draw(venn_plot)

#~~~~~~~~~~~
id = "G933_T"

n = 1
ecDNA_malign_cell_1 <- read.table(paste0(id, ".outputs/", id, "-", n, "_ecDNA_positive_malign_cell.txt"), header = F, stringsAsFactors = F) %>% {colnames(.) = "barcode";.} %>% pull("barcode")

n = 2
ecDNA_malign_cell_2 <- read.table(paste0(id, ".outputs/", id, "-", n, "_ecDNA_positive_malign_cell.txt"), header = F, stringsAsFactors = F) %>% {colnames(.) = "barcode";.} %>% pull("barcode")

n = 3
ecDNA_malign_cell_3 <- read.table(paste0(id, ".outputs/", id, "-", n, "_ecDNA_positive_malign_cell.txt"), header = F, stringsAsFactors = F) %>% {colnames(.) = "barcode";.} %>% pull("barcode")

# hits
hit_12 <- length(intersect(ecDNA_malign_cell_1, ecDNA_malign_cell_2))
hit_13 <- length(intersect(ecDNA_malign_cell_1, ecDNA_malign_cell_3))
hit_23 <- length(intersect(ecDNA_malign_cell_2, ecDNA_malign_cell_3))
hit_123 <- length(intersect(intersect(ecDNA_malign_cell_1, ecDNA_malign_cell_2), ecDNA_malign_cell_3))

venn_plot <- venn.diagram(
  x = list(G933_T_1 = ecDNA_malign_cell_1, G933_T_2 = ecDNA_malign_cell_2, G933_T_3 = ecDNA_malign_cell_3),
  filename = NULL,  
  fill = "white",
  alpha = 0.5,
  cex = 2,
  cat.cex = 1.5
)

grid::grid.draw(venn_plot)






