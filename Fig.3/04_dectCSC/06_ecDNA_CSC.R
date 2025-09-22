#~~~~ step6: identify ecDNA-positive CSC from GBM ~~~#

rm(list = ls())

library(grid)
library(dplyr)  
library(tidyr)
library(stringr)
library(ggplot2) 
library(data.table)
library(VennDiagram)

id =  "G933_T"
n = "3"

#~~~ filtered cells
fileDir <- list.files(paste0("/Users/chupan/Documents/data_su2c/scATAC/filtercell/"), full.names = T)
ids <- gsub(".*/(.*?)\\.filtercell\\.txt$", "\\1", fileDir)

filter_number_list <- lapply(fileDir, function(path) {
  id <- fread(path, header = F, stringsAsFactors = F) %>% as.data.frame() %>% {colnames(.) = "barcode";.}
  length(id$barcode)
})
names(filter_number_list) <- ids

id_number_df <- data.frame(sample_base = names(filter_number_list), size = unlist(filter_number_list))

#~~~ ecDNA-positive cells
fileDir <- list.files(paste0("/Users/chupan/Documents/data_su2c/scATAC/ecDNAcell/"), full.names = T)
ids <- gsub(".*/(.*?)\\.ecDNA_positive_cell\\.txt$", "\\1", fileDir)

ecDNA_cell_list <- lapply(fileDir, function(path) {
  fread(path, header = F, stringsAsFactors = F) %>% as.data.frame() %>% {colnames(.) = "barcode";.}
})
names(ecDNA_cell_list) <- ids

ecDNA_number_list <- lapply(fileDir, function(path) {
  id <- fread(path, header = F, stringsAsFactors = F) %>% as.data.frame() %>% {colnames(.) = "barcode";.}
  length(id$barcode)
})
names(ecDNA_number_list) <- ids
ecDNA_number_df <- data.frame(sample = names(ecDNA_number_list), amount = unlist(ecDNA_number_list))

#~~~ cancer stem cells (CSCs)
CSC_cell_df <- readRDS(file = paste0("outputs/overlap_cells.rds"))
CSC_cell_df <- CSC_cell_df %>% dplyr::filter(type == "T") 
CSC_number_df <- table(CSC_cell_df$sample) %>% as.data.frame() %>% {colnames(.) = c("sample", "count");.}

ecDNA_CSC_number_list <- lapply(seq_along(ecDNA_cell_list), function(i) {
  id <- sub("-\\d+$", "", names(ecDNA_cell_list[i]))
  tmp_df <- CSC_cell_df %>% filter(sample == id) %>% dplyr::select(barcode) %>% mutate(barcode = sub(paste0(id, "#"), "", barcode))
  number <- length(intersect(ecDNA_cell_list[[i]]$barcode, tmp_df$barcode))
  return(number)
})
names(ecDNA_CSC_number_list) <- names(ecDNA_cell_list)

#~~~ ecDNA-positive CSCs
ecDNA_CSC_number_df <- do.call(rbind, ecDNA_CSC_number_list) %>% as.data.frame() %>% tibble::rownames_to_column(var = "sample") %>% {colnames(.) = c("sample", "number");.}
ecDNA_CSC_number_df <- ecDNA_CSC_number_df %>% mutate(sample_base = str_extract(sample, "^[^-]+"))  

#~~~ count dataframe  
df <- left_join(ecDNA_CSC_number_df, CSC_number_df, by = c("sample_base" = "sample")) %>% mutate(percentage = number/count, diff = count - number) 
df <- left_join(df, ecDNA_number_df, by = "sample")
df <- left_join(df, id_number_df, by = "sample_base") 
df <- df %>% dplyr::select(sample, size, amount, count, number, percentage, diff) %>% filter(sample == paste0(id, "-", n))

#~~~ visualization 
total_cells <- df$size
ecDNA_positive <- df$amount
CSCs <- df$count
ecDNA_CSC_cells <- df$number  # ecDNA-positive CSC
ecDNA_negative <- total_cells - ecDNA_positive

#~~~ perform hypergeometric test
p_value <- phyper(ecDNA_CSC_cells-1, ecDNA_positive, ecDNA_negative, CSCs, lower.tail = F)
p_value

#~~~ create Venn diagram
venn.plot <- draw.pairwise.venn(
  area1 = ecDNA_positive,
  area2 = CSCs,
  cross.area = ecDNA_CSC_cells,
  category = c("", ""),
  fill = c("white", "blue"),
  alpha = 1,
  cat.col = NA,
  cat.cex = 1,
  cex = 1,
  lwd = 1
)

grid.newpage() #~~~ draw the Venn diagram with a square
pushViewport(viewport(width = 0.8, height = 0.8))  # adjust viewport size if needed
grid.rect(gp = gpar(lwd = 1, col = "black"))  #~~~ draw a square to represent total cells

#~~~ add text labels for total cells and p-value
grid.text(paste0(id, "-", n), x = 0.12, y = 1.05, gp = gpar(fontsize = 10))
grid.text(paste0(total_cells), x = 0.90, y = 0.05, gp = gpar(fontsize = 10))
grid.text(paste0("Pvalue: ", signif(p_value, 3)), x = 0.70, y = 1.05, gp = gpar(fontsize = 10, col = "black"))
grid.draw(venn.plot) #~~~ draw the Venn diagram inside the square

#~~~ generate a stacked bar
df <- df %>% mutate(csc_pos_pct = number / count, csc_neg_pct = diff / count)

#~~~ reshape for stacked bar plot
df_long <- df %>% dplyr::select(sample, csc_pos_pct, csc_neg_pct) %>% pivot_longer(cols = c(csc_pos_pct, csc_neg_pct), names_to = "type", values_to = "value")
df_long <- df_long %>% group_by(sample) %>% arrange(sample, type) %>% mutate(pos = cumsum(value) - 0.5 *value)

dir.create("outputs.CSC")
pdf(paste0("./outputs.CSC/", id, "-", n, ".stacked.pdf"), width = 1, height = 2.75)

ggplot(df_long, aes(x = sample, y = value, fill = type)) +
  theme_classic() +
  geom_bar(stat = "identity") +
  geom_text(aes(y = pos, label = scales::percent(value, accuracy = 1)), color = "red", size = 3.5) +
  labs( title = "", x = "", y = "Proportion") +
  scale_fill_manual(values = c("csc_pos_pct" = "blue", "csc_neg_pct" = "white")) + 
  geom_hline(yintercept = df_long$value[2], linetype = "solid", color = "red", linewidth = 0.5) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(colour = "black"),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank())

dev.off()
