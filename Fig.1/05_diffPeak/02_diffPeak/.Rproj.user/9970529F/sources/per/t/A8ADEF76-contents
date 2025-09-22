#--- step2: visualize the fold changes ---#

rm(list = ls())

library(tidyr)
library(dplyr)
library(DESeq2)
library(Rsubread)
library(ggplot2)
library(ggpubr)
library(data.table)
"%ni%" <- Negate("%in%")

id = "G958_T"

up <- fread(paste0("./outputs.", id, "/", id, ".up.bed"), header = T, data.table = F, stringsAsFactors = F) %>% mutate(peak = paste0(chr, ":", start, "-", end))
up.within.amplicon <- fread(paste0("./outputs.", id, "/", id, ".up.within.amplicon.bed"), header = F, data.table = F, stringsAsFactors = F) %>% 
  {colnames(.) = c("chr", "start", "end");.} %>% mutate(peak = paste0(chr, ":", start, "-", end))
up.outside.amplicon <- fread(paste0("./outputs.", id, "/", id, ".up.outside.amplicon.bed"), header = F, data.table = F, stringsAsFactors = F) %>% 
  {colnames(.) = c("chr", "start", "end");.} %>% mutate(peak = paste0(chr, ":", start, "-", end))

up.within.amplicon.fc <- up %>% inner_join(up.within.amplicon, by = "peak") %>% mutate(type = "within")
up.outside.amplicon.fc <- up %>% inner_join(up.outside.amplicon, by = "peak") %>% mutate(type = "outside")

peak.matrix <- rbind(up.within.amplicon.fc, up.outside.amplicon.fc) %>% dplyr::select(log2FoldChange, type)
peak.df <- reshape2::melt(peak.matrix) %>% as.data.frame() %>% mutate(variable = as.factor(variable))

counts <- as.data.frame(table(peak.df$type))
colnames(counts) <- c("type", "n")

dir.create("outputs.foldchange")

pal = c("G958_T" = "purple", 
        "G946_T" = "orchid",
        "G933_T" = "orange",
        "A7TK_T" = "salmon", 
        "G900_T" = "yellowgreen",
        "G837_T" = "olivedrab",
        "G797_T" = "navy",
        "AA9S_T" = "darkgreen")


pdf(paste0("./outputs.foldchange/", id, ".foldchange.pdf"), width = 1.35, height = 2.35)

ggplot(peak.df, aes(x = type , y = value, fill = type)) + 
  theme_bw() +    
  geom_boxplot(width = 0.6, size = 0.35,   outlier.shape = NA) +                      
  #geom_jitter(width = 0.2, size = 0.1, alpha = 0.7, color = "black") +    
  scale_fill_manual(values = c("#CCFF00B3", "#FF0099B3")) +
  labs(title = NULL, x = paste0(id), y = "log2FC") +           
  geom_text(data = counts, aes(x = type, y = 6, label = n),         
            inherit.aes = FALSE, size = 3) +
  theme(axis.text.x = element_blank(),  
        legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  stat_compare_means(method = "wilcox.test", label = "p.format")

dev.off()


