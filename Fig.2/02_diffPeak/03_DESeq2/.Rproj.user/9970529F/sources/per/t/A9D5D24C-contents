#~~~ step2: visualize the change folds 

rm(list = ls())

library(tidyr)
library(dplyr)
library(ggpubr)
library(data.table)
library(ggplot2)

id = "G958_T-1"

up_within_amplicon <- readRDS(paste0("./outputs.", id, "/", id, "_up_within_amplicon.rds")) %>% select(log2FoldChange, change) %>% mutate(change = "within")
up_outside_amplicom <- readRDS(paste0("./outputs.", id, "/", id, "_up_outside_amplicon.rds")) %>% select(log2FoldChange, change) %>% mutate(change = "without")

up_amplicon_matrix <- rbind(up_within_amplicon, up_outside_amplicom)
up_amplicon_df <- reshape2::melt(up_amplicon_matrix) %>% as.data.frame() %>% mutate(variable = as.factor(variable))

counts <- as.data.frame(table(up_amplicon_df$change))
colnames(counts) <- c("type", "n")

pdf(file = paste0("outputs.up_amplicon", "/", id, ".up_amplicon.pdf"), height = 2.0, width = 2.5)

p <- ggplot(up_amplicon_df, aes(x = change , y = value, fill = change)) + 
  theme_bw() +    
  geom_boxplot(width = 0.5, size = 0.35,   outlier.shape = NA) +                      
  scale_fill_manual(values = c("#fc4e2a", "#E6AB02")) +
  labs(title = NULL, x = NULL, y = "log2FC") +           
  geom_text(data = counts, aes(x = type, y = 6, label = n), inherit.aes = FALSE, size = 3) +
  theme(axis.text.x = element_blank(), legend.position = "right", 
        axis.text.y = element_text(colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  stat_compare_means(method = "wilcox.test", label = "p.format")
print(p)

dev.off()

#~~~~~~~~~ up amplicon_2
up_within_amplicon_2 <- readRDS(paste0("./outputs.", id, "/", id, "_up_within_amplicon_2.rds")) %>% select(log2FoldChange, change) %>% mutate(change = "within")
up_outside_amplicom_2 <- readRDS(paste0("./outputs.", id, "/", id, "_up_outside_amplicon_2.rds")) %>% select(log2FoldChange, change) %>% mutate(change = "without")

up_amplicon_matrix_2 <- rbind(up_within_amplicon_2, up_outside_amplicom_2)
up_amplicon_df_2 <- reshape2::melt(up_amplicon_matrix_2) %>% as.data.frame() %>% mutate(variable = as.factor(variable))

counts_2 <- as.data.frame(table(up_amplicon_df_2$change))
colnames(counts_2) <- c("type", "n")

pdf(file = paste0("outputs.up_amplicon_2", "/", id, ".up_amplicon_2.pdf"), height = 2.0, width = 2.5)
p <- ggplot(up_amplicon_df_2, aes(x = change , y = value, fill = change)) + 
  theme_bw() +    
  geom_boxplot(width = 0.5, size = 0.35,   outlier.shape = NA) +                      
  scale_fill_manual(values = c("green", "blue")) +
  labs(title = NULL, x = NULL, y = "log2FC") +           
  geom_text(data = counts_2, aes(x = type, y = 6, label = n), inherit.aes = FALSE, size = 3) +
  theme(axis.text.x = element_blank(), legend.position = "right", 
        axis.text.y = element_text(colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  stat_compare_means(method = "wilcox.test", label = "p.format")
print(p)

dev.off()



