#~~~~ step5: count CSC from GBM ~~~#

rm(list = ls())

library(MASS)    
library(ggplot2)
library(dplyr)  
set.seed(1)

# load overlap data
SVD_df <- readRDS(file = paste0("outputs/integrate_SVD.rds"))
overlap_df <- readRDS(file = paste0("outputs/overlap_cells.rds"))

SVD_df <- table(SVD_df$type, SVD_df$sample) %>% as.data.frame() %>% {colnames(.) = c("type", "sample", "overallcount");. }
overlap_df <- table(overlap_df$type, overlap_df$sample) %>% as.data.frame() %>% {colnames(.) = c("type", "sample", "overlapcount");. }

count <- full_join(SVD_df, overlap_df, by = c("type", "sample")) %>% mutate(percentage = overlapcount/overallcount) %>% filter(!is.na(percentage))
count_GBM <- count %>% filter(type == "T")
count_GSC <- count %>% filter(type == "L")
count_GBM$sample <- factor(count_GBM$sample, levels = c("A7TK_T", "G933_T", "G946_T", "G958_T", "AA9S_T", "G797_T", "G837_T", "G900_T"))
count_GSC$sample <- factor(count_GSC$sample, levels = c("G523_L", "G583_L", "G620_L", "G797_L", "G828_L", "G837_L"))

pdf(paste0("./outputs.SVD/count_GBM.pdf"), width = 3.5, height = 3)

ggplot(count_GBM, aes(x = sample, y = percentage, fill = "custom")) +
  theme_bw() +
  geom_bar(stat = "identity", color = "black", width = 0.9) +  
  geom_text(aes(label = paste0(overlapcount, "(", round(percentage *100, 2), "%)")), vjust = 0.35, size = 4, angle = 90, hjust = -0.1) +  
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_manual(values = c("custom" = "blue")) +
  labs(x = NULL, y = "Percentage", fill = "Sample") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position = "none",
        panel.grid.minor = element_blank(), panel.grid.major = element_blank())  

dev.off()

pdf(paste0("./outputs.SVD/count_GSC.pdf"), width = 3.5, height = 3)

ggplot(count_GSC, aes(x = sample, y = percentage, fill = "custom")) +
  theme_bw() +
  geom_bar(stat = "identity", color = "black", width = 0.9) +  
  geom_text(aes(label = paste0(overlapcount, "(", round(percentage *100, 2), "%)")), vjust = 0.35, size = 4, angle = 90, hjust = -0.1) +  
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_manual(values = c("custom" = "blue")) +
  labs(x = NULL, y = "Percentage", fill = "Sample") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position = "none",
        panel.grid.minor = element_blank(), panel.grid.major = element_blank())  

dev.off()

