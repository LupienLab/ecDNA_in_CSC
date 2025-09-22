rm(list = ls())

library(dplyr)
library(readxl)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(magrittr)
library(data.table)

qBF = read_excel("./inputs/GBM5K_qBF.xlsx", col_names = T) %>% as.data.frame() 
qBF <- qBF %>% filter(Gene %in% c("EGFR", "MDM2", "AGAP2", "CDK4", "MYCN")) %>% dplyr::select(Gene, G523, G549, G583, G620) 

qBFdf <- reshape2::melt(qBF, id.vars = "Gene") %>% magrittr::set_colnames(c("Gene", "sample", "qBF"))
qBFdf$genelabels <- ifelse(qBFdf$sample %in% c("G523", "G583"), ifelse(qBFdf$sample %in% "G583", "G583", "G523"), "")
qBFdf$genelabrrr <- ifelse(qBFdf$sample %in% c("G549", "G620"), ifelse(qBFdf$sample %in% "G620", "G620", "G549"), "")
qBFdf$label <- paste(qBFdf$genelabels, qBFdf$genelabrrr, sep = "")
qBFdf <- qBFdf %>% mutate(color = ifelse(qBF > 0, "positive", "negative"))

#~~~
qBFdf$genelabels <- ifelse(qBFdf$sample %in% c("G523", "G583"), ifelse(qBFdf$sample %in% "G583", "G583", "G523"), "")
qBFdf$genelabrrr <- ifelse(qBFdf$sample %in% c("G549", "G620"), ifelse(qBFdf$sample %in% "G620", "G620", "G549"), "")
qBFdf$label <- paste(qBFdf$genelabels, qBFdf$genelabrrr, sep = "")
qBFdf <- qBFdf %>% mutate(color = ifelse(qBF > 0, "positive", "negative"))

pdf(paste0("./outputs/qBF_oncogene.pdf"), width = 6, height = 4) # width = 6

p <- ggplot(qBFdf, aes(x = Gene, y = qBF)) + theme_bw() +
  geom_violin( color = "blue") +
  geom_point(position = position_jitter(seed = 1, width = 0.12), data = qBFdf[-which(qBFdf$label %in% c("G523", "G583", "G549", "G620")), ], aes(x = Gene, y = qBF), colour="black") +
  geom_point(position = position_jitter(seed = 1, width = 0.12), data = qBFdf[which(qBFdf$label %in% c("G523", "G583", "G549", "G620")), ],  aes(x = Gene, y = qBF), colour="red") + 
  geom_text_repel(aes(label = label)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.25) +
  scale_color_manual(values = c("positive" = "red", "negative" = "grey")) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, colour = "black"),
        axis.text.y = element_text(color = "black"),
        panel.grid = element_blank()) + 
  xlab("") + ylab("qBF score of oncogene") 

print(p)

dev.off()



