rm(list = ls())

library(dplyr)
library(readxl)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(magrittr)
library(data.table)

qBF = read_excel("./inputs/GBM5K_qBF.xlsx", col_names = TRUE) %>% as.data.frame() 
SOX <- qBF$Gene[grep("^SOX", qBF$Gene)]
POU <- qBF$Gene[grep("^POU", qBF$Gene)]
bHLH <- c(qBF$Gene[grep("^OLIG", qBF$Gene)], qBF$Gene[(grep("^ASCL", qBF$Gene))])
tfs <- c(SOX, POU, bHLH)

qBF <- qBF %>% filter(Gene %in% tfs) 


#---- highlight the G523, G583, G649 and G620
qBF <- qBF %>% dplyr::select(Gene, G523, G583, G549, G620)
qBFdf <- reshape2::melt(qBF, id.vars = "Gene") %>% magrittr::set_colnames(c("Gene", "cellline", "essentiality"))
qBFdf$TFlabels <- ifelse(qBFdf$essentiality > 0, qBFdf$Gene, "")

pdf(paste0("./outputs/qBF_TF.pdf"), width = 6, height = 4) # width = 6

p1 <- ggplot(qBFdf, aes(x = cellline, y = essentiality, fill = cellline)) + theme_bw() +
  geom_violin(color = "blue") +
  geom_point(position = position_jitter(seed = 1, width = 0.1), data = qBFdf[which(qBFdf$essentiality < 0), ], aes(x = cellline, y = essentiality), colour="black") +
  geom_point(position = position_jitter(seed = 1, width = 0.1), data = qBFdf[which(qBFdf$essentiality > 0), ],  aes(x = cellline, y = essentiality), colour="red") + 
  geom_text_repel(aes(label = TFlabels)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, colour = "black"),
        axis.text.y = element_text(color = "black"),
        panel.grid = element_blank()) + 
  xlab("") + ylab("qBF score of TFs") 

p1

dev.off()



