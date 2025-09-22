#--- step3: genetic annotation of DAPs within ecDNA amplicons for individual ecDNA-harboring tumors ---#

rm(list = ls())

library(tidyr)
library(dplyr)
library(DESeq2)
library(ggplot2)
library(Rsubread)
library(data.table)
library(bedtoolsr)
library(ChIPseeker)
library(GenomicRanges)
library(clusterProfiler)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

ids = c("A7TK_T", "G933_T", "G946_T", "G958_T")

fileGrList <- lapply(ids, function(id){
  fread(paste0("./outputs.", id, "/", id, ".up.within.amplicon.bed")) %>% {colnames(.) = c("chr", "start", "end");.} %>% GRanges()
})
names(fileGrList) <- ids

peakAnnoList <- lapply(fileGrList, annotatePeak, TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb = "org.Hs.eg.db", 
                       tssRegion = c(-3000, 3000), genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"))
p0 <- plotAnnoBar(peakAnnoList) + theme(legend.position = "right") + labs(title = "") 
# Build the plot to access its data
pbuilt <- ggplot_build(p0)

# Extract fill colors used for bars
fill_colors <- unique(pbuilt$data[[1]]$fill)

annoStat <- lapply(peakAnnoList, function(df) df@annoStat)
annostat <- annoStat[[1]] %>% magrittr::set_colnames(c("Feature", names(annoStat)[1]))
for(i in 2:length(annoStat)){
  annostat <- annostat %>% full_join(annoStat[[i]], join_by(Feature))
  colnames(annostat)[dim(annostat)[2]] <- names(annoStat)[i]
}
annostat[is.na(annostat)] <- 0
class(annostat)

df <- annostat
df <- df %>% mutate(Feature = sub("\\s*\\([^)]+\\)", "", Feature))
df <- df %>% mutate(Feature = sub(".*\\b(Exon)\\b.*", "\\1", Feature)) #Exon
df <- df %>% mutate(Feature = sub(".*\\b(Intron)\\b.*", "\\1", Feature)) #Intron

for(i in 1:dim(df)[1]){
  if(df$Feature[i] == "Promoter") df$Feature[i] <- 1
  if(df$Feature[i] == "3' UTR") df$Feature[i] <- 2
  if(df$Feature[i] == "5' UTR") df$Feature[i] <- 3
  if(df$Feature[i] == "Exon") df$Feature[i] <- 4
  if(df$Feature[i] == "Intron") df$Feature[i] <- 5
  if(df$Feature[i] %in% "Distal Intergenic") df$Feature[i] <- 6
}

annostat.df <- df %>% group_by(Feature) %>% summarize(across(everything(), sum, na.rm = TRUE)) %>% ungroup() %>% as.data.frame()
annostat.df <- annostat.df %>% pivot_longer(cols = -Feature, names_to = "Sample", values_to = "Percentage") %>% as.data.frame()


dir.create("outputs.annotation")
pdf(paste0("./outputs.annotation/DAP.annotation.pdf"), width = 3.5, height = 3.25)

p <- ggplot(annostat.df, aes(x = Sample, y = Percentage, fill = as.factor(Feature))) + 
  theme_bw() +
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(values = c("1" = "#3F51B5FF", "2" = "#1f78b4", "3" = '#b2df8a', "4" = "#33a02c", "5" = "#9C27B0FF", "6" = "#e31a1c"), 
                    labels = c("Promoter", "3' UTR", "5' UTR", "Exon", "Intron", "Distal Intergenic"))  +
  labs(title = NULL, x = NULL, y = "Percentage", fill = "Feature") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.text = element_text(size = 10, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key.size = unit(0.5, "cm"),
        legend.key.width = unit(0.5, "cm"))

p

dev.off()

