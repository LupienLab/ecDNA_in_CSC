rm(list = ls())

library(tidyr)
library(dplyr)
library(ggplot2)
library(Rsubread)
library(data.table)
library(bedtoolsr)
library(ChIPseeker)
library(GenomicRanges)
library(clusterProfiler)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
"%ni%" <- Negate("%in%")

id = "G523_L"
control <- c("G564_L|G566_L|G567_L|G789_L|G797_L")

interval <- fread(paste0("/Users/chupan/Documents/data_su2c/WGS/amplicon/", id, "-1_ecDNA_1_intervals.bed"), header = F) %>% as.data.frame() %>% magrittr::set_colnames(c("chr", "start", "end"))
narrowPeak <- fread(paste0("/Users/chupan/Documents/data_su2c/ATAC/narrowPeak/", id, ".peaks.narrowPeak"), header = F) %>% as.data.frame()
peak_interval <- bedtoolsr::bt.intersect(interval, narrowPeak,   wo =T) %>% dplyr::select(c(V4, V5, V6)) %>% magrittr::set_colnames(c("chr", "start", "end")) %>% GRanges()

fileGrList <- list(peak_interval)

fileDir <- list.files(paste0("/Users/chupan/Documents/data_su2c/ATAC/narrowPeak"))
fileDir <- fileDir[grep(control, fileDir)]
sampleName <- unlist(lapply(strsplit(fileDir, split = ".p"), function(x) x[1]))

for(i in 1:length(sampleName)){
  narrowPeak <- fread(paste0("/Users/chupan/Documents/data_su2c/ATAC/narrowPeak/", sampleName[i], ".peaks.narrowPeak"), header = F) %>% as.data.frame() 
  peak_interval <- bedtoolsr::bt.intersect(interval, narrowPeak,   wo =T) %>% dplyr::select(c(V4, V5, V6)) %>% magrittr::set_colnames(c("chr", "start", "end")) %>% GRanges()
  fileGrList <- c(fileGrList, peak_interval)
}

fileList <- setNames(fileGrList, c(id, sampleName))

peakAnnoList <- lapply(fileList, annotatePeak, TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb = "org.Hs.eg.db", 
tssRegion = c(-3000, 3000), genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"))
plotAnnoBar(peakAnnoList) + theme(legend.position = "none") + labs(title = "") 


annoStat <- lapply(peakAnnoList, function(df) df@annoStat)
annostat <- annoStat[[1]] %>% magrittr::set_colnames(c("Feature", names(annoStat)[1]))
for(i in 2:length(annoStat)){
  annostat <- annostat %>% full_join(annoStat[[i]], join_by(Feature))
  colnames(annostat)[dim(annostat)[2]] <- names(annoStat)[i]
}
annostat[is.na(annostat)] <- 0

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

annostat.df <- tidyr::gather(df, key = "sample", value = "percentage", -Feature)

pdf(paste0("./outputs/", id, "_peakAnn.pdf"), width = 3.35, height = 2.75) # width = 6

ggplot(annostat.df, aes(x = sample, y = percentage, fill = Feature)) + theme_bw() + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = c("1" = "#3F51B5FF", "2" = "#1f78b4", "3" = '#b2df8a', "4" = "#33a02c", "5" = "#9C27B0FF", "6" = "#e31a1c"), 
                    labels = c("Promoter", "3' UTR", "Exon", "Intron", "Distal Intergenic")) +
  labs(x = NULL, y = "Percentage (%)") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, colour = "black"), 
        axis.text.y = element_text(colour = "black"),
        legend.text = element_text(size = 10, colour = "black"), 
        legend.key.size = unit(0.35, "cm"),
        legend.key.width = unit(0.35,"cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

dev.off()

