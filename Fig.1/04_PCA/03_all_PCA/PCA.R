#--- workflow of differential peak analysis: union (i.e., consensus) peak -> featureCounts -> DEseq2
rm(list = ls())

library(tidyr)
library(dplyr)
library(DESeq2)
library(ggplot2)
library(Rsubread)
library(data.table)
"%ni%" <- Negate("%in%")

id = "G958_T"
sample = c("A7TK_T", "AA9S_T", "G797_T", "G837_T", "G900_T", "G933_T", "G946_T", "G958_T") 
sample = sample[sample %ni% id]

unifiedPeak <- fread(paste0("/Users/chupan/Documents/data_su2c/scATAC/unifiedPeak/malign.unifiedPeakSet.bed"), header = T, data.table = F) %>% 
  mutate(geneID = 1:dim(.)[1]) %>% mutate(chr = seqnames, start = start, end = end, strand = strand,  sample = GroupReplicate) %>% 
  dplyr::select(geneID, chr, start, end, strand, sample) %>% filter(chr != "chrY")

#~~~ generate peak-read matrix using featureCpunts function 
peakCount <- featureCounts(paste0("/Users/chupan/Documents/data_su2c/scATAC/malignbam/", id, ".maligncell.bam"), annot.ext = unifiedPeak,  nthreads = 4, annot.inbuilt = "hg38", isPairedEnd = T)
readCount <- data.frame(count = as.numeric(peakCount$counts)) %>% magrittr::set_colnames(paste0(id))

for(i in 1:length(sample)){
  peakCount <- featureCounts(paste0("/Users/chupan/Documents/data_su2c/scATAC/malignbam/", sample[i], ".maligncell.bam"), annot.ext = unifiedPeak,  nthreads = 4, annot.inbuilt = "hg38", isPairedEnd = T)
  readCount$tmp <- as.numeric(peakCount$counts)
  colnames(readCount)[dim(readCount)[2]] <- sample[i]
}

readCount <- readCount %>% magrittr::set_rownames(paste0(unifiedPeak$chr, ":", unifiedPeak$start, "-", unifiedPeak$end)) 

#~~~ PCA
group1 <- c("A7TK_T", "G933_T", "G946_T", "G958_T") # ecDNA-harboring group
group2 <- c("AA9S_T", "G797_T", "G837_T", "G900_T") # non-ecDNA-harboring group

count <- readCount[ ,c(group1, group2)]

#~~~ set condition
condition <- factor(ifelse(colnames(count) %in% group1, 'ecDNA', 'control'))

#~~~ dds matrix
colData <- DataFrame(condition)
rownames(colData) <- colnames(count)
dds <- DESeqDataSetFromMatrix(count, colData, design = ~ condition) 
vsd <- vst(dds, blind = F)
plotPCA(vsd)

pca <- plotPCA(vsd, intgroup = c("condition"), returnData = T) 
percentVar <- round(100 * attr(pca, "percentVar"))

custom_colors <- c(
  "G958_T" = "purple", 
  "G946_T" = "orchid",
  "G933_T" = "orange",
  "A7TK_T" = "salmon", 
  "AA9S_T" = "darkgreen",
  "G797_T" = "navy",
  "G837_T" = "olivedrab",
  "G900_T" = "yellowgreen"
)

custom_name_order <- c("AA9S_T", "G797_T", "G837_T", "G900_T", "A7TK_T", "G933_T", "G946_T", "G958_T")
pca$name <- factor(pca$name, levels = custom_name_order)

pdf(paste0("./outputs/PCA_groMalign.pdf"), width = 3.15, height = 4)

plt <- ggplot(pca, aes(x = PC1, y = PC2, color = name)) + 
  geom_point( size = 2, aes(color = name)) + theme_bw() + 
  geom_text(aes(label = name), color = "black", hjust = 1.5, vjust = 1.5) +
  scale_color_manual(values = custom_colors) +
  labs(x = paste0("PC1:", percentVar[1], "% variance"), 
       y = paste0("PC2:", percentVar[2], "% variance")) +
  theme(panel.border = element_rect(colour = "black", fill = NA), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank()) 
plt

dev.off()


