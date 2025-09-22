#~~~ step1: identify DAPs for each ecDNA species 

rm(list = ls())

library(tidyr)
library(dplyr)
library(DESeq2)
library(Rsubread)
library(data.table)
"%ni%" <- Negate("%in%")

id = "G958_T-1"
controls <- c("AA9S_T", "G797_T", "G837_T", "G900_T")

#~~~ load unified peaks
unifiedPeak <- fread(paste0("/Users/chupan/Documents/data_su2c/scATAC/unifiedPeak/", id, ".unifiedPeakSet.bed"), header = T) %>% 
  as.data.frame() %>% mutate(geneID = 1:dim(.)[1]) %>% mutate(chr = seqnames, start = start, end = end, strand = strand,  sample = GroupReplicate) %>% 
  dplyr::select(geneID, chr, start, end, strand, sample) %>% filter(chr != "chrY")

dir.create(paste0("outputs.", id))
write.table(unifiedPeak %>% select(chr, start, end), file = paste0("./outputs.", id, "/", id, "_unifiedPeak.bed"), row.names = F, col.names = F, sep = "\t", quote = F)

#~~~ generate a peak-by-sample read matrix
peakCount <- featureCounts(paste0("/Users/chupan/Documents/data_su2c/scATAC/ecDNA_malign/", id, "/", id, ".ecDNA_positive_malign/", id, ".ecDNA_positive_malign.bam"), annot.ext = unifiedPeak,  nthreads = 4, annot.inbuilt = "hg38", isPairedEnd = T)
readCount <- data.frame(count = as.numeric(peakCount$counts)) %>% magrittr::set_colnames(paste0(id))

for(i in 1:length(controls)){
  peakCount <- featureCounts(paste0("/Users/chupan/Documents/data_su2c/scATAC/malignbam/", controls[i], ".maligncell.bam"), annot.ext = unifiedPeak,  nthreads = 4, annot.inbuilt = "hg38", isPairedEnd = T)
  readCount$tmp <- as.numeric(peakCount$counts)
  colnames(readCount)[dim(readCount)[2]] <- controls[i]
}
readCount <- readCount %>% magrittr::set_rownames(paste0(unifiedPeak$chr, ":", unifiedPeak$start, "-", unifiedPeak$end)) 

save(readCount, file = paste0("outputs.", id, "/", id, ".readCount.RData"))

#~~~ differential accessibility analysis using DEseq2 
load(paste0("./outputs.", id, "/", id, ".readCount.RData"))

group1 <- id
group2 <- controls
count <- readCount[ ,c(group1, group2)]

#~~~ set condition
condition <- factor(ifelse(colnames(count) %in% group1, 'ecDNA', 'control'))

#~~~ peak filter
keepPeak <- rowSums(edgeR::cpm(count) > 0) > 0
count <- count[keepPeak, ]

#~~~ dds matrix
dds <- DESeqDataSetFromMatrix(count, DataFrame(condition), design = ~ condition) 
vsd <- vst(dds, blind = F)
plotPCA(vsd)

#~~~ standardization
dds2 <- DESeq(dds, parallel = T) 

#~~~ threshold for filter DEseq2 results
thres_fd <- 2
thres_pv <- 0.01

diffPeak <- results(dds2, contrast = c("condition", "ecDNA", "control"), cooksCutoff = FALSE) %>% as.data.frame() %>% 
  arrange(-log2FoldChange) %>% drop_na(padj) %>% filter(abs(log2FoldChange) > thres_fd & padj < thres_pv)  

diffPeak$change = ifelse(diffPeak$padj < thres_pv & abs(diffPeak$log2FoldChange) >= thres_fd, ifelse(diffPeak$log2FoldChange > thres_fd,'Up','Down'), 'Stable')

#~~~ up peaks 
up <- diffPeak %>% filter(change == "Up") %>% mutate(peak = rownames(.)) %>% 
  tidyr::separate(peak, sep = ":", into = c("chr","coord")) %>% tidyr::separate(coord, sep = "-", into = c("start", "end")) %>% 
  mutate(start = as.integer(start)) %>% mutate(end = as.integer(end)) %>% select(chr, start, end, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, change) %>% 
  magrittr::set_rownames(1:dim(.)[1])

#~~~ ecDNA amplicon
fileDir <- list.files(paste0("/Users/chupan/Documents/data_su2c/WGS/amplicon/"))
selectedFile <- grep(sub("_.*", "", id), fileDir, value = T)
amplicon <- read.table(paste0("/Users/chupan/Documents/data_su2c/WGS/amplicon/", grep(id, selectedFile, value = T) )) %>% magrittr::set_colnames(c("chr", "start", "end")) %>% mutate(sample = sub("_.*", "", selectedFile[1]))

#~~~ up peak within current ecDNA amplicons 
up_within_amplicon <- bedtoolsr::bt.intersect(up, amplicon, wa = T, wb = F) %>% filter(!duplicated(.)) %>% {colnames(.) = colnames(up);.}
up_outside_amplicon <- anti_join(up, up_within_amplicon, by = c("chr", "start", "end")) %>% filter(!duplicated(.))

write.table(up_within_amplicon %>% select(chr, start, end), file = paste0("./outputs.", id, "/", id, "_up_within_amplicon.bed"), row.names = F, col.names = F, sep = "\t", quote = F)
write.table(up_outside_amplicon %>% select(chr, start, end), file = paste0("./outputs.", id, "/", id, "_up_outside_amplicon.bed"), row.names = F, col.names = F, sep = "\t", quote = F)

saveRDS(up_within_amplicon, file = paste0("./outputs.", id, "/", id, "_up_within_amplicon.rds"))
saveRDS(up_outside_amplicon, file = paste0("./outputs.", id, "/", id, "_up_outside_amplicon.rds"))

#~~~ another ecDNA amplicons
selectedFile_2 <- selectedFile[selectedFile %ni% grep(id, selectedFile, value = T)]
if(length(selectedFile_2) > 0){
  amplicon_2 <- data.frame()
  for(i in 1:length(selectedFile_2)){
    tmp <- read.table(paste0("/Users/chupan/Documents/data_su2c/WGS/amplicon/", selectedFile_2[i])) %>% {colnames(.) = c("chr", "start", "end");.}
    amplicon_2 <- rbind(amplicon_2, tmp)
  }
  
  #~~~ up peak within another ecDNA amplicons
  up_within_amplicon_2 <- bedtoolsr::bt.intersect(up_outside_amplicon, amplicon_2, wa = T, wb = F) %>% filter(!duplicated(.)) %>% {colnames(.) = colnames(up);.}
  write.table(up_within_amplicon_2 %>% select(chr, start, end), file = paste0("./outputs.", id, "/", id, "_up_within_amplicon_2.bed"), row.names = F, col.names = F, sep = "\t", quote = F)
  #~~~ up peak outside any ecDNA amplicons
  up_outside_amplicon_2 <- anti_join(up_outside_amplicon, up_within_amplicon_2, by = c("chr", "start", "end")) %>% filter(!duplicated(.))
  write.table(up_outside_amplicon_2 %>% select(chr, start, end), file = paste0("./outputs.", id, "/", id, "_up_outside_amplicon_2.bed"), row.names = F, col.names = F, sep = "\t", quote = F)
  
  saveRDS(up_within_amplicon_2, file = paste0("./outputs.", id, "/", id, "_up_within_amplicon_2.rds"))
  saveRDS(up_outside_amplicon_2, file = paste0("./outputs.", id, "/", id, "_up_outside_amplicon_2.rds"))
}

