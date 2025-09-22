#--- step2: identify differential accessible peaks for individual ecDNA-harboring tumors ---#

rm(list = ls())

library(tidyr)
library(dplyr)
library(DESeq2)
library(Rsubread)
library(data.table)
"%ni%" <- Negate("%in%")

id = "G946_T"
controls <- c("AA9S_T", "G797_T", "G837_T", "G900_T")

unifiedPeak <- fread(paste0("/Users/chupan/Documents/data_su2c/scATAC/unifiedPeak/", id, ".unifiedPeakSet.bed"), header = T) %>% 
  mutate(geneID = 1:dim(.)[1]) %>% mutate(chr = seqnames, start = start, end = end, strand = strand,  sample = GroupReplicate) %>% 
  dplyr::select(geneID, chr, start, end, strand, sample) %>% filter(chr != "chrY")

#~~~ generate peak-read matrix 
peakCount <- featureCounts(paste0("/Users/chupan/Documents/data_su2c/scATAC/malignbam/", id, ".maligncell.bam"), annot.ext = unifiedPeak,  nthreads = 4, annot.inbuilt = "hg38", isPairedEnd = T)
readCount <- data.frame(count = as.numeric(peakCount$counts)) %>% magrittr::set_colnames(paste0(id))

for(i in 1:length(controls)){
  peakCount <- featureCounts(paste0("/Users/chupan/Documents/data_su2c/scATAC/malignbam/", controls[i], ".maligncell.bam"), annot.ext = unifiedPeak,  nthreads = 4, annot.inbuilt = "hg38", isPairedEnd = T)
  readCount$tmp <- as.numeric(peakCount$counts)
  colnames(readCount)[dim(readCount)[2]] <- controls[i]
}

readCount <- readCount %>% magrittr::set_rownames(paste0(unifiedPeak$chr, ":", unifiedPeak$start, "-", unifiedPeak$end)) 

dir.create(paste0("outputs.", id))
save(readCount, file = paste0("outputs.", id, "/", id, ".readCount.RData"))

#~~~ peak differential analysis using DEseq2 
load(paste0("./outputs.", id, "/", id, ".readCount.RData"))

group1 <- id
group2 <- controls

count <- readCount[ ,c(group1, group2)]
dim(count)

#~~ set condition
condition <- factor(ifelse(colnames(count) %in% group1, 'ecDNA', 'control'))

#--- peak filter
keepPeak <- rowSums(edgeR::cpm(count) > 0) > 0
count <- count[keepPeak, ]
dim(count) 

#--- dds matrix
dds <- DESeqDataSetFromMatrix(count, DataFrame(condition), design = ~ condition) 
vsd <- vst(dds, blind = F)
plotPCA(vsd)
head(dds) 

#~~~ standardization
dds2 <- DESeq(dds, parallel = T) 

#~~~ DEseq2 result with removal of NA 
fd.thres <- 2
pv.thres <- 0.01
diffPeak <- results(dds2, contrast = c("condition", "ecDNA", "control"), cooksCutoff = FALSE) %>% as.data.frame() %>% 
arrange(-log2FoldChange) %>% drop_na(padj) %>% filter(abs(log2FoldChange) > fd.thres & padj < pv.thres)  # select the peaks with significant difference 
diffPeak$change = ifelse(diffPeak$padj < pv.thres & abs(diffPeak$log2FoldChange) >= fd.thres, ifelse(diffPeak$log2FoldChange > fd.thres,'Up','Down'), 'Stable')

#~~~ up peaks 
up <- diffPeak %>% filter(change == "Up") %>% mutate(peak = rownames(.)) %>% tidyr::separate(peak, sep = ":", into = c("chr","coord")) %>% tidyr::separate(coord, sep = "-", into = c("start", "end")) %>% 
  mutate(start = as.integer(start)) %>% mutate(end = as.integer(end)) %>% 
  select(chr, start, end, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, change) %>% 
  magrittr::set_rownames(1:dim(.)[1])

#~~~ ecDNA amplicon
fileDir <- list.files(paste0("/Users/chupan/Documents/data_su2c/WGS/AA/"))
selectedFile <- grep(id, fileDir, value = T)
amplicon <- fread(paste0("/Users/chupan/Documents/data_su2c/WGS/AA/", selectedFile[1])) %>% as.data.frame() %>% magrittr::set_colnames(c("chr", "start", "end")) %>% mutate(sample = sub("_.*", "", selectedFile[1]))
if(length(selectedFile) > 1){
  for(i in 2:length(selectedFile)){
    tmpamplicon <- fread(paste0("/Users/chupan/Documents/data_su2c/WGS/AA/", selectedFile[i])) %>% as.data.frame() %>% magrittr::set_colnames(c("chr", "start", "end")) %>% mutate(sample = sub("_.*", "", selectedFile[i]))
    amplicon <- rbind(amplicon, tmpamplicon)
  }
}

up.within.amplicon <- bedtoolsr::bt.intersect(up %>% select(chr, start, end), amplicon, wa = T, wb = F) %>% filter(!duplicated(.)) %>% magrittr::set_colnames(c("chr", "start", "end"))
up.outside.amplicon <- anti_join(up, up.within.amplicon, by = c("chr", "start", "end")) %>% select(chr, start, end)

write.table(up, file = paste0("./outputs.", id, "/", id, ".up.bed"), row.names = F, col.names = T, sep = "\t", quote = F)
write.table(up.within.amplicon, file = paste0("./outputs.", id, "/", id, ".up.within.amplicon.bed"), row.names = F, col.names = F, sep = "\t", quote = F)
write.table(up.outside.amplicon, file = paste0("./outputs.", id, "/", id, ".up.outside.amplicon.bed"), row.names = F, col.names = F, sep = "\t", quote = F)

