#--- workflow of differential peak analysis: union (i.e., consensus) peak -> featureCounts -> DEseq2
rm(list = ls())

library(tidyr)
library(dplyr)
library(DESeq2)
library(Rsubread)
library(data.table)
"%ni%" <- Negate("%in%")

id = "G523" # sample id
amplicon = "_L_amplicon1_ecDNA_1_"

#--- inputs
ecDNA.interval <- fread(paste0("/Users/cpan/Documents/rawdata/amplicon/", id, "_L_amplicon1_ecDNA_1_intervals.bed")) %>% magrittr::set_colnames(c("chr", "start", "end"))
idPeak <- fread(paste0("/Users/cpan/Documents/rawdata/bulkATAC/narrowPeak/", id, "_peaks.narrowPeak"), header = F) %>% dplyr::select(V1, V2, V3) %>% magrittr::set_colnames(c("chr", "start", "end")) %>% filter(chr != "chrY")
ecDNA.interval.peak <- bedtoolsr::bt.intersect(idPeak, ecDNA.interval, wa = T, wb = T) %>% select(V1, V2, V3) %>% magrittr::set_colnames(c("chr", "start", "end"))

unionPeak <- fread(paste0("/Users/cpan/Documents/rawdata/bulkATAC/unionPeak/", id, ".unionPeakSet.bed"), header = T) %>% 
mutate(geneID = 1:dim(.)[1]) %>% mutate(chr = seqnames, start = start, end = end, strand = strand,  sample = GroupReplicate) %>% 
dplyr::select(geneID, chr, start, end, strand, sample) %>% filter(chr != "chrY")
union.ecDNA.peak <- bedtoolsr::bt.intersect(unionPeak %>% select(chr, start, end), ecDNA.interval, wa = T, wb = T) %>% 
magrittr::set_colnames(c("chr", "start", "end", "chr_ecDNA", "start_ecDNA", "end_ecDNA"))

#---************************** generate peak-read matrix using featureCpunts function *************************---#
peakCount <- featureCounts(paste0("/Users/cpan/Documents/rawdata/bulkATAC/bam/", id, ".rmdup.bam"), annot.ext = unionPeak,  nthreads = 4, annot.inbuilt = "hg38", isPairedEnd = F)
readCount <- data.frame(count = as.numeric(peakCount$counts)) %>% magrittr::set_colnames(paste0(id))

controls = sort(unique(unionPeak$sample)[unique(unionPeak$sample) %ni% id])
for(i in 1:length(controls)){
  peakCount <- featureCounts(paste0("/Users/cpan/Documents/rawdata/bulkATAC/bam/", controls[i], ".rmdup.bam"), annot.ext = unionPeak,  nthreads = 4, annot.inbuilt = "hg38", isPairedEnd = F)
  readCount$tmp <- as.numeric(peakCount$counts)
  colnames(readCount)[dim(readCount)[2]] <- controls[i]
}

readCount <- readCount %>% magrittr::set_rownames(paste0(unionPeak$chr, ":", unionPeak$start, "-", unionPeak$end))

#---******************************** peak differential analysis using DEseq2 **********************************---#
group1 <- id
group2 <- controls

readCount <- readCount[ ,c(group1, group2)]

#--- set condition
condition <- factor(ifelse(colnames(readCount) %in% group1, 'ecDNA', 'control'))
colData <- model.matrix(~-1 + condition) %>% magrittr::set_colnames(levels(condition)) %>% magrittr::set_rownames(colnames(readCount))

#--- peak filter
keepPeak <- rowSums(edgeR::cpm(readCount) > 0) > 0
readCount <- readCount[keepPeak, ]
dim(readCount) 

#--- dds matrix
dds <- DESeqDataSetFromMatrix(countData = readCount, colData = DataFrame(condition), design = ~ condition) 
dds <- dds[rowSums(counts(dds)) > 10, ]
head(dds) 

#--- standardization
dds2 <- DESeq(dds, parallel = T) 

#--- DEseq2 result with removal of NA 
fd.thres <- 1.5
pv.thres <- 0.05
diffPeak <- results(dds2, contrast = c("condition", "ecDNA", "control"), cooksCutoff = FALSE) %>% as.data.frame() %>% 
arrange(-log2FoldChange) %>% drop_na(padj) %>% filter(abs(log2FoldChange) > fd.thres & padj < pv.thres)  # select the peaks with significant difference 
diffPeak$change = ifelse(diffPeak$padj < pv.thres & abs(diffPeak$log2FoldChange) >= fd.thres, ifelse(diffPeak$log2FoldChange > fd.thres,'Up','Down'), 'Stable')
write.table(diffPeak, file = paste0("./outputs/", id, "/", id, "_L.diffPeak.overall.bed"), row.names = T, col.names = T, quote = F)
dim(diffPeak)

#--- differential peaks 
peaks <- diffPeak %>% mutate(peak = rownames(.)) %>% tidyr::separate(peak, sep = ":", into = c("chr","coord")) %>% tidyr::separate(coord, sep = "-", into = c("start", "end")) %>% 
mutate(start = as.integer(start)) %>% mutate(end = as.integer(end)) %>% 
select(chr, start, end, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, change) %>% magrittr::set_rownames(1:dim(.)[1])

#--- differential peaks specific for ecDNA
diffPeak.ecDNA <- bedtoolsr::bt.intersect(peaks, ecDNA.interval.peak,  wa = T, wb = T) %>% 
magrittr::set_colnames(c("chr", "start", "end", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "change", "chr_ecDNA", "start_ecDNA", "end_ecDNA")) %>%
arrange(-log2FoldChange) %>% filter(change == "Up") %>% dplyr::select(chr, start, end) %>% 
magrittr::set_colnames(c("chr", "start", "end")) %>% filter(!duplicated(.))
write.table(diffPeak.ecDNA, file = paste0("./outputs/", id, "/", id, "_L.diffPeak.ecDNA.bed"), row.names = F, col.names = F, quote = F, sep = "\t")

dim(unionPeak)
unionPeak.ecDNA.bg <- anti_join(unionPeak, union.ecDNA.peak, by = c("chr", "start", "end"))
dim(unionPeak.ecDNA.bg)
write.table(unionPeak.ecDNA.bg %>% select(chr, start, end) , file = paste0("./outputs/", id, "/", id, "_L.unionPeak.ecDNA.bg.bed"), row.names = F, col.names = F, quote = F, sep = "\t")

#--- plot 
slices <- c(dim(diffPeak.ecDNA)[1], (length(which(diffPeak$change == "Up")) - dim(diffPeak.ecDNA)[1]))
lbls <- c("peaks within ecDNA", "peaks without ecDNA")
df <- data.frame(category = lbls, count = slices) %>% mutate(percentage = count/sum(count) *100)
pie(slices, labels = paste0(slices, "(", round(slices/sum(slices), 4)*100, "%)"), col = c("#fee0d2", "white"))
library(ggplot2)
ggplot(df, aes(x = 2, y = count, fill = category)) +
  geom_bar(stat = "identity", width = 0.75, color = "white") +
  coord_polar(theta = "y", start = pi/0.85) +
  xlim(0.5, 2.5) + # Adjust the inner radius for the donut chart
  scale_fill_manual(values = c("peaks within ecDNA" = "#d94801", "peaks without ecDNA" = "grey")) +
  theme_void() + # Remove background and axis
  theme(legend.position = "right") + # Adjust legend position
  geom_text(aes(label = paste0(count, "\n", round(percentage, 1), "%")), 
            position = position_stack(vjust = 0.5)) + # Add percentage labels
  theme(plot.title = element_text(hjust = 0.5))

#--- union peaks specific for chrom
diffPeak.chrom <- peaks[which(peaks$change == "Down"), ] %>% select(chr, start, end) %>% filter(!duplicated(.))
write.table(diffPeak.chrom %>% select(chr, start, end) , file = paste0("./outputs/", id, "/", id, "_L.diffPeak.chrom.bed"), row.names = F, col.names = F, quote = F, sep = "\t")

unionPeak.chrom.bg <- anti_join(unionPeak, diffPeak.chrom, by = c("chr", "start", "end"))
write.table(unionPeak.ecDNA.bg %>% select(chr, start, end) , file = paste0("./outputs/", id, "/", id, "_L.unionPeak.chrom.bg.bed"), row.names = F, col.names = F, quote = F, sep = "\t")

