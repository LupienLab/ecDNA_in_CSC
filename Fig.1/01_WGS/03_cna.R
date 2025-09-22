rm(list = ls())

library(data.table)
library(dplyr)
library(magrittr)
library(DNAcopy)
library(ggplot2)
library(scales)
library(BSgenome.Hsapiens.UCSC.hg38)

excludeChr = c("chrY", "chrM")
"%ni%" <- Negate("%in%")

id <- "G958_T" #--- sample name 
genome = BSgenome.Hsapiens.UCSC.hg38 #--- genome 
genomeSeqinfo <- keepStandardChromosomes(Seqinfo(seqnames = names(genome), seqlengths = seqlengths(genome))) #--- chromosome length

#--- read the cna profile
cnvInfo <- read.table(paste0("/Users/cpan/Documents/rawdata/cnvkit/", id, ".rmdup.cns"), header = T, sep = "\t") %>% as.data.frame() %>% mutate(ratio = 2^log2) %>% filter(chromosome %ni% excludeChr)

.getobsposition <- function(x, c1, c2, c3){
    chrIndex <- which(x[c1] == names(seqlengths(c3)))
    return(sum(seqlengths(c3)[1:chrIndex]) - as.numeric(seqlengths(c3)[chrIndex]) + as.numeric(x[c2]))
}

#--- get the absolute position of each focal amplicon
cnvInfo <- cnvInfo %>% mutate(absstart = apply(cnvInfo, 1, .getobsposition, c1 = 'chromosome', c2 = 'start', c3 = genomeSeqinfo)) %>% mutate(absend = absstart + (end - start))
cnvInfo <- cnvInfo %>% mutate(log2 = ifelse(log2 < (-4), -4, log2))

#--- segmentation for cnv profile
CNAobject <- CNA(cnvInfo$log2, cnvInfo$chromosome, cnvInfo$start, data.type = "logratio", sampleid = id)
smoothCNAobject <- smooth.CNA(CNAobject)
segsmoothCNAobject <- segment(smoothCNAobject, verbose = 1)
segmentRatio <- segsmoothCNAobject$output

#--- add segmentation information into the cnv profile
cnvInfo <- cnvInfo %>% arrange(chromosome) %>% mutate(segmean = rep(segmentRatio$seg.mean, segmentRatio$num.mark))

chrrangeStart <- cnvInfo %>% group_by(chromosome) %>% arrange(chromosome, absstart) %>% filter(row_number() == 1) %>% ungroup()   # get the start position of each chromosome
chrrangeEnd  <-  cnvInfo %>% group_by(chromosome) %>% arrange(chromosome, absstart) %>% filter(row_number() == n()) %>% ungroup() # get the end position of each chromosome

#--- create a data frame for chromosomes 
chromrect <- data.frame(chr = chrrangeStart$chromosome, xstart = as.numeric(chrrangeStart$absstart), xend = as.numeric(chrrangeEnd$absend)) %>% 
    arrange(xstart) %>% magrittr::set_rownames(1:dim(.)[1]) %>%       # sorted by xstart and set a rowname
    mutate(xbreaks = rowMeans(matrix(c(xstart, xend), ncol = 2))) %>% # get the middle position of each chromosome 
    mutate(colors = rep_len(c("white", "gray"), nrow(.)))             # define a color for each chromosome 

#--- create a rectangle
ggchrback <- list(geom_rect(data = chromrect, aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf, fill = colors), alpha = .1), scale_fill_identity())

#--- add theme for the rectangle
gtheme <- list(scale_x_continuous(breaks = chromrect$xbreaks, labels = gsub("chr", "", chromrect$chr), position = "top", expand = c(0, 0), #--- define chromosome label
    sec.axis = sec_axis(~., breaks = c(0, 0.5e9, 1e9, 1.5e9, 2e9, 2.5e9, 3e9), #--- define x scale
                        labels = c(0, 0.5, 1, 1.5, 2, 2.5, 3),                 #--- define x label
                        name = "ordered genome position (Gb)")),
    theme_classic(),
    xlab(""),
    ylab("copy number (log2 ratio)"),
    theme(axis.text.x = element_text(angle = 0, vjust = .5, size = 10), axis.text.y = element_text(size = 10), legend.position = "none",
    axis.ticks.x = element_blank(),
    axis.title = element_text(size = 10),
    plot.title = element_text(size = 10),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))
)

#--- colours
solarExtra = c('#3361A5', '#248AF3', '#14B3FF', '#88CEEF', '#C1D5DC', '#EAD397', '#FDB31A', '#E42A2A', '#A31D1D') 

#--- begin to plot ratio starting by rectangle 
p <- ggplot(cnvInfo) + ggchrback + gtheme 
#--- add point 
p <- p + geom_point(aes(absstart, log2, colour = log2), shape = 20, size = 1.5, alpha = 1) + 
     scale_colour_gradientn(colours = solarExtra) + #scale_colour_gradientn(colours = solarExtra, breaks = c(-5, -3, -2, -1, 0, 2, 4, 6, 7))
     theme(legend.position = "right", panel.grid.major = element_blank(), 
           legend.key.size = unit(0.5, 'cm'), legend.title = element_text(size = 8), legend.text = element_text(size = 8)) 
  
#--- add line 
p <- p + geom_line(aes(absstart, segmean), col = "black", linewidth = 0.35)
p <- p + geom_hline(yintercept = log2(5), linetype = "dashed", color = "#A31D1D", linewidth = 0.15)
p

ggsave(paste0("./outputs/", sample, ".wgs.cnv.ratio.pdf"), width = 24, height = 10, units = "cm")



