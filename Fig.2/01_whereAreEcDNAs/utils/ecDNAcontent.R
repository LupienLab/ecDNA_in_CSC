
identecDNAcell <- function(
    obj             = NULL,
    pathToAmplicon  = NULL,
    nth             = NULL
){
  #sample = "G946_T"
  #obj <- readRDS(paste0(sample, ".outputs", "/", sample, ".obj.rds"))
  whereAreEcDNAsObj <- obj
  sampleName <- unique(whereAreEcDNAsObj$sample)
  cellNames <- sub(".*#", "", whereAreEcDNAsObj$cellNames)
  peakScoreMat <- getPeakScoreMat(obj = whereAreEcDNAsObj) # calculate a peak intensity matrix
  
  peakGr <- GRanges(
    seqnames = sub(":.*", "", rownames(peakScoreMat)),
    ranges = IRanges(
      start = as.integer(sub(".*:(\\d+)-\\d+", "\\1", rownames(peakScoreMat))),
      end   = as.integer(sub(".*:(\\d+)-(\\d+)", "\\2", rownames(peakScoreMat)))
    ))
  
  #~~~~ define the outputDir 
  outputDir <- file.path(paste0(sampleName,".outputs"))
  dir.create(outputDir, recursive = TRUE)
  
  #~~~ check the path to amplicon 
  if(is.null(pathToAmplicon)){
    stop("Please enter the path for pathToAmplicon.")
  }
  
  #~~~ read the amplicon data and convert into GRanges format
  amplicon <-  fread(pathToAmplicon) %>% as.data.frame() %>% magrittr::set_colnames(c("chr", "start", "end"))
  ampliconGr <- makeGRangesFromDataFrame(amplicon, keep.extra.columns = TRUE) %>% {names(.) <- 1:length(.); .}
  
  #~~~ peak hits within ecDNA amplicon
  hits <- queryHits(findOverlaps(peakGr, ampliconGr))
  subpeakScoreMat <- peakScoreMat[hits, ] 
  
  #~~~ calculate the U score based on Mann-Whitney U statistic
  uMatrixRank <- matrix(rank(-subpeakScoreMat, ties.method = 'min'), nrow = dim(subpeakScoreMat)[1]) %>% {rownames(.) <- rownames(subpeakScoreMat); .} %>% {colnames(.) <- colnames(subpeakScoreMat); .}
  uStatis <- apply(uMatrixRank, 2, sum) - (nrow(uMatrixRank) * (nrow(uMatrixRank) + 1))/2
  uScore <- 1 - uStatis / (as.numeric(nrow(subpeakScoreMat)) * max(uMatrixRank)) 
  uScore <- sort(uScore)
  
  #~~~ compute density 
  uDensity <- density(uScore)
  dens_df  <- tibble(x = uDensity$x, y = uDensity$y)
  
  #~~~ find inflection point
  delta <- diff(dens_df$y)
  candidates <- which(head(delta, -1) * tail(delta, -1) < 0)
  pick <- candidates[which(delta[candidates] - delta[candidates + 1] < -0.001)][1]
  
  inflection_x <- dens_df$x[pick]
  inflection_y <- dens_df$y[pick]
  
  pdf(paste0("./", sampleName, ".outputs/", sampleName, "_ecDNA_content_density.pdf"), width = 5, height = 3)
  
  p <- ggplot(dens_df, aes(x, y)) + 
    theme_bw() +
    geom_line(linewidth = 0.5) +
    geom_point(aes(x = inflection_x, y = inflection_y), size = 2, color = "blue") +
    annotate("text", x = inflection_x, y = inflection_y, 
             label = sprintf("%.3f", inflection_x),
             vjust = 2, size = 3) +
    labs(title = paste0(sampleName, " ecDNA content density distribution"), x = "ecDNA content", y = "Density") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.position = "none",
          axis.text.x = element_text(size = 10, colour = "black"), 
          axis.text.y = element_text(size = 10, colour = "black"))
   print(p)
  
   dev.off()
  
   #~~~ ecDNA-positive vs. ecDNA-negative cells
   ecDNA_positive <- sort(uScore[which(uScore > inflection_x)], decreasing = T)
   write.table(names(ecDNA_positive), file = paste0(sample, ".outputs", "/", sample, "-", nth, "_ecDNA_positive_cell.txt"), row.names = F, col.names = F, sep = "\t", quote = F)
   
   ecDNA_negative <- sort(uScore[which(uScore <= inflection_x)], decreasing = T)
   write.table(names(ecDNA_negative), file = paste0(sample, ".outputs", "/", sample, "-", nth, "_ecDNA_negative_cell.txt"), row.names = F, col.names = F, sep = "\t", quote = F)
   
   #~~~ ecDNA content score
   write.table(data.frame(barcode = names(uScore), uscore = as.numeric(uScore)), file = paste0(sample, ".outputs", "/", sample, "-", nth, "_ecDNA_content.txt"), row.names = F, col.names = F, sep = "\t", quote = F)
   
   #~~~ save results 
   uScore_df <- data.frame(
     cell   = names(uScore),
     uScore = as.numeric(uScore),
     stringsAsFactors = FALSE
   )
   uScore_df <- uScore_df[match(sub(".*#", "", whereAreEcDNAsObj$cellNames), uScore_df$cell), ]
   uScore_df <- uScore_df %>% mutate(celltype = ifelse(uScore > inflection_x, "ecDNA+", "ecDNA-"))
  
   whereAreEcDNAsObj$ecDNAcontent <- uScore_df$uScore
   whereAreEcDNAsObj$celltype <- uScore_df$celltype
   
   return(whereAreEcDNAsObj)
}

