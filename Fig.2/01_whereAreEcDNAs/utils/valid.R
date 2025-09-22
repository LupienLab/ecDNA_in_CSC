#验证granges
.validGRanges <- function(gr = NULL){
  stopifnot(!is.null(gr))
  if(inherits(gr, "GRanges")){
    return(gr)
  }else{
    stop("Error cannot validate genomic range!")
  }
}
#验证GeneAnnotation返回SimpleList 包含genes exons TSS
.validGeneAnnotation <- function(geneAnnotation = NULL){
  
  if(!inherits(geneAnnotation, "SimpleList")){
    if(inherits(geneAnnotation, "list")){
      geneAnnotation <- as(geneAnnotation, "SimpleList")
    }else{
      stop("geneAnnotation must be a list/SimpleList of 3 GRanges for : Genes GRanges, Exons GRanges and TSS GRanges!")
    }
  }
  if(identical(sort(tolower(names(geneAnnotation))), c("exons", "genes", "tss"))){
    
    gA <- SimpleList()
    gA$genes <- .validGRanges(geneAnnotation[[grep("genes", names(geneAnnotation), ignore.case = TRUE)]])
    gA$exons <- .validGRanges(geneAnnotation[[grep("exons", names(geneAnnotation), ignore.case = TRUE)]])
    gA$TSS <- .validGRanges(geneAnnotation[[grep("TSS", names(geneAnnotation), ignore.case = TRUE)]])
    
  }else{
    stop("geneAnnotation must be a list/SimpleList of 3 GRanges for : Genes GRanges, Exons GRanges and TSS GRanges!")
  }
  return(gA)
}
#验证GenomeAnnotation返回SimpleList 包含 blacklist genome chromSizes其中genome进行validBSgenome检查并导入相关R包
.validGenomeAnnotation <- function(genomeAnnotation = NULL){
  
  if(!inherits(genomeAnnotation, "SimpleList")){
    if(inherits(genomeAnnotation, "list")){
      genomeAnnotation <- as(genomeAnnotation, "SimpleList")
    }else{
      stop("genomeAnnotation must be a list/SimpleList of 3 GRanges for : blacklist GRanges, chromSizes GRanges and genome BSgenome package string (ie hg38 or BSgenome.Hsapiens.UCSC.hg38)!")
    }
  }
  
  if(identical(sort(tolower(names(genomeAnnotation))), c("blacklist", "chromsizes", "genome"))){
    
    gA <- SimpleList()
    gA$blacklist <- .validGRanges(genomeAnnotation[[grep("blacklist", names(genomeAnnotation), ignore.case = TRUE)]])
    if(genomeAnnotation[[grep("genome", names(genomeAnnotation), ignore.case = TRUE)]]=="nullGenome"){
      gA$genome <- "nullGenome"
    }else{
      bsg <- validBSgenome(genomeAnnotation[[grep("genome", names(genomeAnnotation), ignore.case = TRUE)]])
      gA$genome <- bsg@pkgname
    }
    gA$chromSizes <- .validGRanges(genomeAnnotation[[grep("chromsizes", names(genomeAnnotation), ignore.case = TRUE)]])
    
  }else{
    stop("genomeAnnotation must be a list/SimpleList of 3 GRanges for : blacklist GRanges, chromSizes GRanges and genome BSgenome package string (ie hg38 or BSgenome.Hsapiens.UCSC.hg38)!")
  }
  return(gA)
}
#如果 BSgenome 对象有效或genome能加载包返回这个BSgenome对象再取出pkgname作为genome
validBSgenome <- function(genome = NULL, masked = FALSE){
  
  stopifnot(!is.null(genome))
  if(inherits(genome, "BSgenome")){
    return(genome)
  }else if(is.character(genome)){
    genome <- tryCatch({
      .requirePackage(genome)
      bsg <- eval(parse(text = genome))
      if(inherits(bsg, "BSgenome")){
        return(bsg)
      }else{
        stop("genome is not a BSgenome valid class!")
      }
    }, error = function(x){
      BSgenome::getBSgenome(genome, masked = masked)
    })  
    return(genome)
  }else{
    stop("Cannot validate BSgenome options are a valid BSgenome or character for getBSgenome")
  }  
}
#去除掉geneAnnotation中不在genomeAnnotation$chromSizes染色体范围内的chr
.validGeneAnnoByGenomeAnno <- function(geneAnnotation, genomeAnnotation){
  
  allSeqs <- unique(paste0(seqnames(genomeAnnotation$chromSizes)))
  
  geneSeqs <- unique(paste0(seqnames(geneAnnotation$genes)))
  if(!all(geneSeqs %in% allSeqs)){
    geneNotIn <- geneSeqs[which(geneSeqs %ni% allSeqs)]
    message("Found Gene Seqnames not in GenomeAnnotation chromSizes, Removing : ", paste0(geneNotIn, collapse=","))
    geneAnnotation$genes <- .subsetSeqnamesGR(geneAnnotation$genes, names = allSeqs)
  }
  
  exonSeqs <- unique(paste0(seqnames(geneAnnotation$exons)))
  if(!all(exonSeqs %in% allSeqs)){
    exonNotIn <- exonSeqs[which(exonSeqs %ni% allSeqs)]
    message("Found Exon Seqnames not in GenomeAnnotation chromSizes, Removing : ", paste0(exonNotIn, collapse=","))
    geneAnnotation$exons <- .subsetSeqnamesGR(geneAnnotation$exons, names = allSeqs)
  }
  
  TSSSeqs <- unique(paste0(seqnames(geneAnnotation$TSS)))
  if(!all(TSSSeqs %in% allSeqs)){
    TSSNotIn <- TSSSeqs[which(TSSSeqs %ni% allSeqs)]
    message("Found TSS Seqnames not in GenomeAnnotation chromSizes, Removing : ", paste0(TSSNotIn, collapse=","))
    geneAnnotation$TSS <- .subsetSeqnamesGR(geneAnnotation$TSS, names = allSeqs)
  }
  
  return(geneAnnotation)
  
}
#去除不存在chr的函数
.subsetSeqnamesGR <- function(gr = NULL, names = NULL){
  gr <- gr[which(as.character(seqnames(gr)) %in% names),]
  seqlevels(gr) <- as.character(unique(seqnames(gr)))
  return(gr)
}
