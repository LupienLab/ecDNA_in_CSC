
#添加annotation
addWhereAreEcDNAsGenome <- function(genome = NULL){
  
  supportedGenomes <- c("hg19","hg38","mm9","mm10","hg19test")
  
  if(tolower(genome) %ni% supportedGenomes){
    
    message("Genome : ", genome, " is not currently supported by deEcDNA.")
    message("Currently supported genomes : ", paste0(supportedGenomes, collapse = ","))
    
  }else{
    
    genome <- paste(toupper(substr(genome, 1, 1)), substr(genome, 2, nchar(genome)), sep="")
    message("Setting default genome to ", genome, ".")
    options(WhereAreEcDNAs.genome = genome)
    
  }
  invisible(0)
}
#获取添加的注释
getWhereAreEcDNAsGenome <- function(geneAnnotation=FALSE,genomeAnnotation=FALSE){
  supportedGenomes <- c("hg19","hg38","mm9","mm10","hg19test")
  .WhereAreEcDNAsGenome <- options()[["WhereAreEcDNAs.genome"]]
  
  if(!is.null(.WhereAreEcDNAsGenome)){
    ag <- .WhereAreEcDNAsGenome
    if(!is.character(ag)){
      return(NULL)
    }else{
      
      if(tolower(ag) %in% supportedGenomes){
        genome <- paste(toupper(substr(ag, 1, 1)), substr(ag, 2, nchar(ag)), sep="")
        if(geneAnnotation & genomeAnnotation){
          stop("Please request either geneAnnotation or genomeAnnotation, not both!")
        }
        if(geneAnnotation){
          message("Using GeneAnnotation set by addWhereAreEcDNAsGene(",ag,")!")
          geneAnno <- paste0("geneAnno", genome)
          eval(parse(text=paste0("data(geneAnno",genome,")")))
          return(eval(parse(text=geneAnno)))
          
        }else if(genomeAnnotation){
          message("Using GenomeAnnotation set by addWhereAreEcDNAsGenome(",ag,")!")
          genomeAnno <- paste0("genomeAnno", genome)
          eval(parse(text=paste0("data(genomeAnno",genome,")")))
          return(eval(parse(text=genomeAnno)))
          
        }else{
          return(genome)
        }
        
      }else{
        stop("option(WhereAreEcDNAs.genome) : ", ag, " This genomic information is not supported.")
      }
    }
    
  }else{
    return(NULL)
  }
  
}
#返回Gene注释
getGeneAnnotation <- function(WhereAreEcDNAsProj = NULL){
  if(is.null(WhereAreEcDNAsProj)){
    geneAnnotation <- getWhereAreEcDNAsGenome(geneAnnotation = TRUE) 
    if(!is.null(geneAnnotation)){
      return(geneAnnotation)
    }
    stop("getGeneAnnotation : WhereAreEcDNAsPRoj is NULL and there is no genome set with addWhereAreEcDNAsGenome!")
  }
  WhereAreEcDNAsProj@geneAnnotation
}
#返回Genome注释
getGenomeAnnotation <- function(WhereAreEcDNAsProj = NULL){
  if(is.null(WhereAreEcDNAsProj)){
    genomeAnnotation <- getWhereAreEcDNAsGenome(genomeAnnotation = TRUE)
    if(!is.null(genomeAnnotation)){
      return(genomeAnnotation)
    }
    stop("getGenomeAnnotation : WhereAreEcDNAsPRoj is NULL and there is no genome set with addWhereAreEcDNAsGenome!")
  }
  return(WhereAreEcDNAsProj@genomeAnnotation)
}
#返回一个包含空的基因、外显子和转录起始位点（TSS）的列表。创建一个空的基因注释
.nullGeneAnnotation <- function(){
  genes <- GRanges("chr1", IRanges(1,1), symbol = "a")
  genes <- genes[-1]
  exons <- genes
  TSS <- genes
  SimpleList(genes = genes, exons = exons, TSS = TSS)
}
#空的基因组注释
.nullGenomeAnnotation <- function(){
  genome <- "nullGenome"
  chromSizes <- GRanges()
  blacklist <- GRanges()
  SimpleList(blacklist = blacklist, genome = genome, chromSizes = chromSizes)
}
#添加annotation
addWhereAreEcDNAsCancer <- function(cancer = NULL){
  
  supportedcancers <- c("gbm","brca")
  
  if(tolower(cancer) %ni% supportedcancers){
    
    message("cancer : ", cancer, " is not currently supported by deEcDNA.")
    message("Currently supported cancers : ", paste0(supportedcancers, collapse = ","))
    
  }else{
    
    cancer <- paste(toupper(substr(cancer, 1, 1)), substr(cancer, 2, nchar(cancer)), sep="")
    message("Setting default cancer to ", cancer, ".")
    options(WhereAreEcDNAs.cancer = cancer)
    eval(parse(text=paste0("data(cancerSpeciGenes)")))
    #print(cancerSpeciGenes[toupper(cancer)])
    return(unname(unlist(cancerSpeciGenes[toupper(cancer)])))
  }
  invisible(0)
}

getWhereAreEcDNAsCancer <- function(cancer = NULL){
  supportedcancers <- c("gbm","brca")
  .WhereAreEcDNAsCancergenes <- options()[["WhereAreEcDNAs.cancer"]]
  
  if(!is.null(.WhereAreEcDNAsCancergenes)){
    ag <- .WhereAreEcDNAsCancergenes
    if(tolower(ag) %in% supportedcancers){
      genes <- addWhereAreEcDNAsCancer(tolower(ag))
      return(genes)
    }else{
      stop("option(WhereAreEcDNAs.genome) : ", ag, " This genomic information is not supported.")
    }
    
  }else{
    stop("option(WhereAreEcDNAs.genome) : ", ag, " This genomic information is not supported.")
  }
}

#save(cancerSpeciGenes, file = "./data/cancerSpeciGenes.rda")
#x1 <- eval(parse(text=paste0("data(cancerSpeciGenes)")
# getCancerGenes <- function(cancer_type) {
#   cancerSpeciGenes <- list(
#     "ACC" = c("ABCB1", "TP53", "CTNNB1"),
#     "BLCA" = c("FGFR3", "TP53", "RB1"),
#     "BRCA" = c("BRCA1", "BRCA2", "TP53"),
#     "CESC" = c("TP53", "PIK3CA", "PTEN"),
#     "CHOL" = c("TP53", "KRAS", "SMAD4"),
#     "COAD" = c("APC", "KRAS", "TP53"),
#     "DLBC" = c("BCL6", "MYC", "TP53"),
#     "ESCA" = c("TP53", "CDKN2A", "NOTCH1"),
#     "GBM" = c("EGFR", "CDK4", "MDM2"),
#     "HNSC" = c("TP53", "CDKN2A", "EGFR"),
#     "KICH" = c("MET", "TP53", "CDKN2A"),
#     "KIRC" = c("VHL", "PBRM1", "SETD2"),
#     "KIRP" = c("PTEN", "CDKN2A", "TP53"),
#     "LAML" = c("FLT3", "NPM1", "DNMT3A"),
#     "LGG" = c("IDH1", "IDH2", "TP53"),
#     "LIHC" = c("TP53", "CTNNB1", "AXIN1"),
#     "LUAD" = c("EGFR", "KRAS", "TP53"),
#     "LUSC" = c("TP53", "PTEN", "CDKN2A"),
#     "MESO" = c("BAP1", "NF2", "CDKN2A"),
#     "OV" = c("BRCA1", "BRCA2", "TP53"),
#     "PAAD" = c("KRAS", "TP53", "CDKN2A"),
#     "PCPG" = c("SDHB", "SDHD", "RET"),
#     "PRAD" = c("TP53", "PTEN", "AR"),
#     "READ" = c("APC", "KRAS", "TP53"),
#     "SARC" = c("TP53", "CDKN2A", "RB1"),
#     "SKCM" = c("BRAF", "NRAS", "CDKN2A"),
#     "STAD" = c("TP53", "CDH1", "ARID1A"),
#     "TGCT" = c("KIT", "KRAS", "TP53"),
#     "THCA" = c("BRAF", "RAS", "RET"),
#     "THYM" = c("NOTCH1", "TP53", "KMT2D"),
#     "UCEC" = c("PTEN", "TP53", "PIK3CA"),
#     "UCS" = c("TP53", "PTEN", "PIK3CA"),
#     "UVM" = c("GNAQ", "GNA11", "BRAF")
#   )
#   
#   
#   if (cancer_type %in% names(cancer_genes)) {
#     return(cancer_genes[[cancer_type]])
#   } else {
#     return("未找到该癌症类型的致癌基因信息。")
#   }
# }