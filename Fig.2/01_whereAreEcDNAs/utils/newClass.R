Fragment <- setClass(
  Class = "Fragment",
  slots = list(
    path = "character",
    hash = "character",
    cells = "ANY"
  )
)


setClassUnion("characterOrNull", c("character", "NULL"))
setClassUnion("GRangesOrNull", c("GRanges", "NULL"))


setClass("whereAreEcDNAs", 
         representation(
           projectMetadata = "SimpleList",
           sampleColData = "DataFrame",
           sampleMetadata = "SimpleList",
           cellColData = "DataFrame", 
           geneAnnotation = "SimpleList",
           genomeAnnotation = "SimpleList",
           anyMatix = "characterOrNull",
           peak = "GRangesOrNull" 
         )
)


setMethod("show", "whereAreEcDNAs",
          
          function(object) {
            scat <- function(fmt, vals=character(), exdent=2, n = 5, ...){
              vals <- ifelse(nzchar(vals), vals, "''")
              lbls <- paste(S4Vectors:::selectSome(vals, maxToShow = n), collapse=" ")
              txt <- sprintf(fmt, length(vals), lbls)
              cat(strwrap(txt, exdent=exdent, ...), sep="\n")
            }
            cat("class:", class(object), "\n")
            cat("outputDirectory:", object@projectMetadata$outputDirectory, "\n")
            
            o <- tryCatch({
              object@cellColData$sample
            }, error = function(x){
              stop(paste0("\nError accessing sample info from WhereAreEcDNAsProject.",
                          "\nThis is most likely the issue with saving the WhereAreEcDNAsProject as an RDS",
                          "\nand not with save/loadWhereAreEcDNAsProject. This bug has mostly been attributed",
                          "\nto bioconductors DataFrame saving cross-compatability. We added a fix to this.",
                          "\nPlease Try:",
                          "\n\trecoverWhereAreEcDNAsObj(WhereAreEcDNAsObj)",
                          "\n\nIf that does not work please report to Github: https://github.com/LupienLab/WhereAreEcDNAs/issues"
              ))
            })
            
            scat("samples(%d): %s\n", rownames(object@sampleColData))
            #scat("sampleColData names(%d): %s\n", names(object@sampleColData))
            scat("cellColData names(%d): %s\n", names(object@cellColData))
            scat("numberOfCells(%d): %s\n", nrow(object@cellColData))
            
          }
          
)

whereAreEcDNAsObj <- function(
    fragment = NULL,
    sampleName = NULL,
    outputDirectory = NULL, 
    geneAnnotation = getGeneAnnotation(),
    genomeAnnotation = getGenomeAnnotation(),
    threads = 1,
    peaks = NULL
){
  if(is.null(outputDirectory)){
    outputDirectory <- paste0(sampleName, ".outputs")
  }
  geneAnnotation <- .validGeneAnnotation(geneAnnotation)
  genomeAnnotation <- .validGenomeAnnotation(genomeAnnotation)
  geneAnnotation <- .validGeneAnnoByGenomeAnno(geneAnnotation = geneAnnotation, genomeAnnotation = genomeAnnotation)
  
  if(grepl(" ", outputDirectory)){
    stop("outputDirectory cannot have a space in the path! Path : ", outputDirectory)
  }
  dir.create(outputDirectory,showWarnings=FALSE)
  if(grepl(" ", normalizePath(outputDirectory))){
    stop("outputDirectory cannot have a space in the full path! Full path : ", normalizePath(outputDirectory))
  }
  
  message("Getting sampleNames...")
  sampleColData <- DataFrame(row.names = sampleName, fragfile = fragment@path)
  sampleMetadata <- SimpleList(lapply(sampleName, function(x) SimpleList()))
  names(sampleMetadata) <- sampleName
  
  message("Getting Cell Metadata...")
  
  md <- DataFrame(fragment@cells)
  colnames(md) <- "CellNames"
  md$CellNames <- paste0(sampleName,"#",md$CellNames)
  md$sample <- Rle(sampleName, nrow(md))
  md$idx <- 1:nrow(md)
  rownames(md) <- md$CellNames
  md <- md[, -which(colnames(md)=="CellNames")]
  metadataList <- md
  
  message("Merging Cell Metadata...")
  # 获取 metadataList 的列名
  allCols <- unique(c("sample", rev(sort(colnames(metadataList)))))
  cellColData <- metadataList
  avMat <- NULL
  message("Initializing whereAreEcDNAsObj...")
  
  Obj <- new("whereAreEcDNAs",  # 修改为 whereAreEcDNAs 类
               projectMetadata = SimpleList(outputDirectory = normalizePath(outputDirectory)),
               sampleColData = sampleColData,
               sampleMetadata = sampleMetadata,
               cellColData = cellColData,
               geneAnnotation = .validGeneAnnotation(geneAnnotation),
               genomeAnnotation = .validGenomeAnnotation(genomeAnnotation),
               anyMatix = avMat,
               peak = peaks  # 将 peaks 传入并存储在 whereAreEcDNAs 对象中
  )
  
  Obj
}

#Accessor methods adapted from Seurat 

#'Accessing cellColData directly from dollar.sign accessor
#' 
#' This function will allow direct access to cellColData with a `$` accessor.
#'
#' @export
#'
".DollarNames.whereAreEcDNAs" <- function(x, pattern = ''){
  cn <- as.list(c("cellNames",colnames(x@cellColData)))
  names(cn) <- c("cellNames",colnames(x@cellColData))
  return(.DollarNames(x = cn, pattern = pattern))
}

#'Accessing cellColData directly from dollar.sign accessor
#' 
#' This function will allow direct access to cellColData with a `$` accessor.
#'
#' @export
#'
"$.whereAreEcDNAs" <- function(x, i){
  if(i=="cellNames"){
    return(rownames(x@cellColData))
  }else{
    val <- x@cellColData[[i, drop = TRUE]]
    return(as.vector(val))
  }
}

#' Add directly to cellColData directly from dollar.sign accessor
#' 
#' This function will allow adding directly to cellColData with a `$` accessor.
#'
#' @export
#'
"$<-.whereAreEcDNAs" <- function(x, i, value){
  if(i == "sample"){
    stop("sample is a protected column in cellColData. Please do not try to overwrite this column!")
  }
  if(i == "cellNames"){
    stop("cellNames is a protected column in cellColData. Please do not try to overwrite this column!")
  }
  if(i == "nFrags"){
    stop("nFrags is a protected column in cellColData. Please do not try to overwrite this column!")
  }
  if(object.size(Rle(value)) < 2 * object.size(value)){ #Check if Rle is more efficient for storage purposes...
    value <- Rle(value)
  }
  if(!is.null(value)){
    if(length(value)==1){
      value <- Rle(value, lengths = nrow(x@cellColData))
    }
  }
  x@cellColData[[i]] <- value
  return(x)
}



setMethod(
  f = "colnames",
  signature = c("x" = "whereAreEcDNAs"),
  definition = function(x) {
    colnames(x@cellColData)
  }
)

setMethod(
  f = "rownames",
  signature = c("x" = "whereAreEcDNAs"),
  definition = function(x) {
    rownames(x@cellColData)
  }
)

