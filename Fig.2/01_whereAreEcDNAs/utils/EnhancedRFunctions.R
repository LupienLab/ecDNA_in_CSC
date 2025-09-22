"%ni%" <- Negate("%in%")
'%bcin%' <- function(x, table) S4Vectors::match(x, table, nomatch = 0) > 0  #正相匹配
'%bcni%' <- function(x, table) !(S4Vectors::match(x, table, nomatch = 0) > 0) #反向匹配

.requirePackage <- function(x = NULL, load = TRUE, installInfo = NULL, source = NULL){
  if(x %in% rownames(installed.packages())){
    if(load){
      suppressPackageStartupMessages(require(x, character.only = TRUE))
    }else{
      return(0)
    }
  }else{
    if(!is.null(source) & is.null(installInfo)){
      if(tolower(source) == "cran"){
        installInfo <- paste0('install.packages("',x,'")')
      }else if(tolower(source) == "bioc"){
        installInfo <- paste0('BiocManager::install("',x,'")')
      }else{
        stop("Unrecognized package source, available are cran/bioc!")
      }
    }
    if(!is.null(installInfo)){
      stop(paste0("Required package : ", x, " is not installed/found!\n  Package Can Be Installed : ", installInfo))
    }else{
      stop(paste0("Required package : ", x, " is not installed/found!"))
    }
  }
}
###lapply对args运算
.batchlapply <- function(args = NULL, sequential = FALSE){
  
  if(is.null(args$tstart)){
    args$tstart <- Sys.time()
  }
  args$subThreads <- 1
  args$threads <- 1
  
  args <- args[names(args) %ni% c("registryDir", "parallelParam", "subThreading")]
  outlist <- do.call(.safelapply, args)
  return(outlist)
}

.safelapply <- function(..., threads = 1, preschedule = FALSE){
  threads <- 1
  o <- lapply(...)
  return(o)
}
#从x中提取文件扩展名
.fileExtension <- function (x = NULL){
  pos <- regexpr("\\.([[:alnum:]]+)$", x)
  ifelse(pos > -1L, substring(x, pos + 1L), "")
}
#获取索引如果成功返回TRUE，失败尝试 indexTabix(file, format = "bed")
.isTabix <- function(file = NULL){
  tryCatch({
    TabixFile(file)
    return(TRUE)
  }, error = function(x){
    tryCatch({
      if(getWhereAreCNVsVerbose()) message("Attempting to index ", file," as tabix..")
      indexTabix(file, format = "bed")
      return(TRUE)
    }, error = function(y){
      message("Tabix indexing failed for ", file,". Note that deCNV requires bgzipped fragment files which is different from gzip. See samtools bgzip!")
      return(FALSE)
    })
  })
}
#一般返回true，WhereAreEcDNAs的一种保护机制
getWhereAreEcDNAsVerbose <- function(){
  WhereAreEcDNAsVerbose <- options()[["WhereAreCNVs.verbose"]]
  if(!is.logical(WhereAreEcDNAsVerbose)){
    options(WhereAreEcDNAs.verbose = TRUE)
    return(TRUE)
  }
  WhereAreEcDNAsVerbose
}
#隐藏函数执行信息
.suppressAll <- function(expr = NULL){
  suppressPackageStartupMessages(suppressMessages(suppressWarnings(expr)))
}

#重写对稀疏矩阵的验证，该R包特殊的验证
.checkSparseMatrix <- function(x, ncol = NULL){
  isSM <- is(x, 'sparseMatrix')
  if(!isSM){
    if(is.null(ncol)){
      stop("ncol must not be NULL if x is not a matrix!")
    }
    cnames <- tryCatch({
      names(x)
    }, error = function(e){
      colnames(x)
    })
    if(length(cnames) != ncol){
      stop("cnames != ncol!")
    }
    x <- Matrix::Matrix(matrix(x, ncol = ncol), sparse=TRUE)
    colnames(x) <- cnames
  }
  x
}
#对数值进行放缩准备
.scaleDims <- function(x, scaleMax = NULL){
  if(!is.null(scaleMax)){
    .rowZscores(m=x, min=-scaleMax, max = scaleMax, limit = TRUE)
  }else{
    .rowZscores(m=x)
  }
}
#对数值进行放缩
.rowZscores <- function(m = NULL, min = -2, max = 2, limit = FALSE){
  z <- sweep(m - rowMeans(m), 1, matrixStats::rowSds(m),`/`)
  if(limit){
    z[z > max] <- max
    z[z < min] <- min
  }
  return(z)
}

#创建中间文件
.tempfile <- function(pattern = "tmp", tmpdir = "tmp", fileext = "", addDOC = TRUE){
  
  if(!dir.exists(tmpdir)){
    if(file.exists(tmpdir)){
      stop(paste0("Attempted to create temporary directory ", tmpdir," but a file already exists with this name. Please remove this file and try again!"))
    }
  }
  
  dir.create(tmpdir, showWarnings = FALSE)
  
  if(!dir.exists(tmpdir)){
    stop(paste0("Unable to create temporary directory ", tmpdir,". Check file permissions!")) 
  }
  
  if(addDOC){
    doc <- paste0("-Date-", Sys.Date(), "_Time-", gsub(":","-", stringr::str_split(Sys.time(), pattern=" ",simplify=TRUE)[1,2]))
  }else{
    doc <- ""
  }
  
  tempfile(pattern = paste0(pattern, "-"), tmpdir = tmpdir, fileext = paste0(doc, fileext))
  
}

#安全保存RDS
.safeSaveRDS <- function(
    object = NULL, 
    file = "", 
    ascii = FALSE, 
    version = NULL, 
    compress = TRUE, 
    refhook = NULL
){
  #Try to save a test data.frame in location
  testDF <- data.frame(a=1,b=2)
  canSave <- suppressWarnings(tryCatch({
    saveRDS(object = testDF, file = file, ascii = ascii, version = version, compress = compress, refhook = refhook)
    TRUE
  }, error = function(x){
    FALSE
  }))
  if(!canSave){
    dirExists <- dir.exists(dirname(file))
    if(dirExists){
      stop("Cannot saveRDS. File Path : ", file)
    }else{
      stop("Cannot saveRDS because directory does not exist (",dirname(file),"). File Path : ", file)
    }
  }else{
    saveRDS(object = object, file = file, ascii = ascii, version = version, compress = compress, refhook = refhook)
  }
}
#对GRanges对某个基因的上下游进行扩展
extendGR <-  function(gr = NULL, upstream = NULL, downstream = NULL){
  #https://bioinformatics.stackexchange.com/questions/4390/expand-granges-object-different-amounts-upstream-vs-downstream
  isMinus <- BiocGenerics::which(strand(gr) == "-")
  isOther <- BiocGenerics::which(strand(gr) != "-")
  #Forward
  start(gr)[isOther] <- start(gr)[isOther] - upstream
  end(gr)[isOther] <- end(gr)[isOther] + downstream
  #Reverse
  end(gr)[isMinus] <- end(gr)[isMinus] + upstream
  start(gr)[isMinus] <- start(gr)[isMinus] - downstream
  return(gr)
}

#获取与该片段文件配套的 Tabix 索引文件（通常是 .tbi），并返回其路径。
GetIndexFile <- function(fragment, verbose = TRUE) {
  # 构建可能的索引文件名：.tbi 和 .csi
  index.paths <- c(paste0(fragment, ".tbi"),
                   paste0(fragment, ".csi"))
  # 仅保留真实存在的文件
  existing <- index.paths[file.exists(index.paths)]
  
  # 如果都不存在，则报错
  if (length(existing) == 0) {
    stop("Fragment 文件未建立索引，请先用 samtools 或 tabix 创建 .tbi 或 .csi 索引。")
  }
  
  # 如果两种索引都存在，优先使用 .tbi
  if (length(existing) == 2) {
    if (verbose) {
      message("同时检测到 .tbi 和 .csi 索引，优先使用 .tbi。")
    }
    return(existing[1])
  }
  
  # 只有一种索引时，直接返回它
  return(existing)
}

"Cells<-.Fragment" <- function(x, ..., value) {
  if (is.null(x = names(x = value))) {
    stop("Cells must be a named vector")
  }
  slot(object = x, name = "cells") <- value
  if (!ValidateCells(object = x, verbose = FALSE, ...)) {
    stop("Cells not present in fragment file")
  } else {
    return(x)
  }
}

Cells.Fragment <- function(x, ...) {
  cells <- slot(object = x, name = "cells")
  return(names(x = cells))
}

"Cells<-" <- function(x, ..., value) {
  UseMethod(generic = "Cells<-", object = x)
}

ChunkGRanges <- function(granges, nchunk) {
  if (length(x = granges) < nchunk) {
    nchunk <- length(x = granges)
  }
  chunksize <- as.integer(x = (length(granges) / nchunk))
  range.list <- sapply(X = seq_len(length.out = nchunk), FUN = function(x) {
    chunkupper <- (x * chunksize)
    if (x == 1) {
      chunklower <- 1
    } else {
      chunklower <- ((x - 1) * chunksize) + 1
    }
    if (x == nchunk) {
      chunkupper <- length(x = granges)
    }
    return(granges[chunklower:chunkupper])
  })
  return(range.list)
}


GRangesToString <- function(grange, sep = c("-", "-")) {
  regions <- paste0(
    as.character(x = seqnames(x = grange)),
    sep[[1]],
    start(x = grange),
    sep[[2]],
    end(x = grange)
  )
  return(regions)
}

ChunkGRanges <- function(granges, nchunk) {
  if (length(x = granges) < nchunk) {
    nchunk <- length(x = granges)
  }
  chunksize <- as.integer(x = (length(granges) / nchunk))
  range.list <- sapply(X = seq_len(length.out = nchunk), FUN = function(x) {
    chunkupper <- (x * chunksize)
    if (x == 1) {
      chunklower <- 1
    } else {
      chunklower <- ((x - 1) * chunksize) + 1
    }
    if (x == nchunk) {
      chunkupper <- length(x = granges)
    }
    return(granges[chunklower:chunkupper])
  })
  return(range.list)
}
#对GRanges对某个基因的上下游进行扩展
extendGR <-  function(gr = NULL, upstream = NULL, downstream = NULL){
  #https://bioinformatics.stackexchange.com/questions/4390/expand-granges-object-different-amounts-upstream-vs-downstream
  isMinus <- BiocGenerics::which(strand(gr) == "-")
  isOther <- BiocGenerics::which(strand(gr) != "-")
  #Forward
  start(gr)[isOther] <- start(gr)[isOther] - upstream
  end(gr)[isOther] <- end(gr)[isOther] + downstream
  #Reverse
  end(gr)[isMinus] <- end(gr)[isMinus] + upstream
  start(gr)[isMinus] <- start(gr)[isMinus] - downstream
  return(gr)
}
