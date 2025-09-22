
#创建fragment对象
createFragmentObject <- function(
    path,
    cells = NULL,
    validate.fragments = TRUE, #是否检查 cells 中的 barcode 是否都存在于文件里（默认 TRUE）
    verbose = TRUE,  # 是否打印进度消息（默认 TRUE）
    ...
) {
  validate.fragments <- FALSE
  #如果是远程文件，则自动关闭后续对片段内容的完整性验证（因为无法快速随机访问远程文件）。
  # is.remote <- isRemote(x = path)
  # if (is.remote) {
  #   validate.fragments <- FALSE
  # }
  if (!file.exists(path)) {
    stop("Fragment file does not exist.")
  }
  #获取与该片段文件配套的 Tabix 索引文件（通常是 .tbi），并返回其路径。
  index.file <- GetIndexFile(fragment = path, verbose = verbose)
  # if (is.remote) {
  #   con <- gzcon(con = url(description = path))
  # } else {
  #   con <- path
  # }
  con <- path
  df <- readLines(con = con, n = 10000)
  for (i in df) {
    if (grepl(pattern = '^#', x = i)) {
      next
    } else {
      if (length(x = strsplit(x = i, split = "\t")[[1]]) != 5) {
        stop("Incorrect number of columns found in fragment file")
      } else {
        break
      }
    }
  }
  if (!is.null(x = cells)) {
    if (is.null(names(x = cells))) {
      # assume cells are as they appear in the assay
      names(x = cells) <- cells
    }
  }
  # compute hash of the file and index
  if (verbose) {
    message("Computing hash")
  }
  # if (!is.remote) {
  #   path <- normalizePath(path = path, mustWork = TRUE)
  # }
  path <- normalizePath(path = path, mustWork = TRUE)
  # will be NA if file remote
  hashes <- md5sum(files = c(path, index.file))
  # create object
  frags <- new(
    Class = "Fragment",
    path = path,
    hash = unname(obj = hashes),
    cells = cells
  )
  return(frags)
  # # validate cells
  # if (!is.null(x = cells) & validate.fragments) {
  #   if (ValidateCells(object = frags, verbose = verbose, ...)) {
  #     return(frags)
  #   } else {
  #     stop("Not all cells requested could be found in the fragment file.")
  #   }
  # } else {
  #   return(frags)
  # }
}

AddMissing <- function(x, cells, features = NULL) {
  # add columns with zeros for cells or features not in matrix
  missing.cells <- setdiff(x = cells, y = colnames(x = x))
  if (!(length(x = missing.cells) == 0)) {
    null.mat <- sparseMatrix(
      i = c(),
      j = c(),
      dims = c(nrow(x = x), length(x = missing.cells))
    )
    rownames(x = null.mat) <- rownames(x = x)
    colnames(x = null.mat) <- missing.cells
    x <- cbind(x, null.mat)
  }
  x <- x[, cells, drop = FALSE]
  if (!is.null(x = features)) {
    missing.features <- setdiff(x = features, y = rownames(x = x))
    if (!(length(x = missing.features) == 0)) {
      null.mat <- sparseMatrix(
        i = c(),
        j = c(),
        dims = c(length(x = missing.features), ncol(x = x))
      )
      rownames(x = null.mat) <- missing.features
      colnames(x = null.mat) <- colnames(x = x)
      x <- rbind(x, null.mat)
    }
    x <- x[features, , drop = FALSE]
  }
  return(x)
}
