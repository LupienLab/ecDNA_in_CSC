CallPeaks <- function(object, ...) {
  print(object)
  UseMethod(generic = "CallPeaks", object = object)
}

GetFragmentData <- function(object, slot = "path") {
  return(methods::slot(object = object, name = slot))
}

CallPeaks.Fragment <- function(
    object,
    macs2.path = NULL,
    outdir = tempdir(),
    broad = FALSE,
    format = "BED",
    effective.genome.size = 2.7e9,
    extsize = 200,
    shift = -extsize/2,
    additional.args = NULL,
    name = "macs2",
    cleanup = TRUE,
    verbose = TRUE,
    ...
) {
  fragpath <- GetFragmentData(object = object, slot = "path")
  gr <- CallPeaks(
    object = fragpath,
    macs2.path = macs2.path,
    outdir = outdir,
    broad = broad,
    format = format,
    effective.genome.size = effective.genome.size,
    extsize = extsize,
    shift = shift,
    additional.args = additional.args,
    name = name,
    cleanup = cleanup,
    verbose = verbose,
    ...
  )
  return(gr)
}

SetIfNull <- function(x, y) {
  if (is.null(x = x)) {
    return(y)
  } else {
    return(x)
  }
}

CallPeaks.default <- function(
    object,
    macs2.path = NULL,
    outdir = tempdir(),
    broad = FALSE,
    format = "BED",
    effective.genome.size = 2.7e9,
    extsize = 200,
    shift = -extsize/2,
    additional.args = NULL,
    name = "macs2",
    cleanup = TRUE,
    verbose = TRUE,
    ...
) {
  if (!dir.exists(paths = outdir)) {
    stop("Requested output directory does not exist")
  }
  # find macs2
  macs2.path <- SetIfNull(
    x = macs2.path,
    y = unname(obj = Sys.which(names = "macs2"))
  )
  if (nchar(x = macs2.path) == 0) {
    stop("MACS2 not found. Please install MACS:",
         "https://macs3-project.github.io/MACS/")
  }
  name <- gsub(pattern = " ", replacement = "_", x = name)
  name <- gsub(pattern = .Platform$file.sep, replacement = "_", x = name)
  
  # if list of paths given, collapse to a single space-separated string
  if (length(x = object) > 1) {
    object <- sapply(
      X = object, FUN = function(x) paste0("'", x, "'"), USE.NAMES = FALSE
    )
    object <- Reduce(f = paste, x = object)
  } else {
    object <- paste0("'", object, "'")
    if (object == "''") {
      stop("No fragment files supplied")
    }
  }
  
  broadstring <- ifelse(test = broad, yes = " --broad ", no = "")
  nomod_str <- ifelse(
    test = format == "BED",
    yes = paste0(" --nomodel --extsize ",
                 as.character(x = extsize),
                 " --shift ",
                 as.character(x = shift)
    ),
    no = ""
  )
  
  cmd <- paste0(
    macs2.path,
    " callpeak -t ",
    object,
    " -g ",
    as.character(x = effective.genome.size),
    broadstring,
    " -f ",
    format,
    nomod_str,
    " -n ",
    "'",
    as.character(x = name),
    "'",
    " --outdir ",
    outdir,
    " ",
    additional.args
  )
  
  # call macs2
  system(
    command = cmd,
    wait = TRUE,
    ignore.stderr = !verbose,
    ignore.stdout = !verbose
  )
  
  if (broad) {
    # read in broadpeak
    df <- read.table(
      file = paste0(outdir, .Platform$file.sep, name, "_peaks.broadPeak"),
      col.names = c("chr", "start", "end", "name",
                    "score", "strand", "fold_change",
                    "neg_log10pvalue_summit", "neg_log10qvalue_summit")
    )
    files.to.remove <- paste0(
      name,
      c("_peaks.broadPeak", "_peaks.xls", "_peaks.gappedPeak")
    )
  } else {
    # read in narrowpeak file
    df <- read.table(
      file = paste0(outdir, .Platform$file.sep, name, "_peaks.narrowPeak"),
      col.names = c("chr", "start", "end", "name",
                    "score", "strand", "fold_change",
                    "neg_log10pvalue_summit", "neg_log10qvalue_summit",
                    "relative_summit_position")
    )
    files.to.remove <- paste0(
      name,
      c("_peaks.narrowPeak", "_peaks.xls", "_summits.bed")
    )
  }
  
  gr <- GenomicRanges::makeGRangesFromDataFrame(df = df, keep.extra.columns = TRUE, starts.in.df.are.0based = TRUE)
  if (cleanup) {
    files.to.remove <- paste0(outdir, .Platform$file.sep, files.to.remove)
    for (i in files.to.remove) {
      if (file.exists(i)) {
        file.remove(i)
      }
    }
  }
  gr <- keepStandardChromosomes(gr, pruning.mode = "coarse")
  return(gr)
}