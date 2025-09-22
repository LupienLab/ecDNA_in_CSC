rm(list = ls())

source("./utils/librarys.R")
addWhereAreEcDNAsGenome("hg38")

sample <- "G946_T"
nth <- 1

fragPath <- paste0("/Users/chupan/Documents/data_su2c/scATAC/fragment/", sample, ".fragments.tsv.gz")
cellBarcode <- readLines(paste0("/Users/chupan/Documents/ArchR/outputs.",  sample, "/", sample, ".filtercell.txt"))
ampliconDir <- list.files("/Users/chupan/Documents/data_su2c/WGS/amplicon/", full.names = T)
ampliconPath <- ampliconDir[grep(sample, ampliconDir)][nth]

#~~~ prepare the inputs
fragment <- createFragmentObject(path = fragPath, cells = cellBarcode)
narrowPeak <- RunMACS2(fragment, sampleName = sample, pathToMacs2 = "/Users/chupan/miniconda3/envs/cpan/bin/macs2")

#~~~ create a whereAreEcDNAs object 
object <- whereAreEcDNAsObj(
  fragment = fragment,
  peak = narrowPeak,
  sampleName = sample,
  geneAnnotation = getGeneAnnotation(),
  genomeAnnotation = getGenomeAnnotation()
) %>% identecDNAcell(., pathToAmplicon = ampliconPath, nth = nth)

saveRDS(object, file = paste0(sample, ".outputs", "/", sample, "-", nth, ".obj.rds"))


