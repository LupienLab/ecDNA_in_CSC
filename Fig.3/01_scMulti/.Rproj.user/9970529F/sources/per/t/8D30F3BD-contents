#******** step1: create Seurat object for individual GSC cultures ********#

rm(list = ls())

library(Seurat) #~~~ 5.2.1
library(Signac) #~~~ 1.14.0
library(dplyr)
library(ggplot2)
library(patchwork)
library(GenomicRanges)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

id = "G523_L" ## G523_L  G583_L  G620_L  G797_L  G828_L  G837_L

#~~~ load scRNA-seq counts
rna_counts <- Read10X(data.dir = paste0("/Users/chupan/Documents/data_su2c/scMulti/cellranger-arc/", id, "/filtered_feature_bc_matrix"))
seurat_obj <- CreateSeuratObject(counts = rna_counts$`Gene Expression`, assay = "RNA")

#~~~ load scATAC-seq fragments
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- "UCSC"  
genome(annotations) <- "hg38"

fragments <- paste0("/Users/chupan/Documents/data_su2c/scMulti/cellranger-arc/", id, "/atac_fragments.tsv.gz")
chrom_assay <- CreateChromatinAssay(counts = rna_counts$Peaks, sep = c(":", "-"), genome = 'hg38', fragments = fragments, annotation = annotations, min.cells = 10, min.features = 200)
seurat_obj[["ATAC"]] <- chrom_assay

#~~~ QC for the RNA side
DefaultAssay(seurat_obj) <- "RNA"
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-", assay = "RNA")
length(seurat_obj$orig.ident)
#~~~ visualize QC metrics
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 500 & nFeature_RNA < 8000 & percent.mt < 20) #subset = nFeature_RNA > 500 & nFeature_RNA < 8000 & percent.mt < 20
length(seurat_obj$orig.ident)

#~~~ QC for the ATAC side
DefaultAssay(seurat_obj) <- "ATAC"
seurat_obj <- NucleosomeSignal(seurat_obj)
seurat_obj <- TSSEnrichment(seurat_obj, fast = FALSE)
length(seurat_obj$orig.ident)

# visualize
VlnPlot(seurat_obj, features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"), ncol = 4, pt.size = 0)
seurat_obj <- subset(seurat_obj, subset = nCount_ATAC > 4000 & nCount_ATAC < 80000 & nucleosome_signal < 2 & TSS.enrichment > 2)
seurat_obj$id <- id

# save results 
dir.create(paste0("outputs.", id))
write.table(names(seurat_obj$orig.ident), file = paste0("./outputs.", id, "/", id, ".filtercell.txt"), row.names = F, col.names = F, sep = "\t", quote = F)
saveRDS(seurat_obj, file = paste0("./outputs.", id, "/", id, ".rds"))

# G523: 4820
# G583: 6935
# G620: 6958
# G797: 5007
# G828: 6111
# G837: 6168




