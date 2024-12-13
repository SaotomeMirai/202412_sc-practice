###########################
# 2024-12-13 
# pbmc2700
###########################

# 设置seurat对象
library(dplyr)
library(Seurat)
library(patchwork)
library(R.utils)

pbmc.data <- Read10X(data.dir = "D:\\Rwork\\202312_PBMC2700\\PBMC2700\\filtered_gene_bc_matrices\\hg19")
pbmc <- CreateSeuratObject(counts = pbmc.data,project = "pbmc2700",
                           min.cells = 3,min.features = 2000)

# 识别线粒体基因
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
# Show QC metrics for the first 5 cells
head(pbmc@meta.data, 5)
