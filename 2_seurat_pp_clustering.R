library(Seurat)
packageVersion('Seurat')
library(dplyr)
library(ggplot2)
library(SeuratData)
library(SeuratDisk)
library(patchwork)
library(dplyr)
library(anndata)
reticulate::use_python("/home/chrissy1/.conda/envs/dcats/bin/python")
out = '/home/chrissy1/spatial/stomics/ovary_froz/redo/seurat/bin20_processed'
data_rep = '/home/chrissy1/spatial/stomics/spDiff/datasets/stereoseq/ov_froz'  # /binsize_[1, 10, 20, 50, 100].h5ad

bin50 = paste0(data_rep, '/binsize_50.h5ad')
Convert(bin50, dest = "h5seurat", overwrite = TRUE)
bin50 = LoadH5Seurat(paste0(data_rep, '/binsize_50.h5seurat'), meta.data = FALSE, misc = FALSE)
meta = read.csv(paste0(data_rep, '/binsize_50_obs.csv'))
bin50 = CreateSeuratObject(counts = bin50@assays$RNA@counts, meta.data = meta)
my_vector = row.names(bin50@meta.data)
last_character <- substr(my_vector, nchar(my_vector), nchar(my_vector))
bin50@meta.data$sample = last_character
saveRDS(bin50, file = paste0(out, '/bin50.rds'))


# bin20 = paste0(data_rep, '/binsize_20.h5ad')
# Convert(bin20, dest = "h5seurat", overwrite = TRUE)
# bin20 = LoadH5Seurat(paste0(data_rep, '/binsize_20.h5seurat'), meta.data = FALSE, misc = FALSE)
# meta = read.csv(paste0(data_rep, '/binsize_20_obs.csv'))
# bin20 = CreateSeuratObject(counts = bin20@assays$RNA@counts, meta.data = meta)
# my_vector = row.names(bin20@meta.data)
# last_character <- substr(my_vector, nchar(my_vector), nchar(my_vector))
# bin20@meta.data$sample = last_character
# saveRDS(bin20, file = paste0(out, '/bin20.rds'))



# ## Preprocessing and clustering
out = '/home/chrissy1/spatial/stomics/ovary_froz/redo/seurat'
bin50 = readRDS(paste0(out, '/bin50.rds'))

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
bin50[["percent.mt"]] <- PercentageFeatureSet(bin50, pattern = "^MT-")


# Visualize QC metrics as a violin plot
png(paste0(out, '/bin50_vlnplot_qc.png'), width = 10, height = 10, units = 'in', res = 300)
VlnPlot(bin50, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "sample", pt.size = 0)
dev.off()

VlnPlot(bin50, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "sample", pt.size = 0)

plot1 <- FeatureScatter(bin50, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(bin50, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
bin50 <- subset(bin50, subset = nFeature_RNA > 100 & nFeature_RNA < 2500 & percent.mt < 15)
# split datasets and process without integration
bin50[["RNA"]] <- split(bin50[["RNA"]], f = bin50$sample)
bin50 <- SCTransform(bin50)
bin50 <- RunPCA(bin50, npcs = 30, verbose = F)
bin50 <- RunUMAP(bin50, dims = 1:30)

png(paste0(out, '/bin50_umap_sct_no_integration.png'), width = 10, height = 10, units = 'in', res = 300)
DimPlot(bin50, reduction = "umap", group.by = c("sample"))
dev.off()

DimPlot(bin50, reduction = "umap", group.by = c("sample"))
saveRDS(bin50, file = paste0(out, '/bin50_processed.rds'))
bin50

# integrate datasets
integrate datasets
bin50 <- IntegrateLayers(object = bin50, 
                        method = HarmonyIntegration, 
                        normalization.method = "SCT", 
                        verbose = F)
saveRDS(bin50, file = paste0(out, '/bin50_processed.rds'))
bin50

bin50 = readRDS("/home/chrissy1/spatial/stomics/ovary_froz/redo/seurat/bin50_processed.rds")
# bin50 <- FindNeighbors(bin50, reduction = "harmony", dims = 1:30)
# saveRDS(bin50, file = paste0(out, '/bin50_processed.rds'))
bin50 <- FindClusters(bin50, resolution = 1, algorithm = 2) 
bin50 <- RunUMAP(bin50, dims = 1:30, reduction = "harmony")

# Visualization
png(paste0(out, '/bin50_umap_sct_after_integration.png'), width = 10, height = 10, units = 'in', res = 300)
DimPlot(bin50, reduction = "umap", group.by = c("sample"))
dev.off()
DimPlot(bin50, reduction = "umap", group.by = c("sample"))

png(paste0(out, '/bin50_umap_sct_after_integration_seurat_clusters.png'), width = 10, height = 10, units = 'in', res = 300)
DimPlot(bin50, reduction = "umap", group.by = c("seurat_clusters"))
dev.off()
DimPlot(bin50, reduction = "umap", group.by = c("seurat_clusters"))

# sctransform v2 results for deg analysis: https://satijalab.org/seurat/archive/v4.3/sctransform_v2_vignette
bin50 = readRDS("/home/chrissy1/spatial/stomics/ovary_froz/redo/seurat/bin50_processed.rds")
bin50 = PrepSCTFindMarkers(bin50)
de_markers <- FindAllMarkers(bin50, assay = "SCT", only.pos = TRUE)

saveRDS(bin50, file = paste0(out, '/bin50_processed.rds'))
write.csv(de_markers, file = paste0(out, '/bin50_de_markers_res_1.csv'))


#### bin20 ####
# DEBUG for sctransform: 
# install.packages("https://cran.r-project.org/src/contrib/Archive/matrixStats/matrixStats_1.1.0.tar.gz", repos = NULL, type = "source")

# bin20 = readRDS(paste0(out, '/bin20.rds'))

# # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
# bin20[["percent.mt"]] <- PercentageFeatureSet(bin20, pattern = "^MT-")

# # Visualize QC metrics as a violin plot
# png(paste0(out, '/figures/bin20_vlnplot_qc.png'), width = 10, height = 10, units = 'in', res = 300)
# VlnPlot(bin20, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "sample", pt.size = 0)
# dev.off()


# bin20 <- subset(bin20, subset = nFeature_RNA > 30 & nFeature_RNA < 2500 & percent.mt < 15 & nCount_RNA < 2000)
# # split datasets and process without integration
# bin20[["RNA"]] <- split(bin20[["RNA"]], f = bin20$sample)
# bin20 <- SCTransform(bin20)
# bin20 <- RunPCA(bin20, npcs = 30, verbose = F)
# bin20 <- RunUMAP(bin20, dims = 1:30)

# png(paste0(out, '/figures/bin20_umap_sct_no_integration.png'), width = 10, height = 10, units = 'in', res = 300)
# DimPlot(bin20, reduction = "umap", group.by = c("sample"))
# dev.off()
# DimPlot(bin20, reduction = "umap", group.by = c("sample"))

# saveRDS(bin20, file = paste0(out, '/bin20_processed.rds'))

# # integrate datasets
# bin20 <- IntegrateLayers(object = bin20, 
#                         method = HarmonyIntegration, 
#                         normalization.method = "SCT", 
#                         verbose = F)
# saveRDS(bin20, file = paste0(out, '/bin20_processed.rds'))

# bin20 <- FindNeighbors(bin20, reduction = "harmony", dims = 1:30)
# bin20 <- FindClusters(bin20, resolution = 1, algorithm = 2)
# bin20 <- RunUMAP(bin20, dims = 1:30, reduction = "harmony")

# # Visualization
# png(paste0(out, '/figures/bin20_umap_sct_after_integration.png'), width = 10, height = 10, units = 'in', res = 300)
# DimPlot(bin20, reduction = "umap", group.by = c("sample"))
# dev.off()

# DimPlot(bin20, reduction = "umap", group.by = c("sample"))

# png(paste0(out, '/figures/bin20_umap_sct_after_integration_seurat_clusters.png'), width = 10, height = 10, units = 'in', res = 300)
# DimPlot(bin20, reduction = "umap", group.by = c("seurat_clusters"))
# dev.off()
# DimPlot(bin20, reduction = "umap", group.by = c("seurat_clusters"))

# bin20 = PrepSCTFindMarkers(bin20)
# de_markers <- FindAllMarkers(bin20, assay = "SCT", only.pos = TRUE)

# write.csv(de_markers, file = paste0(out, '/bin20_de_markers_res_1.csv'))
# saveRDS(bin20, file = paste0(out, '/bin20_processed.rds'))
