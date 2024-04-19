library(utils, pos = "package:base", verbose = FALSE)
library(Matrix)

writeCountData <- function(seo, filePath, assay = "RNA", slots=c("counts", "data", "scale.data"), reductions=c("pca", "tsne", "umap")) {
    # Create the directory if it does not exist
    dir.create(filePath, recursive = TRUE, showWarnings = FALSE) 
    
    for (s in slots) {
        if (s == 'scale.data') {
            countsPath <- paste0(filePath, '/', assay, '_', s, '.csv')
            if (!file.exists(countsPath)) {
                write.csv(slot(seo@assays[[assay]], s), countsPath)
            }
        } else {
            countsPath <- paste0(filePath, '/', assay, '_', s, '.mtx')
            if (!file.exists(countsPath)) {
                writeMM(slot(seo@assays[[assay]], s), countsPath)
            }
        }
        # countsPath <- paste0(filePath, '/', assay, '_', s, '.mtx')
        # if (!file.exists(countsPath)) {
        #     if (s == 'scale.data') {
        #         writeMM(as(slot(seo@assays[[assay]], s), 'sparseMatrix'), countsPath)
        #     }
        #     writeMM(slot(seo@assays[[assay]], s), countsPath)
        #     if (s == 'scale.data') {
        #         sctgenesPath = paste0(filePath, '/', assay, '_sctgenes.csv')
        #         sctGenes = data.frame(rownames(slot(bin20@assays[['SCT']], 'scale.data')))
        #         colnames(sctGenes) <- 'Gene'
        #         write.csv(sctGenes, sctgenesPath, quote = FALSE, row.names = FALSE)
        # }
        # }
    }
    
    barcodesPath <- paste0(filePath, '/', assay, '_barcodes.csv')
    genesPath <- paste0(filePath, '/', assay, '_genes.csv')
    cellMetaPath <- paste0(filePath, '/', assay, '_cellMeta.csv')

    if (!file.exists(barcodesPath)) {
        barcodes <- data.frame(colnames(seo))
        colnames(barcodes) <- 'Barcode'
        write.csv(barcodes, barcodesPath, quote = FALSE, row.names = FALSE)
    }
    
    if (!file.exists(genesPath)) {    
        genes <- data.frame(rownames(seo@assays[[assay]]))
        colnames(genes) <- 'Gene'
        write.csv(genes, genesPath, quote = FALSE, row.names = FALSE)
    }

    if (!file.exists(cellMetaPath)) {
        cellMeta <- seo@meta.data
        for (reduction in reductions) {
            if (reduction %in% names(seo@reductions)) {
                cellMeta <- cbind(cellMeta, seo@reductions[[reduction]]@cell.embeddings)
            }
        }
        write.csv(cellMeta, cellMetaPath, quote = FALSE, row.names = TRUE)
    }
}

bin20 = readRDS("/home/chrissy1/spatial/stomics/ovary_froz/redo/seurat/bin50_processed/bin50_processed.rds")

writeCountData(bin20, "/home/chrissy1/spatial/stomics/ovary_froz/redo/seurat/bin50_processed", assay='SCT', 
                slots=c("counts", "data", "scale.data"), reductions=c("pca", "umap", "harmony"))
