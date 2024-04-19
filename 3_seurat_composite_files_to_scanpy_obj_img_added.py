import os
import scanpy as sc
import pandas as pd
from scipy import io, sparse
import importlib
import numpy as np
from PIL import Image

def loadDataSavedFromSeurat(path: str, prefix: str = "RNA", X_slot='data', slots=["counts", "data", "scale.data"], 
                            reductions=["pca", "umap", "harmony"]):
    """
    Load data from a specified path.

    Args:
        path (str): Path to the directory containing the input files.
        prefix (str, optional): Prefix for the input file names (default: 'RNA').

    Returns:
        Anndata object with loaded counts, barcodes, genes, cell meta data and umap coordinates.
    """
    
    # Load the count matrix from the provided path
    if X_slot not in slots:
        raise ValueError("X_slot not found in slots to load. Please check the slots and X_slot.")
    
    for slot in slots:
        try:
            slot = X_slot
            counts = io.mmread(f"{path}/{prefix}_{slot}.mtx")
            counts = sparse.csr_matrix(counts)
            adata = sc.AnnData(counts.T)
            for slot in slots:
                if slot != X_slot:
                    if slot=='scale.data':
                        counts = pd.read_csv(f"{path}/{prefix}_{slot}.csv")
                        # counts = sparse.csr_matrix(counts)
                        adata.uns[f'{prefix}_{slot}'] = counts.T
                    else:
                        counts = io.mmread(f"{path}/{prefix}_{slot}.mtx")
                        counts = sparse.csr_matrix(counts)
                        adata.raw = sc.AnnData(counts.T)
        except FileNotFoundError:
            raise FileNotFoundError("Count matrix file not found. Please check the path and prefix.")
    del counts
    # Load the barcodes and genes information
    try:
        barcodes = pd.read_csv(f"{path}/{prefix}_barcodes.csv")
        genes = pd.read_csv(f"{path}/{prefix}_genes.csv")
    except FileNotFoundError:
        raise FileNotFoundError("Barcodes or genes file not found. Please check the path and prefix.")
                
    # Set the observation names to the barcodes
    adata.obs_names = barcodes['Barcode'].values
    
    # Set the variable names to the gene names
    adata.var_names = genes['Gene'].values
    
    # Load the cell metadata and set it to the observation annotations in adata object
    try:
        cellMeta = pd.read_csv(f"{path}/{prefix}_cellMeta.csv")
        adata.obs = cellMeta
    except FileNotFoundError:
        raise FileNotFoundError("Cell metadata file not found. Please check the path and prefix.")
        
    # Get UMAP coordinates from cell metadata and set them to the observation embeddings in adata object
    for reduction in reductions:
        relevant_cols = [colname for colname in cellMeta.columns.values if reduction in colname]
        if len(relevant_cols)>0:
            adata.obsm['X_' + reduction] = cellMeta.loc[:, relevant_cols].values
        else:
            print(f"{reduction} coordinates not found. Skipping {reduction} coordinates...")

    if all(pd.Series(['x', 'y', 'z']).isin(adata.obs.columns)):
        adata.obsm['spatial'] = cellMeta.loc[:, ['x', 'y', 'z']].values
    elif all(pd.Series(['x', 'y']).isin(adata.obs.columns)):
        adata.obsm['spatial'] = cellMeta.loc[:, ['x', 'y']].values

    # Delete the variables that are no longer needed to free up memory
    del barcodes, genes, cellMeta
    return adata

binsize = 50

adata = loadDataSavedFromSeurat(path = f'/home/chrissy1/spatial/stomics/ovary_froz/redo/seurat/bin{binsize}_processed/for_export', prefix='SCT', 
                                X_slot = 'data',
                                slots = ["data", "counts", "scale.data"], # 
                                reductions=["umap", "harmony"]) # "pca", 

adata.write(f'/home/chrissy1/spatial/stomics/ovary_froz/redo/seurat/bin{binsize}_processed/SCT.h5ad')




#### Read in images (absent when converting from seurat to scanpy)
Image.MAX_IMAGE_PIXELS = 2000000000
highres_scale = 0.7
lowres_scale = 0.1
batch = 'batch'
binsize = 50

adata = sc.read_h5ad(f'/home/chrissy1/spatial/stomics/ovary_froz/redo/seurat/bin{binsize}_processed/SCT.h5ad')
data_rep = '/home/chrissy1/spatial/stomics/spDiff/datasets/stereoseq/ov_froz'  # /binsize_[1, 10, 20, 50, 100].h5ad
bin50_ori = sc.read_h5ad(os.path.join(data_rep, 'binsize_50.h5ad'))  # original bin50 data converted by stereopy

file_paths = ['/home/chrissy1/spatial/stomics/ovary_froz/rawdata/E-SCT20221109001_SZD20230526112_1/Result/Fresh/01.StandardWorkflow_Result/Register/D01154B5_regist.tif',
                '/home/chrissy1/spatial/stomics/ovary_froz/rawdata/E-SCT20221109001_SZD20230526112_2/Result/Slow/01.StandardWorkflow_Result/Register/D01970B5_regist.tif',
                '/home/chrissy1/spatial/stomics/ovary_froz/rawdata/E-SCT20221109001_SZD20230526112_2/Result/Vitri/01.StandardWorkflow_Result/Register/B01320B6_regist.tif']

for i,impath in enumerate(file_paths):
    im = Image.open(impath)
    sample = impath.split('/')[-4]
    mask = [im_tis_mask[sample][int(row[1]), int(row[0])] for row in adata.obsm['spatial'][adata.obs[batch]==sample]]
    bin50_ori.obs.loc[bin50_ori.obs[batch]==sample, 'in_tissue'] = mask

bin50_ori.obs.to_csv(os.path.join('/home/chrissy1/spatial/stomics/ovary_froz/redo/seurat/bin50_processed/for_export', 'binsize_50_ori_obs.csv'))

obs = pd.read_csv('/home/chrissy1/spatial/stomics/ovary_froz/redo/seurat/bin50_processed/for_export/binsize_50_ori_obs.csv')

# Load the spatial_utils module
file_path = '/home/chrissy1/spatial/stomics/ovary_froz/spatial_utils.py'
spec = importlib.util.spec_from_file_location("spatial_utils", file_path)
spatial_utils = importlib.util.module_from_spec(spec)
spec.loader.exec_module(spatial_utils)

file_paths = ['/home/chrissy1/spatial/stomics/ovary_froz/rawdata/E-SCT20221109001_SZD20230526112_1/Result/Fresh/01.StandardWorkflow_Result/Register/D01154B5_regist.tif',
                '/home/chrissy1/spatial/stomics/ovary_froz/rawdata/E-SCT20221109001_SZD20230526112_2/Result/Slow/01.StandardWorkflow_Result/Register/D01970B5_regist.tif',
                '/home/chrissy1/spatial/stomics/ovary_froz/rawdata/E-SCT20221109001_SZD20230526112_2/Result/Vitri/01.StandardWorkflow_Result/Register/B01320B6_regist.tif']

tissue_cut_paths = ['/home/chrissy1/spatial/stomics/ovary_froz/rawdata/E-SCT20221109001_SZD20230526112_1/Result/Fresh/01.StandardWorkflow_Result/Register/D01154B5_tissue_cut.tif',
                '/home/chrissy1/spatial/stomics/ovary_froz/rawdata/E-SCT20221109001_SZD20230526112_2/Result/Slow/01.StandardWorkflow_Result/Register/D01970B5_tissue_cut.tif',
                '/home/chrissy1/spatial/stomics/ovary_froz/rawdata/E-SCT20221109001_SZD20230526112_2/Result/Vitri/01.StandardWorkflow_Result/Register/B01320B6_tissue_cut.tif']

for i,impath in enumerate(file_paths):
    im = Image.open(impath)
    sample = impath.split('/')[-4]
    hi_dim = (int(im.size[0]//(1/highres_scale)), int(im.size[1]//(1/highres_scale)))
    low_dim = (int(im.size[0]//(1/lowres_scale)), int(im.size[1]//(1/lowres_scale)))
    im_hires = {sample: im.resize(hi_dim)}
    im_lowres = {sample: im.resize(low_dim)}
    im_tis_mask = {sample: np.array(Image.open(tissue_cut_paths[i]))}
    spatial_utils.add_spatial_img(adata=adata, bins=binsize, highres_scale=highres_scale, lowres_scale=lowres_scale, 
                                im_lowres=im_lowres, im_hires=im_hires)
    mask = [im_tis_mask[sample][int(row[1]), int(row[0])] for row in adata.obsm['spatial'][adata.obs[batch]==sample]]
    adata.obs.loc[adata.obs[batch]==sample, 'in_tissue'] = mask
    mask = [im_tis_mask[sample][int(row[1]), int(row[0])] for row in bin50_ori.obsm['spatial'][bin50_ori.obs[batch]==sample]]

adata.write(f'/home/chrissy1/spatial/stomics/ovary_froz/redo/seurat/bin{binsize}_processed/SCT.h5ad')
