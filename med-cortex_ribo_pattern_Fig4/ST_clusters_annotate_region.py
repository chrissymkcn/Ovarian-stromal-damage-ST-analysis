import pandas as pd
import numpy as np
import scanpy as sc
import os
import importlib
from PIL import Image
import pickle as pkl
import json
import gseapy as gp
import matplotlib.pyplot as plt
from gseapy import barplot, dotplot

file_path = '/home/chrissy1/spatial/stomics/ovary_froz/spatial_utils.py'
spec = importlib.util.spec_from_file_location("spatial_utils", file_path)
spatial_utils = importlib.util.module_from_spec(spec)
spec.loader.exec_module(spatial_utils)
binsize = 50
savedir = f'/home/chrissy1/spatial/stomics/ovary_froz/redo/seurat/bin{binsize}_processed'
os.chdir(savedir)
adata = sc.read_h5ad('SCT.h5ad')
sep_dict = {'Fresh': [14000, 0],
            'Slow': [0, 14000], 
            'Vitri': [15000, 0]}

for sample in sep_dict.keys():
    cf = sep_dict[sample]
    if cf[0]!=0:
        cell_filt = (adata[adata.obs['batch']==sample].obsm['spatial'][:,0]>cf[0])
    elif cf[1]!=0: 
        cell_filt = (adata[adata.obs['batch']==sample].obsm['spatial'][:,1]>cf[1])
    slice1 = adata[adata.obs['batch']==sample][cell_filt,].obs_names
    slice2 = adata[adata.obs['batch']==sample][~cell_filt,].obs_names
    adata.obs.loc[slice1, 'slice'] = f'{sample}_slice1'
    adata.obs.loc[slice2, 'slice'] = f'{sample}_slice2'
adata.obs['slice'].value_counts()
annot_mean_counts = pd.read_csv('/home/chrissy1/spatial/stomics/ovary_froz/redo/seurat/bin50_processed/other_result/annot_mean.csv', index_col=0)
sample_sortbys = {
    'Fresh': ['x',True],
    'Slow': ['y', True],
    'Vitri': ['x', False]
}
distance = {}
for batch in adata.obs['slice'].unique(): 
    ad = adata[adata.obs['slice']==batch]
    sortby = sample_sortbys[batch.split('_')[0]]
    sp_loc = pd.DataFrame(ad.obsm['spatial'][:,:2], columns=['x', 'y'], index=ad.obs.index)
    sp_loc['annot'] = ad.obs['annot']
    annot_spatial_loc_mean = sp_loc.groupby('annot').mean()
    annot_spatial_loc_mean['counts'] = ad.obs['annot'].value_counts()[annot_spatial_loc_mean.index]
    annot_spatial_loc_mean['mean_nCount_RNA'] = ad.obs.groupby('annot').apply(lambda x: x['nCount_RNA'].mean())[annot_spatial_loc_mean.index]
    distance[batch] = annot_spatial_loc_mean.sort_values(sortby[0], ascending=sortby[1])
before_mid_cortex = {}
after_mid_cortex = {}
for k,df in distance.items():
    before = [i for i in df.loc[:'ST(KRT19)', :].index]
    after = [i for i in df.loc['ST(KRT19)':, :].index]
    before_mid_cortex[k] = before
    after_mid_cortex[k] = after
before_sets = [set(x) for x in before_mid_cortex.values()]
after_sets = [set(x) for x in after_mid_cortex.values()]
most_common_before = list(set.intersection(*before_sets))
most_common_after = list(set.intersection(*after_sets)-set('ST(KRT19)'))
adata.obs['region'] = '0'
adata.obs.loc[adata.obs['annot'].isin(most_common_before), 'region'] = 'medulla'
adata.obs.loc[adata.obs['annot'].isin(most_common_after), 'region'] = 'cortex'
adata.obs.loc[adata.obs['annot'].isin(['ST(KRT19)', 'ST(SECTM1)']), 'region'] = 'mid-cortex'
adata.obs['region'].replace({'0':'mid-cortex'}, inplace=True)
adata.write(f'{savedir}/SCT.h5ad')
