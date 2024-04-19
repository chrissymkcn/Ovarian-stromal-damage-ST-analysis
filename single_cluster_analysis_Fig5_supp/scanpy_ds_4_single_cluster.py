import scanpy as sc
import pandas as pd
import numpy as np
import os
import importlib
from PIL import Image
import pickle as pkl
import gseapy as gp
import matplotlib.pyplot as plt

file_path = '/home/chrissy1/spatial/stomics/ovary_froz/spatial_utils.py'
spec = importlib.util.spec_from_file_location("spatial_utils", file_path)
spatial_utils = importlib.util.module_from_spec(spec)
spec.loader.exec_module(spatial_utils)
binsize = 50

annot_col = 'annot'
savedir = f'/home/chrissy1/spatial/stomics/ovary_froz/redo/seurat/bin{binsize}_processed'
os.chdir(savedir)
resdir = f'{savedir}/annot_celltype_ind_analysis'

gsea_col = 'leiden'

# for ct in ['ENDO_SM(CLDN5)', 'ST(BEX2)', 'O(PDCD5)', 'ST(KRT19)', 'ST(RPL41)']:
#     ct_dir = os.path.join(resdir, ct)
#     os.makedirs(ct_dir) if not os.path.exists(ct_dir) else None
#     os.chdir(ct_dir)
#     figdir = 'figures'
#     if not os.path.exists(figdir):
#         os.makedirs(figdir)    
#     subsetted = sc.read_h5ad('SCT.h5ad')
#     leiden_cls = subsetted.obs[gsea_col].unique()
#     for leiden_cl in leiden_cls:
#         if f'{leiden_cl}_res.pkl' not in os.listdir():
#             rest_gs = [x for x in leiden_cls if x!=leiden_cl]
#             matrix = subsetted.to_df().T
#             design = ['pos' if x == leiden_cl else 'neg' for x in subsetted.obs[gsea_col]]
#             SELECTED_GS = ['GO_Biological_Process_2023', 'MSigDB_Hallmark_2020', 'Reactome_2022']
#             result = gp.gsea(data=matrix, # or data='./P53_resampling_data.txt'
#                                 gene_sets=SELECTED_GS,
#                                 cls=design,
#                                 permutation_type='phenotype',
#                                 pheno_pos='pos',
#                                 pheno_neg='neg',
#                                 permutation_num=1000, # reduce number to speed up test
#                                 outdir=None,  # do not write output to disk
#                                 method='signal_to_noise',
#                                 threads=40, seed=7)
#             with open(f'{leiden_cl}_res.pkl', 'wb') as f: 
#                 pkl.dump(result, f, protocol=pkl.HIGHEST_PROTOCOL)
#         with open(f'{leiden_cl}_res.pkl', 'rb') as f:
#             res = pkl.load(f).res2d
#         res = res[~res.Term.str.startswith('KEGG')]  # remove KEGG pathways 
#         res['Term'] = res['Term'].str.split('__').str[1]
#         res['Term'] = res['Term'].str.split('\(GO').str[0]
#         res['Term'] = res['Term'].str.split(' R-').str[0]
#         top_lgs = [','.join(lgs.split(';')[:5]) for lgs in res['Lead_genes'].values]
#         res['Term'] = [f"{t} ({top_lgs[i]})" for i,t in enumerate(res['Term'].values)]
#         posres = res[(res['FDR q-val'] < 0.05) & (res['NES']>0)].sort_values(by='NES', ascending=False)
#         # for positive or negative NES create the directory to save the dotplot
#         n = posres.shape[0]
#         try:
#             gp.dotplot(posres,	
#                         column="FDR q-val",
#                         cmap=plt.cm.viridis,
#                         size=5,
#                         top_term=min(80,n),
#                         figsize=(8,20), cutoff=0.05,
#                         ofname=f'{figdir}/{leiden_cl}_dotplot.pdf')
#             gp.dotplot(posres,	
#                         column="FDR q-val",
#                         cmap=plt.cm.viridis,
#                         size=5,
#                         top_term=min(15,n),
#                         figsize=(5,6), cutoff=0.05,
#                         ofname=f'{figdir}/{leiden_cl}_dotplot_top.pdf')
#         except ValueError:
#             print("ValueError: Warning: No enrich terms when cutoff = 0.05 for leiden: ", leiden_cl)

# import pickle as pkl
# import pandas as pd
# import numpy as np
# import requests
# import json
# import time
# import pickle
# # import stereo as st
# import os
# import importlib

# file_path = '/home/chrissy1/spatial/stomics/ovary_froz/spatial_utils.py'
# spec = importlib.util.spec_from_file_location("spatial_utils", file_path)
# spatial_utils = importlib.util.module_from_spec(spec)
# spec.loader.exec_module(spatial_utils)

# gsea_col = 'leiden'

# for ct in ['O(PDCD5)', 'ENDO_SM(CLDN5)', 'ST(BEX2)', 'ST(KRT19)', 'ST(RPL41)']:
#     ct_dir = os.path.join(resdir,ct)
#     os.chdir(ct_dir)
#     csv = 'leiden_markers.csv'
#     markers = pd.read_csv(csv)
#     markers = markers[(markers['pvals_adj']<0.05) & (markers['logfoldchanges']>0)]
#     markers = markers[~markers['names'].str.contains('MT-') & ~markers['names'].str.contains('RPL') & ~markers['names'].str.contains('RPS')]
#     markers = markers.groupby('group').apply(lambda x: x.sort_values('logfoldchanges', ascending=False).head(100))
#     spatial_utils.scrape_mks(markers, 
#                 df_gene_col = 'names', 
#                 extra_cols=['Molecular function', 'single cell', 'biological process'],
#                 outfile = csv.replace('.csv', '_scraped.csv'))    
    
import scanpy as sc
import pandas as pd
import numpy as np
import os
import importlib
from PIL import Image
import matplotlib.pyplot as plt
file_path = '/home/chrissy1/spatial/stomics/ovary_froz/spatial_utils.py'
spec = importlib.util.spec_from_file_location("spatial_utils", file_path)
spatial_utils = importlib.util.module_from_spec(spec)
spec.loader.exec_module(spatial_utils)

binsize = 50
resdir = f'/home/chrissy1/spatial/stomics/ovary_froz/redo/seurat/bin{binsize}_processed/annot_celltype_ind_analysis'
for savedir in ['ST(BEX2)']: # ['ST(BEX2)', 'ENDO_SM(CLDN5)', 'O(PDCD5)', 'ST(KRT19)', 'ST(RPL41)']: # os.listdir(resdir):    
    os.chdir(os.path.join(resdir,savedir))   
    adata = sc.read_h5ad(f'SCT.h5ad')

    for batch in adata.obs.batch.unique():
        batch_figdir = f'figures/show/{batch}'
        os.makedirs(batch_figdir) if not os.path.exists(batch_figdir) else None
        annot_col = 'leiden'
        top = 5
        # p = dict(zip(adata.obs[annot_col].unique(), adata.uns[f'{annot_col}_colors']))
        p = 'tab20'
        adata.obs['batch'] = adata.obs['batch'].astype(str)
        adata.obs[annot_col] = adata.obs[annot_col].astype(str)

        markers = pd.read_csv(f'/home/chrissy1/spatial/stomics/ovary_froz/redo/seurat/bin{binsize}_processed/bin{binsize}_annot.csv').dropna()
        marker_dict = {k: markers[markers['annot']==k]['gene'].unique()[0] for k in markers['annot'].unique()}

        for cluster in adata.obs[annot_col].unique(): # ['O(PDCD5)'] ['ST_LIP', 'ST_SM', 'B_Plasma']:
            print('CELL TYPE', cluster)
            tiles = spatial_utils.find_tiles(adata[adata.obs['batch']==batch], annot_col=annot_col, main_cli=cluster, top_n=top, thres=0, tile_size=[1500, 1500])
            for i, row in tiles.iterrows():
                subsetted = adata[adata.obs['batch']==batch]
                subsetted = subsetted[(subsetted.obsm['spatial'][:, 0] >= row['coord_xmin']) & (subsetted.obsm['spatial'][:, 0] <= row['coord_xmax']) & (subsetted.obsm['spatial'][:, 1] >= row['coord_ymin']) & (subsetted.obsm['spatial'][:, 1] <= row['coord_ymax'])]
                groups = list(subsetted.obs[annot_col].value_counts().index[:5].astype(str).values)
                groups = groups + [cluster] if cluster not in groups else groups
                f1 = '-'.join([cluster, str(i)])
                f2 = '-'.join([cluster, 'annot', str(i)])
                f3 = '-'.join([cluster, 'ref', str(i)])
                sc.pl.spatial(subsetted, color=annot_col, library_id=batch, 
                                    crop_coord=row[['coord_xmin', 'coord_xmax', 'coord_ymin', 'coord_ymax']], 
                                    palette = p,
                                    save=f'/{batch}/{f1}', 
                                    groups=groups
                                    )
                plt.close()
                sc.pl.spatial(subsetted, color=annot_col, library_id=batch, 
                                    crop_coord=row[['coord_xmin', 'coord_xmax', 'coord_ymin', 'coord_ymax']], 
                                    palette = p,
                                    alpha=0,
                                    save=f'/{batch}/{f3}'
                                    )
                plt.close()