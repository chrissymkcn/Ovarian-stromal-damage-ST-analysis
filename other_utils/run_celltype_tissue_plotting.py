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
savedir = f'/home/chrissy1/spatial/stomics/ovary_froz/redo/seurat/bin{binsize}_processed'

os.chdir(savedir)   
adata = sc.read_h5ad(f'SCT.h5ad')
annot_col = 'annot'
top = 10


for batch in adata.obs.batch.unique():
    batch_figdir = f'figures/show/{batch}'
    os.makedirs(batch_figdir) if not os.path.exists(batch_figdir) else None
    # p = dict(zip(adata.obs[annot_col].unique(), adata.uns[f'{annot_col}_colors']))
    p = 'tab20'
    markers = pd.read_csv(f'/home/chrissy1/spatial/stomics/ovary_froz/redo/seurat/bin{binsize}_processed/bin{binsize}_annot.csv').dropna()
    marker_dict = {k: markers[markers['annot']==k]['gene'].unique()[0] for k in markers['annot'].unique()}
    for cluster in ['ERY_ENDO(HBB)']: # adata.obs[annot_col].unique() ['O(PDCD5)'] ['ST_LIP', 'ST_SM', 'B_Plasma']:
        print('CELL TYPE', cluster)
        tiles = spatial_utils.find_tiles(adata[adata.obs['batch']==batch], annot_col=annot_col, main_cli=cluster, top_n=top, thres=0, tile_size=[1000, 1000])
        for i, row in tiles.iterrows():
            subsetted = adata[adata.obs['batch']==batch]
            subsetted = subsetted[(subsetted.obsm['spatial'][:, 0] >= row['coord_xmin']) & (subsetted.obsm['spatial'][:, 0] <= row['coord_xmax']) & (subsetted.obsm['spatial'][:, 1] >= row['coord_ymin']) & (subsetted.obsm['spatial'][:, 1] <= row['coord_ymax'])]
            groups = list(subsetted.obs[annot_col].value_counts().index[:5].astype(str).values)
            groups = groups + [cluster] if cluster not in groups else groups
            f1 = '-'.join([cluster, str(i)])
            f2 = '-'.join([cluster, 'no_marker', str(i)])
            f3 = '-'.join([cluster, 'ref', str(i)])
            f4 = '-'.join([cluster, marker_dict[cluster], str(i)])
            sc.pl.spatial(subsetted, color=annot_col, library_id=batch, 
                            crop_coord=row[['coord_xmin', 'coord_xmax', 'coord_ymin', 'coord_ymax']], 
                            palette = p,
                            save=f'/{batch}/{f1}', 
                            groups=groups,
                            # na_action="ignore"
                            )
            plt.close()
            # sc.pl.spatial(subsetted, color='annot_no_marker', library_id=batch, 
            #                     crop_coord=row[['coord_xmin', 'coord_xmax', 'coord_ymin', 'coord_ymax']], 
            #                     palette = p,
            #                     save=f'/{batch}/{f2}'
            #                     )
            # plt.close()
            sc.pl.spatial(subsetted, color=annot_col, library_id=batch, 
                                crop_coord=row[['coord_xmin', 'coord_xmax', 'coord_ymin', 'coord_ymax']], 
                                palette = p,
                                alpha=0,
                                save=f'/{batch}/{f3}',
                                # na_action="ignore"
                                )
            plt.close()
            sc.pl.spatial(subsetted, color=marker_dict[cluster], library_id=batch, 
                        crop_coord=row[['coord_xmin', 'coord_xmax', 'coord_ymin', 'coord_ymax']], 
                        save=f'/{batch}/{f4}', use_raw=False,
                        )
            plt.close()