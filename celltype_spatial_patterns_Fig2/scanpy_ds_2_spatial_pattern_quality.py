import squidpy as sq
import matplotlib.pyplot as plt
import scanpy as sc
import pandas as pd
import numpy as np
from scipy.cluster import hierarchy as sch
from copy import deepcopy
import os
binsize = 50
savedir = f'/home/chrissy1/spatial/stomics/ovary_froz/redo/seurat/bin{binsize}_processed'
os.chdir(savedir)
# adata_full = sc.read_h5ad(f'SCT.h5ad')
output_dir = 'spatial_patterns'
os.makedirs(output_dir) if not os.path.exists(output_dir) else None
os.chdir(output_dir)
adata_full = sc.read_h5ad(f'{savedir}/SCT.h5ad')

for sample in ['Vitri']:  # 'Fresh', 'Slow', 
    # adata_full = sc.read_h5ad(f'SCT.h5ad')
    output_dir = f'{savedir}/spatial_patterns'
    os.makedirs(output_dir) if not os.path.exists(output_dir) else None
    os.chdir(output_dir)

    figdir = f'{output_dir}/{sample}'
    if not os.path.exists(figdir):
        os.makedirs(figdir)

    cluster_key = 'annot'
    adata = adata_full[adata_full.obs['batch']==sample] if sample != 'all' else adata_full
    sq.gr.spatial_neighbors(adata, coord_type="generic", spatial_key="spatial")
    sq.gr.nhood_enrichment(adata, cluster_key=cluster_key)

    # from scipy.cluster import hierarchy as sch

    # leiden_clusters = adata.obs[cluster_key].cat.categories
    # df_nhood_enr = pd.DataFrame(
    #     adata.uns[f"{cluster_key}_nhood_enrichment"]["zscore"],
    #     columns=leiden_clusters,
    #     index=leiden_clusters,
    # )

    # crude_ct = 'ST_LIP'
    # st_cl = [x for x in adata.obs['annot'].unique() if crude_ct in x]
    # sq.pl.spatial_scatter(
    #     adata,
    #     groups=st_cl,
    #     color='annot',
    #     # size=15,
    #     img=False,
    #     palette='tab20',
    #     library_id=sample,
    #     save=f'{sample}/{crude_ct}.png'
    # )
    # crude_ct = 'ST_SM'
    # st_cl = [x for x in adata.obs['annot'].unique() if crude_ct in x]
    # sq.pl.spatial_scatter(
    #     adata,
    #     groups=st_cl,
    #     color='annot',
    #     # size=15,
    #     img=False,
    #     palette='tab20',
    #     library_id=sample,
    #     save=f'{sample}/{crude_ct}.png'
    # )
    # crude_ct = ['ENDO_SM', 'FIB']
    # name = '-'.join(crude_ct)
    # st_cl = [x for x in adata.obs['annot'].unique() if any([True for c in crude_ct if c in x])]
    # sq.pl.spatial_scatter(
    #     adata,
    #     groups=st_cl,
    #     color='annot',
    #     # size=15,
    #     img=False,
    #     palette='tab20',
    #     library_id=sample,
    #     save=f'{sample}/{name}.png'
    # )
    # crude_ct = ['BEX2', 'RPL', 'SEC']
    # st_cl = [x for x in adata.obs['annot'].unique() if any([True for c in crude_ct if c in x])]
    # name = '-'.join(st_cl)
    # sq.pl.spatial_scatter(
    #     adata,
    #     groups=st_cl,
    #     color='annot',
    #     # size=15,
    #     img=False,
    #     palette='tab20',
    #     library_id=sample,
    #     save=f'{sample}/{name}.png'
    # )

    # crude_ct = ['RAD', 'RPK', 'BNC', 'EXT', 'MAG', 'FNDC', 'MAML', 'PRKG1', 'RBF', 'SLCO', 'ZFPM']
    # st_cl = [x for x in adata.obs['annot'].unique() if any([True for c in crude_ct if c in x])]
    # name = '-'.join(st_cl)
    # sq.pl.spatial_scatter(
    #     adata,
    #     groups=st_cl,
    #     color='annot',
    #     # size=15,
    #     img=False,
    #     palette='tab20',
    #     library_id=sample,
    #     save=f'{sample}/inner_ST.png'
    # )
    # crude_ct = ['KRT']
    # st_cl = [x for x in adata.obs['annot'].unique() if any([True for c in crude_ct if c in x])]
    # name = '-'.join(st_cl)
    # sq.pl.spatial_scatter(
    #     adata,
    #     groups=st_cl,
    #     color='annot',
    #     # size=15,
    #     img=False,
    #     palette='tab20',
    #     library_id=sample,
    #     save=f'{sample}/{name}.png'
    # )

    #### Autocorrelation: Moran’s I Score

    sq.gr.spatial_autocorr(adata, mode="moran")
    num_view = 20
    top_autocorr = (
        adata.uns["moranI"]["I"].sort_values(ascending=False).head(num_view).index.tolist()
    )
    bot_autocorr = (
        adata.uns["moranI"]["I"].sort_values(ascending=True).head(num_view).index.tolist()
    )
    # print statistics of all moran's I values
    adata.uns["moranI"][adata.uns['moranI']['pval_norm'] < 0.05].shape, adata.uns["moranI"]['I'].describe()
    ser_counts = adata.obs[cluster_key].value_counts()
    ser_counts.name = "cell counts"
    meta_leiden = pd.DataFrame(ser_counts)

    cat_name = cluster_key
    sig_leiden = pd.DataFrame(
        columns=adata.var_names, index=adata.obs[cat_name].cat.categories
    )
    for clust in adata.obs[cat_name].cat.categories:
        sig_leiden.loc[clust] = adata[adata.obs[cat_name].isin([clust]), :].X.mean(0)
    sig_leiden = sig_leiden.transpose()
    leiden_clusters = sig_leiden.columns.tolist()
    sig_leiden.columns = leiden_clusters
    meta_leiden.index = sig_leiden.columns.tolist()
    meta_leiden[cluster_key] = pd.Series(
        meta_leiden.index.tolist(), index=meta_leiden.index.tolist()
    )
    sq.pl.spatial_scatter(
        adata, color=top_autocorr, cmap="Reds", img=False, figsize=(5, 5), library_id=sample,
        save=f"{sample}/top_MI_{sample}.jpg", use_raw=False
    )
    # top cell types based on average expression of top_autocorr genes
    sig_leiden.loc[top_autocorr].mean(axis=0).sort_values(ascending=False).index.tolist()[
        :5
    ]
    sq.pl.spatial_scatter(
        adata, color=bot_autocorr, cmap="Reds", img=False, figsize=(5, 5), library_id=sample,
        save=f"{sample}/bot_MI_{sample}.jpg", use_raw=False
    )
    # top cell types based on average expression of bot_autocorr genes
    sig_leiden.loc[bot_autocorr].mean(axis=0).sort_values(ascending=False).index.tolist()[
        :5
    ]
    # adata.uns['moranI'].to_csv(f'moranI_{sample}.csv')
    # print('n genes significantly spatially correlated', adata.uns['moranI'][adata.uns['moranI']['pval_norm']<0.05].shape[0])
    # correlated = adata.uns['moranI'][adata.uns['moranI']['pval_norm']<0.05]
    # plt.hist(correlated['I'], bins=30)
    # correlated['I'].describe()
    