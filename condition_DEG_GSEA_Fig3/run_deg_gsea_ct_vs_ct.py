# import scanpy as sc
# import h5py
# import numpy as np
# import pandas as pd
# from PIL import Image
# import matplotlib.pyplot as plt
# import seaborn as sns
# import gseapy as gp
# import os

# binsize = 50
# savedir = f'/home/chrissy1/spatial/stomics/ovary_froz/redo/seurat/bin{binsize}_processed'
# os.chdir(savedir)
# res_adata_path = f'SCT.h5ad'
# adata = sc.read_h5ad(res_adata_path)
# ## Sample wise deg
# # One vs the rest
# if 'deg' not in os.listdir():
#     os.makedirs('deg')
# os.chdir('deg')

# sc.tl.rank_genes_groups(adata, 'batch', method='wilcoxon', n_genes=1000, pts=True)
# df = sc.get.rank_genes_groups_df(adata, group=None)
# df['names'] = adata.var_names[df['names'].astype(int)]
# df.to_csv(f'sample_markers.csv', sep=',', index=False) 

# df = pd.read_csv('sample_markers.csv')
# top_n_degs = 10
# df['names'] = df['names'].astype(str)
# df = df[~df.names.str.contains('MT') & ~df.names.str.contains('RPS') & ~df.names.str.contains('RPL')]
# df = df[(df['pvals_adj'] < 0.05) & (df['logfoldchanges'] > 0)]
# top = df.groupby('group').apply(lambda x: x.sort_values('scores', ascending=False).head(top_n_degs)).reset_index(drop=True)  # get top 20 de genes
# markers = {k: v['names'].tolist() for k,v in top.groupby('group')}  # get top 20 de gene names
# # flip the dotplot axes
# adata_batch_key = 'batch'
# pdf_name = f'batch_deg_dotplot{top_n_degs}.png'
# plt.rcParams["figure.figsize"] = (15, 12)
# sc.pl.dotplot(adata, var_names=markers, groupby=adata_batch_key,
#             use_raw=False, dendrogram=False,
#             cmap='Blues',
#             dot_max=0.2, dot_min=0, 
#             swap_axes=True, 
#             save=pdf_name)

#     ## Reference wise deg
#     for batch in adata.obs.batch.unique():
#         sc.tl.rank_genes_groups(adata, 'batch', method='wilcoxon', n_genes=1000, pts=True, reference=batch)
#         df = sc.get.rank_genes_groups_df(adata, group=None)
#         df.to_csv(f'ref_{batch}_sample_deg.csv', sep=',', index=False) 

#     csvs = [x for x in os.listdir() if x.endswith('.csv') and x.startswith('ref')]

#     for csv in csvs:
#         df = pd.read_csv(csv)
#         top_n_degs = 10
#         df = df[~df.names.str.contains('MT') & ~df.names.str.contains('RPS') & ~df.names.str.contains('RPL')]
#         df = df[(df['pvals_adj'] < 0.05) & (df['logfoldchanges'] > 0)]
#         top = df.groupby('group').apply(lambda x: x.sort_values('logfoldchanges', ascending=False).head(top_n_degs)).reset_index(drop=True)  # get top 20 de genes
#         markers = {k: v['names'].tolist() for k,v in top.groupby('group')}  # get top 20 de gene names
#         # flip the dotplot axes
#         adata_batch_key = 'batch'
#         pdf_name = csv.replace('.csv', '_dotplot.png')
#         plt.rcParams["figure.figsize"] = (5, 12)
#         sc.pl.dotplot(adata, var_names=markers, groupby=adata_batch_key,
#                     use_raw=False, dendrogram=True,
#                     cmap='Blues',
#                     dot_max=0.3, dot_min=0, 
#                     swap_axes=True, 
#                     save=pdf_name)

# GSEA across all cells / specific cell types
# gsea
import scanpy as sc
import h5py
import numpy as np
import pandas as pd
from PIL import Image
import matplotlib.pyplot as plt
import seaborn as sns
import gseapy as gp
import os
import anndata
import pickle as pkl

binsize = 50
savedir = f'/home/chrissy1/spatial/stomics/ovary_froz/redo/seurat/bin{binsize}_processed'
os.chdir(savedir)
gsea_res_dir = os.path.join(savedir, 'gsea_region')
other_res_dir = os.path.join(savedir, 'other_result')
if not os.path.exists(gsea_res_dir):
    os.makedirs(gsea_res_dir)
res_adata_path = f'SCT.h5ad'
adata_full = sc.read_h5ad(res_adata_path)
# adata_full = adata_full[adata_full.obs['annot']!='Unknown', :]
perc = 1

subset_col = 'annot_no_marker'
crude_celltypes = [['ST', 'ST_LIP', 'ST_SM'], ['B_Plasma']]
gsea_col = 'annot'

for crude_celltype in crude_celltypes:
    adata = adata_full[adata_full.obs[subset_col].isin(crude_celltype), :] if crude_celltype!='all' else adata_full
    gsea_crude_celltype_dir = os.path.join(gsea_res_dir, '-'.join(crude_celltype))
    if not os.path.exists(gsea_crude_celltype_dir):
        os.makedirs(gsea_crude_celltype_dir)
    os.chdir(gsea_crude_celltype_dir)
    
    # remove cell types with less than 10 cells in any sample
    celltypes = list(adata.obs[gsea_col].unique())
    cnts = adata.obs[['batch', gsea_col]].value_counts()
    to_remove = cnts[cnts<10].index.get_level_values(1).unique()
    adata = adata[~adata.obs[gsea_col].isin(to_remove), :]
    celltypes = list(adata.obs[gsea_col].unique())
    
    for ct in celltypes:
        #### subsetting the data randomly to speed up the gsea
        cells = np.random.choice(adata_full.obs.index, int(perc*adata_full.n_obs), replace=False)
        adata = adata_full[cells, :]
        rest_gs = [x for x in adata.obs[gsea_col].unique() if x!=ct]
        print('gsea', ct)
        matrix = adata.to_df().T
        design = ['pos' if x == ct else 'neg' for x in adata.obs[gsea_col]]
        SELECTED_GS = ['GO_Biological_Process_2023', 'MSigDB_Hallmark_2020', 'Reactome_2022']
        result = gp.gsea(data=matrix, # or data='./P53_resampling_data.txt'
                            gene_sets=SELECTED_GS,
                            cls=design,
                            permutation_type='phenotype',
                            pheno_pos='pos',
                            pheno_neg='neg',
                            permutation_num=1000, # reduce number to speed up test
                            outdir=None,  # do not write output to disk
                            method='signal_to_noise',
                            threads=40, seed=7)
        with open(f'{ct}_res.pkl', 'wb') as f: 
            pkl.dump(result, f, protocol=pkl.HIGHEST_PROTOCOL)
    figdir = 'figures'
    if not os.path.exists(figdir):
        os.makedirs(figdir)    
    for ct in celltypes:
        pklFile = f'{ct}_res.pkl'
        with open(pklFile, 'rb') as f:
            res = pkl.load(f).res2d
        res = res[~res.Term.str.startswith('KEGG')]  # remove KEGG pathways 
        res['Term'] = res['Term'].str.split('__').str[1]
        res['Term'] = res['Term'].str.split('\(GO').str[0]
        res['Term'] = res['Term'].str.split(' R-').str[0]
        posres = res[(res['FDR q-val'] < 0.05) & (res['NES']>0)].sort_values(by='NES', ascending=False)
        # for positive or negative NES create the directory to save the dotplot
        n = posres.shape[0]
        try:
            gp.dotplot(posres,	
                        column="FDR q-val",
                        cmap=plt.cm.viridis,
                        size=5,
                        top_term=min(80,n),
                        figsize=(8,20), cutoff=0.05,
                        ofname=f'{figdir}/{ct}_dotplot.pdf')
            gp.dotplot(posres,	
                        column="FDR q-val",
                        cmap=plt.cm.viridis,
                        size=5,
                        top_term=min(15,n),
                        figsize=(5,6), cutoff=0.05,
                        ofname=f'{figdir}/{ct}_dotplot_top.pdf')
        except ValueError:
            print("ValueError: Warning: No enrich terms when cutoff = 0.05 for celltype: ", ct)
