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
res_adata_path = f'SCT.h5ad'
adata_full = sc.read_h5ad(res_adata_path)

perc = 1
gsea_col = 'batch'
celltype_annot = 'annot_no_marker'
gsea_res_dir = os.path.join(savedir, 'gsea')
if not os.path.exists(gsea_res_dir):
    os.makedirs(gsea_res_dir)
os.chdir(gsea_res_dir)

# remove cell types with less than 10 cells in any sample
celltypes = list(adata_full.obs[celltype_annot].unique())
cnts = adata_full.obs[['batch', celltype_annot]].value_counts()
to_remove = cnts[cnts<10].index.get_level_values(1).unique()
celltypes = ['ST_LIP', 'ENDO_SM']
# celltypes = ['all'] + [x for x in celltypes if x not in to_remove]
# celltypes = [x for x in celltypes if x not in to_remove]
# celltypes = ['FIB', 'ST_FIB', 'ST_SM','ENDO_SM', 'ST_LIP']

for ct in celltypes:
    gsea_ct_res_dir = os.path.join(gsea_res_dir, ct)
    if not os.path.exists(gsea_ct_res_dir):
        os.makedirs(gsea_ct_res_dir)
    adata = adata_full[adata_full.obs[celltype_annot]==ct, :] if ct != 'all' else adata_full
    if ct == 'all':
        #### subsetting the data randomly to speed up the gsea
        cells = np.random.choice(adata.obs.index, int(perc*adata.n_obs), replace=False)
        adata = adata[cells, :]
    gs = adata.obs[gsea_col].unique()
    for g in gs:
        rest_gs = [x for x in adata.obs[gsea_col].unique() if x!=g]
        for rest_g in rest_gs:
            combo = '-'.join([g, rest_g]) 
            print('gsea', combo, ct)
            sub_ad = adata[adata.obs[gsea_col].isin([g,rest_g])]
            matrix = sub_ad.to_df().T
            design = ['pos' if x == g else 'neg' for x in sub_ad.obs[gsea_col]]
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
            with open(f'{gsea_ct_res_dir}/{combo}.pkl', 'wb') as f: 
                pkl.dump(result, f, protocol=pkl.HIGHEST_PROTOCOL)    

for ct in celltypes:
    gsea_ct_res_dir = os.path.join(gsea_res_dir, ct)
    os.chdir(gsea_ct_res_dir)
    figdir = os.path.join(gsea_ct_res_dir, 'figures')
    # if os.path.exists(figdir):
    #     continue
    if not os.path.exists(figdir):
        os.makedirs(figdir)
    pklFiles = [x for x in os.listdir() if x.endswith('.pkl')]
    dfs = []
    for pklFile in pklFiles:
        with open(pklFile, 'rb') as f:
            res = pkl.load(f).res2d
        res['batch'] = pklFile.split('.')[0]  # add batch column to denote which batch comparison the results are from
        dfs.append(res)
    res = pd.concat(dfs, axis=0, ignore_index=True)
    res = res[~res.Term.str.startswith('KEGG')]  # remove KEGG pathways 
    res['Term'] = res['Term'].str.split('__').str[1]
    res['Term'] = res['Term'].str.split('\(GO').str[0]
    res['Term'] = res['Term'].str.split(' R-').str[0]
    posres = res[(res['FDR q-val'] < 0.05) & (res['NES']>0)].sort_values(by='NES', ascending=False)
    # for positive or negative NES create the directory to save the dotplot
    for combo in posres.batch.unique():
        n = posres[posres.batch == combo].shape[0]
        dotplot_name = combo
        combo_res = posres[(posres.batch == combo)]
        try:
            gp.dotplot(combo_res,	
                        column="FDR q-val",
                        cmap=plt.cm.viridis,
                        size=5,
                        top_term=min(80,n),
                        figsize=(8,20), cutoff=0.05,
                        ofname=os.path.join(figdir, f'{dotplot_name}_dotplot.pdf'))
            gp.dotplot(combo_res,	
                        column="FDR q-val",
                        cmap=plt.cm.viridis,
                        size=5,
                        top_term=min(15,n),
                        figsize=(5,6), cutoff=0.05,
                        ofname=os.path.join(figdir, f'{dotplot_name}_dotplot_top.pdf'))
        except ValueError:
            print("ValueError: Warning: No enrich terms when cutoff = 0.05 for combo: ", combo)
    # allTerms = [['collagen', 'extracellular', 'catabolic','wnt', 'mRNA','ecm', 'tnf', 'epithelial', 'Angiogenesis'], 
    #             ['mito','oxid', 'mTORC1', 'apop'],
    #             [' stabil','telomer', 'antigen', 'immune']]
    # sortbys = [['slow'], ['fresh'], ['vitri']]  # 
    # lgDir = os.path.join(gsea_ct_res_dir, 'leading_genes')
    # if not os.path.exists(lgDir):
    #     os.makedirs(lgDir)
    # adata = adata_full[adata_full.obs[celltype_annot]==ct, :] if ct != 'all' else adata_full
    # gs = adata.obs[gsea_col].unique()
    # for g in gs:
    #     rest_gs = [x for x in adata.obs[gsea_col].unique() if x!=g]
    #     for rest_g in rest_gs:
    #         combo = '-'.join([g, rest_g]) 
    #         sc.tl.rank_genes_groups(adata, groupby = gsea_col, 
    #                                 method='wilcoxon', n_genes=1000, pts=True,
    #                                 reference=rest_g, groups=[g], key_added=combo)
    # for i,terms in enumerate(allTerms):
    #     for term in terms:
    #         # extract result terms containing a given keyword and compile the leading genes by batch
    #         sortby = sortbys[i]
    #         term_res = posres[posres.Term.str.contains(term, case=False)]
    #         if term_res.shape[0]>1:
    #             term_LGs = term_res.groupby('batch').apply(lambda d: 
    #                 (pd.Series([v for x in d['Lead_genes'].values for v in x.split(';')]).value_counts()))  # get the leading genes for each combo
    #             term_LGs = pd.DataFrame(term_LGs)
    #             for combo in term_LGs.index.get_level_values(0).unique():  # append terms for each combo and each leading gene
    #                 comboRes = term_res[term_res['batch']==combo]
    #                 combo_deg = sc.get.rank_genes_groups_df(adata, key=combo, group=None)
    #                 combo_deg.to_csv(os.path.join(lgDir, f'{combo}_deg.csv'), sep=',')
    #                 if term_LGs.index.nlevels > 1: # if there is more than one combo
    #                     for gene in term_LGs.index.get_level_values(1):
    #                         t = '_'.join(comboRes[comboRes.Lead_genes.str.contains(gene)].Term.values)
    #                         term_LGs.loc[(combo, gene), 'terms'] = t
    #                         if gene in combo_deg.names.values:
    #                             for col in ['pvals_adj', 'scores', 'logfoldchanges']:
    #                                 term_LGs.loc[(combo, gene), col] = combo_deg[combo_deg.names==gene].logfoldchanges.values[0]
    #                 elif term_LGs.index.nlevels == 1:   
    #                     term_LGs = term_LGs.T
    #                     for gene in term_LGs.index:
    #                         t = '_'.join(comboRes[comboRes.Lead_genes.str.contains(gene)].Term.values)
    #                         term_LGs.loc[gene, 'terms'] = t
    #                         if gene in combo_deg.names.values:
    #                             for col in ['pvals_adj', 'scores', 'logfoldchanges']:
    #                                 term_LGs.loc[gene, col] = combo_deg[combo_deg.names==gene].logfoldchanges.values[0]
    #                 term_res.to_csv(os.path.join(lgDir, f'{term}ResTerms.csv'), sep=',')
    #             term_LGs.to_csv(os.path.join(lgDir, f'ComboLGCountsTerms_{term}.csv'), sep=',')
    
