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
adata_full = sc.read_h5ad('SCT.h5ad')
# adata_full.obs = adata_full.obs[['Unnamed: 0', 'orig.ident', 'nCount_RNA', 'nFeature_RNA', 'X', 'batch', 'x', 'y', 'z', 'in_tissue', 'sample', 'percent.mt', 'nCount_SCT', 'nFeature_SCT', 'SCT_snn_res.1', 'seurat_clusters', 'PC_1', 'PC_2', 'PC_3', 'PC_4', 'PC_5', 'PC_6', 'PC_7', 'PC_8', 'PC_9', 'PC_10', 'PC_11', 'PC_12', 'PC_13', 'PC_14', 'PC_15', 'PC_16', 'PC_17', 'PC_18', 'PC_19', 'PC_20', 'umap_1', 'umap_2', 'harmony_1', 'harmony_2', 'harmony_3', 'harmony_4', 'harmony_5', 'harmony_6', 'harmony_7', 'harmony_8', 'harmony_9', 'harmony_10', 'harmony_11', 'harmony_12', 'harmony_13', 'harmony_14', 'harmony_15', 'harmony_16', 'harmony_17', 'harmony_18', 'harmony_19', 'harmony_20', 'harmony_21', 'harmony_22', 'harmony_23', 'harmony_24', 'harmony_25', 'harmony_26', 'harmony_27', 'harmony_28', 'harmony_29','harmony_30', 'annot', 'annot_no_marker']]

def run(adata_full, cl, savedir, top_gsea_terms, annot_col, groupby_col, comparison):
    adata = adata_full[adata_full.obs[annot_col]==cl, :] if cl != 'all' else adata_full  # subset adata to cell type
    batch = comparison.split('-')[0]
    # adata = adata_full
    gseadir = f'{savedir}/gsea/{cl}'
    pkl_files = [f for f in os.listdir(gseadir) if f.endswith('.pkl')]
    mod_score_dir = f'{gseadir}/module_scores'
    figdir = f'{mod_score_dir}/figures'
    if not os.path.exists(figdir):
        os.makedirs(figdir)
    os.chdir(mod_score_dir)
    # load gsea results
    if 'chosen_posres.pkl' in os.listdir(mod_score_dir):
        posres_comparisons = pkl.load(open(f'{mod_score_dir}/chosen_posres.pkl', 'rb'))
    else:
        lg_comparisons = {}
        posres_comparisons = {}
        ranking_comparisons = {}
        for one_comp in pkl_files:
            # one_comp = pkl_files[0]
            store_key = one_comp.split('.')[0]
            gs = pkl.load(open(os.path.join(gseadir, one_comp), 'rb'))
            # gs = pkl.load(open('/home/chrissy1/spatial/stomics/ovary_froz/redo/seurat/bin20_processed/gsea/all/Slow-Fresh.pkl', 'rb'))
            res = gs.res2d
            posres = res[(res['FDR q-val'] < 0.05) & (res['NES']>0)].sort_values(by='NES', ascending=False).head(top_gsea_terms)
            posres_comparisons[store_key] = posres
            lg = (';'.join(posres['Lead_genes'].tolist())).split(';')
            lg = pd.DataFrame(np.unique(lg, return_counts=True), index=['gene', 'count']).T.sort_values('count', ascending=False)
            lg_comparisons[store_key] = lg
            ranking_comparisons[store_key] = gs.ranking
        for k, posres in posres_comparisons.items():
            posres = posres[~posres.Term.str.startswith('KEGG')]  # remove KEGG pathways 
            posres['Term'] = posres['Term'].values
            posres['Term'] = posres['Term'].str.split('__').str[1]
            posres['Term'] = posres['Term'].str.split('\(GO').str[0]
            posres['Term'] = posres['Term'].str.split(' R-').str[0]
            posres['v1_Term'] = [row['Term']+' ('+','.join(row['Lead_genes'].split(';')[:5])+')' for i, row in posres.iterrows()]
            posres = posres.sort_values('Term')
            posres_comparisons[k] = posres
        with open(f'{mod_score_dir}/chosen_posres.pkl', 'wb') as f:
            pkl.dump(posres_comparisons, f)
    for k, posres in posres_comparisons.items():
        posres_comparisons[k]['Comparison'] = k
        posres_comparisons[k]['Term'] = [x.split(r' (')[0] for x in posres_comparisons[k]['Term'].values]
        posres_comparisons[k] = posres_comparisons[k].sort_values('FDR q-val')
    # cross sample gsea dotplot
    total_results = pd.concat([posres_comparisons[k] for k in posres_comparisons.keys()])
    ax = dotplot(total_results,
                column="NES",
                x='Comparison', # set x axis, so you could do a multi-sample/library comparsion
                size=3,
                top_term=10,
                figsize=(5,15),
                xticklabels_rot=45, # rotate xtick labels
                show_ring=True, # set to False to revmove outer ring
                marker='o', 
                ofname=f'{figdir}/gseapy_dotplot_cross_sample.pdf')
    if f"module_scores_key{cl}.csv" in os.listdir(mod_score_dir):  # load module scores if already computed, and check if all have been computed
        mod_scores = pd.read_csv(f"{mod_score_dir}/module_scores_key{cl}.csv", index_col=0)
        mod_scores.drop(columns='batch_annot', inplace=True, errors='ignore')
        mod_scores.index = adata.obs.index.values
        for col in mod_scores.columns:
            adata.obs[col] = mod_scores[col].values
        if comparison is not None:
            posres = posres_comparisons[comparison]
            print(posres.shape)
            for i,row in posres.iterrows():
                if row['v1_Term'] not in adata.obs.columns:
                    sc.tl.score_genes(adata, gene_list=row['Lead_genes'].split(';'), score_name=row['v1_Term'], use_raw = False)
                    mod_scores[row['v1_Term']] = adata.obs[row['v1_Term']].values
                    print('new term added', row['v1_Term'])
        else:
            for k, posres in posres_comparisons.items():
                for i,row in posres.iterrows():
                    if row['v1_Term'] not in adata.obs.columns:
                        sc.tl.score_genes(adata, gene_list=row['Lead_genes'].split(';'), score_name=row['v1_Term'], use_raw = False)
                        mod_scores[row['v1_Term']] = adata.obs[row['v1_Term']].values
        mod_scores.to_csv(f"{mod_score_dir}/module_scores_key{cl}.csv")
    else:
        if comparison is not None:
            posres = posres_comparisons[comparison]
            for i,row in posres.iterrows():
                if row['v1_Term'] not in adata.obs.columns:
                    sc.tl.score_genes(adata, gene_list=row['Lead_genes'].split(';'), score_name=row['v1_Term'], use_raw = False)
        else:
            for k, posres in posres_comparisons.items():
                for i,row in posres.iterrows():
                    if row['v1_Term'] not in adata.obs.columns:
                        sc.tl.score_genes(adata, gene_list=row['Lead_genes'].split(';'), score_name=row['v1_Term'], use_raw = False)
        adata.obs.loc[:, posres['v1_Term']].to_csv(f"{mod_score_dir}/module_scores_key{cl}.csv")
    # grouping = pd.read_csv('/home/chrissy1/spatial/stomics/ovary_froz/redo/seurat/Book.csv')
    # grouping = {row['Functional Category']:row['Gene Set Terms'].split(', ') for i,row in grouping.iterrows()}
    grouping = {
        'ECM organization': ['Nuclear', 'Tubulin',
                            'Extracellular', 'NCAM', 'Integrin', 'ECM'],
        'EMT': ['Epithelial'],
        'GAG': ['GAG', 'Glycosaminoglycan', 'Glycosulation'],
        'Collagen organization': ['Collagen'],
        'Angiogenesis': ['Angio'],
        # 'Vascular process': ['Fibrinolysis', 'Coagulation', 'Plasminogen'],
        'Ion-related process': ['Ion', 'Cation'],
        'Inflammation': ['TNF', 'IL-', 'Interleukin', 'Inflamm', 'Cytokine'],
        'Heat related': ['Heat', 'HSP'],	
        'Stress response': ['Stress', 'Oxidative'],
        'Myc targets': ['Myc'],
        'mRNA stability': ['mRNA Stabil'],
        'mRNA processes': ['Translation', 'Elongation', 'Splicing', 'mRNA Activation'],
        'Respiratory process': ['Mitochondrial', 'Oxydative', 'NADH', 'Respiration'],
        # 'Cell cycle': ['Ubiquitin', 'Synthesis', 'Phase'],
        'Other Signaling': ['Hedgehog', 'Reproductive'],
        'B cell signaling': ['Fc', 'B Cell'],
        'GF Signaling': ['EGF', 'FGF', 'PDGF', 'TGF', 'HGF', 'IGF', 'Insulin', 'Growth factor'],
    }
    posres_grouping = {}
    for k in grouping.keys():
        allterms = grouping[k]
        if comparison is not None:
            cols = posres_comparisons[comparison]['v1_Term'].values
        else:
            cols = [list(posres['v1_Term'].values) for k, posres in posres_comparisons.items()]
            cols = [item for sublist in cols for item in sublist]
        for term in allterms:
            terms_found = [x for x in cols if term in x]
            if k not in posres_grouping.keys():
                posres_grouping[k] = terms_found
            elif k in posres_grouping.keys():
                posres_grouping[k] = posres_grouping[k]+terms_found
    # topic_module_means = {}
    # for k, v in posres_grouping.items():
    #     adata.obs = adata.obs.T.drop_duplicates().T
    #     intersected_columns = [x for x in v if x in adata.obs.columns.values]
    #     print('unique cols only', len(adata.obs.columns.unique())==len(adata.obs.columns))
    #     print('levels', adata.obs.columns.nlevels)  
    #     topic_module_means[k] = adata.obs.groupby('batch_annot').apply(lambda x: x[intersected_columns].mean().mean())
    # res = pd.DataFrame(topic_module_means)
    # max_val_indices = res.apply(lambda s: pd.Series(s.nlargest(3).index))
    # print(res, max_val_indices)
    # with open(f'{mod_score_dir}/module_scores_means_top_ranking_{cl}_groupby_{groupby_col}_from{comparison}.pkl', 'wb') as f:
    #     d = {'res':res, 'max_moduleScore_celltype': max_val_indices}
    #     pkl.dump(d, f)
    ## plot leading genes as module scores
    # for k,posres in posres_comparisons.items():
    #     # try:
    #     # sc.pl.dotplot(adata[adata.obs['batch'].isin(ks),:], posres['Term'].values, groupby='batch_annot', dendrogram=False, save=f'{k}_modulescore_annot_top.pdf', swap_axes=True)
    #     posres = posres.sort_values('Lead_genes')
    #     sc.pl.dotplot(adata, posres['v1_Term'], groupby='batch_annot', dendrogram=False, swap_axes=False, save=f'_{k}_modulescore_annot_top')
    #     plt.savefig(f'figures/withTotal_{k}_modulescore_annot_top.pdf')
    lg_dict = {}
    for k,terms in posres_grouping.items():
        terms = [x for x in terms if x in adata.obs.columns]
        term_lg_dict = {term:posres_comparisons[comparison].loc[posres_comparisons[comparison]['v1_Term']==term, 'Lead_genes'].values[0].split(';') for term in terms}
        unique_lgs = np.unique([item for sublist in term_lg_dict.values() for item in sublist])
        lg_dict[k] = unique_lgs
        if len(terms)>0:
        #     adata.obs['batch_annot'] = [str(x) + '_' + str(adata.obs[groupby_col].values[i]) for i,x in enumerate(adata.obs['batch'])]
        #     sc.pl.dotplot(adata, terms, groupby='batch_annot', dendrogram=False, swap_axes=True, save=f'{k}_modulescore_{annot_col}_from{comparison}.pdf')
            sc.pl.dotplot(adata, terms, groupby='batch', dendrogram=False, swap_axes=True, use_raw=False, save=f'{k}_modulescore_{annot_col}_from{comparison}.pdf')
            # sc.pl.stacked_violin(adata, terms, groupby='batch', swap_axes=True, save=f'{k}_modulescore_by_{annot_col}_from{comparison}.pdf')
    if 'rank_genes_groups_batch.csv' not in os.listdir(mod_score_dir):
        sc.tl.rank_genes_groups(adata, groupby='batch', method='wilcoxon', use_raw=False, n_genes=1000, pts=True)
        markers = sc.get.rank_genes_groups_df(adata, group=None)
        markers.to_csv(f'{mod_score_dir}/rank_genes_groups_batch.csv')
    else: 
        markers = pd.read_csv(f'{mod_score_dir}/rank_genes_groups_batch.csv')
    markers = markers[(markers['pvals_adj']<0.05)&(markers['logfoldchanges']>0)].sort_values('scores', ascending=False)
    # for batch in markers['group'].unique():
    batch_markers = markers[markers['group']==batch]
    filtered_lgs = {k:[x for x in batch_markers['names'] if x in v][:6] for k,v in lg_dict.items()}
    k_to_del = [k for k,v in filtered_lgs.items() if len(v)==0]
    for k in k_to_del:
        del filtered_lgs[k]
    sc.pl.dotplot(adata, var_names=filtered_lgs, groupby='batch', dendrogram=False, swap_axes=True, save=f'lg_{annot_col}_from{comparison}_batchwise_{batch}sample.pdf', use_raw=False)


top_gsea_terms = 50
annot_col = 'annot_no_marker'
groupby_col = 'annot_no_marker'
comparison = 'Slow-Fresh'
cl = 'ST_SM'
cls = [x for x in os.listdir(f'{savedir}/gsea') if '(' in x]
for cl in ['ST_SM']:  # 'all', 'ST', 
    print(f'PLOTTING FOR {cl}')
    run(adata_full, cl, savedir, top_gsea_terms=top_gsea_terms, annot_col=annot_col, groupby_col=groupby_col, comparison=comparison)