# import necessary modules
import requests
import os
import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
try:  # unnecessary if the environment does not have it 
    import gseapy as gp
    import spatialdm as sdm
    import spatialdm.plottings as pl
    from PIL import Image
    from wordcloud import WordCloud
    import anndata2ri
    from rpy2.robjects.packages import importr
    from rpy2.robjects import r, pandas2ri
    anndata2ri.activate()
    pandas2ri.activate()
except:
    pass
import itertools
import pickle
import scipy.sparse as sp
import numpy as np
import time
from PIL import Image
import re
import json


def plot_tiles(v, k, annot_col, main_cli, bin_figure_dir='', top_n=10, thres=20, 
               res='hires', alpha=[0, 0.6], subset_to_cli=False, tile_size=[1000, 1000], show=True, **kwargs):
    """
    Plot tiles based on the given data and parameters.

    Parameters:
        v (adata): adata containing spot locations.
        k (str): Name of the sample.
        annot_col (str): Name of the annotation column.
        main_cli (str): Name of the main cluster of interest.
        bin_figure_dir (str): Path to the directory where the figures will be saved. If not specified, default to scanpy's default current_dir+'figures' directory.
        top_n (int): Number of tiles to plot.
        thres (int): Threshold for the number of cells in a tile.
        res (str): Resolution of the image. Default is 'hires'.
        a_list (list): List of alpha values for the plots. Default is [0, 0.6].
        main_cli (list): List of main CLI. Default is [].
        tile_size (list): List of integers for the size of each tile for cropping the image into tiles.
        show (bool): Whether to show the plot or not. Default is True.
        kwargs: Additional arguments to be passed to scanpy's pl.spatial() function.

    Returns:
        None
    """
    default_figdir = os.getcwd()+'/figures'
    showbin_figure_dir = bin_figure_dir if bin_figure_dir[:4]=='show' else 'show'+bin_figure_dir
    full_figdir = os.path.join(default_figdir, showbin_figure_dir)
    # Create the directory if it doesn't exist
    if not os.path.exists(full_figdir):
        os.makedirs(full_figdir)
    # Select the cells that belong to the main CLI
    X = v.obsm['spatial'][v.obs[annot_col]==main_cli]
    # Count the occurrences of each tuple
    tile_w, tile_h = tile_size[0], tile_size[1]
    entry_counts = Counter(list(zip(X[:, 0]//tile_w, X[:, 1]//tile_h)))
    # Get the unique entries and their counts
    unique_entries = list(entry_counts.keys())
    entry_counts = list(entry_counts.values())
    # Keep only the tiles with more than the specified number of cells
    l = [[entry[0], entry[1], c] for entry, c in zip(unique_entries, entry_counts) if c > thres]
    l = pd.DataFrame(l)
    l.columns = ['entry_x', 'entry_y', 'count']
    l = l.sort_values('count', ascending=False)
    l.loc[:,['coord_xmin', 'coord_xmax', 'coord_ymin', 'coord_ymax']] = None
    l = l.iloc[:top_n, :]
    print('tile cell counts', '\n', l)
    # Plot the top n tiles
    for i, r in l.iterrows():
        entry_x, entry_y, count = r['entry_x'], r['entry_y'], r['count']
        entry = [entry_x, entry_y]
        coord = [entry[0]*tile_w, (entry[0]+1)*tile_w, entry[1]*tile_h, (entry[1]+1)*tile_h]
        coord_str = '_'.join([str(s) for s in entry])
        alpha = [alpha] if isinstance(alpha, float) else alpha
        x_filt = (v.obsm['spatial'][:,0]>=coord[0]) & (v.obsm['spatial'][:,0]<=coord[1])
        y_filt = (v.obsm['spatial'][:,1]>=coord[2]) & (v.obsm['spatial'][:,1]<=coord[3])
        sub_v = v[x_filt & y_filt, ]
        if not show and bin_figure_dir!='':
            for a in alpha:
                sc.pl.spatial(sub_v, img_key=res, 
                                color=annot_col, 
                                crop_coord=coord, 
                                alpha=a, 
                                save=os.path.join(bin_figure_dir, f'{k}_{coord_str}_a{a}.png'), 
                                **kwargs)
        elif show or bin_figure_dir=='':
            for a in alpha:
                sc.pl.spatial(sub_v, img_key=res, 
                                color=annot_col, 
                                crop_coord=coord, 
                                alpha=a, 
                                **kwargs)
        if subset_to_cli and bin_figure_dir!='':
            a = [x for x in alpha if x>0][0]
            sc.pl.spatial(sub_v, img_key=res, 
                                color=annot_col, 
                                group=main_cli,
                                crop_coord=coord, 
                                alpha=a, 
                                save=os.path.join(bin_figure_dir, f'{k}_{coord_str}_a{a}_maincli.png'), 
                                **kwargs)
        l.loc[i,['coord_xmin', 'coord_xmax', 'coord_ymin', 'coord_ymax']] = coord
    return l

def find_tiles(v, annot_col, main_cli, top_n=10, thres=20, 
               tile_size=[1000, 1000]):
    """
    Find tiles of adata based on the given annotation column and parameters.

    Parameters:
        v (adata): adata containing spot locations.
        annot_col (str): Name of the annotation column.
        main_cli (str): Name of the main cluster of interest.
        top_n (int): Number of tiles to plot.
        thres (int): Threshold for the number of cells in a tile.
        tile_size (list): List of integers for the size of each tile for cropping the image into tiles.

    Returns:
        The list of tiles with columns ['entry_x', 'entry_y', 'count', 'coord_xmin', 'coord_xmax', 'coord_ymin', 'coord_ymax'].
    """
    # Select the cells that belong to the main CLI
    X = v.obsm['spatial'][v.obs[annot_col]==main_cli]

    # Count the occurrences of each tuple
    tile_w, tile_h = tile_size[0], tile_size[1]
    entry_counts = Counter(list(zip(X[:, 0]//tile_w, X[:, 1]//tile_h)))

    # Get the unique entries and their counts
    unique_entries = list(entry_counts.keys())
    entry_counts = list(entry_counts.values())

    # Keep only the tiles with more than the specified number of cells
    l = [[entry[0], entry[1], c] for entry, c in zip(unique_entries, entry_counts) if c > thres]
    l = pd.DataFrame(l)
    if l.shape[0]==0:
        return l
    l.columns = ['entry_x', 'entry_y', 'count']
    l = l.sort_values('count', ascending=False)
    l.loc[:,['coord_xmin', 'coord_xmax', 'coord_ymin', 'coord_ymax']] = None
    l = l.iloc[:top_n, :]

    # Plot the top n tiles
    for i, r in l.iterrows():
        entry_x, entry_y, count = r['entry_x'], r['entry_y'], r['count']
        entry = [entry_x, entry_y]
        coord = [entry[0]*tile_w, (entry[0]+1)*tile_w, entry[1]*tile_h, (entry[1]+1)*tile_h]
        l.loc[i,['coord_xmin', 'coord_xmax', 'coord_ymin', 'coord_ymax']] = coord
    return l
            
def get_hpa_columns(q, fuzzy=True):
    column_dict = json.load(open('/home/chrissy1/spatial/stomics/ovary_froz/hpa_cols.json', 'r'))
    if q in column_dict.keys():
        return [column_dict[q]]
    elif fuzzy:
        return [column_dict[k] for k in column_dict.keys() if re.search(q, k, re.IGNORECASE)]
    else:
        return [column_dict[k] for k in column_dict.keys() if q in k]
    
def run_sctransform(adata, layer=None, kwargs_dict=dict(), slot_for_X='data'):
    if layer:
        mat = adata.layers[layer]
    else:
        mat = adata.X
    # Set names for the input matrix
    cell_names = adata.obs_names
    gene_names = adata.var_names
    r.assign('mat', mat.T)
    r.assign('cell_names', cell_names)
    r.assign('gene_names', gene_names)
    r('colnames(mat) <- cell_names')
    r('rownames(mat) <- gene_names')
    seurat = importr('Seurat')
    r('seurat_obj <- CreateSeuratObject(mat)')
    # Run
    for k, v in kwargs_dict.items():
        r.assign(k, v)
    kwargs_str = ','.join([f'{k}={v}' for k,v in kwargs_dict.items()])
    r(f'seurat_obj <- SCTransform(seurat_obj, {kwargs_str})') # return.only.var.genes=FALSE, 
    
    # save data and scale.data
    if 'return.only.var.genes' not in kwargs_dict.keys():
        kwargs_dict['return.only.var.genes'] = 'FALSE'
    else:
        kwargs_dict['return.only.var.genes'] = 'TRUE'
        
    for slot in ['counts', 'data']:
        sct_data = np.asarray(r['as.matrix'](r(f'seurat_obj@assays$SCT@{slot}'))).T
        sct_data = sp.csr_matrix(sct_data)
        adata.uns[f'SCT_{slot}'] = sct_data

    adata.uns['SCT_var_genes'] = r('rownames(seurat_obj@assays$SCT@scale.data)')
    adata.uns['all_genes'] = r(f'rownames(seurat_obj@assays$SCT@counts)')
    inds = [i for i,x in enumerate(adata.uns['all_genes']) if x in adata.uns['SCT_var_genes']]
    adata = adata[:, inds]  # overwrite the raw counts with the SCT (log transformed) data
    adata.layers['SCT_scale.data'] = sp.csr_matrix(np.asarray(r['as.matrix'](r(f'seurat_obj@assays$SCT@scale.data'))).T)

    del sct_data
    adata.X = adata.uns[f'SCT_{slot_for_X}'] if slot_for_X!='scale.data' else adata.layers['SCT_scale.data']
    
    return adata

def add_spatial_img(adata, bins, highres_scale, lowres_scale, im_lowres, im_hires):
    """
    Customizes the data_list object by modifying specific attributes based on the provided parameters.
    Parameters:
    - adata (object): The adata object to be customized.
    - bins (int): The value for the 'spot_diameter_fullres' scale factor.
    - highres_scale (float): The value for the 'tissue_hires_scalef' scale factor.
    - lowres_scale (float): The value for the 'tissue_lowres_scalef' scale factor.
    - im_lowres (dict): A dictionary containing low-resolution images.
    - im_hires (dict): A dictionary containing high-resolution images.
    Returns:
    - None
    """
    Image.MAX_IMAGE_PIXELS = 2000000000
    if 'spatial' not in adata.uns.keys():
        adata.uns['spatial'] = dict()
    if len(im_lowres.keys()) > 1:
        print('Multiple images found. Use the first key accessed.')
    if type(im_lowres) != dict:
        print('The images are not in a dictionary format. Will use 0 as key.')
        key = 0
    else:
        key = list(im_lowres.keys())[0]
    adata.uns['spatial'][key] = dict()
    adata.uns['spatial'][key]['images'] = dict()
    adata.uns['spatial'][key]['scalefactors'] = dict()
    adata.uns['spatial'][key]['scalefactors'] = {'spot_diameter_fullres': bins,
                                                          'tissue_hires_scalef': highres_scale,
                                                          'fiducial_diameter_fullres': 100,
                                                          'tissue_lowres_scalef': lowres_scale}
    adata.uns['spatial'][key]['images']['lowres'] = np.array(im_lowres[key])
    adata.uns['spatial'][key]['images']['hires'] = np.array(im_hires[key])

def scrape_mks(df, df_gene_col='names', std_cols='g', extra_cols=['single cell'], fuzzy=True, outfile=None):
    """
    Scrape the Human Protein Atlas for information on the genes in the dataframe.
    df: pandas dataframe
        The dataframe containing the gene names to be scraped
    df_gene_col: str
        The name of the column in the dataframe containing the gene names
    std_cols: str
        The standard columns to be scraped from the HPA
    extra_cols: list
        Additional columns to be scraped from the HPA. Default is ['single cell'] which fetch single cell related column names from HPA
    fuzzy: bool
        Whether to use fuzzy matching for extra columns.
    outfile: str
        The path to the file to save the scraped data to. If None, the data is not saved.        
    """
    if extra_cols:
        celltype_query_columns = []
        for celltype in extra_cols:
            query_key = get_hpa_columns(celltype, fuzzy=fuzzy)
            celltype_query_columns += query_key
        celltype_query_columns = ','.join(celltype_query_columns)
    else:
        celltype_query_columns = ''
    full_query_columns = ','.join([std_cols, celltype_query_columns])
    url = 'https://www.proteinatlas.org/api/search_download.php'
    for (i, gene) in enumerate(df[df_gene_col].unique()):
        print(i)
        params = {
            'search': gene,
            'format': 'json',
            'columns': full_query_columns, 
            'compress': 'no'
        }
        result = requests.get(url.format(gene=gene), params=params)
        data = result.json()
        if len(data) == 0:
            print(f'No data found for {gene}')
        for item in data:
            retrieved_gene = item['Gene']
            if retrieved_gene==gene:  # check if the gene name matches the query completely
                for k in item.keys():
                    if type(item[k]) == list:
                        item[k] = '|'.join(item[k])
                    elif type(item[k]) == dict:
                        kv_pairs = ['-'.join([k, str(v)]) for k, v in item[k].items()]
                        item[k] = '|'.join(kv_pairs)
                    if k not in df.columns:
                        if type(item[k]) == str:
                            df[k] = ''
                        elif type(item[k]) == float:
                            df[k] = 0.0
                        elif type(item[k]) == int:
                            df[k] = 0
                    df.loc[df[df_gene_col] == gene, k] = item[k]                
        time.sleep(0.2) 
    if outfile:
        df.to_csv(outfile, sep=',', index=False)
    return df

def gsea_np2_combo(celltypes, adata, celltype_annot_col, gsea_col, gsea_dir, extra_info_abt_dataset, norm='sct'):
    for celltype in celltypes:
        if len(celltype)==0:
            celltype = list(adata.obs[celltype_annot_col].unique())
        # Filter to celltype of interest
        ad = adata[[True if x in celltype else False for x in adata.obs[celltype_annot_col]]].copy() if len(celltype)>0 else adata.copy()

        ## Gsea
        gene_set_libraries = gp.get_library_name()
        SELECTED_GS = ['GO_Biological_Process_2023', 'MSigDB_Hallmark_2020', 'KEGG_2021_Human', 'Reactome_2022']

        for grouping in gsea_col: # 
            print('grouping is: ', grouping)
            gs = ad.obs[grouping].unique()
            # gs = list(itertools.combinations(gs, 2))
            for g in gs:
                # g = list(g)
                rest_gs = [x for x in ad.obs[grouping].unique() if x!=g]
                for rest_g in rest_gs:
                    sub_ad = ad[ad.obs[grouping].isin([g,rest_g])]
                    matrix = sub_ad.to_df().T
                    design = ['pos' if x == g else 'neg' for x in sub_ad.obs[grouping]]
                    result = gp.gsea(data=matrix, # or data='./P53_resampling_data.txt'
                                        gene_sets=SELECTED_GS,
                                        cls=design,
                                        permutation_type='phenotype',
                                        pheno_pos='pos',
                                        pheno_neg='neg',
                                        permutation_num=1000, # reduce number to speed up test
                                        outdir=None,  # do not write output to disk
                                        method='signal_to_noise',
                                        threads=10, seed=7)
                    combo = '-'.join([g, rest_g]) 
                    celltype_pickle_name = '-'.join(celltype) if isinstance(celltype, list) else celltype
                    dataset_gsea_dir = os.path.join(gsea_dir, extra_info_abt_dataset) if not gsea_dir.endswith(extra_info_abt_dataset) else gsea_dir
                    if not os.path.exists(dataset_gsea_dir):
                        os.mkdir(dataset_gsea_dir)
                    names = [x.replace('_', '-') for x in [grouping,celltype_pickle_name,norm]]    
                    fullPklName = '_'.join(names)
                    with open(os.path.join(dataset_gsea_dir, f'{fullPklName}.pkl'), 'wb') as f: 
                        pickle.dump(result, f, protocol=pickle.HIGHEST_PROTOCOL)    
                            
                            
def gsea_default(celltypes, adata, celltype_annot_col, gsea_col, gsea_dir, extra_info_abt_dataset, norm='sct'):
    for celltype in celltypes:
        # Filter to celltype of interest
        ad = adata[[True if x in celltype else False for x in adata.obs[celltype_annot_col]]].copy() if len(celltype)>0 else adata.copy()

        ## Gsea
        gene_set_libraries = gp.get_library_name()
        SELECTED_GS = ['GO_Biological_Process_2023', 'MSigDB_Hallmark_2020', 'KEGG_2021_Human', 'Reactome_2022']

        for grouping in gsea_col: # 
            print('grouping is: ', grouping)
            gs = ad.obs[grouping].unique()
            # gs = list(itertools.combinations(gs, 2))
            for g in gs:
                matrix = ad.to_df().T
                design = ['pos' if x == g else 'neg' for x in ad.obs[grouping]]
                result = gp.gsea(data=matrix, # or data='./P53_resampling_data.txt'
                                    gene_sets=SELECTED_GS,
                                    cls=design,
                                    permutation_type='phenotype',
                                    pheno_pos='pos',
                                    pheno_neg='neg',
                                    permutation_num=1000, # reduce number to speed up test
                                    outdir=None,  # do not write output to disk
                                    method='signal_to_noise',
                                    threads=10, seed=7)
                group_name = '-'.join([grouping, g, 'vsRest']) 
                celltype_pickle_name = '-'.join(celltype) if isinstance(celltype, list) else celltype
                dataset_gsea_dir = os.path.join(gsea_dir, extra_info_abt_dataset) if not gsea_dir.endswith(extra_info_abt_dataset) else gsea_dir
                if not os.path.exists(dataset_gsea_dir):
                    os.mkdir(dataset_gsea_dir)
                names = [x.replace('_', '-') for x in [group_name,celltype_pickle_name,norm]]    
                fullPklName = '_'.join(names)
                with open(os.path.join(dataset_gsea_dir, f'{fullPklName}.pkl'), 'wb') as f: 
                    pickle.dump(result, f, protocol=pickle.HIGHEST_PROTOCOL)    


def find_files_with_suffix(directory, suffix):
    file_paths = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith(suffix):
                file_path = os.path.join(root, file)
                file_paths.append(file_path)
    return file_paths

def generate_word_cloud(text_list, exclude_words, phrase_length, **kwargs):
    # Join the list of strings into a single text
    d = {}
    for t in text_list:
        for w in exclude_words:
            t = t.replace(w, '')
        t = ' '.join(t.split(' ')[:phrase_length])    
        d[t] = d[t] + 1 if t in d.keys() else 1
    # Generate the word cloud
    cloud = WordCloud(collocations=False, background_color='white').generate_from_frequencies(d)    
    # Display the word cloud using matplotlib
    # plot the WordCloud image                        
    plt.figure(**kwargs, facecolor = None)
    plt.imshow(cloud, interpolation='bilinear')
    plt.axis('off')
    plt.show()


