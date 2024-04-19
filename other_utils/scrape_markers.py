
# Scrape markers from HPA
import pandas as pd
import numpy as np
import requests
import json
import time
import pickle
# import stereo as st
import os
import importlib

file_path = '/home/chrissy1/spatial/stomics/ovary_froz/spatial_utils.py'
spec = importlib.util.spec_from_file_location("spatial_utils", file_path)
spatial_utils = importlib.util.module_from_spec(spec)
spec.loader.exec_module(spatial_utils)

out = '/home/chrissy1/spatial/stomics/ovary_froz/redo/seurat/bin50_processed'

# markers = pd.read_csv(out+'/bin50_de_markers_res_1.csv')
# markers = markers[(markers['p_val']<0.05) & (markers['pct.1']>0.1) & (markers['avg_log2FC']>0)]
# markers = markers[~markers['gene'].str.contains('MT-') & ~markers['gene'].str.contains('RPL') & ~markers['gene'].str.contains('RPS')]
# markers = markers.groupby('cluster').apply(lambda x: x.sort_values('avg_log2FC', ascending=False).head(100))
# scrape_mks(markers, 
#             df_gene_col = 'gene', 
#             outfile = out+'/bin50_de_markers_res_1_scraped.csv')

# markers = pd.read_csv(out+'/bin20_de_markers_res_1.csv')
# markers = markers[(markers['p_val']<0.05) & (markers['pct.1']>0.1) & (markers['avg_log2FC']>0)]
# markers = markers[~markers['gene'].str.contains('MT-') & ~markers['gene'].str.contains('RPL') & ~markers['gene'].str.contains('RPS')]
# markers = markers.groupby('cluster').apply(lambda x: x.sort_values('avg_log2FC', ascending=False).head(2)).head(10)
# markers = spatial_utils.scrape_mks(markers, 
#             df_gene_col='gene', 
#             outfile=out+'/bin20_de_markers_res_1_scraped.csv')


## DEG markers
savedir = f'{out}/deg'
os.chdir(savedir)
csv = 'sample_markers.csv'
markers = pd.read_csv(csv)
markers = markers[(markers['pvals_adj']<0.05) & (markers['logfoldchanges']>0)]
markers = markers[~markers['names'].str.contains('MT-') & ~markers['names'].str.contains('RPL') & ~markers['names'].str.contains('RPS')]
markers = markers.groupby('group').apply(lambda x: x.sort_values('logfoldchanges', ascending=False).head(100))
spatial_utils.scrape_mks(markers, 
            df_gene_col = 'names', 
            extra_cols=['Molecular function', 'Disease', 'Evidence'],
            outfile = csv.replace('.csv', '_scraped.csv'))