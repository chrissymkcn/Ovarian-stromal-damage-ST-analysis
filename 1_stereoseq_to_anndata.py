import stereo as st
import warnings
warnings.filterwarnings('ignore')
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pickle
import os
from scipy import stats
import copy
import matplotlib.pyplot as plt
from natsort import natsorted
import sys
import dill
from stereo.core.ms_data import MSData
from stereo.core.ms_pipeline import slice_generator
from copy import deepcopy
import dill
import importlib.util
from PIL import Image
    
def saveGef(obj, filename='obj'):
    with open(f'{filename}.pkl', 'wb') as f:
        dill.dump(obj, f)

def loadGef(filename):
    with open(filename+'.pkl', 'rb') as f:
        d = dill.load(f)
    return d


def GetMarkers(msdata, group=None, topkey=None, key='marker_genes', ref='rest'):
    if topkey:
        result_toplevel_key = topkey
    else:
        result_toplevel_key = list(msdata.tl.result.keys())
        if len(result_toplevel_key) == 1:
            result_toplevel_key = result_toplevel_key[0]  
        else: 
            result_toplevel_key = result_toplevel_key[0]
            # raise('Multiple toplevel keys found in data.tl.result')
    if key not in msdata.tl.result[result_toplevel_key].keys():
        raise KeyError(f"Key '{key}' not found in data.tl.result")
    if group != None:
        if type(group) != str:
            group = str(group)
        cluster = group+'.vs.'+ref
        markers = msdata.tl.result[result_toplevel_key][key][cluster]
        markers['group'] = group
    elif group == None:
        mks_dict = msdata.tl.result[result_toplevel_key][key].copy()
        mks_dict.pop('parameters', None)
        for k in mks_dict:
            mks_dict[k]['group'] = k.split('.')[0]
        markers = pd.concat(mks_dict.values(), axis=0)
    return markers


def find_files(dir_path, file_ext):
    result = []
    for root, dirs, files in os.walk(dir_path):
        for file in files:
            if file.endswith(file_ext):
                result.append(os.path.join(root, file))
    return result


def find_bin_files(directory, key):
    return find_files(directory, key+'.gef')

def mergeGefs(gef_list):
    merged = copy.deepcopy(gef_list[0])    
    # for i in range(1, len(gef_list)):
    #     merged = st.utils.data_helper.merge(merged, gef_list[i])
    merged = st.utils.data_helper.merge(merged, *gef_list[1:])
    return merged

datasetdir = '/home/chrissy1/spatial/stomics/ovary_froz/datasets/'
output_dir = '/home/chrissy1/spatial/stomics/ovary_froz/redo/bidcell'

files = find_bin_files(datasetdir, 'tissue')

binsizes = [20, 50]  # 1, 10, 20, 50, 100, 150, 200
for binsize in binsizes:
    data_list = {path.split('/')[-2]: st.io.read_gef(file_path=path, bin_size=binsize) for path in files}

    ms_data = MSData(_relationship='other', _var_type='intersect')
    for k,sample in data_list.items():
        ms_data += sample
    ms_data.integrate()
    ms_data.tl.raw_checkpoint()
    data = deepcopy(ms_data.merged_data)

    adata = st.io.stereo_to_anndata(data,flavor='seurat',output=f'{output_dir}/binsize_{binsize}.h5ad', split_batches=False)
    adata.obs['batch'] = adata.obs['batch'].replace(dict(adata.uns['sn'].values))
    
    # Load the spatial_utils module
    file_path = '/home/chrissy1/spatial/stomics/ovary_froz/spatial_utils.py'
    spec = importlib.util.spec_from_file_location("spatial_utils", file_path)
    spatial_utils = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(spatial_utils)
    
    # load the images
    file_paths = ['/home/chrissy1/spatial/stomics/ovary_froz/rawdata/E-SCT20221109001_SZD20230526112_1/Result/Fresh/01.StandardWorkflow_Result/Register/D01154B5_regist.tif',
                    '/home/chrissy1/spatial/stomics/ovary_froz/rawdata/E-SCT20221109001_SZD20230526112_2/Result/Slow/01.StandardWorkflow_Result/Register/D01970B5_regist.tif',
                    '/home/chrissy1/spatial/stomics/ovary_froz/rawdata/E-SCT20221109001_SZD20230526112_2/Result/Vitri/01.StandardWorkflow_Result/Register/B01320B6_regist.tif']

    tissue_cut_paths = ['/home/chrissy1/spatial/stomics/ovary_froz/rawdata/E-SCT20221109001_SZD20230526112_1/Result/Fresh/01.StandardWorkflow_Result/Register/D01154B5_tissue_cut.tif',
                    '/home/chrissy1/spatial/stomics/ovary_froz/rawdata/E-SCT20221109001_SZD20230526112_2/Result/Slow/01.StandardWorkflow_Result/Register/D01970B5_tissue_cut.tif',
                    '/home/chrissy1/spatial/stomics/ovary_froz/rawdata/E-SCT20221109001_SZD20230526112_2/Result/Vitri/01.StandardWorkflow_Result/Register/B01320B6_tissue_cut.tif']

    Image.MAX_IMAGE_PIXELS = 2000000000
    highres_scale = 0.2
    lowres_scale = 0.1
    batch = 'batch'
    im_hires = {}
    im_lowres = {}
    im_tis_mask = {}
    for i,impath in enumerate(file_paths):
        im = Image.open(impath)
        sample = impath.split('/')[-4]
        hi_dim = (int(im.size[0]//(1/highres_scale)), int(im.size[1]//(1/highres_scale)))
        low_dim = (int(im.size[0]//(1/lowres_scale)), int(im.size[1]//(1/lowres_scale)))
        im_hires[sample] = im.resize(hi_dim)
        im_lowres[sample] = im.resize(low_dim)
        im_tis_mask[sample] = np.array(Image.open(tissue_cut_paths[i]))
    for k in adata.obs[batch].unique():
        spatial_utils.add_spatial_img(adata=adata, key=k, bins=binsize, highres_scale=highres_scale, lowres_scale=lowres_scale, im_lowres=im_lowres, im_hires=im_hires)
        mask = [im_tis_mask[row[0], row[1]] for row in adata.obsm['spatial']]
        adata.obs['in_tissue'] = mask
    adata.write_h5ad(f'{output_dir}/binsize_{binsize}.h5ad')