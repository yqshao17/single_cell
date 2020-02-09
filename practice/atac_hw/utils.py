import pandas as pd
from scipy.io import mmread, mmwrite
from scipy.sparse import csr_matrix
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
import os
import time
import scanpy as sc

def load_data(path, transpose=False):
    print("Loading  data ...")
    t0 = time.time()
    if os.path.isdir(path):
        df = read_mtx(path)
    elif os.path.isfile(path):
        df = read_csv(path)
    else:
        raise ValueError("File {} not exists".format(path))     
    if transpose:
        df = df.T
    print('Original data contains {} features x {} cells'.format(*df.shape))
    print("Finished loading takes {:.2f} min".format((time.time()-t0)/60))
    return df

def read_anndata(path):
    for filename in glob(path+'/*'):
        basename = os.path.basename(filename)
        if (('count' in basename) or ('matrix' in basename)) and ('mtx' in basename):
            count = mmread(filename)
        elif 'barcode' in basename or 'cell' in basename:
            cell_id = pd.read_csv(filename, sep='\t', header=None)[0].values
        elif 'gene' in basename or 'peak' in basename:
            feature = pd.read_csv(filename, sep='\t', header=None)[0].values
    adata = sc.AnnData(count, cell_id, feature)
    return adata

def read_mtx(path):
    for filename in glob(path+'/*'):
        basename = os.path.basename(filename)
        if (('count' in basename) or ('matrix' in basename)) and ('mtx' in basename):
            count = mmread(filename).todense()
        elif 'barcode' in basename or 'cell' in basename:
            cell_id = pd.read_csv(filename, sep='\t', header=None)[0].values
        elif 'gene' in basename or 'peak' in basename:
            feature = pd.read_csv(filename, sep='\t', header=None)[0].values
    df = pd.DataFrame(count, feature, cell_id)
    return df

def read_csv(path):
    if ('.txt' in path) or ('tsv' in path):
        sep = '\t'
    elif '.csv' in path:
        sep = ','
    else:
        raise ValueError("File {} not in format txt or csv".format(path))
    df = pd.read_csv(path, sep=sep, index_col=0).astype('float32')
    return df


