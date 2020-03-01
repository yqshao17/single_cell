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
        if 'barcode' in basename or 'cell' in basename:
            cell_id = pd.read_csv(filename, sep='\t', header=None)[0].values
        elif ('gene' in basename) or ('peak' in basename) or ('feature' in basename):
            feature = pd.read_csv(filename, sep='\t', header=None)[0].values
        elif (('count' in basename) or ('matrix' in basename)) and ('mtx' in basename):
            count = mmread(filename)
    adata = sc.AnnData(count.T.tocsr())
    adata.obs_names=cell_id
    adata.var_names=feature
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


def write_mtx(path, adata):
    mmwrite(path+'/count.mtx', csr_matrix(adata.X.T, dtype=int))
    np.savetxt(path+'/barcodes.txt',adata.var_names,fmt="%s")
    np.savetxt(path+'/peaks.txt',adata.obs_names,fmt="%s")
    
def modi_peak(infile, outfile, inft, outft):
    # ':-': chr:start-end
    # '__' or '_': chr_start_end
    with open(infile) as input, open(outfile,'w') as output:
        for line in input:
            chr = line.split(inft[0])[0]
            s, e = line.split(inft[-1])[1]
            output.write(chr+outft[0]+s+outft[-1]+e)
            
#Creat new dir
#************************************************#
def Mkdir(path):
    path=path.strip()
    path=path.rstrip('\\')
    isExists=os.path.exists(path)
    if not isExists:
        os.makedirs(path)
        return True
    else:
        return False

def generate_cm(df, c1, c2):
    # generate confusion matrix
    # c1,c2 is column names of df
    # return a dataframe with index as c1 and columns as c2
    cm=df.groupby([c1,c2]).size()
    cm_df=pd.DataFrame(cm)
    cm_df.reset_index(inplace=True)
    new_df=cm_df.pivot(index=c1, columns=c2, values=0)
    new_df=new_df.fillna(0)
    return new_df