import pandas as pd
import scanpy as sc
import numpy as np

def normalize_rna(rna_adata,n_top_genes):
### normalize and filter data
        sc.pp.filter_cells(rna_adata, min_genes=100)
        
        sc.pp.normalize_per_cell(rna_adata, counts_per_cell_after=1e4)
        sc.pp.log1p(rna_adata)
        #rna_adata.raw = rna_adata
        sc.pp.highly_variable_genes(rna_adata, min_mean=0.0125, max_mean=3, min_disp=0.5, n_top_genes=n_top_genes)
        rna_adata = rna_adata[:, rna_adata.var['highly_variable']]
        print('Highly variable genes: ', rna_adata.shape[1])
        return rna_adata

def normalize_atac(atac_adata):
        sc.pp.filter_cells(atac_adata, min_genes=100)
        sc.pp.filter_genes(atac_adata, min_cells=10)
        #sc.pp.normalize_per_cell(atac_adata, counts_per_cell_after=1e4)
        #sc.pp.log1p(atac_adata)       
        return atac_adata

    
def overlap_adata(data1, data2):
    barcodes=list(set(data1.obs_names)&set(data2.obs_names))       
    data1=data1[barcodes, :]
    data2=data2[barcodes, :]
    return data1, data2

def overlap_feature(data1, data2):
    features=list(set(data1.var_names)&set(data2.var_names))       
    data1=data1[:, features]
    data2=data2[:, features]
    return data1, data2
    
def group_cells(data1, data2):   
    meta_cells=data1.obs['ClusterID'].unique()
    meta1=[]
    meta2=[]
    weight=[]
    for cluster in meta_cells:
        idx=np.where(data1.obs['ClusterID']==cluster)
        bc_set=data1.obs['ClusterID'].index[idx]
        try:
            meta1.append(data1.layers['norm_data'][idx].mean(axis=0))
        except:
            meta1.append(data1.X[idx].mean(axis=0))
        meta2.append(data2[bc_set,].X.mean(axis=0))
        weight.append(len(idx[0]))
    df1=pd.DataFrame(np.array(meta1), columns=data1.var_names, index=meta_cells)
    df2=pd.DataFrame(np.array(meta2), columns=data2.var_names, index=meta_cells)
    adata1=sc.AnnData(df1)
    adata2=sc.AnnData(df2)
    adata1.obs['Weights']=weight
    adata2.obs['Weights']=weight
    return adata1, adata2


def pseudo_pair_adata(data1, data2, CellID):
    adata1=data1[data1.obs['ClusterID']==CellID]
    adata2=data2[data2.obs['ClusterID']==CellID]
    return adata1, adata2

