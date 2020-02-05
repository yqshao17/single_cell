# group cells by RNA
# get regulatory matrix for each pseudocell
# merge matrix of all pseudocells

import pickle
import sys
sys.path.insert(1, '/home/xionglei/yanqiu')
from utils import load_data, write_mtx
from gene_peak_fit import get_gene_peak_fit, pseudo_regulatory_matrix
#from plot import plot_info,plot_gene_peak
import os
import pandas as pd
import scanpy as sc
from plot import plot_info, plot_gene_peak
from process_data import normalize_rna,normalize_atac,pseudo_pair_adata,overlap_adata
from multiprocessing import Process

import argparse
parser = argparse.ArgumentParser()

parser.add_argument('--outdir','-o', help='output dir')
parser.add_argument('--nthreads','-t', help='nthreads', default=20)
parser.add_argument('--rnadir', '-r', help='rna directory or h5ad file')
parser.add_argument('--atacdir', '-a', help='atac directory or h5ad file')
parser.add_argument('--method', '-m', help='regress, pearson or spearman', default='regress')
parser.add_argument('--species', '-s', help='mm, hg or mm_hg')
parser.add_argument('--bin_size', '-bs', help='bin size: *k', default=250)
parser.add_argument('--n_top_genes', '-gn', help='number of filtered variable genes', default=None)
parser.add_argument('--group_by','-g', help='group cells by RNA, ATAC or None', default='RNA')


args = parser.parse_args()

outdir=args.outdir
nthreads=int(args.nthreads)
rna_dir=args.rnadir
atac_dir=args.atacdir
species=args.species
method=args.method
bin_size=int(args.bin_size)
n_top_genes=args.n_top_genes

if n_top_genes:
    n_top_genes=int(n_top_genes)
    
group_by=args.group_by
    
    
if __name__=='__main__':
    
    if not os.path.exists(outdir):
        os.makedirs(outdir)   
    
    if rna_dir.endswith('h5ad'):
        rn=sc.read_h5ad(rna_dir)    
    elif not rna_dir.endswith('loom'):
        rna=load_data(rna_dir)
        rna_adata=sc.AnnData(rna.T)
        rn=normalize_rna(rna_adata, n_top_genes)    
    
    if atac_dir.endswith('h5ad'):
        an=sc.read_h5ad(atac_dir)
    else:
        atac=load_data(atac_dir)   
        atac_adata=sc.AnnData(atac.T)
        an=normalize_atac(atac_adata)
        
    '''
    if rna_dir.endswith('loom'):
        pass
    else:
        os.system('Rscript group_cells.R %s %s'%(rna_dir, outdir))
    meta_rn=sc.read_loom(outdir+'/group_labeled.loom', sparse=False, obs_names='CellID', var_names='Gene', dtype='float32')
    
    rn, an=overlap_adata(meta_rn, an)
    an.obs['ClusterID']=rn.obs['ClusterID']
    print('Overlap cells: %d'%(rn.shape[0]))
    print('RNA genes: %d'%(rn.shape[1]))
    print('ATAC peaks: %d'%(an.shape[1]))    
    rn.write(outdir+'/rn.h5ad')
    an.write(outdir+'/an.h5ad')
    '''
    
    an.obs['ClusterID']=rn.obs['ClusterID']
    gene_info = pickle.load(open('/home/xionglei/yanqiu/lib/gene_info_%s.pkl'%species, 'rb'))
            
    
    meta_cells=rn.obs['ClusterID'].unique()

    reg_df=pd.DataFrame()
    
    
    def pseudo_reg(rn, an, cluster, outdir, gene_info):
        outdir=outdir+'/pseudo_cells'
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        adata1, adata2 = pseudo_pair_adata(rn, an, cluster)
        sc.pp.filter_genes(adata1, min_cells=3)
        sc.pp.filter_genes(adata2, min_cells=3)
        # calculate gene-peak regression    
        gene_peak_fit = get_gene_peak_fit(gene_info, adata1, adata2, outdir+'/c'+str(cluster)+'_', bin_size, method)
        # gene_peak_fit=pickle.load(open(outdir+'/gene_peak_regress_250k.pkl','rb')) 
        # get regulatory matrix
        pseudo_reg_df=pseudo_regulatory_matrix(gene_peak_fit['filtered_peaks'], str(cluster))
        pseudo_reg_df.to_csv(outdir+'/c'+str(cluster)+'_reg_df.txt', header=True, index=True, sep='\t')
        
        #if not os.path.exists(outdir+'/reg_df'):
        #    os.makedirs(outdir+'/reg_df')
        #return pseudo_reg_df
    
    Pros=[]
    i=0
    for cluster in meta_cells:
        i+=1
        print('Starting processing %d' %i)
        p = Process(target=pseudo_reg, args=(rn, an, cluster, outdir, gene_info))
        Pros.append(p)
        p.start()
    for t in Pros:
        t.join()
    
    print('Meta cells: %d'%(i))
