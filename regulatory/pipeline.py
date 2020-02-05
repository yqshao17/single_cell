import pickle
import sys
sys.path.insert(1, '/home/xionglei/yanqiu')
from utils import load_data, write_mtx
from gene_peak_fit import get_gene_peak_fit, get_regulatory_matrix
#from plot import plot_info,plot_gene_peak
import os
import pandas as pd
import scanpy as sc
from plot import plot_info, plot_gene_peak
from process_data import normalize_rna,normalize_atac,group_cells,overlap_adata


import argparse
parser = argparse.ArgumentParser()

parser.add_argument('--outdir','-o', help='output dir')
parser.add_argument('--rnadir', '-r', help='rna directory or h5ad file')
parser.add_argument('--atacdir', '-a', help='atac directory or h5ad file')
parser.add_argument('--method', '-m', help='regress, pearson or spearman', default='regress')
parser.add_argument('--species', '-s', help='mm, hg or mm_hg')
parser.add_argument('--bin_size', '-bs', help='bin size: *k', default=250)
parser.add_argument('--n_top_genes', '-gn', help='number of filtered variable genes', default=None)
parser.add_argument('--group_by','-g', help='group cells by RNA, ATAC or None', default=None)


args = parser.parse_args()

outdir=args.outdir
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

    if not group_by:
            rn, an=overlap_adata(rn, an)
            weight=[1]*rn.shape[0]
            rn.obs['Weights']=weight
            an.obs['Weights']=weight
            print('Overlap cells: %d'%(rn.shape[0]))
            print('RNA genes: %d'%(rn.shape[1]))
            print('ATAC peaks: %d'%(an.shape[1]))
        
    if group_by=='RNA':
        if rna_dir.endswith('loom'):
            pass
        else:
            os.system('Rscript group_cells.R %s %s'%(rna_dir, outdir))
        meta_rn=sc.read_loom(outdir+'/group_labeled.loom', sparse=False, obs_names='CellID', var_names='Gene', dtype='float32')
        rn, an=overlap_adata(meta_rn, an)
        print('Overlap cells: %d'%(rn.shape[0]))
        rn, an=group_cells(rn, an) #rn is from seurate generated loom
        print('Meta cells: %d'%(rn.shape[0]))
        print('RNA genes: %d'%(rn.shape[1]))
        print('ATAC peaks: %d'%(an.shape[1]))
        
    
    rn.write(outdir+'/rn.h5ad')
    an.write(outdir+'/an.h5ad')
    
    
    #rn=sc.read_h5ad(outdir+'/rn.h5ad')
    #an=sc.read_h5ad(outdir+'/an.h5ad')

    gene_info = pickle.load(open('/home/xionglei/yanqiu/regulatory/lib/gene_info_%s.pkl'%species, 'rb'))
    # calculate gene-peak regression
    print('Calculating gene-peak fit')        
    gene_peak_fit = get_gene_peak_fit(gene_info, rn, an, outdir,bin_size, method)
    #gene_peak_fit=pickle.load(open(outdir+'/gene_peak_regress_250k.pkl','rb'))   
         
    # get regulatory matrix
    print('Getting gene-peak regulatory matrix')
    reg_df=get_regulatory_matrix(gene_peak_fit['filtered_peaks'],rn,an)
    if not os.path.exists(outdir+'/reg_df'):
        os.makedirs(outdir+'/reg_df')
    write_mtx(reg_df, outdir+'/reg_df')
    
    # plot information
    if method=='regress':
        plot_info(gene_peak_fit, save=outdir+'/gene_peaks_info.pdf')
        plot_gene_peak(gene_peak_fit, rn, an, save=outdir+'/gene_peaks_regression.pdf')
        

