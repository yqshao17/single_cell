# test different methods for getting regulatory matrix
import sys
indir=sys.argv[1]
outdir=sys.argv[2]

import os
if not os.path.exists(outdir):
    os.makedirs(outdir)

import pickle
import scanpy as sc
gene_peak_fit=pickle.load(open(indir+'/gene_peak_regress_250k.pkl','rb'))
rn=sc.read_h5ad(indir+'/rn.h5ad')
an=sc.read_h5ad(indir+'/an.h5ad')

from gene_peak_fit import get_regulatory_matrix
reg_df=get_regulatory_matrix(gene_peak_fit['filtered_peaks'],rn,an)

#sys.path.insert(1, '/home/xionglei/yanqiu')
#from utils import write_mtx
#write_mtx(reg_df, outdir)
import pandas as pd
reg_df.to_csv(outdir+'/count.txt', sep='\t', header=True, index=True)