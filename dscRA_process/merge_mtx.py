# merge RNA mtx
from utils import get_prefix, read_mtx, write_mtx
import sys
import pandas as pd

indir=sys.argv[1]
outdir=sys.argv[2]

prefix_list=get_prefix(indir+'/sample.list')

df_merge=pd.DataFrame()

for prefix in prefix_list:
    df=read_mtx(indir+'/DGE_filtered/'+prefix+'_RNA')
    df_merge=pd.merge(df_merge, df, left_index=True, right_index=True, how='outer')

df_merge=df_merge.fillna(0)
df_merge=df_merge.astype('int')

write_mtx(df_merge, outdir)
