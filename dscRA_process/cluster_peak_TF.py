import pandas as pd
import sys


def split_bed(bedfile, out_prefix, cluster_file):
    
    cluster=pd.read_csv(cluster_file,sep='\t', header=0, index_col=0)['clusters']
    classes=cluster.unique()
    cluster_dict=cluster.to_dict()
    
    inbed=pd.read_csv(bedfile, sep='\t', header=None, index_col=None)
    output={}
    for c in classes:
        output[c]=open(out_prefix+'_shift_'+c+'.bed','w')
    
    
    
    for c in classes:
        output.close()