# stats of public data
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import pandas as pd
import numpy as np
import seaborn as sns

sys.path.insert(1, '/home/xionglei/yanqiu/')
from utils import plot_umi, load_data, read_mtx
import config

data_config=config.data_config

umi_atac=pd.DataFrame()
for dataset in data_config:
    df=load_data(data_config[dataset]['atac_path'])
    umi=df.sum(axis=0)
    umi=pd.DataFrame(umi, columns=['n_umi'])
    umi['sample']=dataset
    umi_atac=umi_atac.append(umi)

umi_atac.to_csv('atac_umi.txt',sep='\t')

plt.figure(figsize=(8,5))
sns.boxplot(x="sample",y="n_umi",data=gene_rna,palette="Set1",linewidth=1,fliersize=1)
plt.ylim(0,10000)
plt.savefig('atac_umi.png')
