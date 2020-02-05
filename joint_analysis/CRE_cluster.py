### CRE of clusters
import pandas as pd
import scanpy as sc
import sys
import os
import pickle
from multiprocessing import Process

script_path='/home/xionglei/yanqiu/regulatory'
sys.path.insert(1,script_path)
sys.path.insert(1,'/home/xionglei/yanqiu')
sys.path.insert(1, '/home/xionglei/yanqiu/joint_analysis/')
from process_data import overlap_adata
from identify_CRE import identify_cre
from utils import load_data, write_mtx


outdir='/home/xionglei/yanqiu/joint_analysis/WYF191116'
gene_info = pickle.load(open('/home/xionglei/yanqiu/lib/gene_info_hg.pkl', 'rb'))

'''
cluster_df=pd.read_csv(outdir+'/WYF191116_cluster.txt', sep='\t', index_col=0, header=0)
clusters=cluster_df['clusters'].unique()
cluster_dict=cluster_df.to_dict()
'''

genedata=sc.read_h5ad(outdir+'/gene_expression.h5ad')
prodata=sc.read_h5ad(outdir+'/promoter_accessibility.h5ad')
genedata.obs['clusters']=genedata.obs['seurat_clusters']
atac_adata=sc.read_h5ad(outdir+'/atac.h5ad')

'''
atac_dir='/home/xionglei/data/joint_ATAC_RNA/WYF191116/ATAC/peak'
atac=load_data(atac_dir)   
atac_adata=sc.AnnData(atac.T)
sc.pp.filter_cells(atac_adata, min_genes=1)
sc.pp.filter_genes(atac_adata, min_cells=1)
atac_adata.write(outdir+'/atac.h5ad')
'''

rn, an=overlap_adata(genedata, atac_adata)
rn, pro=overlap_adata(genedata, prodata)

#print(rn.shape, an.shape, pro.shape)

def cre_cluster(rn, an, pro, cluster, outdir):
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    cells=rn[rn.obs['clusters']==cluster].obs_names
    #print(len(cells))
    rn_f=rn[cells, ]
    an_f=an[cells, ]
    pro_f=pro[cells, ]
    #print(rn_f.shape, an_f.shape, pro_f.shape)
    CRE_df=identify_cre(rn_f, an_f, pro_f, gene_info, outdir)


clusters=rn.obs['clusters'].unique()
# cannot run together??
Pros=[]
i=0
for c in ['2','4']:
        i+=1
        print('Starting processing %d' %i)
        p = Process(target=cre_cluster, args=(rn, an, pro, c, outdir+'/CRE/'+str(c)))
        Pros.append(p)
        p.start()
for t in Pros:
        t.join()



#cluster=sys.argv[1]
#for c in clusters:
#    cre_cluster(rn, an, pro, c, outdir+'/CRE/'+str(c))