import sys
import os
import numpy as np
import pandas as pd
import scanpy as sc
sys.path.insert(1, '/home/xionglei/yanqiu/')
sys.path.insert(1, '/home/xionglei/yanqiu/joint_analysis/')

import matplotlib
import matplotlib.pyplot as plt
from utils import read_mtx,write_mtx, overlap_adata, plot_heatmap
from gene_promoter import gene_promoter_heatmap
%matplotlib inline

import pickle
gene_info = pickle.load(open('/home/xionglei/yanqiu/lib/gene_info_hg.pkl', 'rb'))
cell_label=pd.read_csv('/home/xionglei/data/joint_ATAC_RNA/WYF191116/cell_label/cell_label2.txt', sep='\t', header=None, index_col=0)
cell_label=cell_label[1].to_dict()

#==========================
# use scanpy to process rna data
#==========================

results_file='/home/xionglei/yanqiu/WYF191116/rna.h5ad'
# read data and basic info
rna=read_mtx('/home/xionglei/data/joint_ATAC_RNA/WYF191116/RNA/')
rna_adata_raw=sc.AnnData(rna.T)
sc.pp.filter_cells(rna_adata_raw, min_genes=1)
sc.pp.filter_genes(rna_adata_raw, min_cells=1)
rna_adata_raw.obs['n_counts'] = rna_adata_raw.X.sum(axis=1)
sc.pl.violin(rna_adata_raw, ['n_genes', 'n_counts'],
                     jitter=0.4, multi_panel=True)
sc.pl.scatter(rna_adata_raw, x='n_counts', y='n_genes')
rna_median_counts=np.median(rna_adata_raw.obs['n_counts'])

# anno gene and cell label
rna_adata_raw.obs['cell_label']=[cell_label[bc] for bc in rna_adata_raw.obs.index.values]
rna_adata_raw.var.index=[gene_info['gene_id_to_name'][gene] for gene in rna_adata_raw.var.index.values]
rna_adata_raw.var_names_make_unique()

# filter and normalize
sc.pp.filter_cells(rna_adata_raw, min_genes=50)
sc.pp.filter_genes(rna_adata_raw, min_cells=5)
sc.pl.highest_expr_genes(rna_adata_raw, n_top=20)
rna_adata=rna_adata_raw
sc.pp.normalize_per_cell(rna_adata, counts_per_cell_after=1e4)
sc.pp.log1p(rna_adata)
rna_adata.raw = rna_adata

# get highly variable genes
sc.pp.highly_variable_genes(rna_adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(rna_adata)
rna_adata = rna_adata[:, rna_adata.var['highly_variable']]
sc.pp.scale(rna_adata, max_value=10)

# cluster by pca, knn and lovain
sc.tl.pca(rna_adata, svd_solver='arpack')
sc.pl.pca(rna_adata, color=['cell_label','PTPRC'])
sc.pl.pca_variance_ratio(rna_adata, log=True)
sc.pp.neighbors(rna_adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(rna_adata)
sc.tl.louvain(rna_adata)
sc.pl.umap(rna_adata, color=['cell_label','PTPRC', 'louvain'])

# get marker genes
sc.tl.rank_genes_groups(rna_adata, 'louvain', method='wilcoxon')
sc.pl.rank_genes_groups(rna_adata, n_genes=25, sharey=False)
marker_genes = pd.DataFrame(rna_adata.uns['rank_genes_groups']['names']).loc[0:2,:].values.flatten()
ax = sc.pl.dotplot(rna_adata, marker_genes, groupby='louvain')
ax = sc.pl.stacked_violin(rna_adata, marker_genes, groupby='louvain', rotation=90)

# write h5ad file
rna_adata.write(results_file)

#==========================
# use seurat for clustering
#==========================
'''
seurat.ipynb
'''

#==========================
#import clusters from seurat
#==========================

import scanpy as sc
seurat = sc.read_loom('seurate_cluster_RNA.loom', sparse=False, obs_names='CellID', var_names='Gene', dtype='float32')
seurat.obs['clusters']=seurat.obs['seurat_clusters']
seurat_raw=seurat

seurat.var.index=[gene_info['gene_id_to_name'][gene] for gene in seurat.var.index.values]
seurat.var_names_make_unique()

sc.tl.rank_genes_groups(seurat, 'clusters', method='wilcoxon')
sc.pl.rank_genes_groups(seurat, n_genes=25, sharey=False)

# plot heatmap of specific genes
specific_genes=pd.DataFrame(seurat.uns['rank_genes_groups']['names']).loc[0:20,:].T.values.flatten()
X=pd.DataFrame(seurat[:,specific_genes].layers['norm_data']).T
X.index=specific_genes
X.columns=seurat.obs['clusters']
from scale_plot_feature import plot_heatmap
plot_heatmap(X, y=seurat.obs['clusters'], #row_labels=specific_genes, 
                     ncol=3, cmap='Reds',vmax=1, row_cluster=False, legend_font=6, cax_title='Peak Value',
                    figsize=(8, 10), bbox_to_anchor=(0.4, 1.2), position=(0.8, 0.76, 0.1, 0.015),
                     save='test_specific_feature.png')

marker_genes={
        'CD45':['PTPRC'],
        'CD8_T':['CD8A', 'CD8B'],
        'exhausted_CD8':['LAG3', 'CD244', 'EOMES', 'PTGER4'],
        'T_cell':['CD6', 'CD3D', 'CD3E', 'SH2D1A', 'TRAT1', 'CD3G'],
        'B_cell':['BLK', 'CD19', 'FCRL2', 'MS4A1', 'KIAA0125', 'TNFRSF17', 'TCL1A', 'SPIB', 'PNOC'],
        'NK':['XCL1', 'XCL2', 'NCR1'],
        'NK_CD56':['KIR2DL3', 'KIR3DL1', 'KIR3DL2', 'IL21R'],
        'DC':['CCL13', 'CD209', 'HSD11B1'],
        'macro':['CD68','CD84', 'CD163', 'MS4A4A'],     
        'mast':['TPSB2', 'TPSAB1', 'CPA3', 'MS4A2', 'HDC'],
        'neutrophil':['FPR1', 'SIGLEC5', 'CSF3R', 'FCAR', 'FCGR3B', 'CEACAM3', 'S100A12'],
        'Th1':['TBX21'],   
        'T_reg':['FOXP3'],
        'cytotoxic':['PRF1', 'GZMA', 'GZMB', 'NKG7', 'GZMH', 'KLRK1', 'KLRB1', 'KLRD1', 'CTSW', 'GNLY'],
    
        'macro2':['CD14', 'AIF1', 'FCER1G', 'FCGR3A', 'TYROBP', 'CSF1R'],
        'T_cell2': ['CD2', 'CD3D', 'CD3E', 'CD3G'],
    
        'oligodendrocytes': ['MBP', 'TF', 'PLP1', 'MAG', 'MOG', 'CLDN11'],
    
        'Low grade astrocytomas':['CHI3L2','NMB','VCAM1'],
        'Anaplastic astrocytomas (grade III)':['CHI3L2','SAA1','SMOC1','ADORA3','NMB','CRF','CSRP2','VCAM1'],
        'Glioblastomas (grade IV)':['AQP1','TYMS','TOP2A','ABCC3','SAA1','CHI3L2','NMB','MGP']
}
kept_markers=[]
for group in marker_genes:
    for gene in marker_genes[group]:
        if gene in seurat.var.index.values:
            kept_markers.append(gene)
            
kept_markers=list(set(kept_markers))
kept_markers.remove('PTPRC')
kept_markers.remove('MBP')
kept_markers.remove('SMOC1')

ax = sc.pl.dotplot(seurat, kept_markers, groupby='seurat_clusters')


#==========================
# match promoter accessibility
#==========================
from utils import overlap_adata
from gene_promoter import proc_prodata, gene_promoter_heatmap
promoter=read_mtx('/home/xionglei/data/joint_ATAC_RNA/WYF191116/promoter/')
pro_adata=sc.AnnData(promoter.T)

pro_adata_join=proc_prodata(pro_adata)

df=pd.DataFrame(pro_adata_join.X.astype('int')).T
df.index=pro_adata_join.var_names
df.columns=pro_adata_join.obs_names
path='/home/xionglei/data/joint_ATAC_RNA/WYF191116/promoter_new'
write_mtx(df, path)

# overlap obs and var
genedata, prodata=overlap_adata(seurat, pro_adata_join)

cell_label=pd.read_csv('/home/xionglei/data/joint_ATAC_RNA/WYF191116/cell_label/cell_label2.txt', sep='\t', header=None, index_col=0)
cell_label=cell_label[1].to_dict()
genedata.obs['cell_label']=[cell_label[bc] for bc in genedata.obs.index.values]
#genedata.obsm['X_umap']=genedata.obsm['umap_cell_embeddings']
prodata.obs['cell_label']=[cell_label[bc] for bc in prodata.obs.index.values]


genedata.obs['clusters']=genedata.obs['seurat_clusters']
genedata.obs[['clusters']].to_csv('WYF191116_cluster.txt',sep='\t')

genedata.write('gene_expression.h5ad')
prodata.write('promoter_accessibility.h5ad')


gene_promoter_heatmap(genedata, prodata, 50, group=False, filter_cor=True)
gene_promoter_heatmap(genedata, prodata, 50, group=False, filter_cor=False)
gene_promoter_heatmap(genedata, prodata, 50, group=True, filter_cor=False)
gene_promoter_heatmap(genedata, prodata, 50, group=True, filter_cor=True)

#==========================
# use chromVAR for TF motif enrichment
#==========================
'''
Rscript chromVAR.R -d ~/data/joint_ATAC_RNA/WYF191116/ATAC/peak \
    -p ~/data/joint_ATAC_RNA/WYF191116/ATAC/peaks_f.bed \
    -g hg38 \
    -o .
'''
genedata=sc.read_h5ad('gene_expression.h5ad')
prodata=sc.read_h5ad('promoter_accessibility.h5ad')

tfdata=pd.read_csv('dev.txt', sep='\t', header=0, index_col=0)
tfdata.index=[el.split('_')[-1] for el in tfdata.index]
tfdata=sc.AnnData(tfdata.T)
tfdata_f=tfdata[genedata.obs_names,]
tfdata_f.obsm['X_umap']=genedata.obsm['umap_cell_embeddings']
tfdata_f.obs['clusters']=genedata.obs['clusters']

cell_label=pd.read_csv('/home/xionglei/data/joint_ATAC_RNA/WYF191116/cell_label/cell_label2.txt', sep='\t', header=None, index_col=0)
cell_label=cell_label[1].to_dict()
tfdata_f.obs['cell_label']=[cell_label[bc] for bc in tfdata_f.obs.index.values]
sc.pl.umap(tfdata_f, color=['cell_label','clusters','CUX1',	'HOXD13','MIXL1','DMRT3','EMX1','BATF3','ETV2','XBP1'])
sc.pp.scale(tfdata_f, max_value=10)

X=pd.DataFrame(tfdata_f.X).T
X.index=tfdata_f.var_names
X.columns=tfdata_f.obs['clusters']
X_group=X.T
X_group['clusters']=X_group.index
X_group=X_group.groupby('clusters').mean().T
plot_heatmap(X_group, X_group.columns, ncol=3, cmap='inferno',vmax=1, row_cluster=True, legend_font=6, cax_title='TF deviation',
                figsize=(8, 10), bbox_to_anchor=(0.4, 1.2), position=(0.8, 0.76, 0.1, 0.015))


#==========================
# identify CRE
#==========================
import sys
script_path='/home/xionglei/yanqiu/regulatory'
sys.path.insert(1,script_path)
sys.path.insert(1,'/home/xionglei/yanqiu')
sys.path.insert(1, '/home/xionglei/yanqiu/joint_analysis/')
import os
import pandas as pd
import scanpy as sc
import pickle
import matplotlib.pyplot as plt
%matplotlib inline
from utils import load_data, write_mtx

#from plot import plot_info,plot_gene_peak
from process_data import normalize_rna,normalize_atac,group_cells,overlap_adata
from identify_CRE import identify_cre

rna_dir='/home/xionglei/data/joint_ATAC_RNA/WYF191116/RNA'
atac_dir='/home/xionglei/data/joint_ATAC_RNA/WYF191116/ATAC/peak'
outdir='/home/xionglei/yanqiu/joint_analysis/WYF191116'

gene_info = pickle.load(open('/home/xionglei/yanqiu/lib/gene_info_hg.pkl', 'rb'))

# input atac and rna data
atac=load_data(atac_dir)   
atac_adata=sc.AnnData(atac.T)
atac_adata=normalize_atac(atac_adata)
# group cells by seurat
os.system('Rscript %s/group_cells.R %s %s'%(script_path, rna_dir, outdir))
meta_rn=sc.read_loom(outdir+'/group_labeled.loom', sparse=False, obs_names='CellID', var_names='Gene', dtype='float32')
rn_o, an_o=overlap_adata(meta_rn, atac_adata)
print('Overlap cells: %d'%(rn_o.shape[0]))
rn, an=group_cells(rn_o, an_o) #rn is from seurate generated loom
print('Meta cells: %d'%(rn.shape[0]))
print('RNA genes: %d'%(rn.shape[1]))
print('ATAC peaks: %d'%(an.shape[1]))

# input promoter accessibility data, var_name is gene id
pro_adata=sc.read_h5ad('promoter_accessibility.h5ad')
pro_adata_o=pro_adata[rn_o.obs_names,]
pro_adata_o.obs['ClusterID']=rn_o.obs['ClusterID']
pro, an_pro=group_cells(pro_adata_o, an_o) # an=an_pro

CRE_df=identify_cre(rn, an, pro, gene_info, outdir)


### CRE of clusters
outdir='/home/xionglei/yanqiu/joint_analysis/WYF191116'
atac_dir='/home/xionglei/data/joint_ATAC_RNA/WYF191116/ATAC/peak'
gene_info = pickle.load(open('/home/xionglei/yanqiu/lib/gene_info_hg.pkl', 'rb'))

cluster_df=pd.read_csv(outdir+'/WYF191116_cluster.txt', sep='\t', index_col=0, header=0)
clusters=cluster_df['clusters'].unique()
cluster_dict=cluster_df.to_dict()

genedata=sc.read_h5ad(outdir+'/gene_expression.h5ad')
prodata=sc.read_h5ad(outdir+'/promoter_accessibility.h5ad')
genedata.obs['clusters']=genedata.obs['seurat_clusters']

atac=load_data(atac_dir)   
atac_adata=sc.AnnData(atac.T)
sc.pp.filter_cells(atac_adata, min_genes=1)
sc.pp.filter_genes(atac_adata, min_cells=1)
atac_adata.write(outdir+'/atac.h5ad')
rn, an=overlap_adata(genedata, atac_adata)
rn, pro=overlap_adata(genedata, prodata)

def cre_cluster(rn, an, pro, cluster, outdir):
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    cells=rn[rn.obs['clusters']==cluster].obs_names
    rn=rn[cells, ]
    an=an[cells, ]
    pro=pro[cells, ]
    CRE_df=identify_cre(rn, an, pro, gene_info, outdir)

clusters=rn.obs['clusters'].unique()
# cannot run together??
Pros=[]
i=0
for c in clusters[2:7]:
        i+=1
        print('Starting processing %d' %i)
        p = Process(target=cre_cluster, args=(rn, an, pro, c, outdir+'/CRE/'+str(c)))
        Pros.append(p)
        p.start()
for t in Pros:
        t.join()
        
        
#==========================
# snp enrichment
#==========================

### snp.py
from glob import glob
import os

oddsratio=[]
samples=['genome_wide','gene_around','around_peak']
for sample in samples:
    oddsratio.append(pd.read_csv(outdir+'/snps_enrich/%s_snp_enrichment.txt'%sample,sep='\t',header=0, index_col=None)['oddsratio'])
oddsratio=np.array(oddsratio)
plt.figure(figsize=(2,4))
plt.bar(range(len(oddsratio)), 
              np.mean(oddsratio, axis=1), width=0.6, 
              #yerr=np.std(oddsratio, axis=1), 
              align='center', alpha=0.7, ecolor='black', capsize=2)
plt.xticks(range(len(oddsratio)), samples, rotation=90)
plt.ylabel('Odds ratio')
plt.savefig(outdir+'/snp_enrichment.png',bbox_inches='tight')

### snp enrichment for each cluster
## snp.py

outdir='/home/xionglei/yanqiu/joint_analysis/WYF191116'
cluster_df=pd.read_csv(outdir+'/WYF191116_cluster.txt', sep='\t', index_col=0, header=0)
clusters=cluster_df['clusters'].unique()
clusters.sort()
oddsratio=[]
for c in clusters:
    data=pd.read_csv(outdir+'/snps_enrich/%s/around_peak_snp_enrichment.txt'%c,sep='\t',header=0, index_col=None)['oddsratio']
    data[np.isnan(data)]=0
    data[np.isinf(data)]=max(data[~np.isinf(data)])
    oddsratio.append(data)
oddsratio=np.array(oddsratio)
#oddsratio[np.isnan(oddsratio)]=0
#oddsratio[np.isinf(oddsratio)] = 100
plt.figure(figsize=(2,4))
plt.bar(range(len(oddsratio)), 
              np.mean(oddsratio, axis=1), width=0.6, 
              #yerr=np.std(oddsratio, axis=1), 
              align='center', alpha=0.7, ecolor='black', capsize=2)
plt.xticks(range(len(oddsratio)), clusters) #, rotation=90
plt.ylabel('Odds ratio')
plt.ylim(0,20)
plt.xlabel('clusters')
plt.savefig(outdir+'/cluster_snp_enrichment.png',bbox_inches='tight')


#==========================
# eQTL enrichment
#==========================

### eqtl.py

from glob import glob
import os

files=glob('/home/xionglei/yanqiu/lib/eQTL/GTEx_Analysis_v8_eQTL/*signif_variant_gene_pairs.txt.gz')
eqtl_data={}
files=sorted(files)
for filename in files:
    basename = os.path.basename(filename)
    sample=basename.split('.')[0]
    eqtl_data[sample]=filename

oddsratio=[]
for sample in eqtl_data.keys():
    oddsratio.append(pd.read_csv(outdir+'/eqtl_enrich/%s_eqtl_enrichment.txt'%sample,sep='\t',header=0, index_col=None)['oddsratio'])
oddsratio=np.array(oddsratio)
plt.figure(figsize=(12,4))
plt.bar(range(len(oddsratio)), 
              np.mean(oddsratio, axis=1), width=0.6, 
              yerr=np.std(oddsratio, axis=1), 
              align='center', alpha=0.7, ecolor='black', capsize=2)
plt.xticks(range(len(oddsratio)), list(eqtl_data.keys()), rotation=90)
plt.savefig(outdir+'/eqtl_enrichment.png',bbox_inches='tight')



#==========================
# trajectory analysis
#==========================

