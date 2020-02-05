import pandas as pd
import scipy.stats as stats
import scanpy as sc
from utils import plot_heatmap
import sys
sys.path.insert(1, '/home/xionglei/yanqiu/')
sys.path.insert(1, '/home/xionglei/yanqiu/joint_analysis/')

def overlap_adata(data1, data2):
    features=list(set(data1.var_names)&set(data2.var_names)) 
    barcodes=list(set(data1.obs_names)&set(data2.obs_names))       
    data1=data1[barcodes, :]
    data2=data2[barcodes, :]
    data1=data1[:, features]
    data2=data2[:, features]
    return data1, data2


def proc_genedata(genedata):
    sc.pp.filter_cells(genedata, min_genes=100)
    sc.pp.filter_genes(genedata, min_cells=10)
    sc.pp.normalize_per_cell(genedata, counts_per_cell_after=1e4)
    sc.pp.log1p(genedata)
    return genedata


def proc_prodata(prodata):
    # change promoter region to gene name
    def pro2gene(pro_adata):
        gene_promoter=pd.read_csv('/home/xionglei/yanqiu/lib/hg38_Gencode_promoter.bed', sep='\t', header=None)
        promoters=gene_promoter.apply(lambda x: x[0]+':'+str(x[1])+'-'+str(x[2]), axis=1)
        genes=gene_promoter.apply(lambda x: x[4].split(',')[0], axis=1)
        pairs={}
        for p,g in zip(promoters, genes):
            pairs[p]=g
        pro_adata.var.index=[pairs[promoter] for promoter in pro_adata.var.index.values]
        pro_adata.var_names_make_unique()
    return pro_adata

    # merge promoters for gene
    def join_promoters(pro_adata):
        genes=pro_adata.var_names
        genes_rm_dup={}
        matrix=[]
        for gene in genes:
            gene_uniq=gene.split('-')[0]
            genes_rm_dup.setdefault(gene_uniq, [])
            genes_rm_dup[gene_uniq].append(pro_adata.obs_vector(gene))  
        for gene in genes_rm_dup:
            obs_vector=np.array(genes_rm_dup[gene]).sum(axis=0)
            matrix.append(obs_vector)
        new_data=pd.DataFrame(np.array(matrix), index=genes_rm_dup.keys(), columns=pro_adata.obs.index.values)
        pro_adata_new=sc.AnnData(new_data.T)
    return pro_adata_new

    prodata=pro2gene(prodata)
    pro_adata_join=join_promoters(pro_adata)
    # normalize prodata
    sc.pp.filter_cells(pro_adata_join, min_genes=1)
    sc.pp.filter_genes(pro_adata_join, min_cells=1)
    #pro_adata_join.obs['n_counts'] = pro_adata_join.X.sum(axis=1)
    #sc.pl.violin(pro_adata_join, ['n_genes', 'n_counts'], jitter=0.4, multi_panel=True)
    sc.pp.normalize_per_cell(pro_adata_join, counts_per_cell_after=1e4)
    sc.pp.log1p(pro_adata_join)
    return pro_adata_join


def gene_promoter_heatmap(genedata, prodata, n_genes, outdir='.', group=True, filter_cor=False):  
    
    # plot cluster-specific genes and promoter accessibility
    # genedata is clusterd by Seurat or other methods, 
    # and saved and read as adata, with 'clusters' as obs key
    
    # whether normalize genedata or prodata?
    save_prefix = outdir+'/'
    genedata2, prodata2=overlap_adata(genedata, prodata)

    sc.tl.rank_genes_groups(genedata2, 'clusters', method='wilcoxon')
    specific_genes=pd.DataFrame(genedata2.uns['rank_genes_groups']['names']).loc[0:n_genes,:].T.values.flatten()
    
    if filter_cor:
        save_prefix+='filtered_'
        # filter specific genes that have high correlation 
        # between gene expression and gene promoter accessibility
        specific_genes_filtered=[]
        for gene in specific_genes:
            cor,pval=stats.pearsonr(genedata2.obs_vector(gene),prodata2.obs_vector(gene))
            if cor>0 and pval<0.05:
                specific_genes_filtered.append(gene)
        specific_genes=specific_genes_filtered
    
    X1=pd.DataFrame(genedata2[:,specific_genes].layers['norm_data']).T
    #X1=pd.DataFrame(genedata2[:,specific_genes].X).T
    X1.index=specific_genes
    X1.columns=genedata2.obs['clusters']
    X2=pd.DataFrame(prodata2[:,specific_genes].X).T
    X2.index=specific_genes
    X2.columns=genedata2.obs['clusters']
    
    
    if group:
        # group each cluster by mean expression as one cell
        X_group=X1.T
        X_group['clusters']=X_group.index
        X_group=X_group.groupby('clusters').mean().T
        X1=X_group
        X2_group=X2.T
        X2_group['clusters']=X2_group.index
        X2_group=X2_group.groupby('clusters').mean().T
        X2=X2_group
        save_prefix+='grouped_'
    #print(X1.iloc[:5,:5])
    
    plot_heatmap(X1, y=X1.columns, #row_labels=specific_genes, 
                 ncol=3, cmap='Reds',vmax=1, row_cluster=False, legend_font=6, cax_title='Gene Value',
                figsize=(8, 10), bbox_to_anchor=(0.4, 1.2), position=(0.8, 0.76, 0.1, 0.015),
                save=save_prefix+'genes_expression.png')

    
    plot_heatmap(X2, y=X2.columns, #row_labels=specific_genes, 
                 ncol=3, cmap='Reds',vmax=1, row_cluster=False, legend_font=6, cax_title='Promoter Value',
                figsize=(8, 10), bbox_to_anchor=(0.4, 1.2), position=(0.8, 0.76, 0.1, 0.015),
                save=save_prefix+'promoter_accessibility.png')

