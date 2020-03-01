from plot import plot_embedding,plot_clustercounts,plot_heatmap

## Fig2.a
# clusters is sorted cluster names
plot_embedding(np.array(hema[['UMAP1','UMAP2']]), hema['ClusterName'], classes=clusters,
               cmap='nipy_spectral', figsize=(8,8),markersize=1) # legend sorted by classes


cluster_counts=hema['ClusterName'].value_counts()
plot_clustercounts(cluster_counts,clusters=clusters, cmap='nipy_spectral',width=0.7,figsize=(2,8),yticks=None)


## Fig2.b
adata=sc.read_h5ad('data/Hematopoiesis-All/adata.h5ad')
specific_genes=pd.DataFrame(adata.uns['rank_genes_groups']['names']).loc[0:20,:].T.values.flatten()
X2=pd.DataFrame(adata[:,specific_genes].X.todense()).T
X2.index=specific_genes
X2.columns=adata.obs['clusters'].values
X2_group=X2.T
X2_group['clusters']=X2_group.index
X2_group=X2_group.groupby('clusters').mean().T
X2=X2_group
y_cluster=np.array([el.split('-')[0] for el in X2.columns])
X2.columns=y_cluster
plot_heatmap(X2, y=y_cluster, ncol=3, vmax=1, row_cluster=True, col_cluster=True, 
             z_score=1, cmap='RdYlBu_r',
             colormap=sns.color_palette('nipy_spectral', n_colors=len(y_cluster)),
            legend_font=6, cax_title='zscore',show_legend=False,
            figsize=(8, 10), bbox_to_anchor=(0.4, 1.2), position=(1, 0.6, .04, .1), 
             save='specific_peaks.png')

# chromVar.R: make sure peak and barcode format is correct

# cicero.R