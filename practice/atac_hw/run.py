from plot import plot_embedding

# clusters is sorted cluster names
plot_embedding(np.array(hema[['UMAP1','UMAP2']]), hema['ClusterName'], classes=clusters,
               cmap='nipy_spectral', figsize=(8,8),markersize=1) # legend sorted by classes

cluster_counts=hema['ClusterName'].value_counts()
plot_clustercounts(cluster_counts,clusters=clusters, cmap='nipy_spectral',width=0.7,figsize=(2,8),yticks=None):