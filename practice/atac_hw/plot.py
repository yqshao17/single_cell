import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def plot_embedding(X, labels, classes=None, method='tSNE', cmap='tab20', figsize=(4, 4), markersize=4, marker=None,
                   return_emb=False, save=False, save_emb=False, show_legend=True, show_axis_label=True, **legend_params):
    if marker is not None:
        X = np.concatenate([X, marker], axis=0)
    N = len(labels)
    if X.shape[1] != 2:
        if method == 'tSNE':
            from sklearn.manifold import TSNE
            X = TSNE(n_components=2, random_state=124).fit_transform(X)
        if method == 'UMAP':
            from umap import UMAP
            X = UMAP(n_neighbors=30, min_dist=0.3, metric='correlation').fit_transform(X)
        if method == 'PCA':
            from sklearn.decomposition import PCA
            X = PCA(n_components=2, random_state=124).fit_transform(X)
        
    plt.figure(figsize=figsize)
    if classes is None:
        classes = np.unique(labels)

    if cmap is not None:
        cmap = cmap
    elif len(classes) <= 10:
        cmap = 'tab10'
    elif len(classes) <= 20:
        cmap = 'tab20'
    else:
        cmap = 'husl'
    colors = sns.color_palette(cmap, n_colors=len(classes))
    #np.random.shuffle(colors)
    
    class_scs=[]
    for i, c in enumerate(classes):
        class_sc=plt.scatter(X[:N][labels==c, 0], X[:N][labels==c, 1], s=markersize, color=colors[i], label=c)
        class_scs.append(class_sc)

    legend_params_ = {'loc': 'center left',
                     'bbox_to_anchor':(1.0, 0.45),
                     'fontsize': 10,
                     'ncol': 1,
                     'frameon': False,
                     'markerscale': 5,
                      'handles':class_scs
                    }
    legend_params_.update(**legend_params)
    if show_legend:
        plt.legend(**legend_params_)
    sns.despine(offset=10, trim=True)
    if show_axis_label:
        plt.xlabel(method+' dim 1', fontsize=12)
        plt.ylabel(method+' dim 2', fontsize=12)

    if save:
        plt.savefig(save, format='pdf', bbox_inches='tight')
    else:
        plt.show()
        
    if save_emb:
        np.savetxt(save_emb, X)
    if return_emb:
        return X
    
def plot_clustercounts(cluster_counts,clusters=None, cmap='nipy_spectral',width=0.7,figsize=(2,8),yticks=None,save=None):
    plt.figure(figsize=figsize)
    if clusters: # sort by clusters
        cluster_counts=cluster_counts.loc[clusters]
    inds=np.arange(len(cluster_counts)) - width/2
    plt.barh(inds, cluster_counts.loc[::-1], 
         width, color=sns.color_palette(cmap, n_colors=len(cluster_counts))[::-1])
    if yticks:
        plt.yticks(inds, clusters[::-1])
    else:
        plt.yticks([])
    plt.ylim(-1,len(cluster_counts)-1)
    if save:
        plt.savefig(save, format='pdf', bbox_inches='tight')
    else:
        plt.show()