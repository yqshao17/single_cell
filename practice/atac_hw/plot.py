import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def plot_embedding(X, labels, classes=None, method='tSNE', cmap=None, custom_cmap=None, figsize=(4, 4), markersize=4,       marker=None,return_emb=False, save=False, save_emb=False, show_legend=True, show_axis_label=True, **legend_params):
    # input array
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
    
    if custom_cmap is not None:
        colors = custom_cmap # list of colors
    else:
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
    # dataframe or series
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
        
def sort_by_classes(X, y, classes):
    if classes is None:
        classes = np.unique(y)
    index = []
    for c in classes:
        ind = np.where(y==c)[0]
        index.append(ind)
    index = np.concatenate(index)
    X = X.iloc[:, index]
    y = y[index]
    return X, y, classes, index

def plot_heatmap(X, y, classes=None, y_pred=None, row_labels=None, colormap=None, col_cluster=False, 
                 row_cluster=False,
                 cax_title='', xlabel='', ylabel='', yticklabels=False,xticklabels=True, legend_font=10, 
                 show_legend=True, show_cax=True, tick_color='black', ncol=3,
                 bbox_to_anchor=(0.5, 1.3), position=(1.1, 0.78, .04, .1), return_grid=False,
                 # position=(0.8, 0.78, .1, .04),
                 save=None, **kw):
    """
    plot hidden code heatmap with labels

    Params:
        X: fxn array, n is sample number, f is feature
        y: a array of labels for n elements or a list of array
    """

    import matplotlib.patches as mpatches  # add legend

    X, y, classes, index = sort_by_classes(X, y, classes)


    if colormap is None:
            colormap = plt.cm.tab20
            colors = {c:colormap(i) for i, c in enumerate(classes)}
    else:
            colors = {c:colormap[i] for i, c in enumerate(classes)}
    col_colors = [ colors[c] for c in y ]
        
    legend_TN = [mpatches.Patch(color=color, label=c) for c, color in colors.items()]

    if row_labels is not None:
        row_colors = [ colors[c] for c in row_labels ]
        kw.update({'row_colors':row_colors})

    kw.update({'col_colors':col_colors})

    cbar_kws={"orientation": "vertical"}
    grid = sns.clustermap(X, yticklabels=True, 
            col_cluster=col_cluster,
            row_cluster=row_cluster,
            cbar_kws=cbar_kws, **kw)
    if show_cax:
        grid.cax.set_position(position)
        grid.cax.tick_params(length=1, labelsize=8, rotation=0)
        grid.cax.set_title(cax_title, fontsize=10, y=0.35)

    if show_legend:
        grid.ax_heatmap.legend(loc='upper center', 
                               bbox_to_anchor=bbox_to_anchor, 
                               handles=legend_TN, 
                               fontsize=legend_font, 
                               frameon=False, 
                               ncol=ncol)
        grid.ax_col_colors.tick_params(labelsize=6, length=0, labelcolor='orange')
        grid.cax.set_position(position)
        grid.cax.tick_params(length=1, labelsize=4, rotation=0)
        grid.cax.set_title(cax_title, fontsize=6, y=0.35)

    if (row_cluster==True) and (yticklabels):
        yticklabels = yticklabels[grid.dendrogram_row.reordered_ind]

    grid.ax_heatmap.set_xlabel(xlabel)
    grid.ax_heatmap.set_ylabel(ylabel, fontsize=8)
    if not xticklabels:
        grid.ax_heatmap.set_xticklabels('')
    #grid.ax_heatmap.set_yticklabels(yticklabels, color=tick_color)
    if not yticklabels:
        grid.ax_heatmap.set_yticklabels('')
    grid.ax_heatmap.yaxis.set_label_position('left')
    grid.ax_heatmap.tick_params(axis='x', length=0)
    grid.ax_heatmap.tick_params(axis='y', labelsize=6, length=0, rotation=0, labelleft=True, labelright=False)
    grid.ax_row_dendrogram.set_visible(False)
    grid.ax_col_dendrogram.set_visible(False)
    grid.cax.set_visible(show_cax)
    grid.row_color_labels = classes

    if save:
        plt.savefig(save, bbox_inches='tight')
    else:
        plt.show()
    if return_grid:
        return grid

def SmoothMat(mat,k):
    return np.convolve(mat,np.ones(k),'same')/k

def heatmap_smooth(X, k):
    # input dataframe
    X_smooth= X.apply(lambda x: SmoothMat(x,k),axis=1)
    X_smooth=pd.DataFrame(np.array(list(X_smooth)))
    X_norm=X_smooth.divide(X_smooth.max(axis=1), axis=0)
    plt.figure(figsize=(6,8))
    sns.heatmap(X_norm,cmap='RdYlBu_r',xticklabels=False, yticklabels=False,vmax=1,vmin=0)
    
def scatter_traj(adata, var_name):
    plt.figure(figsize=(8,2))
    x=range(adata.shape[0])
    y=adata.obs_vector(var_name)
    plt.scatter(x,y,s=1, c=x, cmap='viridis')
    plt.plot(x,SmoothMat(y, 500),c='k',linewidth=1)
    
def scatter_traj(adata, var_names):
    # input adta
    x=range(tf_traj_s.shape[0])
    if type(var_names)==type([]):
        n_rows = len(var_names)
        fig, ax = plt.subplots(nrows=n_rows, ncols=1, figsize=(8,2.5*n_rows))
        for i in range(n_rows):
            y=tf_traj_s.obs_vector(var_names[i])
            ax[i].scatter(x,y,s=1, c=range(tf_traj_s.shape[0]), cmap='viridis')
            ax[i].plot(x,SmoothMat(y, 500),c='k',linewidth=1)
            ax[i].set_title(var_names[i])
            plt.tight_layout() # avoid subplots overlap
    else:
        fig, ax = plt.subplots(figsize=(8,2.5))
        y=tf_traj_s.obs_vector(var_names)
        ax.scatter(x,y,s=1, c=range(tf_traj_s.shape[0]), cmap='viridis')
        ax.plot(x,SmoothMat(y, 500),c='k',linewidth=1)
        ax.set_title(var_names)
        plt.tight_layout()
        
def plot_confusion_matrix(cm, x_classes=None, y_classes=None,
                          normalize=False,
                          title='',
                          cmap=plt.cm.Blues,
                          figsize=(4,4),
                          mark=True,
                          save=None,
                          rotation=45,
                          show_cbar=True,
                          show_xticks=True,
                          show_yticks=True,
                        ):
    """
    This function prints and plots the confusion matrix.
    Normalization can be applied by setting `normalize=True`.

    Params:
        cm: confusion matrix, MxN 
        x_classes: N
        y_classes: M
    """
    import itertools
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes


    if normalize:
        cm = cm.astype('float') / cm.sum(axis=0)[np.newaxis, :]

    fig = plt.figure(figsize=figsize)
    plt.imshow(cm, interpolation='nearest', cmap=cmap)

    plt.title(title)
    x_tick_marks = np.arange(len(x_classes))
    y_tick_marks = np.arange(len(y_classes))
    plt.xticks(x_tick_marks, x_classes, rotation=rotation, ha='right')
    plt.yticks(y_tick_marks, y_classes)
    
    ax=plt.gca()
    if not show_xticks:
        ax.axes.get_xaxis().set_ticks([])
        ax.axes.get_xaxis().set_ticklabels([])
    if not show_yticks:
        ax.axes.get_yaxis().set_ticks([])
        ax.axes.get_yaxis().set_ticklabels([])
    else:
        plt.ylabel('Predicted Cluster')


    fmt = '.2f' if normalize else 'd'
    thresh = cm.max() / 2.
    if mark:
        for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
            if cm[i, j] > 0.1:
                plt.text(j, i, format(cm[i, j], fmt),
                     horizontalalignment="center",
                     color="white" if cm[i, j] > thresh else "black")

    plt.tight_layout()
    if show_cbar:
        plt.colorbar(shrink=0.8) 
    if save:
        plt.savefig(save, format='pdf', bbox_inches='tight')
    plt.show()
