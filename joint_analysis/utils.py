# read mtx
import pandas as pd
from scipy.io import mmread, mmwrite
from scipy.sparse import csr_matrix
import numpy as np

from glob import glob
import os
import time

def load_data(path, transpose=False):
    print("Loading  data ...")
    t0 = time.time()
    if os.path.isdir(path):
        df = read_mtx(path)
    elif os.path.isfile(path):
        df = read_csv(path)
    else:
        raise ValueError("File {} not exists".format(path))     
    if transpose:
        df = df.T
    print('Original data contains {} features x {} cells'.format(*df.shape))
    print("Finished loading takes {:.2f} min".format((time.time()-t0)/60))
    return df

def read_mtx(path):
    for filename in glob(path+'/*'):
        basename = os.path.basename(filename)
        if (('count' in basename) or ('matrix' in basename)) and ('mtx' in basename):
            count = mmread(filename).todense()
        elif ('barcode' in basename) or ('cell' in basename):
            cell_id = pd.read_csv(filename, sep='\t', header=None)[0].values
        elif ('gene' in basename) or ('peak' in basename) or ('feature' in basename):
            feature = pd.read_csv(filename, sep='\t', header=None)[0].values
    df = pd.DataFrame(count, feature, cell_id)
    return df

def read_csv(path):
    if ('.txt' in path) or ('tsv' in path):
        sep = '\t'
    elif '.csv' in path:
        sep = ','
    else:
        raise ValueError("File {} not in format txt or csv".format(path))
    df = pd.read_csv(path, sep=sep, index_col=0).astype('float32')
    return df

def write_mtx(df, path):
    if not os.path.exists(path):
        os.makedirs(path)
    mmwrite(path+'/count.mtx', csr_matrix(df))
    np.savetxt(path+'/barcodes.txt',df.columns.values,fmt="%s")
    np.savetxt(path+'/peaks.txt',df.index.values,fmt="%s")
    

def overlap_adata(data1, data2):
    features=list(set(data1.var_names)&set(data2.var_names)) 
    barcodes=list(set(data1.obs_names)&set(data2.obs_names))       
    data1=data1[barcodes, :]
    data2=data2[barcodes, :]
    data1=data1[:, features]
    data2=data2[:, features]
    return data1, data2
    
# plot umi
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import seaborn as sns

def plot_umi(df_list, samples=None, save=None):
    # input df dataframe or dataframe list
    # n_feature*n_sample
    if isinstance(df_list, pd.DataFrame):
        df=df_list
        umi=df.sum(axis=0)
        plt.boxplot(umi)
        #sns.kdeplot(np.array(umi), cumulative=False)
    else:
        umi_merge=pd.DataFrame()
        if not samples:
            samples=['sample'+str(i) for i in range(len(df_list))]
        for df, sample in zip(df_list,samples):
            umi=df.sum(axis=0)
            umi=pd.DataFrame(umi, columns=['umi'])
            umi['sample']=sample
            umi_merge=umi_merge.append(umi)
        sns.violinplot(x="sample",y="umi",data=umi_merge,palette="Set3")
    if save:
        plt.savefig(save)
    else:
        plt.show()
    return umi_merge
    
# plot dropout (non-zero peak rate in cells)
def plot_dropout(df_list, samples=None, save=None):
    # n_feature*n_sample
    if isinstance(df_list, pd.DataFrame):
        df=df_list
        droprate=1-df[df>0].count(axis=0)/df.shape[0]
        plt.boxplot(droprate)
    else:
        drop_merge=pd.DataFrame()
        if not samples:
            samples=['sample'+str(i) for i in range(len(df_list))]
        for df,sample in zip(df_list, samples):
            droprate=1-df[df>0].count(axis=0)/df.shape[0]
            droprate=pd.DataFrame(droprate, columns=['sparsity'])
            droprate['sample']=sample
            drop_merge=drop_merge.append(droprate)
        sns.boxplot(x="sample",y="sparsity",data=drop_merge,palette="Set3")
    if save:
        plt.savefig(save)
    else:
        plt.show()
        
# plot non-zero peak num
def plot_dropout_n(df_list, samples=None, save=None):
    # n_feature*n_sample
    if isinstance(df_list, pd.DataFrame):
        df=df_list
        peak_num=df[df>0].count(axis=0)
        plt.boxplot(droprate)
    else:
        drop_merge=pd.DataFrame()
        if not samples:
            samples=['sample'+str(i) for i in range(len(df_list))]
        for df,sample in zip(df_list, samples):
            droprate=1-df[df>0].count(axis=0)/df.shape[0]
            droprate=pd.DataFrame(droprate, columns=['true_peak_num'])
            droprate['sample']=sample
            drop_merge=drop_merge.append(droprate)
        sns.boxplot(x="sample",y="true_peak_num",data=drop_merge,palette="Set3")
    if save:
        plt.savefig(save)
    else:
        plt.show()
        
# plot dropout (cell num of peaks)
def plot_dropout_peak(df_list, samples=None, save=None):
    # n_feature*n_sample
    if isinstance(df_list, pd.DataFrame):
        df=df_list
        droprate=1-df[df>0].count(axis=1)/df.shape[1]
        plt.boxplot(droprate)
    else:
        drop_merge=pd.DataFrame()
        if not samples:
            samples=['sample'+str(i) for i in range(len(df_list))]
        for df,sample in zip(df_list, samples):
            droprate=1-df[df>0].count(axis=1)/df.shape[1]
            droprate=pd.DataFrame(droprate, columns=['sparsity'])
            droprate['sample']=sample
            drop_merge=drop_merge.append(droprate)
        sns.boxplot(x="sample",y="sparsity",data=drop_merge,palette="Set3")
    if save:
        plt.savefig(save)
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

def plot_heatmap(X, y, classes=None, y_pred=None, row_labels=None, colormap=None, row_cluster=False,
                 cax_title='', xlabel='', ylabel='', yticklabels='', legend_font=10, 
                 show_legend=True, show_cax=True, tick_color='black', ncol=3,
                 bbox_to_anchor=(0.5, 1.3), position=(0.8, 0.78, .1, .04), return_grid=False,
                 save=None, **kw):
    """
    plot hidden code heatmap with labels

    Params:
        X: fxn array, n is sample number, f is feature
        y: a array of labels for n elements or a list of array
    """

    import matplotlib.patches as mpatches  # add legend
    # if classes is not None:
    X, y, classes, index = sort_by_classes(X, y, classes)
    # else:
        # classes = np.unique(y)

    if y_pred is not None:
        y_pred = y_pred[index]
        classes = list(classes) + list(np.unique(y_pred)) 
        if colormap is None:
            colormap = plt.cm.tab20
            colors = {c:colormap(i) for i, c in enumerate(classes)}
        else:
            colors = {c:colormap[i] for i, c in enumerate(classes)}
        col_colors = []
        col_colors.append([colors[c] for c in y])
        col_colors.append([colors[c] for c in y_pred])
    else:
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

    cbar_kws={"orientation": "horizontal"}
    grid = sns.clustermap(X, yticklabels=True, 
            col_cluster=False,
            row_cluster=row_cluster,
            cbar_kws=cbar_kws, **kw)
    if show_cax:
        grid.cax.set_position(position)
        grid.cax.tick_params(length=1, labelsize=4, rotation=0)
        grid.cax.set_title(cax_title, fontsize=6, y=0.35)

    if show_legend:
        grid.ax_heatmap.legend(loc='upper center', 
                               bbox_to_anchor=bbox_to_anchor, 
                               handles=legend_TN, 
                               fontsize=legend_font, 
                               frameon=False, 
                               ncol=ncol)
        grid.ax_col_colors.tick_params(labelsize=6, length=0, labelcolor='orange')

    if (row_cluster==True) and (yticklabels is not ''):
        yticklabels = yticklabels[grid.dendrogram_row.reordered_ind]

    grid.ax_heatmap.set_xlabel(xlabel)
    grid.ax_heatmap.set_ylabel(ylabel, fontsize=8)
    grid.ax_heatmap.set_xticklabels('')
    grid.ax_heatmap.set_yticklabels(yticklabels, color=tick_color)
    grid.ax_heatmap.yaxis.set_label_position('left')
    grid.ax_heatmap.tick_params(axis='x', length=0)
    grid.ax_heatmap.tick_params(axis='y', labelsize=6, length=0, rotation=0, labelleft=True, labelright=False)
    grid.ax_row_dendrogram.set_visible(False)
    grid.cax.set_visible(show_cax)
    grid.row_color_labels = classes

    if save:
        plt.savefig(save, bbox_inches='tight')
    else:
        plt.show()
    if return_grid:
        return grid
