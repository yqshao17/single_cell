import pandas as pd
from scipy.io import mmread, mmwrite
from scipy.sparse import csr_matrix
import numpy as np
import matplotlib.pyplot as plt
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
        elif 'barcode' in basename or 'cell' in basename:
            cell_id = pd.read_csv(filename, sep='\t', header=None)[0].values
        elif 'gene' in basename or 'peak' in basename:
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

def read_mtxR(prefix):
    count = mmread(prefix+'_DGE.mtx').todense()
    cell_id = pd.read_csv(prefix+'_cell_metadata.csv', sep=',')['cell_barcode'].values
    feature = pd.read_csv(prefix+'_genes.csv', sep=',')['gene_id'].values
    df = pd.DataFrame(count, cell_id, feature).T
    return df

def write_mtx(df, path):
    if not os.path.exists(path):
        os.makedirs(path)
    mmwrite(path+'/count.mtx', csr_matrix(df))
    np.savetxt(path+'/barcodes.txt',df.columns.values,fmt="%s")
    np.savetxt(path+'/peaks.txt',df.index.values,fmt="%s")

from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.cluster import KMeans
from sklearn.mixture import GaussianMixture
import scipy.cluster.hierarchy as sch
from sklearn.cluster import AgglomerativeClustering

def kmeans_cluster(X,k):
    # X: n_sample*n_feature
    y_pred = KMeans(n_clusters=k, random_state=124).fit_predict(X)
    return y_pred
def gmm_cluster(X,k):
    y_pred=GaussianMixture(n_components=k).fit_predict(X)
    return y_pred
def hc_cluster(X,k):
    hc = AgglomerativeClustering(n_clusters=k, affinity = 'euclidean', linkage = 'ward')
    y_pred = hc.fit_predict(X)
    return y_pred

def cluster_tsne_pca(df,k, method='kmeans'):
    # df: n_sample*n_feature
    pca = PCA(n_components=30)
    pca.fit(df)
    X = pca.transform(df)    
    X2 = TSNE(n_components=2).fit_transform(X)
    if method=='kmeans':
        y_pred = kmeans_cluster(X2,k)
    elif method=='gmm':
        y_pred = gmm_cluster(X2,k)
    elif method=='hc':
        y_pred = hc_cluster(X2,k)
    plt.scatter(X2[:,0],X2[:,1],alpha=0.5,s=10,c=y_pred)
    return X, X2, y_pred    




def get_prefix(prefix_list):
    # items split by space or tab
    # no space in each item name
    names=[]
    with open(prefix_list) as input:
        for line in input:
            line=line.strip()
            if line:
                names.append(line.split()[0])
    return names

def get_sample_bc(sample_bc_list):
	sample_bc_dic={}
	with open(sample_bc_list) as input:
		for line in input:
			if line.strip():
				[sample, bc]=line.strip().split()
				sample_bc_dic[bc]=sample
	return sample_bc_dic

def plot_cor(x1,x2,tl1,tl2, method='pearson',normalize=None,save=None):
    if normalize:
        x1=x1/sum(x1)*normalize
        x2=x2/sum(x2)*normalize
    l1=np.log2(np.array(x1)+1)
    l2=np.log2(np.array(x2)+1)
    if method=='pearson':
        cor=scipy.stats.pearsonr(l1, l2)
    elif method=='spearman':
        cor=scipy.stats.spearmanr(l1, l2)
    print('Correlation',cor)
    plt.figure(figsize=(6,6))
    plt.scatter(l1, l2, s=4, alpha=0.6)
    plt.text(min(l1), max(l2)*0.9,'cor=%.2f\np=%.4f\nN=%d'%(cor[0],cor[1],len(l1)))
    plt.xlabel(tl1)
    plt.ylabel(tl2)
    if save:
        plt.savefig(save)
    else:
        plt.show()
    return cor


def peak_cor(df1,df2,tl1,tl2, method='pearson',normalize=None,save=None):
    # n_features*n_barcodes
    # normalize by barcodes
    if normalize:
        df1=df1/df1.sum(axis=0)*normalize
        df2=df2/df2.sum(axis=0)*normalize
    peak_count1=df1.sum(axis=1) 
    peak_count2=df2.sum(axis=1)
    peak_count_dict1=peak_count1.to_dict()
    peak_count_dict2=peak_count2.to_dict()
    count1=[]
    count2=[]
    for peak in peak_count_dict1.keys():
        if (peak in peak_count_dict1) and (peak in peak_count_dict2):
            # get peaks in both intersect_bed
            count1.append(peak_count_dict1[peak])
            count2.append(peak_count_dict2[peak])
    l1=np.log2(np.array(count1)+1)
    l2=np.log2(np.array(count2)+1)
    if method=='pearson':
        cor=scipy.stats.pearsonr(l1, l2)
    elif method=='spearman':
        cor=scipy.stats.spearmanr(l1, l2)
    print('Correlation of peak counts',cor)
    plt.figure(figsize=(6,6))
    plt.scatter(l1, l2, s=4, alpha=0.6)
    plt.text(min(l1), max(l2)*0.9,'cor=%.2f\np=%.4f\nN=%d'%(cor[0],cor[1],len(l1)))
    plt.xlabel(tl1)
    plt.ylabel(tl2)
    if save:
        plt.savefig(save)
    else:
        plt.show()
    return cor
