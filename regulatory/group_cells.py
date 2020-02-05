# group cells for sparse data

import numpy as np
from sklearn.decomposition import PCA
from sklearn.decomposition import TruncatedSVD
#from sklearn.cluster import KMeans
from sklearn.manifold import TSNE



# Dimension reduction 
def reduce_dim(X, method, dim):
    # X: n_sample*n_feature
    if method=='PCA':
        pca = PCA(n_components=dim)
        pca.fit(df)
        X = pca.transform(df)
    elif method=='LSA':
        svd = TruncatedSVD(n_components=dim)
        svd.fit(X)  
        X=svd.fit_transform(X)
    elif method=='TSNE':
        X = TSNE(n_components=dim).fit_transform(X)
    return X



# K nearest neighbour



# Lovain cluster

# group cells

