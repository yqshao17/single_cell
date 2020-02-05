import numpy as np
import sys
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats
from scipy.stats import gaussian_kde

def get_outlier_indices(data, max_deviation=200):
        """
        The method is based on the median absolute deviation. See
        Boris Iglewicz and David Hoaglin (1993),
        "Volume 16: How to Detect and Handle Outliers",
        The ASQC Basic References in Quality Control:
        Statistical Techniques, Edward F. Mykytka, Ph.D., Editor.

        returns the list, without the outliers

        The max_deviation=200 is like selecting a z-score
        larger than 200, just that it is based on the median
        and the median absolute deviation instead of the
        mean and the standard deviation.
        """
        median = np.median(data)
        b_value = 1.4826  # value set for a normal distribution
        mad = b_value * np.median(np.abs(data))
        outliers = []
        if mad > 0:
            deviation = abs(data - median) / mad
            """
            outliers = data[deviation > max_deviation]
            print "outliers removed {}".format(len(outliers))
            print outliers
            """
            outliers = np.flatnonzero(deviation > max_deviation)
        return outliers

def remove_outliers(matrix, verbose=True):
        """
        get the outliers *per column* using the median absolute
        deviation method

        Returns the filtered matrix
        """

        unfiltered = len(matrix)
        to_remove = None
        for col in matrix.T:
            outliers = get_outlier_indices(col)
            if to_remove is None:
                to_remove = set(outliers)
                #to_remove = outliers
            else:
                # only set to remove those bins in which
                # the outliers are present in all cases (colums)
                # that's why the intersection is used
                to_remove = to_remove.intersection(outliers)
                #to_remove = np.union1d(to_remove, outliers)
        if len(to_remove):
            to_keep = [x for x in range(matrix.shape[0])
                       if x not in to_remove]
            matrix = matrix[to_keep, :]
            if verbose:
                sys.stderr.write(
                    "total/filtered/left: "
                    "{}/{}/{}\n".format(unfiltered,
                                        unfiltered - len(to_keep),
                                        len(to_keep)))

        return matrix
    
def remove_rows_of_zeros(matrix):
        # remove rows containing all zeros or all nans
        _mat = np.nan_to_num(matrix)
        to_keep = _mat.sum(1) != 0
        matrix = matrix[to_keep, :]
        return matrix


def plot_cor(x1,x2,tl1,tl2, method='pearson', removeOutliers=True,
             normalize=None,save=None,with_density=False):
    x=np.array(x1)
    y=np.array(x2)
    matrix=np.array([x,y]).T
    matrix=remove_rows_of_zeros(matrix)
    if removeOutliers:
        matrix=remove_outliers(matrix)
    x=matrix.T[0]
    y=matrix.T[1]
    if normalize:
        x=x/sum(x)*normalize
        y=y/sum(y)*normalize
    x=np.log2(np.array(x)+1)
    y=np.log2(np.array(y)+1)
    if method=='pearson':
        cor=scipy.stats.pearsonr(x, y)
    elif method=='spearman':
        cor=scipy.stats.spearmanr(x,y)
    print('Correlation',cor)
    
    
    # to slow
    if with_density:
        plt.figure(figsize=(8,6))
        xy = np.vstack([x,y])
        z = gaussian_kde(xy)(xy)
        plt.scatter(x,y,c=z,s=4,alpha=0.5)
        plt.colorbar()
        '''
        if ax is None :
            fig , ax = plt.subplots()
        data , x_e, y_e = np.histogram2d( x, y, bins = 20)
        z = interpn( ( 0.5*(x_e[1:] + x_e[:-1]) , 0.5*(y_e[1:]+y_e[:-1]) ) , data , np.vstack([x,y]).T , method = "splinef2d", bounds_error = False )
        idx = z.argsort()
        x, y, z = x[idx], y[idx], z[idx]
        plt.scatter(x,y,c=z,s=4,alpha=0.5)

        '''
    else:
        plt.figure(figsize=(6,6))
        plt.scatter(x, y, s=4, alpha=0.6)
    lim = max(max(x), max(y))*1.1
    plt.xlim(0, lim)
    plt.ylim(0, lim)
    plt.text(0.5, max(y)*0.9,'cor=%.2f\np=%.4f\nN=%d'%(cor[0],cor[1],len(x)))
    plt.xlabel(tl1)
    plt.ylabel(tl2)
    
    if save:
        plt.savefig(save)
    else:
        plt.show()
    return cor