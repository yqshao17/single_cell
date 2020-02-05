
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn import preprocessing
import numpy as np
import statsmodels.api as sm

def plot_info(gene_peak_regress, save=None):
    fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(12,9))
    #plt.title('Total counts of potential gene-peak pairs: %d'%sum([df.shape[0] for df in gene_peak_regress['filtered_peaks'].values()]))
    bin_num=20
    ax[0,0].hist([len(peaks) for peaks in gene_peak_regress['peaks'].values()], bin_num, facecolor='green', alpha=0.5)
    ax[0,0].set_xlabel('peak number around gene')
    ax[1,0].hist([df.shape[0] for df in gene_peak_regress['filtered_peaks'].values()], bin_num, facecolor='green', alpha=0.5)
    ax[1,0].set_xlabel('filtered peak number')
    ax[0,1].hist(gene_peak_regress['R2'].values(), bin_num, facecolor='green', alpha=0.5)
    ax[0,1].set_xlabel('R2')
    #ax[1,1].hist(gene_peak_regress['R2_adj'].values(), bin_num, facecolor='green', alpha=0.5)
    #ax[1,1].set_xlabel('R2 adjusted')
    if save:
        plt.savefig(save,bbox_inches='tight')
    else:
        plt.show()

def plot_gene_peak(gene_peak_regress, rn, an, save=None):
    
    scaler = preprocessing.StandardScaler()
    # select four top best regressed genes
    R2_sorted=sorted(gene_peak_regress['R2_adj'].items(), key=lambda x: x[1], reverse=True)
    selected_genes=[geneR[0] for geneR in R2_sorted[:4]]
    selected_peaks=[gene_peak_regress['filtered_peaks'][gene].index[0] for gene in selected_genes]
    fig, ax = plt.subplots(nrows=4, ncols=3, figsize=(12,20))
    for gene, peak, i in zip(selected_genes, selected_peaks, range(4)):
        # scatter plot
        ax[i,0].scatter(rn.obs_vector(gene), an.obs_vector(peak), alpha=0.5)
        ax[i,0].set_title(gene+' '+peak)
        # predictions and true y
        predictions=gene_peak_regress['predictions'][gene]
        ax[i,1].scatter(rn.obs_vector(gene), predictions, alpha=0.5)
        ax[i,1].set_title('R2 adjusted: %f'%R2_sorted[i][1])
        # residual distribution
        resid=gene_peak_regress['resid'][gene]
        scaled = scaler.fit_transform(np.array(resid).reshape(-1,1))
        sm.qqplot(scaled.flatten(),line='45', ax=ax[i,2])
        ax[i,2].set_title('resid distribution')
        
    if save:
        plt.savefig(save,bbox_inches='tight')
    else:
        plt.show()
                

if __name__ == '__main__':
    import scanpy as sc
    import pickle
    sys.path.insert(1, '/home/xionglei/yanqiu/')
    from utils import write_mtx
    gene_peak_regress=pickle.load(open('SNARE-seq/ad/gene_peak_regress_250k.pkl','rb'))
    gene_peak_regress.keys()
    rn=sc.read_h5ad('SNARE-seq/ad/rn.h5ad')
    an=sc.read_h5ad('SNARE-seq/ad/an.h5ad')
    reg_df=get_regulatory_matrix(gene_peak_regress['filtered_peaks'],rn,an) 
    write_mtx(reg_df, 'SNARE-seq/ad/reg_df')
