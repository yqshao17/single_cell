import sys
sys.path.insert(1, '/Share2/home/zhangqf6/yanqiu/bin/SCALE/scale')
from plot import plot_heatmap,feature_specifity
from specifity import cluster_specific, mat_specificity_score

import pandas as pd
from scipy.io import mmread, mmwrite
from scipy.sparse import csr_matrix
import numpy as np

from glob import glob
import os

def read_mtx(path):
    for filename in glob(path+'/*'):
        basename = os.path.basename(filename)
        if (('count' in basename) or ('matrix' in basename)) and ('mtx' in basename):
            count = mmread(filename).todense()
        elif 'barcode' in basename:
            cell_id = pd.read_csv(filename, sep='\t', header=None)[0].values
        elif 'gene' in basename or 'peak' in basename:
            feature = pd.read_csv(filename, sep='\t', header=None)[0].values
    df = pd.DataFrame(count, feature, cell_id)
    return df
    

# plot feature.txt
def plot_feature(output_dir):
	y = pd.read_csv(output_dir+'/cluster_assignments.txt', sep='\t', index_col=0, header=None)[1]
	feature = pd.read_csv(output_dir+'/feature.txt', sep='\t', index_col=0, header=None)
	plot_heatmap(feature.T, y,
             figsize=(8, 3), cmap='RdBu_r', vmax=8, vmin=-8, center=0,
             ylabel='Feature dimension', yticklabels=np.arange(10)+1,
             cax_title='Feature value', legend_font=6, ncol=1,
             bbox_to_anchor=(1.1, 1.1), position=(0.92, 0.15, .08, .04),
             save=output_dir+'/feature.png')

# plot specific feature
def plot_specific_feature(df, y, save):
        score_mat = mat_specificity_score(df, y)
        peak_index, peak_labels = cluster_specific(score_mat, np.unique(y), top=200)
        plot_heatmap(df.iloc[peak_index], y=y, row_labels=peak_labels, ncol=3, cmap='Reds',
                         vmax=1, row_cluster=False, legend_font=6, cax_title='Peak Value',
                         figsize=(8, 10), bbox_to_anchor=(0.4, 1.2), position=(0.8, 0.76, 0.1, 0.015),
                         save=save)

if __name__ == '__main__':
	input_dir=sys.argv[1]
	output_dir=sys.argv[2]
	type=sys.argv[3]
	y = pd.read_csv(output_dir+'/cluster_assignments.txt', sep='\t', index_col=0, header=None)[1]
	if type=='feature':
		plot_feature(output_dir)
	elif type=='raw_specific':
		raw=read_mtx(input_dir)
		imputed = read_mtx(output_dir+'/binary_imputed')
		#imputed = pd.read_csv(output_dir+'/binary_imputed_data.txt', sep='\t', index_col=0)
		raw_filtered=raw.loc[imputed.index,imputed.columns]
		print(raw.shape, imputed.shape, raw_filtered.shape,y.shape)
		plot_specific_feature(raw_filtered, y, output_dir+'/raw_specific_feature.png')
	elif type=='imputed_specific':
		#imputed = pd.read_csv(output_dir+'/binary_imputed_data.txt', sep='\t', index_col=0)
		imputed = read_mtx(output_dir+'/binary_imputed')
		plot_specific_feature(imputed, y, output_dir+'/imputed_specific_feature.png')

