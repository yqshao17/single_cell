# get correlation of peak counts
# input: reference peak, bed1, bed2
# reference peak: chr, start, end, peakname
import matplotlib as mpl 
mpl.use("Agg")
import matplotlib.pyplot as plt
import subprocess
import pandas as pd
import scipy.stats
import numpy as np
import os

def get_peak_count(refpeak, bed, outdir):
	#intersect=subprocess.check_output(['bedtools', 'intersect', '-a', refpeak, '-b', bed, '-wa','-wb'])
	refpeak_name=refpeak.split('/')[-1].split('.')[0]
	bed_name=bed.split('/')[-1].split('.')[0]
	out_intersect=outdir+'/'+refpeak_name+'_'+bed_name+'_intersect.bed'
	#if os.path.exists(out_intersect):
	#	out_intersect=outdir+'/'+refpeak_name+'_'+bed_name+'_intersect2.bed'
	#with open(out_intersect,'wb') as output:
	#	output.write(intersect)
	peak_df=pd.read_table(out_intersect)
	n_col=peak_df.shape[1]
	peak_df.columns=['chr','start','end']+['']*(n_col-3)
	peak_df['peak']=peak_df.apply(lambda x: x['chr']+':'+str(x['start'])+'-'+str(x['end']), axis=1)
	peak_count=peak_df.groupby('peak').size()
	peak_count_dict=peak_count.to_dict()
	return peak_count_dict


def correlation(l1,l2,tl1,tl2,figname):
	cor=scipy.stats.pearsonr(l1, l2)
	print('Correlation of peak counts',cor)
	plt.figure(figsize=(6,6))
	plt.scatter(l1, l2, s=4, alpha=0.6)
	plt.text(min(l1), max(l2)*0.9,'cor=%.2f\np=%.4f\nN=%d'%(cor[0],cor[1],len(l1)))
	plt.xlabel(tl1)
	plt.ylabel(tl2)
	plt.savefig(figname)


def main(refpeak, bed1, bed2, outdiri, normalize=2000000): #peak length is 500, normalize=10**9/500
	bed1_name=bed1.split('/')[-1].split('.')[0]
	bed2_name=bed2.split('/')[-1].split('.')[0]
	peak_count_dict1=get_peak_count(refpeak, bed1, outdir)
	peak_count_dict2=get_peak_count(refpeak, bed2, outdir)
	count1=[]
	count2=[]
	for peak in peak_count_dict1.keys():
		if (peak in peak_count_dict1) and (peak in peak_count_dict2):
			# get peaks in both intersect_bed
			count1.append(peak_count_dict1[peak])
			count2.append(peak_count_dict2[peak])
	if normalize:
		count1=np.array(count1)/sum(count1)*normalize
		count2=np.array(count2)/sum(count2)*normalize
	l1=np.log2(np.array(count1)+1)
	l2=np.log2(np.array(count2)+1)
	correlation(l1,l2, 'log2 peak counts '+bed1_name, 'log2 peak counts '+bed2_name, outdir+'/'+bed1_name+'_'+bed2_name+'_cor.pdf')

if __name__ == '__main__':
	import sys
	refpeak=sys.argv[1]
	bed1=sys.argv[2]
	bed2=sys.argv[3]
	outdir=sys.argv[4]
	main(refpeak, bed1, bed2, outdir)

