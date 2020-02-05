import matplotlib as mpl 
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np 
import pysam
import os
import pandas as pd 
from scipy.stats import gaussian_kde
#import seaborn as sns

def get_UmiCount(samplebam):
	bamfile = pysam.Samfile(samplebam, "rb")
	UmiCount={}
	for p2_rds in bamfile:
		if p2_rds.is_unmapped: continue
		if p2_rds.mapq<30:continue
		if p2_rds.is_reverse:continue
		else:
			bc=p2_rds.query_name.split(':')[0]
			UmiCount.setdefault(bc,0)
			UmiCount[bc]+=1
	return UmiCount

def get_InsertMat(tssbed, samplebam):
	
	InsertMat={}

	def asn_mat(bc,pos,InsertMat,s_int,e_int,ilen,weight):
		InsertMat.setdefault(bc,np.zeros(4000))
		if float(pos)>=s_int and float(pos)<e_int-1 and ilen<1000:
			InsertMat[bc][pos-s_int] += weight
		return InsertMat

	bamfile = pysam.Samfile(samplebam, "rb")
	for line in open(tssbed):
		items=line.strip().split('\t')
		chrx_tss,t1,t2=items[:3]
		tss=int((int(t1)+int(t2))/2)
		s_int=tss-2000
		e_int=tss+2000
		for p2_rds in bamfile.fetch(chrx_tss, max(0,s_int-2000), e_int+2000):
			if p2_rds.is_unmapped: continue
			if p2_rds.mapq<30:continue
			if p2_rds.is_reverse:continue
			else:
				bc=p2_rds.query_name.split(':')[0]
				l_pos = p2_rds.pos+4
				ilen = abs(p2_rds.tlen)-9
				r_pos=p2_rds.pos+abs(p2_rds.tlen)-5
				InsertMat = asn_mat(bc,l_pos,InsertMat,s_int,e_int,ilen,1)
				InsertMat = asn_mat(bc,r_pos,InsertMat,s_int,e_int,ilen,1)
	return InsertMat

def get_TssEnrichmentMat(insert_mat):
	return insert_mat/np.mean(np.hstack((insert_mat[:100],insert_mat[-100:])))

def plot_tss_enrich(tssbed,samplebam, prefix):
	UmiCount=get_UmiCount(samplebam)
	InsertMat=get_InsertMat(tssbed,samplebam)	
	TssEnrich=[] # barcode: float
	UMI=[]
	for bc in InsertMat.keys():
		TssEnrich.append(get_TssEnrichmentMat(InsertMat[bc])[2000])
		UMI.append(UmiCount[bc])

	OverallMat=sum(InsertMat.values())
	OverallTssEnrichMat=get_TssEnrichmentMat(OverallMat)
	TotalUmi=sum(UMI)

	# plot TSS enrichment score vs umi counts
	plt.figure(figsize=(6,6))
	x=np.log(np.array(UMI)+1)
	y=np.array(TssEnrich)
	plt.scatter(x,y,alpha=0.5,s=4)
	plt.xlabel('Log10 UMI counts')
	plt.ylabel('TSS enrichment score')
	plt.title('TSS enrichment score vs UMI counts')
	plt.savefig(prefix+'TSSEnrichmentScore_UMI.pdf')

	# plot TSS enrichment score with color map
	ix=np.where(~np.isnan(y))[0]
	x=x[ix]
	y=y[ix]
	plt.figure()
	xy = np.vstack([x,y])
	z = gaussian_kde(xy)(xy)
	plt.scatter(x,y,c=z,s=4,alpha=0.5)
	plt.colorbar()
	y=np.sort(y)
	ylim=y[int(len(y)*0.9)] # show 90 percent cell barcodes
	plt.ylim(0,ylim)
	plt.xlabel("Log10 UMI")
	plt.ylabel("TSS enrichment")
	plt.savefig(prefix+'TSSEnrichmentScore_colormap.pdf')

	# plot overalll tss enrichment score
	plt.figure(figsize=(6,4))
	plt.scatter(range(-2000,2000),OverallTssEnrichMat,alpha=0.5,s=4,c='k')
	smooth_mat=np.convolve(OverallMat,np.ones(50),'same')/50
	plt.plot(range(-2000,2000),get_TssEnrichmentMat(smooth_mat),'r')
	plt.xlabel('Distance to TSS center')
	plt.ylabel('TSS enrichment score')
	plt.savefig(prefix+'Overall_TssEnrichment%.2f_%d.pdf'%(OverallTssEnrichMat[2000], TotalUmi))

	outdf=pd.DataFrame()
	outdf['barcode']=InsertMat.keys()
	outdf['UMI']=UMI
	outdf['TssEnrichScore']=TssEnrich
	outdf.to_csv(prefix+'TssEnrichUMI.txt', sep='\t', header=True, index=None)


def is_sorted_coordinate(header):                                                                                                                                                                                                                                            
	"""
	Check if bam fiel is sorted by read name.
	"""
	if("HD" in header):
		if("SO" in header["HD"]):
			if(header["HD"]["SO"] == "coordinate"):
				return True
	else:
		return False

def TSS_ENRICH(tssbed,samplebam,prefix):
	num_threads=4
	# check sorted
	sortedbam = None
	samfile = pysam.AlignmentFile(samplebam, "rb")
	if not is_sorted_coordinate(samfile.header):
		sortedbam = samplebam[:-4]+'.sorted.bam'
		os.system('samtools sort %s -o %s -@ %d'%(samplebam, sortedbam, num_threads))
	samfile.close()
	
	if sortedbam:
		os.system('samtools index '+sortedbam)
		plot_tss_enrich(tssbed,sortedbam,prefix)
	else:
		os.system('samtools index '+samplebam)
		plot_tss_enrich(tssbed,samplebam,prefix)


if __name__ == '__main__':

	import sys
	prefix=sys.argv[1]
	samplebam=sys.argv[2]
	tssbed=sys.argv[3] 
	TSS_ENRICH(tssbed,samplebam,prefix)
	#TSS_ENRICH('human_epdnew_TGfph.bed','/Share2/home/zhangqf5/yanqiu/scAR/output/20190823_test/atac_qc/test_filtered.bam','test')
	
	#samtools merge test_filtered.bam *bam
