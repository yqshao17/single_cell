# snp enrichment
import sys
import os
sys.path.insert(1, '/home/xionglei/yanqiu/joint_analysis/')
import pyranges as pr
from pyranges import PyRanges
import random
import scipy.stats as stats
from identify_CRE import get_rand_region, count_snp_region
import pickle
import numpy as np
import pandas as pd
import scipy.stats as stats
from multiprocessing import Process


#def snp_enrich(CRE_bed, snps_file, gene_info, outdir):
def snp_enrich(CRE_bed, snps_file, gene_info, outdir, method, gene_around_peaks=None, atac_peak=None):
    # gene_around_peaks is needed when method is around_peak
    # atac_peak is needed when method is atac_peak
    gene_info = pickle.load(open(gene_info, 'rb'))

    CRE_df=pd.read_csv(CRE_bed, header=None, sep='\t', index_col=None)
    CRE_df.columns=['Chromosome', 'Start' , 'End']
    length_dist=CRE_df.apply(lambda x: x['End']-x['Start'], axis=1)
    length=int(np.mean(length_dist))
    count=CRE_df.shape[0]
    CRE_pr=pr.PyRanges(CRE_df)

    snps=pd.read_csv(snps_file, sep='\t', header=None)
    snps.columns=['CHRM','pos','ic','ref','alt','qual','filter','info']
    snps['chrm']=snps['CHRM'].apply(lambda x: 'chr'+str(x))
    snps.head()
    total=snps.shape[0]

    snps_in_CRE=count_snp_region(snps[['chrm','pos']].values, CRE_pr)
    print('snps_in_CRE',snps_in_CRE)
    
    print (gene_around_peaks)
    ### random sampling
    def enrich(length, count, method, outdir, gene_around_peaks=None):
        print (gene_around_peaks)
        with open('%s/%s_snp_enrichment.txt'%(outdir,method),'w') as output:

            output.write('snps_in_CRE\tsnps_in_random\tsnps_total\toddsratio\tpval')

            for i in range(10):
                if method=='gene_around':
                    random_pr=get_rand_region(length, count, 'gene_around', bin_size=250, gene_info=gene_info) 
                elif method=='genome_wide':
                    random_pr = get_rand_region(length, count, 'genome_wide')
                elif method=='atac_peak':
                    random_pr=get_rand_region(length, count, 'atac_peak',atac_peak=atac_peak)
                elif method=='around_peak':
                    gene_around_peaks_dict=pickle.load(open(gene_around_peaks,'rb'))
                    random_pr=get_rand_region(length, count, 'around_peak', gene_around_peaks=gene_around_peaks_dict)

                snps_in_random=count_snp_region(snps[['chrm','pos']].values, random_pr)
                oddsratio,pval=stats.fisher_exact([[snps_in_CRE, total-snps_in_CRE], [snps_in_random, total-snps_in_random]])
                output.write('\n%d\t%d\t%d\t%f\t%f'%(snps_in_CRE, snps_in_random, total, oddsratio, pval))

    if type(method)==type([]):
        Pros=[]
        i=0
        for m in method:
            i+=1
            print('Starting processing %d: %s' %(i,method))
            p = Process(target=enrich, args=(length, count, m, outdir,gene_around_peaks))
            Pros.append(p)
            p.start()
        for t in Pros:
            t.join()
    else:
        enrich(length, count, method, outdir,gene_around_peaks)
        
        
if __name__ == '__main__':
    gene_info='/home/xionglei/yanqiu/lib/gene_info_hg.pkl'
    snps_file='/home/xionglei/yanqiu/lib/ClinVar/clinvar_20191125_cancer.vcf'
    '''
    CRE_bed='CRE.bed'    
    gene_around_peaks='gene_peak_cor/gene_around_peaks_250k.pkl'
    outdir='snps_enrich'
    method=['genome_wide','gene_around','around_peak']
    snp_enrich(CRE_bed, snps_file, gene_info, outdir, method, gene_around_peaks=gene_around_peaks)
    '''    
    
    path='/home/xionglei/yanqiu/joint_analysis/WYF191116'
    method='around_peak'
    
    cluster_df=pd.read_csv(path+'/WYF191116_cluster.txt', sep='\t', index_col=0, header=0)
    clusters=cluster_df['clusters'].unique()
    
    # cannot run together??
    Pros=[]
    i=0
    for c in ['3']:
            i+=1
            print('Starting processing %d' %i)
            outdir=path+'/snps_enrich/%s'%c
            if not os.path.exists(outdir):
                os.makedirs(outdir)
            CRE_bed=path+'/CRE/%s/CRE.bed'%c   
            gene_around_peaks=path+'/CRE/%s/gene_peak_cor/gene_around_peaks_250k.pkl'%c  
            p = Process(target=snp_enrich, args=(CRE_bed, snps_file, gene_info, outdir, method, gene_around_peaks))
            Pros.append(p)
            p.start()
    for t in Pros:
            t.join()
            
    #cluster=sys.argv[1]
    #snp_enrich(CRE_bed, snps_file, gene_info, outdir, method, gene_around_peaks)