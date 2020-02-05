import sys
sys.path.insert(1, '/home/xionglei/yanqiu/joint_analysis/')
import pyranges as pr
from pyranges import PyRanges
from identify_CRE import count_snp_region,get_rand_link
import numpy as np
import pandas as pd
import pickle
import scipy.stats as stats

gene_cre_pairs=pd.read_csv('gene_CRE_pairs.txt',sep='\t',header=None, index_col=None)
gene_cre_pairs.columns=['gene','CRE']
gene_cre_pairs['Chromosome']=gene_cre_pairs['CRE'].apply(lambda x: x.split(':')[0])
gene_cre_pairs['Start']=gene_cre_pairs['CRE'].apply(lambda x: x.split(':')[1].split('-')[0])
gene_cre_pairs['End']=gene_cre_pairs['CRE'].apply(lambda x: x.split(':')[1].split('-')[1])
cre_grouped=gene_cre_pairs.groupby('gene')

gene_around_peaks=pickle.load(open('gene_peak_cor/gene_around_peaks_250k.pkl','rb'))

def enrich(f, name, outdir):
    eqtl=pd.read_csv(f,sep='\t',header=0, index_col=None)
    eqtl['gene']=eqtl['gene_id']
    eqtl['chrom']=eqtl['variant_id'].apply(lambda x: x.split('_')[0])
    eqtl['pos']=eqtl['variant_id'].apply(lambda x: x.split('_')[1])
    eqtl_grouped=eqtl.groupby('gene')
    total=eqtl.shape[0]
    
    genes=list(set(gene_cre_pairs['gene'].unique())&set(eqtl['gene'].unique()))
    eqtl_in_cre=0
    for gene in genes:
        cre_gene=cre_grouped.get_group(gene)[['Chromosome','Start','End']]
        cre_gene=pr.PyRanges(cre_gene)
        eqtl_gene=eqtl_grouped.get_group(gene)[['chrom','pos']].values
        eqtl_in_cre+=count_snp_region(eqtl_gene, cre_gene)
    
    # randomly assign gene-around peaks to gene
    #eqtls_random=[]
    with open('%s/%s_eqtl_enrichment.txt'%(outdir,name),'w') as output:
        
        output.write('eqtl_in_cre\teqtl_in_random\teqtl_total\toddsratio\tpval')
        
        for i in range(10):
            random_df=get_rand_link(gene_around_peaks, gene_cre_pairs.shape[0])
            random_grouped=random_df.groupby('gene')
            genes=list(set(random_df['gene'].unique())&set(eqtl['gene'].unique()))
            eqtl_in_random=0
            for gene in genes:
                cre_gene=random_grouped.get_group(gene)[['Chromosome','Start','End']]
                cre_gene=pr.PyRanges(cre_gene)
                eqtl_gene=eqtl_grouped.get_group(gene)[['chrom','pos']].values
                eqtl_in_random+=count_snp_region(eqtl_gene, cre_gene)

            oddsratio,pval=stats.fisher_exact([[eqtl_in_cre, total-eqtl_in_cre], [eqtl_in_random, total-eqtl_in_random]])
    
            output.write('\n%d\t%d\t%d\t%f\t%f'%(eqtl_in_cre, eqtl_in_random, total, oddsratio, pval))
        
        #eqtls_random.append(eqtl_in_random)

    #eqtl_in_random=np.mean(eqtls_random)    
    
    

from glob import glob
import os
from multiprocessing import Process

eqtl_data={}
outdir='eqtl_enrich'
files=glob('/home/xionglei/yanqiu/lib/eQTL/GTEx_Analysis_v8_eQTL/*signif_variant_gene_pairs.txt.gz')
for filename in files:
    basename = os.path.basename(filename)
    sample=basename.split('.')[0]
    eqtl_data[sample]=filename

nthreads=len(list(eqtl_data.keys()))
Pros=[]
i=0
for sample in eqtl_data.keys():
    i+=1
    print('Starting processing %d: %s' %(i,sample))
    p = Process(target=enrich, args=(eqtl_data[sample], sample, outdir))
    Pros.append(p)
    p.start()
for t in Pros:
    t.join()
    


    
