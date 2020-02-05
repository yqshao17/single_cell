import pyranges as pr
from pyranges import PyRanges
import random
import scipy.stats as stats

import sys
script_path='/home/xionglei/yanqiu/regulatory'
sys.path.insert(1,script_path)
sys.path.insert(1,'/home/xionglei/yanqiu')
import os
import numpy as np
import pandas as pd
import scanpy as sc
import pickle
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

from gene_peak_fit import get_gene_peak_fit, get_peak


def get_link(gene_info, rn, an, outdir, bin_size, save=None):
    # get peaks around gene high correlated with gene expression
    # or
    # get peaks around gene high correlated with gene promoter accessibility
    #gene_peak_fit = get_gene_peak_fit(gene_info, rn, an, outdir,250,'pearson')
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    gene_peak_fit = get_gene_peak_fit(gene_info, rn, an, outdir,250,'spearman')
    n_peaks=[len(peaks) for peaks in gene_peak_fit['peaks'].values()]
    n_peaks_f=[df.shape[0] for df in gene_peak_fit['filtered_peaks'].values()]
    bin_num=20
    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(8,4))
    ax[0].hist(n_peaks, bin_num, facecolor='green')
    ax[0].set_xlabel('peak number around gene')
    ax[1].hist(n_peaks_f, bin_num, facecolor='green')
    ax[1].set_xlabel('filtered peak number')
    if save:
        plt.savefig(outdir+'/gene_peak_fit_basic.png',bbox_inches='tight')
    else:
        plt.show()
    return gene_peak_fit

def identify_cre(rn, an, pro, gene_info, outdir):
    # identify atac peaks correlated to gene expression or promoter accessibility or both
    # rn, an, pro should have the same obs_names
    assert list(rn.obs_names)==list(an.obs_names)==list(pro.obs_names), "Error: input adata do not have the same obs_names"
    # filter all-zero features
    sc.pp.filter_genes(rn, min_cells=1)
    sc.pp.filter_genes(an, min_cells=1)
    sc.pp.filter_genes(pro, min_cells=1)
    
    gene_peak_fit_gene = get_link(gene_info, rn, an, outdir+'/gene_peak_cor/', 250, save=True)
    gene_peak_fit_pro = get_link(gene_info, pro, an, outdir+'/promoter_peak_cor/', 250, save=True)

    overlap_link_f=[]
    for gene in gene_peak_fit_gene['filtered_peaks']:
        if gene in gene_peak_fit_pro['filtered_peaks']:
            for peak in gene_peak_fit_gene['filtered_peaks'][gene].index.values:
                if peak in gene_peak_fit_pro['filtered_peaks'][gene].index.values:
                    overlap_link_f.append((gene, peak))

    df=pd.DataFrame(overlap_link_f)
    df.to_csv(outdir+'/gene_CRE_pairs.txt', header=None, sep='\t', index=None)
    CRE_peaks=df[1].unique()
    CRE_peaks_=[get_peak(p) for p in CRE_peaks]
    CRE_df=pd.DataFrame(CRE_peaks_, columns=['Chromosome','Start','End'])
    CRE_df.to_csv(outdir+'/CRE.bed', header=None, sep='\t', index=None)

    # plot venn
    n_peaks_gene=[df.shape[0] for df in gene_peak_fit_gene['filtered_peaks'].values()]
    n_peaks_pro=[df.shape[0] for df in gene_peak_fit_pro['filtered_peaks'].values()]
    link_gene_peak=sum(n_peaks_gene)
    link_pro_peak=sum(n_peaks_pro)
    link_overlap=len(overlap_link_f)
    plt.figure()
    venn2(subsets = (link_gene_peak-link_overlap, link_pro_peak-link_overlap, link_overlap), 
          set_labels = ('gene-CRE links', 'promoter-CRE links'))
    plt.savefig(outdir+'/gene_promoter_CRE_overlap.png',bbox_inches='tight')

    # plot distance
    distance=[]
    for gene, cre in overlap_link_f:
        gene_center=(gene_info['gene_starts'][gene]+gene_info['gene_ends'][gene])/2
        strand=gene_info['gene_id_to_strand'][gene]
        cre_start, cre_end=cre.split(':')[1].split('-')
        cre_center=(int(cre_start)+int(cre_end))/2
        if strand=='+':
            dist=cre_center-gene_center
        else:
            dist=gene_center-cre_center
        distance.append(dist)
    plt.figure()
    x=plt.hist(distance, 100, facecolor='brown')
    plt.title('Distance between CRE and linked genes')
    plt.savefig(outdir+'/Distance_CRE_genes.png',bbox_inches='tight')

    # plot n_genes_per_CRE and n_CRE_per_genes
    n_peaks_per_gene=df[0].value_counts()
    n_genes_per_peak=df[1].value_counts()
    bin_num=20
    plt.figure()
    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(8,4))
    ax[0].hist(n_peaks_per_gene, bin_num, facecolor='green')
    ax[0].set_xlabel('Number of CREs linked to each gene')
    ax[1].hist(n_genes_per_peak, bin_num, facecolor='darkblue', alpha=0.8)
    ax[1].set_xlabel('Number of genes linked to each CRE')
    plt.savefig(outdir+'/identified_CREs_basic.png',bbox_inches='tight')

    # plot length dist
    length_dist=CRE_df.apply(lambda x: x['End']-x['Start'], axis=1)
    length=int(np.mean(length_dist))
    plt.figure()
    x=plt.hist(length_dist, 20, facecolor='green')
    plt.title('length_dist: %d'%length)
    plt.savefig(outdir+'/CRE_length_dist.png',bbox_inches='tight')
    
    return CRE_df
    



#===============
# SNP enrichment
#===============

def count_snp_region(snps, pr):
    # count the number of snps in pr regions
    # snps: [(chrm, pos),(),...]
    # regions: pr
    num_in=0
    for [chrm, pos] in snps:
        if len(pr[chrm,int(pos):int(pos)+1])>0:
            num_in+=1
    return num_in


def get_rand_region(length, count, method, bin_size=250, gene_info=None, atac_peak=None):
    
    # get random regions around genes
    def get_random_cre(gene_info, length, bin_size, count):
        bin_size=bin_size*1000
        random_cres=[]
        for i in range(count):
            gene=random.choice(list(gene_info['gene_id_to_name'].keys()))       
            chrm=gene_info['gene_id_to_chrom'][gene].split('_')[-1]
            random_start=random.randrange(gene_info['gene_starts'][gene]-bin_size, gene_info['gene_ends'][gene]+bin_size-length)
            end=random_start+length
            random_cres.append((chrm, random_start, end))
        random_df=pd.DataFrame(random_cres)
        random_df.columns=['Chromosome','Start','End']
        return random_df

    # get random genome regions
    def get_random_region(length, count):
        genome_size='/home/xionglei/yanqiu/lib/hg38.chrom.sizes'
        genome_size=pd.read_csv(genome_size,sep='\t',header=None, index_col=0)[1]
        choose_prob=(genome_size.cumsum()/genome_size[1].sum()).to_dict()
        random_regions=[]
        for i in range(count):
            p=random.random()
            if 0<=p<=list(choose_prob.values())[0]:
                chrm=list(choose_prob.keys())[0]
                random_start=random.randrange(0, genome_size[0]-length)
            else:
                for i in range(1, genome_size.shape[0]):
                    if list(choose_prob.values())[i-1]<p<=list(choose_prob.values())[i]:
                        chrm=list(choose_prob.keys())[i]
                        random_start=random.randrange(0, genome_size[1][i]-length)
            end=random_start+length        
            random_regions.append((chrm, random_start,end))
        random_df=pd.DataFrame(random_regions)
        random_df.columns=['Chromosome','Start','End']
        return random_df

    def get_random_atac(atac_peak,count):
        n_peaks=atac_peak.shape[0]
        random_index = []
        for j in range(count): 
            random_index.append(random.randint(0, n_peaks-1)) 
        random_df=atac_peak.iloc[random_index, :]
        return random_df
        

    if method=='gene_around':
        random_df = get_random_cre(gene_info, length, bin_size, count)
    elif method=='genome_wide':
        random_df = get_random_region(length, count)
    elif method=='atac_peak':
        random_df = get_random_atac(atac_peak,count)
    
    random_pr=pr.PyRanges(random_df)

    return random_pr


def get_rand_region(length, count, method, bin_size=250, gene_info=None, atac_peak=None, gene_around_peaks=None):
    
    # get random regions around genes
    def get_random_cre(gene_info, length, bin_size, count):
        bin_size=bin_size*1000
        random_cres=[]
        for i in range(count):
            gene=random.choice(list(gene_info['gene_id_to_name'].keys()))       
            chrm=gene_info['gene_id_to_chrom'][gene].split('_')[-1]
            random_start=random.randrange(gene_info['gene_starts'][gene]-bin_size, gene_info['gene_ends'][gene]+bin_size-length)
            end=random_start+length
            random_cres.append((chrm, random_start, end))
        random_df=pd.DataFrame(random_cres)
        random_df.columns=['Chromosome','Start','End']
        return random_df

    # get random genome regions
    def get_random_region(length, count):
        genome_size='/home/xionglei/yanqiu/lib/hg38.chrom.sizes'
        genome_size=pd.read_csv(genome_size,sep='\t',header=None, index_col=0)[1]
        choose_prob=(genome_size.cumsum()/genome_size[1].sum()).to_dict()
        random_regions=[]
        for i in range(count):
            p=random.random()
            if 0<=p<=list(choose_prob.values())[0]:
                chrm=list(choose_prob.keys())[0]
                random_start=random.randrange(0, genome_size[0]-length)
            else:
                for i in range(1, genome_size.shape[0]):
                    if list(choose_prob.values())[i-1]<p<=list(choose_prob.values())[i]:
                        chrm=list(choose_prob.keys())[i]
                        random_start=random.randrange(0, genome_size[1][i]-length)
            end=random_start+length        
            random_regions.append((chrm, random_start,end))
        random_df=pd.DataFrame(random_regions)
        random_df.columns=['Chromosome','Start','End']
        return random_df

    def get_random_atac(atac_peak,count):
        n_peaks=atac_peak.shape[0]
        random_index = []
        for j in range(count): 
            random_index.append(random.randint(0, n_peaks-1)) 
        random_df=atac_peak.iloc[random_index, :]
        return random_df
    
    def get_random_around_peak(gene_around_peaks, count):        
        total_peaks=[]
        for gene in gene_around_peaks:
            total_peaks+=gene_around_peaks[gene]
        random_peaks=[]
        for i in range(count):
            peak=random.choice(total_peaks)
            random_peaks.append(get_peak(peak))
        random_df=pd.DataFrame(random_peaks)
        random_df.columns=['Chromosome','Start','End']
        return random_df
        
    if method=='gene_around':
        random_df = get_random_cre(gene_info, length, bin_size, count)
    elif method=='genome_wide':
        random_df = get_random_region(length, count)
    elif method=='atac_peak':
        random_df = get_random_atac(atac_peak,count)
    elif method=='around_peak':
        random_df = get_random_around_peak(gene_around_peaks, count)
            
    random_pr=pr.PyRanges(random_df)

    return random_pr


def get_rand_link(gene_around_peaks, count):
    total_pairs=[]
    for gene in gene_around_peaks:
        peaks=gene_around_peaks[gene]
        for peak in peaks:
            total_pairs.append((gene, peak))
    random_pairs=[]
    for i in range(count):
        random_pairs.append(random.choice(total_pairs))
    random_df=pd.DataFrame(random_pairs)
    random_df.columns=['gene','peak']
    random_df['Chromosome']=random_df['peak'].apply(lambda x: x.split(':')[0])
    random_df['Start']=random_df['peak'].apply(lambda x: x.split(':')[1].split('-')[0])
    random_df['End']=random_df['peak'].apply(lambda x: x.split(':')[1].split('-')[1])   
    return random_df