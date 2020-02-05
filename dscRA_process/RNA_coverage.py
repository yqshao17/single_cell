import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pyBigWig
from numpy import mean
import pandas as pd
import numpy as np

def rpkm(geneinfo, genecounts):
    # genecounts: n_features*n_samples
    gene_len={}
    for gene in geneinfo['gene_starts']:
        gene_len[gene]=geneinfo['gene_ends'][gene]-geneinfo['gene_starts'][gene]+1 
    genelen=np.array([gene_len[gene] for gene in genecounts.index])
    fpkm=pd.DataFrame()
    for bc in genecounts.columns:
        fpkm[bc]=genecounts[bc]/sum(genecounts[bc])/genelen*10**9
    return fpkm

def bin_coverage(gene_info, genes, bw, bin_num):
    gene_id_to_chrom=gene_info['gene_id_to_chrom']
    gene_starts=gene_info['gene_starts']
    gene_ends=gene_info['gene_ends']
    gene_strands=gene_info['gene_id_to_strand']
    def gene_cov(gene):
        bin_cover=np.zeros(bin_num)
        chrom=gene_id_to_chrom[gene]
        end=gene_ends[gene]
        start=gene_starts[gene]
        strand=gene_strands[gene]
        mean_gene_cov=mean(bw.values(chrom, start, end))
        if end-start>bin_num and mean_gene_cov>0:         
            bin_len=(end-start)/bin_num
            bins=[int(start+ bin_len*i) for i in range(bin_num+1)]
            bin_cover=[mean(bw.values(chrom, bins[i], bins[i+1]))/mean_gene_cov for i in range(bin_num)]
            if strand=='-':
                bin_cover=bin_cover[::-1]
        return np.array(bin_cover)
    '''
    coverage_overall=np.zeros(bin_num)
    for gene in genes:
        bin_cover=gene_coverage(gene, bw, bin_num)
        coverage_overall=coverage_overall+bin_cover
    '''
    coverage_overall=[]
    for gene in genes:
        bin_cover=gene_cov(gene)
        coverage_overall.append(bin_cover)
    coverage_overall=pd.DataFrame(coverage_overall, index=genes)
    return coverage_overall

if __name__=='__main__':
    import pickle
    import sys
    
    bwfile=sys.argv[1]
    gene_file=sys.argv[2]
    bin_num=int(sys.argv[3])
    outprefix=sys.argv[4]
    
    with open('/Share2/home/zhangqf5/yanqiu/scAR_old/split_seq/split-seq-pipeline/INDEX_mm_hg/gene_info.pkl', 'rb') as f:
        gene_info = pickle.load(f)
    bw = pyBigWig.open(bwfile)
    genes=pd.read_csv(gene_file)['gene_id']
    hgenes=[]
    mgenes=[]
    for gene in genes:
        if gene.startswith('ENSG'):
            hgenes.append(gene)
        elif gene.startswith('ENSMUSG'):
            mgenes.append(gene)
    if hgenes:
        hg_coverage_overall=bin_coverage(gene_info, hgenes, bw, bin_num)
        hg_coverage_overall.to_csv(outprefix+'_hg_coverage.txt', sep='\t')
        plt.figure()
        plt.scatter(range(bin_num), hg_coverage_overall.sum(axis=0))
        plt.savefig(outprefix+'_hg_coverage.pdf')
    if mgenes:
        mm_coverage_overall=bin_coverage(gene_info, mgenes, bw, bin_num)
        mm_coverage_overall.to_csv(outprefix+'_mm_coverage.txt', sep='\t')
        plt.figure()
        plt.scatter(range(bin_num), mm_coverage_overall.sum(axis=0))
        plt.savefig(outprefix+'_mm_coverage.pdf')

