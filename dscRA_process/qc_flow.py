from utils import get_prefix
import re
import pandas as pd

def qc_prefix(indir, prefix, seq_type, trim_tso):
    trim_log=open('%s/trimmed/%s.log'%(indir, prefix)).read()
    raw=re.search("Total read pairs processed:\s+([\d,]+)\n", trim_log).group(1).replace(',','')
    trimmed=re.search("Pairs written.*\s+([\d,]+)\s[(].*\n", trim_log).group(1).replace(',','')
    
    parse_log=open('%s/parse_bc/%s_pipeline_stats.txt'%(indir, prefix)).read()
    bc_parsed=re.search("fastq_valid_barcode_reads\s+([\d,]+)\n", parse_log).group(1)
    splited=re.search("--%s_valid_barcode_reads\s+([\d,]+)\n"%seq_type, parse_log).group(1)
    
    if trim_tso:
        tmp_log=open('%s/trimmed/%s_tmp.log'%(indir, prefix)).read()
        raw=re.search("Total read pairs processed:\s+([\d,]+)\n", tmp_log).group(1).replace(',','')
        with_tso=re.search("Pairs written.*\s+([\d,]+)\s[(].*\n", tmp_log).group(1).replace(',','')
        return raw, with_tso, trimmed, bc_parsed, splited
    else:
        return raw, trimmed, bc_parsed, splited
    
def qc_sample(indir, sample, seq_type):
    split_log=open('%s/split_sample/%s_split_pipeline_stats.txt'%(indir, sample)).read()
    reads_ATAC=re.search("fastq_reads_ATAC\s+([\d,]+)\n", split_log).group(1)
    reads_RNA=re.search("fastq_reads_RNA\s+([\d,]+)\n", split_log).group(1)
    reads_UNKNOWN=re.search("fastq_reads_UNKNOWN\s+([\d,]+)", split_log).group(1)
    total_reads=int(reads_ATAC)+int(reads_RNA)+int(reads_UNKNOWN)
    #total_reads=re.search("fastq_reads_total\s+([\d,]+)\n", split_log).group(1)
    splited_reads=re.search("fastq_reads_%s\s+([\d,]+)\n"%seq_type, split_log).group(1)
    
    map_log=open('%s/mapping/%s_%s.Log.final.out'%(indir, sample, seq_type)).read()   
    mapped=re.search("Uniquely mapped reads number \|\s+([\d,]+)\n", map_log).group(1)
    
    molecule_log=open('%s/molecule_info/%s_%s_pipeline_stats.txt'%(indir, sample, seq_type)).read()
    total_umis=re.search("total_umis\s+([\d,]+)", molecule_log).group(1)
    usable_reads=re.search("total_read_count\s+([\d,]+)\n", molecule_log).group(1)
    
    summary=open('%s/analysis/%s_%s_analysis_summary.csv'%(indir, sample, seq_type)).read()
    fraction_in_cell=re.search('Fraction Reads in Cells,([\d.]+)\n',summary).group(1)
    cell_filtered_umis=int(float(fraction_in_cell)*int(total_umis))  
    
    if seq_type=='RNA':
        retrim_log=open('%s/split_sample/%s_%s.log'%(indir, sample, seq_type)).read()
        retrimmed=re.search("Reads written.*\s+([\d,]+)\s[(].*\n", retrim_log).group(1).replace(',','')      
        gene_rate=re.search("Reads Mapped to gene,([\d.]+)\n", summary).group(1)
        mapped_to_gene=int(float(gene_rate)*int(cell_filtered_umis))
        return total_reads, splited_reads, retrimmed, mapped, usable_reads, total_umis, cell_filtered_umis, mapped_to_gene
    else:        
        multi_mapped=re.search("Number of reads mapped to multiple loci \|\s+([\d,]+)\n", map_log).group(1)
        mapped=int(mapped)+int(multi_mapped)
        bam_stats=open('%s/molecule_info/%s_%s_bam_stats.txt'%(indir, sample, seq_type)).read()        
        rm_chrM=re.search("number of proper pairs\s+([\d,]+)\n", bam_stats).group(1)
        #rm_chrM=int(int(rm_chrM)/2)    
        return total_reads, splited_reads, mapped, rm_chrM, usable_reads, total_umis, cell_filtered_umis
    
def qc_flow(indir, seq_type, trim_tso):
    stats_prefix=[]
    prefix_list=get_prefix(indir+'/preproc.list')
    sample_list=get_prefix(indir+'/sample.list')
    for prefix in prefix_list:
        stats = qc_prefix(indir, prefix, seq_type, trim_tso)
        stats_prefix.append(map(int, list(stats)))
    stats_prefix=pd.DataFrame(stats_prefix, index=prefix_list)
    if trim_tso:
        stats_prefix.columns=['raw', 'with_tso', 'trimmed', 'bc_parsed', 'splited']       
    else:
        stats_prefix.columns=['raw', 'trimmed', 'bc_parsed', 'splited']
    
    stats_sample=[]
    for sample in sample_list:
        stats = qc_sample(indir, sample, seq_type)
        stats_sample.append(map(int, list(stats)))
    stats_sample=pd.DataFrame(stats_sample, index=sample_list)
    if seq_type=='RNA':
        stats_sample.columns=['total_reads', 'splited_reads', 'retrimmed', 'mapped', 'usable_reads', 'total_umis', 'cell_filtered_umis', 'mapped_to_gene']
    else:
        stats_sample.columns=['total_reads', 'splited_reads', 'mapped', 'rm_chrM', 'usable_reads', 'total_umis', 'cell_filtered_umis']
    frac_prefix=pd.DataFrame([stats_prefix[col]/stats_prefix['raw'] for col in stats_prefix.columns], index=stats_prefix.columns).T
    frac_prefix=frac_prefix.round(3)
    frac_sample_cols=stats_sample.columns[:5]
    frac_sample=pd.DataFrame([stats_sample[col]/stats_sample['total_reads'] for col in frac_sample_cols], index=frac_sample_cols).T
    frac_sample=frac_sample.round(3)
    return stats_prefix, stats_sample, frac_prefix, frac_sample
        
if __name__=='__main__':
    import sys
    import os
    indir=sys.argv[1]
    seq_type=sys.argv[2]
    trim_tso=sys.argv[3]
    if trim_tso=='trim_tso':
        trim_tso=True
    else:
        trim_tso=False
    stats_prefix, stats_sample, frac_prefix, frac_sample=qc_flow(indir, seq_type, trim_tso)
    if not os.path.exists(indir+'/qc_flow'):
        os.makedirs(indir+'/qc_flow')
    stats_prefix.to_csv(indir+'/qc_flow/stats_prefix.txt', sep='\t')
    stats_sample.to_csv(indir+'/qc_flow/stats_sample.txt', sep='\t')
    frac_prefix.to_csv(indir+'/qc_flow/frac_prefix.txt', sep='\t')
    frac_sample.to_csv(indir+'/qc_flow/frac_sample.txt', sep='\t')
    
    
    
