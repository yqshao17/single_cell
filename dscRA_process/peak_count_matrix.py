import pandas as pd
import subprocess
from scipy.io import mmread, mmwrite
from scipy.sparse import csr_matrix
import numpy as np
import os

def summit_extend(summit_file,out_file,extend_size):
    with open(summit_file) as input:
        output=open(out_file,'w')
        for line in input:
            items=line.rstrip('\n').split('\t')
            [chr,start,end,name,val]=items
            name=name.split('/')[-1]
            start = int(int(start) - extend_size/2)
            if start<0:
                start=0
            end = int(int(end) + extend_size/2-1)
            output.write('\t'.join(map(str, [chr,start,end,name,val,'.']))+'\n')


def peak_count_matrix(atac_inter_bed, out_prefix):
    # first three columns is 'chr','start','end'
    # last column is barcode
    inter_peak=pd.read_csv(atac_inter_bed,header=None,sep='\t')
    inter_peak.columns=['chr','start','end']+['']*(inter_peak.shape[1]-4)+['barcode']
    #inter_peak['peak']=inter_peak['name'].apply(lambda x: x.split('/')[-1])
    inter_peak['peak']=inter_peak.apply(lambda x: x['chr']+':'+str(x['start'])+'-'+str(x['end']), axis=1)
    bc_peak_counts=inter_peak.groupby(['barcode','peak']).size()
    df=pd.DataFrame(bc_peak_counts)
    df.reset_index(inplace=True)
    new_df=df.pivot(index='peak', columns='barcode', values=0)
    new_df.columns.name = None
    new_df.index.name = None
    new_df=new_df.fillna(0)
    new_df=new_df.astype('int')
    mmwrite(out_prefix+'count.mtx', csr_matrix(new_df))
    np.savetxt(out_prefix+'peaks.txt',new_df.index.values,fmt="%s")
    np.savetxt(out_prefix+'barcodes.txt',new_df.columns.values,fmt="%s")

    
def intersect(refpeak, bed, outdir):
    intersect=subprocess.check_output(['bedtools', 'intersect', '-a', refpeak, '-b', bed, '-wa','-wb'])
    refpeak_name=refpeak.split('/')[-1].split('.')[0]
    bed_name=bed.split('/')[-1].split('.')[0]
    out_intersect=outdir+'/'+refpeak_name+'_'+bed_name+'_intersect.bed'
    #if os.path.exists(out_intersect):
    #    out_intersect=outdir+'/'+refpeak_name+'_'+bed_name+'_intersect2.bed'
    with open(out_intersect,'wb') as output: # py3
        output.write(intersect)
    return out_intersect

def get_count_matrix(refpeak, bed, outdir, sample):
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    if not os.path.exists(outdir+'/'+sample):
        os.makedirs(outdir+'/'+sample)
    out_intersect=intersect(refpeak, bed, outdir)
    peak_count_matrix(out_intersect, outdir+'/'+sample+'/')


    
if __name__ == '__main__':
    import sys
    refpeak=sys.argv[1]
    bed=sys.argv[2]
    outdir=sys.argv[3]
    sample=sys.argv[4]
    get_count_matrix(refpeak, bed, outdir, sample)

    
    
'''
indir='/Share2/home/zhangqf5/yanqiu/scAR/output/ZJSMix0928/peak_cor'
peak_count_matrix(indir+'/K562_ATAC_summits_K562_ATAC_shift_intersect2.bed', 
                  indir+'/coA_K562/')

peak_count_matrix(indir+'/K562_ATAC_summits_K562_ATAC_shift_intersect.bed', 
                  indir+'/onlyA_K562/')

peak_count_matrix(indir+'/3T3_ATAC_summits_3T3_ATAC_shift_intersect2.bed', 
                  indir+'/coA_3T3/')

peak_count_matrix(indir+'/3T3_ATAC_summits_3T3_ATAC_shift_intersect.bed', 
                  indir+'/onlyA_3T3/')


peak_count_matrix(indir+'/Mix_ATAC_summits_Mix_ATAC_shift_intersect2.bed', 
                  indir+'/coA_Mix/')

peak_count_matrix(indir+'/Mix_ATAC_summits_Mix_ATAC_shift_intersect.bed', 
                  indir+'/onlyA_Mix/')
'''
