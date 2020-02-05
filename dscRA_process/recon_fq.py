# reconstruct RNA and ATAC

import gzip
import re


def reconstruct_RNA(indir, sample, outdir, lane):
    fq1 = gzip.open(indir+'/'+sample+'_'+lane+'_R1_001.fastq.gz','rb')
    fq2 = gzip.open(indir+'/'+sample+'_'+lane+'_R3_001.fastq.gz','rb')
    bcf = gzip.open(indir+'/'+sample+'_'+lane+'_R2_001.fastq.gz','rb')
    fout1 = open(outdir+'/'+sample+'_recon_'+lane+'_R1_001.fastq','w')
    fout2 = open(outdir+'/'+sample+'_recon_'+lane+'_R2_001.fastq','w')
    total = 0
    count = 0
    while True:
        seqname1 = fq1.readline().decode("utf-8")
        if not seqname1:
            break
        seq1 = fq1.readline().decode("utf-8")
        strand1 = fq1.readline().decode("utf-8")
        qual1 = fq1.readline().decode("utf-8")
        seqname2 = fq2.readline().decode("utf-8")
        seq2 = fq2.readline().decode("utf-8")
        strand2 = fq2.readline().decode("utf-8")
        qual2 = fq2.readline().decode("utf-8")
        seqname3 = bcf.readline().decode("utf-8")
        seq3 = bcf.readline().decode("utf-8")
        strand3 = bcf.readline().decode("utf-8")
        qual3 = bcf.readline().decode("utf-8")
        
        total += 1 
        assert seqname1.strip().split()[0]==seqname2.strip().split()[0]==seqname3.strip().split()[0], "Unpaired reads"
        if re.search('TTTTTTT',seq1):
            count += 1
            seq11 = seq3.strip()+seq1[:12]+'\n'
            qual11 = qual3.strip()+qual1[:12]+'\n'
            fout1.write(seqname1+seq11+strand1+qual11)
            fout2.write(seqname2+seq2+strand2+qual2)
    print('Total reads: %d'%total)
    print('Filtered reads with TTT: %d'%count)     
    

def reconstruct_ATAC(indir, sample, outdir, lane):
    fq1 = gzip.open(indir+'/'+sample+'_'+lane+'_R1_001.fastq.gz','rb')
    fq2 = gzip.open(indir+'/'+sample+'_'+lane+'_R3_001.fastq.gz','rb')
    bcf = gzip.open(indir+'/'+sample+'_'+lane+'_R2_001.fastq.gz','rb')
    fout1 = open(outdir+'/'+sample+'_recon_'+lane+ '_R1_001.fastq','w')
    fout2 = open(outdir+'/'+sample+'_recon_'+lane+'_R2_001.fastq','w')
    total = 0
    while True:
        seqname1 = fq1.readline().decode("utf-8")
        if not seqname1:
            break
        seq1 = fq1.readline().decode("utf-8")
        strand1 = fq1.readline().decode("utf-8")
        qual1 = fq1.readline().decode("utf-8")
        seqname2 = fq2.readline().decode("utf-8")
        seq2 = fq2.readline().decode("utf-8")
        strand2 = fq2.readline().decode("utf-8")
        qual2 = fq2.readline().decode("utf-8")
        seqname3 = bcf.readline().decode("utf-8")
        seq3 = bcf.readline().decode("utf-8")
        strand3 = bcf.readline().decode("utf-8")
        qual3 = bcf.readline().decode("utf-8")
        
        total += 1 
        assert seqname1.strip().split()[0]==seqname2.strip().split()[0]==seqname3.strip().split()[0], "Unpaired reads"
        seq11 = seq3.strip()+seq1
        qual11 = qual3.strip()+qual1
        fout1.write(seqname1+seq11+strand1+qual11)
        fout2.write(seqname2+seq2+strand2+qual2)
    fout1.close()
    fout2.close()
    print('Total reads: %d'%total)     


if __name__=="__main__":
    import argparse
    from multiprocessing import Process
    import os
    parser = argparse.ArgumentParser()
    parser.add_argument('--datatype','-d', help='ATAC or RNA')
    parser.add_argument('--indir','-i', help='input path')
    parser.add_argument('--sample','-s', help='sample name')
    parser.add_argument('--outdir','-o', help='output path')
    parser.add_argument('--lanes','-l', help='L001,L002')
    args = parser.parse_args()
    
    indir=args.indir
    outdir=args.outdir
    sample=args.sample
    
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    lanes = args.lanes.split(',')
    Pros = []
    for lane in lanes:
        if args.datatype=='RNA':
            p = Process(target=reconstruct_RNA, args=(indir, sample, outdir, lane))
        elif args.datatype=='ATAC':
            p = Process(target=reconstruct_ATAC, args=(indir, sample, outdir, lane))
        Pros.append(p)
        p.start()
    for t in Pros:
        t.join()
    
    if len(lanes)>1:
        os.system("""cat {0}/{1}_recon_*_R1_001.fastq > {0}/{1}_recon_R1.fastq""".format(outdir, sample))
        os.system("""cat {0}/{1}_recon_*_R2_001.fastq > {0}/{1}_recon_R2.fastq""".format(outdir, sample))
        