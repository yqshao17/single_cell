import sys
import re
import os
import datetime
import pickle
import gzip
import numpy as np
import pandas as pd 
from multiprocessing import Process, Manager

from generage_bc_dicts import free_divergence,write_pkl


def get_min_edit_dists(bc,edit_dict,max_d=3):
    """Returns a list of nearest edit dist seqs
    Input 8nt barcode, edit_dist_dictionary
    Output <list of nearest edit distance seqs>, <edit dist>"""
    bc_matches = edit_dict[0][bc]
    edit_dist = 0
    if (len(bc_matches)==0) and (max_d>=1):
        edit_dist+=1
        bc_matches = edit_dict[1][bc]
    if (len(bc_matches)==0) and (max_d>=2):
        edit_dist+=1
        bc_matches = edit_dict[2][bc]
    if (len(bc_matches)==0) and (max_d>=3):
        edit_dist+=1
        bc_matches = edit_dict[3][bc]
    return bc_matches,edit_dist


def split_fq(output_dir, prefix, nthreads):
    fastq1=output_dir+'/trimmed/'+prefix+'_R1.fq.gz'
    fastq2=output_dir+'/trimmed/'+prefix+'_R2.fq.gz'
    fq1 = gzip.open(fastq1,'rb')
    fq2 = gzip.open(fastq2,'rb')
    fq1_chunks = {}
    fq2_chunks = {}
    for i in range(nthreads):
        fq1_chunks[i] = open(output_dir + '/trimmed/'+prefix+'_R1.chunk%d.fastq' %(i+1),'w')
        fq2_chunks[i] = open(output_dir + '/trimmed/'+prefix+'_R2.chunk%d.fastq' %(i+1),'w')
        
    # Get the total number of aligned reads:
    with open(output_dir +'/trimmed/'+prefix+'.log') as f:
        log=f.read()
        try:
            num_reads = re.search("Pairs written.*\s+([\d,]+)\s[(].*\n", log).group(1)
        except:
            sys.exit("Invalid trim log file!")
        num_reads = int(num_reads.replace(',',''))

    # Number of reads per file. Round up to ensure we don't write (nthreads +1) files.
    reads_per_chunk = int(np.ceil(num_reads/nthreads))

    c = 0
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
        
        chunk = int(np.floor(c/reads_per_chunk))
        fq1_chunks[chunk].write(seqname1+seq1+strand1+qual1)
        fq2_chunks[chunk].write(seqname2+seq2+strand2+qual2)
        c+=1
            
    for i in range(nthreads):
        fq1_chunks[i].close()
        fq2_chunks[i].close()
    fq1.close()
    fq2.close()
    
def join_fq(output_dir, prefix, nthreads):
    files1 = [output_dir +'/parse_bc/'+prefix+ '_tag_barcode_R1.chunk%d.fastq' %i for i in range(1,nthreads+1)]
    files2 = [output_dir +'/parse_bc/'+prefix+ '_tag_barcode_R2.chunk%d.fastq' %i for i in range(1,nthreads+1)]
    fout1=open(output_dir + '/parse_bc/'+prefix+ '_tag_barcode_R1.fastq', "w")
    fout2=open(output_dir + '/parse_bc/'+prefix+ '_tag_barcode_R2.fastq', "w")
    for f in files1:
        fin=open(f)
        fout1.write(fin.read())
        fin.close()
    fout1.close()
    for f in files2:
        fin=open(f)
        fout2.write(fin.read())
        fin.close()
    fout2.close()
    

def preprocess_fastq_chunk(output_dir, prefix, bcset, bc4, counts, config, return_dict, chunk=None):
    
    """
    get barcode
    """   
    parse_bc_dir=output_dir+'/parse_bc'
    
    if chunk is None: #not used
        fastq1=output_dir+'/trimmed/'+prefix+'_R1.fq.gz'
        fastq2=output_dir+'/trimmed/'+prefix+'_R2.fq.gz'
        fastq1_out = parse_bc_dir +'/'+prefix+ '_tag_barcode_R1.fastq'
        fastq2_out = parse_bc_dir +'/'+prefix+ '_tag_barcode_R2.fastq'
    else:
        fastq1=output_dir+'/trimmed/'+prefix+'_R1.chunk%d.fastq'%chunk
        fastq2=output_dir+'/trimmed/'+prefix+'_R2.chunk%d.fastq'%chunk
        fastq1_out = parse_bc_dir +'/'+prefix+ '_tag_barcode_R1.chunk%d.fastq'%chunk
        fastq2_out = parse_bc_dir +'/'+prefix+ '_tag_barcode_R2.chunk%d.fastq'%chunk
    
    umi_bc_len = config.umi_bc_len
    umi_bc_starts = config.umi_bc_starts
    tag_seq=config.tag_seq
    tag_start=config.tag_start
    tag_length=config.tag_length
    method=config.method
    bc_edit_dist=config.bc_edit_dist
    tag_edit_dist=config.tag_edit_dist

    bc_edit_dist1 = int(bc_edit_dist[0])
    bc_edit_dist2 = int(bc_edit_dist[1])
    bc_edit_dist3 = int(bc_edit_dist[2])

    tag_edit_dist_atac = int(tag_edit_dist[0])
    tag_edit_dist_rna = int(tag_edit_dist[1])


    def check_pkl(barcode_n):
        if not os.path.exists(bcset+'/barcode%s.txt'%barcode_n):
            sys.exit('Error: No barcode%s.txt'%barcode_n)
        elif not os.path.exists(bcset+'/bc_dict_%s_%s.pkl'%(barcode_n,method)):
            write_pkl(method,bcset)
    check_pkl('1')
    check_pkl('2')
    check_pkl('3')

    with open(bcset+'/bc_dict_1_%s.pkl'%method, 'rb') as f:
        edit_dict_v1 = pickle.load(f)
    with open(bcset+'/bc_dict_2_%s.pkl'%method, 'rb') as f:
        edit_dict_v2 = pickle.load(f)
    with open(bcset+'/bc_dict_3_%s.pkl'%method, 'rb') as f:
        edit_dict_v3 = pickle.load(f)
    
    parse_bc_dir=output_dir+'/parse_bc'
    if not os.path.exists(parse_bc_dir):
        os.makedirs(parse_bc_dir)
    
    
    # Read in barcode sequences
    bc1 = pd.read_csv(bcset+'/barcode1.txt',names=['barcode'],index_col=0, sep='\t').barcode.values
    bc2 = pd.read_csv(bcset+'/barcode2.txt',names=['barcode'],index_col=0, sep='\t').barcode.values
    bc3 = pd.read_csv(bcset+'/barcode3.txt',names=['barcode'],index_col=0, sep='\t').barcode.values

    bc1_edit_dict = edit_dict_v1
    bc2_edit_dict = edit_dict_v2
    bc3_edit_dict = edit_dict_v3
    

    
    def correct_barcodes(bc1,bc2,bc3, counts, bc1_dict=bc1_edit_dict,bc2_dict=bc2_edit_dict,bc3_dict=bc3_edit_dict):
        bc1_matches,edit_dist1 = get_min_edit_dists(bc1,bc1_dict,max_d=bc_edit_dist1)
        bc2_matches,edit_dist2  = get_min_edit_dists(bc2,bc2_dict,max_d=bc_edit_dist2)
        bc3_matches,edit_dist3  = get_min_edit_dists(bc3,bc3_dict,max_d=bc_edit_dist3)

        if 0==edit_dist1==edit_dist2==edit_dist3:
            bc1_m=bc1_matches[0]
            bc2_m=bc2_matches[0]
            bc3_m=bc3_matches[0]
            return 'perfect',bc1_m,bc2_m,bc3_m
        else:
            matches = 0
            for bc1_m in bc1_matches:
                for bc2_m in bc2_matches:
                    for bc3_m in bc3_matches:
                        try:
                            cur_counts = counts[(bc1_m,bc2_m,bc3_m)]
                        except:
                            cur_counts = 0
                        if cur_counts>0:
                            bc1_fixed = bc1_m
                            bc2_fixed = bc2_m
                            bc3_fixed = bc3_m
                            matches += 1
            #print('Fixed:',matches,\
            #      'Nearest_bcs',(len(bc1_matches),len(bc2_matches),len(bc3_matches)),\
            #      'Edit_dist',(edit_dist1,edit_dist2,edit_dist3))
            if matches==1:               
                return 'correct',bc1_fixed,bc2_fixed,bc3_fixed
            elif matches>1:
                return 'multi_match','','',''
            else:
                return 'no_match','','',''

    def check_atac_rna(tag_readcut):  
        distance = free_divergence(tag_seq, tag_readcut)
        if distance <= tag_edit_dist_atac:
            return 'A' # atac reads
        elif distance >= tag_edit_dist_rna:
            return 'R' # rna read
        else: 
            return 'W' # waste

    fastq_reads = 0
    fastq_valid_barcode_reads = 0
    bc1_Q30_sum = 0
    bc2_Q30_sum = 0
    bc3_Q30_sum = 0
    umi_Q30_sum = 0
    R1_cDNA_Q30_sum = 0
    R2_cDNA_Q30_sum = 0
    perfect_bc_n = 0
    correct_bc_n = 0
    multimatch_bc_n = 0
    nomatch_bc_n = 0
    atac_reads = 0
    rna_reads = 0
    waste_reads = 0
    with open(fastq1,'r') as f1, open(fastq2,'r') as f2, \
        open(fastq1_out,'w') as fout1,\
        open(fastq2_out,'w') as fout2:
        
        while True:
            header2 = f2.readline()
            #bc4 = header.strip().split(':')[-1]
            if len(header2)==0:
                break
            seq2 = f2.readline()
            bc1 = seq2[umi_bc_starts[1]:umi_bc_starts[1]+umi_bc_len[1]]
            bc2 = seq2[umi_bc_starts[2]:umi_bc_starts[2]+umi_bc_len[2]]
            bc3 = seq2[umi_bc_starts[3]:umi_bc_starts[3]+umi_bc_len[3]]
            umi = seq2[umi_bc_starts[0]:umi_bc_starts[0]+umi_bc_len[0]]
            
            strand2 = f2.readline()
            qual2 = f2.readline()
            
            map_tag, bc1,bc2,bc3 = correct_barcodes(bc1,bc2,bc3, counts)
            if map_tag=='perfect': perfect_bc_n+=1
            elif map_tag=='correct': correct_bc_n+=1
            elif map_tag=='multi_match': multimatch_bc_n+=1
            elif map_tag=='no_match': nomatch_bc_n+=1
            
            cellbc_umi = bc1 + bc2 + bc3 +':' + umi

            header1 = f1.readline()
            seq1 = f1.readline()
            strand1 = f1.readline()
            qual1 = f1.readline()
            
            if len(cellbc_umi)==sum(umi_bc_len)+1:
                tag_readcut = seq2[tag_start:(tag_start+tag_length)]
                tag = check_atac_rna(tag_readcut)
                
                header1 = '@' + tag + '_' + bc4 + cellbc_umi + ':'+header1[1:].split()[0] +'\n'
                header2 = '@' + tag + '_' + bc4 + cellbc_umi + ':'+header2[1:].split()[0] +'\n'
                
                fout1.write(header1)
                fout1.write(seq1)
                fout1.write(strand1)
                fout1.write(qual1)
                fout2.write(header2)
                
                fastq_valid_barcode_reads += 1

                if tag=='A':
                    atac_reads += 1
                    fout2.write(seq2[(tag_start+tag_length):])
                    fout2.write(strand2)
                    fout2.write(qual2[(tag_start+tag_length):])
                elif tag=='R':
                    rna_reads += 1
                    fout2.write(seq2[tag_start:])
                    fout2.write(strand2)
                    fout2.write(qual2[tag_start:])
                elif tag=='W':
                    waste_reads +=1
                    fout2.write(seq2[(tag_start+tag_length):])
                    fout2.write(strand2)
                    fout2.write(qual2[(tag_start+tag_length):])
            fastq_reads += 1
            # bc1 refers to the first barcode seen in the sequencing read, but is actually
            # bc3 in terms of rounds
            bc1_Q30_sum += np.mean([ord(c)>62 for c in qual2[umi_bc_starts[1]:umi_bc_starts[1]+umi_bc_len[1]]])
            bc2_Q30_sum += np.mean([ord(c)>62 for c in qual2[umi_bc_starts[2]:umi_bc_starts[2]+umi_bc_len[2]]])
            bc3_Q30_sum += np.mean([ord(c)>62 for c in qual2[umi_bc_starts[3]:umi_bc_starts[3]+umi_bc_len[3]]])
            umi_Q30_sum += np.mean([ord(c)>62 for c in qual2[umi_bc_starts[0]:umi_bc_starts[0]+umi_bc_len[0]]])
            # Make sure R2 length > (bc_starts[3]+umi_bc_len[3])
            R2_cDNA_Q30_sum += np.mean([ord(c)>62 for c in qual2[umi_bc_starts[3]+umi_bc_len[3]:-1]])
            R1_cDNA_Q30_sum += np.mean([ord(c)>62 for c in qual1[:-1]])
            
    return_dict[chunk]=np.array([fastq_reads, fastq_valid_barcode_reads, bc1_Q30_sum, bc2_Q30_sum, bc3_Q30_sum, umi_Q30_sum, R1_cDNA_Q30_sum, R2_cDNA_Q30_sum, perfect_bc_n, correct_bc_n, multimatch_bc_n, nomatch_bc_n, atac_reads, rna_reads,waste_reads])
     

def preprocess_fastq(output_dir, prefix, bcset, bc4, config, nthreads):
    
    parse_bc_dir=output_dir+'/parse_bc'
    
    umi_bc_len = config.umi_bc_len
    umi_bc_starts = config.umi_bc_starts
    nthreads = int(nthreads)
    # Read in barcode sequences
    bc1 = pd.read_csv(bcset+'/barcode1.txt',names=['barcode'],index_col=0, sep='\t').barcode.values
    bc2 = pd.read_csv(bcset+'/barcode2.txt',names=['barcode'],index_col=0, sep='\t').barcode.values
    bc3 = pd.read_csv(bcset+'/barcode3.txt',names=['barcode'],index_col=0, sep='\t').barcode.values
    
    def get_perfect_bc_counts(fastq2, n_reads=1000000):
        #randomly choose n_reads
        lines=0
        f=gzip.open(fastq2,'rb')
        for line in f:
            lines+=1
        total_reads = lines/4
        n_bin = round(int(total_reads)/n_reads-0.5)

        quality_scores = []
        seqs = []
        n = 0
        with gzip.open(fastq2,'rb') as f:
            while n<n_reads:
                seqname=f.readline().decode("utf-8")
                if len(seqname)==0:
                    break
                seq = f.readline().decode("utf-8")
                f.readline().decode("utf-8")
                qual = f.readline().decode("utf-8")
                seqs.append(seq)
                quality_scores.append(qual)
                if n %100000==0:
                    print(datetime.datetime.now(),n)
                n+=1
                # pass n_bin reads
                for i in range(n_bin):
                    for j in range(4):
                        f.readline()
        seqs = pd.Series(seqs)
        bc_df = pd.DataFrame()
        bc_df['bc1'] = seqs.str.slice(umi_bc_starts[1],umi_bc_starts[1]+umi_bc_len[1])
        bc_df['bc2'] = seqs.str.slice(umi_bc_starts[2],umi_bc_starts[2]+umi_bc_len[2])
        bc_df['bc3'] = seqs.str.slice(umi_bc_starts[3],umi_bc_starts[3]+umi_bc_len[3])
        bc_df['bc1_valid'] = bc_df['bc1'].apply(lambda s: s in bc1)
        bc_df['bc2_valid'] = bc_df['bc2'].apply(lambda s: s in bc2)
        bc_df['bc3_valid'] = bc_df['bc3'].apply(lambda s: s in bc3)
        counts = bc_df.query('bc1_valid & bc2_valid & bc3_valid').groupby(['bc1','bc2','bc3']).size().sort_values(ascending=False)       
        perfect_rate= [bc_df.query('bc1_valid & bc2_valid & bc3_valid').shape[0]*1.0/n,
          bc_df.query('bc1_valid').shape[0]*1.0/n,
          bc_df.query('bc2_valid').shape[0]*1.0/n,
          bc_df.query('bc3_valid').shape[0]*1.0/n]
        return counts.to_dict(), perfect_rate
    
    print(datetime.datetime.now(), 'Getting perfect barcodes from fastq2.')
    fastq2=output_dir+'/trimmed/'+prefix+'_R2.fq.gz'
    counts, perfect_rate= get_perfect_bc_counts(fastq2)
    
    print(datetime.datetime.now(), 'Reading and writing fastq and correcting barcodes.')
    print('Method: ', config.method,', Barcode distance: ', config.bc_edit_dist)
    print('Tag distance for ATAC: <=', config.tag_edit_dist[0], ', Tag distance for RNA: >=', config.tag_edit_dist[1])
    split_fq(output_dir, prefix, nthreads)
     
    manager = Manager()
    return_dict = manager.dict()
    Pros = []
    for i in range(1,nthreads+1):
        print('Starting thread %d' %i)
        p = Process(target=preprocess_fastq_chunk, args=(output_dir, prefix, bcset, bc4, counts, config, return_dict, i))
        Pros.append(p)
        p.start()
    for t in Pros:
        t.join()

    join_fq(output_dir, prefix, nthreads)

    stats_array=np.zeros(15)
    for i in range(1,nthreads+1):
        stats_array+=return_dict[i]
    
    #print(stats_array)
    
    [fastq_reads, fastq_valid_barcode_reads, bc1_Q30_sum, bc2_Q30_sum, bc3_Q30_sum, umi_Q30_sum, R1_cDNA_Q30_sum, R2_cDNA_Q30_sum, perfect_bc_n, correct_bc_n, multimatch_bc_n, nomatch_bc_n, atac_reads, rna_reads, waste_reads]=list(stats_array)
    
    with open(parse_bc_dir +'/'+prefix+ '_sequencing_stats.txt', 'w') as f:
        f.write('umi_Q30\t%0.4f\n' %(umi_Q30_sum/fastq_reads))
        f.write('bc1_Q30\t%0.4f\n' %(bc1_Q30_sum/fastq_reads))
        f.write('bc2_Q30\t%0.4f\n' %(bc2_Q30_sum/fastq_reads))
        f.write('bc3_Q30\t%0.4f\n' %(bc3_Q30_sum/fastq_reads)) 
        f.write('R2_tag+cDNA_Q30\t%0.4f\n' %(R2_cDNA_Q30_sum/fastq_reads))
        f.write('R1_cDNA_Q30\t%0.4f\n' %(R1_cDNA_Q30_sum/fastq_reads))
        
    with open(parse_bc_dir +'/'+prefix+ '_pipeline_stats.txt', 'w') as f:
        f.write('fastq_reads\t%d\n' %fastq_reads)
        f.write('fastq_valid_barcode_reads\t%d\n' %fastq_valid_barcode_reads)
        f.write('--ATAC_valid_barcode_reads\t%d\n' %atac_reads)
        f.write('--RNA_valid_barcode_reads\t%d\n' %rna_reads)
        f.write('--Undetermined_valid_barcode_reads\t%d\n' %waste_reads)
        f.write('fastq_valid_rate\t%0.4f\n' %(fastq_valid_barcode_reads*1.0/fastq_reads))     
        f.write('--ATAC_valid\t%0.4f\n' %(atac_reads*1.0/fastq_reads))      
        f.write('--RNA_valid\t%0.4f\n' %(rna_reads*1.0/fastq_reads))
        f.write('--Undetermined_valid\t%0.4f\n' %(waste_reads*1.0/fastq_reads))
        f.write('--perfect_match\t%0.4f\n' %(perfect_bc_n*1.0/fastq_reads))
        f.write('--correct_match\t%0.4f\n' %(correct_bc_n*1.0/fastq_reads))
        f.write('multi_match\t%0.4f\n' %(multimatch_bc_n*1.0/fastq_reads))
        f.write('no_match\t%0.4f\n' %(nomatch_bc_n*1.0/fastq_reads))
        f.write('total_barcode_num\t%d\n' %len(counts.keys()))
        f.write('perfect_barcodes_all\t%0.4f\n'%perfect_rate[0])
        f.write('perfect_barcode1\t%0.4f\n'%perfect_rate[1])
        f.write('perfect_barcode2\t%0.4f\n'%perfect_rate[2])
        f.write('perfect_barcode3\t%0.4f\n'%perfect_rate[3])
        
    for i in range(1,int(nthreads)+1):
        os.remove(output_dir +'/trimmed/'+prefix+ '_R1.chunk%d.fastq' %i)
        os.remove(output_dir +'/trimmed/'+prefix+ '_R2.chunk%d.fastq' %i)
        os.remove(output_dir +'/parse_bc/'+prefix+ '_tag_barcode_R1.chunk%d.fastq' %i)
        os.remove(output_dir +'/parse_bc/'+prefix+ '_tag_barcode_R2.chunk%d.fastq' %i)

        
if __name__=='__main__':
    import sys
    output_dir='/Share2/home/zhangqf5/yanqiu/scAR/output/test_set'
    prefix='test'
    bcset='barcodes/bcset0g'
    bc4='testxxx.'
    nthreads=10
    
    #preprocess_fastq(output_dir, prefix, bcset, bc4, nthreads)
    join_fq(output_dir, prefix, nthreads)
    for i in range(1,int(nthreads)+1):
        os.remove(output_dir +'/trimmed/'+prefix+ '_R1.chunk%d.fastq' %i)
        os.remove(output_dir +'/trimmed/'+prefix+ '_R2.chunk%d.fastq' %i)
        os.remove(output_dir +'/parse_bc/'+prefix+ '_tag_barcode_R1.chunk%d.fastq' %i)
        os.remove(output_dir +'/parse_bc/'+prefix+ '_tag_barcode_R2.chunk%d.fastq' %i)
    
