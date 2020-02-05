import sys
import re
import os
import datetime
import pickle
import gzip
import numpy as np
import pandas as pd 
from itertools import product
from multiprocessing import Process, Manager

from generage_bc_dicts import free_divergence,write_pkl

#import pdb

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
        fq1_chunks[i] = open(output_dir + '/trimmed/'+prefix+'_R1.chunk%d.fq' %(i+1),'w')
        fq2_chunks[i] = open(output_dir + '/trimmed/'+prefix+'_R2.chunk%d.fq' %(i+1),'w')
        
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
    
def join_fq(output_dir, prefix, read_type, nthreads):
    if 'R1' in read_type:
        files1 = [output_dir +'/parse_bc/'+prefix+ '_barcode_R1.chunk%d.fq' %i for i in range(1,nthreads+1)]
        fout1=open(output_dir + '/parse_bc/'+prefix+ '_barcode_R1.fq', "w")
        for f in files1:
            fin=open(f)
            fout1.write(fin.read())
            fin.close()
        fout1.close()
       
    if 'R2' in read_type:
        files2 = [output_dir +'/parse_bc/'+prefix+ '_barcode_R2.chunk%d.fq' %i for i in range(1,nthreads+1)]
        fout2=open(output_dir + '/parse_bc/'+prefix+ '_barcode_R2.fq', "w")
        for f in files2:
            fin=open(f)
            fout2.write(fin.read())
            fin.close()
        fout2.close()

def preprocess_fastq_chunk(output_dir, prefix, bcset, bc4, counts, config, return_dict=None, chunk=None):
    
    """
    get barcode
    """   
    parse_bc_dir=output_dir+'/parse_bc'
    
    if chunk is None: #not used
        fastq1=output_dir+'/trimmed/'+prefix+'_R1.fq.gz'
        fastq2=output_dir+'/trimmed/'+prefix+'_R2.fq.gz'
        fastq1_out = parse_bc_dir +'/'+prefix+ '_barcode_R1.fq'
        fastq2_out = parse_bc_dir +'/'+prefix+ '_barcode_R2.fq'
    else:
        fastq1=output_dir+'/trimmed/'+prefix+'_R1.chunk%d.fq'%chunk
        fastq2=output_dir+'/trimmed/'+prefix+'_R2.chunk%d.fq'%chunk
        fastq1_out = parse_bc_dir +'/'+prefix+ '_barcode_R1.chunk%d.fq'%chunk
        fastq2_out = parse_bc_dir +'/'+prefix+ '_barcode_R2.chunk%d.fq'%chunk
    
    data_type = config.data_type
    read_type = config.read_type
    
    bc_type = config.bc_type
    bc_starts = config.bc_starts
    bc_len = config.bc_len
    bc_edit_dist=config.bc_edit_dist # list of edit distance for barcodes
    
    umi_type = config.umi_type
    umi_start = config.umi_start
    umi_len = config.umi_len
    
    read_start = config.read_start
    method=config.method
    barcode_set_num = len(bc_edit_dist)
    
    # check_config
    if data_type == 'ATAC_RNA' and (not tag_type):
        sys.exit('Error: No tag for spliting ATAC and RNA. Check data type or tag config.')
    if len(read_start)!=2:
        sys.exit('Error: read starts should be a list with 2 elements. Give any number if you do not want to keep R1 or R2.')
    # check barcode num
    if type(bc_starts)==type(0):
        bc_starts = [bc_starts]
    if type(bc_len)==type(0):
        bc_len = [bc_len]
    if type(bc_edit_dist)==type(0):
        bc_edit_dist = [bc_edit_dist]
            
    if not len(bc_starts)==len(bc_len)==len(bc_edit_dist):
        sys.exit('Error: Not consistent barcode number. Check barcode config.')

    def check_pkl(barcode_n):
        if not os.path.exists(bcset+'/barcode%s.txt'%barcode_n):
            sys.exit('Error: No barcode%s.txt'%barcode_n)
        elif not os.path.exists(bcset+'/bc_dict_%s_%s.pkl'%(barcode_n,method)):
            #write_pkl(method,bcset)
            sys.exit('Error: No barcode%s.pkl'%barcode_n)
    
    for i in range(barcode_set_num):
        check_pkl(str(i+1))


    parse_bc_dir=output_dir+'/parse_bc'
    if not os.path.exists(parse_bc_dir):
        os.makedirs(parse_bc_dir)
    
    bc_edit_dict = [] # list of barcode.pkl for barcode sets
    for i in range(barcode_set_num):
        with open(bcset+'/bc_dict_%d_%s.pkl'%(i+1,method), 'rb') as f:
            edit_dict_v = pickle.load(f)
            bc_edit_dict.append(edit_dict_v)
    
    # Read in barcode sequences
    
    def correct_barcodes(bcsets, counts, bc_edit_dict): 
        bc_matches_list = []
        edit_dist_list = []
        
        for i in range(barcode_set_num):
            bc_matches,edit_dist = get_min_edit_dists(bcsets[i],bc_edit_dict[i],max_d=bc_edit_dist[i])
            #print(bc_matches,edit_dist)
            bc_matches_list.append(bc_matches)
            edit_dist_list.append(edit_dist)
        #print(bc_matches_list,edit_dist_list)
            
        if sum(edit_dist_list)==0: # perfect match
            bc_m = [bc_matches[0] for bc_matches in bc_matches_list]
            return 'perfect', bc_m
        else:
            matches = 0
            if len(bcsets)==1:
                combinations=bc_matches_list[0]
            else:
                combinations=list(product(*bc_matches_list))
                # bc_comb is a tuple of combinations of matched barcode sets
            for bc_comb in combinations: 
                try:
                    cur_counts = counts[bc_comb]
                except:
                    cur_counts = 0
                if cur_counts>0:
                    bc_fixed = bc_comb
                    matches += 1
            if matches==1:
                if len(bcsets)==1:
                    return 'correct', list([bc_fixed])
                else:
                    return 'correct', list(bc_fixed)
            elif matches>1:
                return 'multi_match',['']*len(bcsets)
            else:
                return 'no_match',['']*len(bcsets)

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
    bc_Q30_sum = [0]*barcode_set_num
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
    with gzip.open(fastq1,'rb') as f1, gzip.open(fastq2,'rb') as f2:
        if 'R1' in read_type:
            fout1=open(fastq1_out,'w')
        if 'R2' in read_type:
            fout2=open(fastq2_out,'w')
            
        while True:
            header2 = f2.readline().decode("utf-8")
            if len(header2)==0:
                break
            seq2 = f2.readline().decode("utf-8")
            strand2 = f2.readline().decode("utf-8")
            qual2 = f2.readline().decode("utf-8")      
            header1 = f1.readline().decode("utf-8")
            seq1 = f1.readline().decode("utf-8")
            strand1 = f1.readline().decode("utf-8")
            qual1 = f1.readline().decode("utf-8")
            
            reads ={'R1':{'header':header1,'seq':seq1,'strand':strand1,'qual':qual1},
                    'R2':{'header':header2,'seq':seq2,'strand':strand2,'qual':qual2}}
            
            bcsets = []
            for i in range(barcode_set_num):
                bc=reads[bc_type]['seq'][bc_starts[i]:bc_starts[i]+bc_len[i]]
                bcsets.append(bc)
            
            umi = reads[umi_type]['seq'][umi_start:umi_start+umi_len]
            
            map_tag, bcsets = correct_barcodes(bcsets, counts, bc_edit_dict)
            
            if map_tag=='perfect': 
                perfect_bc_n+=1
                fastq_valid_barcode_reads+=1
            elif map_tag=='correct': 
                correct_bc_n+=1
                fastq_valid_barcode_reads+=1
            elif map_tag=='multi_match': multimatch_bc_n+=1
            elif map_tag=='no_match': nomatch_bc_n+=1
            
            cellbc_umi = ''.join(bcsets) +':' + umi
            
            if len(cellbc_umi)==sum(bc_len+[umi_len])+1:
                if 'R1' in read_type:
                    header1 = '@' +  bc4 + cellbc_umi + ':'+header1[1:].split()[0] +'\n'
                    fout1.write(header1+seq1[read_start[0]:]+strand1+qual1[read_start[0]:])
                if 'R2' in read_type:
                    header2 = '@' +  bc4 + cellbc_umi + ':'+header2[1:].split()[0] +'\n'
                    fout2.write(header2+seq2[read_start[1]:]+strand2+qual2[read_start[1]:])
                
            fastq_reads += 1
            # bc1 refers to the first barcode seen in the sequencing read, but is actually
            # bc3 in terms of rounds
            for i in range(barcode_set_num):
                bc=reads[bc_type]['seq'][bc_starts[i]:bc_starts[i]+bc_len[i]]
                bc_Q30_sum[i]+=np.mean([ord(c)>62 for c in reads[bc_type]['qual'][bc_starts[i]:bc_starts[i]+bc_len[i]]])
            
            if umi_len:
                umi_Q30_sum += np.mean([ord(c)>62 for c in reads[umi_type]['qual'][umi_start:umi_start+umi_len]])
            # Make sure R2 length > (bc_starts[3]+umi_bc_len[3])
            cDNA_Q30_sum = []
            if 'R1' in read_type:
                R1_cDNA_Q30_sum += np.mean([ord(c)>62 for c in qual1[read_start[0]:]])
                cDNA_Q30_sum.append(R1_cDNA_Q30_sum)
            if 'R2' in read_type:
                R2_cDNA_Q30_sum += np.mean([ord(c)>62 for c in qual2[read_start[1]:]])
                cDNA_Q30_sum.append(R2_cDNA_Q30_sum)

    if 'R1' in read_type:
        fout1.close()
    if 'R2' in read_type:
        fout2.close()
    if chunk:
        return_dict[chunk]=np.array([fastq_reads, fastq_valid_barcode_reads, *bc_Q30_sum, umi_Q30_sum, *cDNA_Q30_sum, perfect_bc_n, correct_bc_n, multimatch_bc_n, nomatch_bc_n])
    else:
        return np.array([fastq_reads, fastq_valid_barcode_reads, *bc_Q30_sum, umi_Q30_sum, *cDNA_Q30_sum, perfect_bc_n, correct_bc_n, multimatch_bc_n, nomatch_bc_n])

def preprocess_fastq(output_dir, prefix, bcset, bc4, config, nthreads):
    
    parse_bc_dir=output_dir+'/parse_bc'
    
    bc_len = config.bc_len
    bc_starts = config.bc_starts
    bc_type = config.bc_type
    read_type = config.read_type
    barcode_set_num = len(bc_len)
    nthreads = int(nthreads)
    # Read in barcode sequences
    bcsets_list = [] # list of barcode sets
    for i in range(barcode_set_num):
        bc = pd.read_csv('%s/barcode%d.txt'%(bcset,i+1),names=['barcode'],index_col=False, sep='\t').barcode.values
        # set dictionary for faster search
        bc_dict = {}
        for el in bc:
            bc_dict[el]=''
        bcsets_list.append(bc_dict)

    
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
                #if n %100000==0:
                #    print(datetime.datetime.now(),n)
                n+=1
                # pass n_bin reads
                for i in range(n_bin):
                    for j in range(4):
                        f.readline()
        seqs = pd.Series(seqs)
        bc_df = pd.DataFrame()
        query = []
        group = []
        perfect_rate = []
        for i in range(barcode_set_num):
            #print (i, seqs.head(), bc_starts, bc_len)
            bc_df['bc%d'%(i+1)] = seqs.str.slice(bc_starts[i],bc_starts[i]+bc_len[i]) # slice barcode seq
            bc_df['bc%d_valid'%(i+1)] = bc_df['bc%d'%(i+1)].apply(lambda s: s in bcsets_list[i]) # check if bc in bcset
            query.append('bc%d_valid'%(i+1)) 
            group.append('bc%d'%(i+1))
            perfect_rate.append(bc_df.query('bc%d_valid'%(i+1)).shape[0]*1.0/n) # bc perfect match rate
        query='&'.join(query)
        counts = bc_df.query(query).groupby(group).size().sort_values(ascending=False)       
        perfect_rate.append(bc_df.query(query).shape[0]*1.0/n)
        #print(perfect_rate)
        return counts.to_dict(), perfect_rate
    
    print(datetime.datetime.now(), 'Getting perfect barcodes from %s.' % bc_type)
    fastq=output_dir+'/trimmed/'+prefix+'_%s.fq.gz'%bc_type
    counts, perfect_rate= get_perfect_bc_counts(fastq)
    print(datetime.datetime.now(), 'Reading and writing fastq and correcting barcodes.')
    
    stats_array = preprocess_fastq_chunk(output_dir, prefix, bcset, bc4, counts, config)
    '''
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
    
    #pdb.set_trace()
    
    stats_array=np.zeros(len(return_dict[1]))
    
    for i in range(1,nthreads+1):
        stats_array+=return_dict[i]
    '''
    #[fastq_reads, fastq_valid_barcode_reads, bc_Q30_sum, umi_Q30_sum, cDNA_Q30_sum, perfect_bc_n, correct_bc_n, multimatch_bc_n, nomatch_bc_n]=list(stats_array)
    fastq_reads, fastq_valid_barcode_reads = stats_array[:2]
    bc_Q30_sum = stats_array[2:2+barcode_set_num]
    umi_Q30_sum = stats_array[2+barcode_set_num]
    cDNA_Q30_sum = stats_array[3+barcode_set_num:-4]
    perfect_bc_n, correct_bc_n, multimatch_bc_n, nomatch_bc_n = stats_array[-4:]
    
    with open(parse_bc_dir +'/'+prefix+ '_sequencing_stats.txt', 'w') as f:
        f.write('umi_Q30\t%0.4f\n' %(umi_Q30_sum/fastq_reads))
        for i in range(barcode_set_num):
            f.write('bc%d_Q30\t%0.4f\n' %(i+1, bc_Q30_sum[i]/fastq_reads))
        
        f.write('cDNA_Q30')
        for el in cDNA_Q30_sum:
            f.write('\t%0.4f'%(el/fastq_reads))
        f.write('\n')

    with open(parse_bc_dir +'/'+prefix+ '_pipeline_stats.txt', 'w') as f:
        f.write('fastq_reads\t%d\n' %fastq_reads)
        f.write('fastq_valid_barcode_reads\t%d\n' %fastq_valid_barcode_reads)
        f.write('fastq_valid_rate\t%0.4f\n' %(fastq_valid_barcode_reads*1.0/fastq_reads))
        f.write('--perfect_match\t%0.4f\n' %(perfect_bc_n*1.0/fastq_reads))
        f.write('--correct_match\t%0.4f\n' %(correct_bc_n*1.0/fastq_reads))
        f.write('multi_match\t%0.4f\n' %(multimatch_bc_n*1.0/fastq_reads))
        f.write('no_match\t%0.4f\n' %(nomatch_bc_n*1.0/fastq_reads))
        f.write('total_barcode_num\t%d\n' %len(counts.keys()))
        for i in range(barcode_set_num):
            f.write('perfect_barcode%d\t%0.4f\n'%(i+1,perfect_rate[i]))
        f.write('perfect_barcodes_all\t%0.4f\n'%perfect_rate[-1])
    
    '''
    join_fq(output_dir, prefix, read_type, nthreads)
    for i in range(1,int(nthreads)+1):
        #os.remove(output_dir +'/trimmed/'+prefix+ '_R1.chunk%d.fq' %i)
        #os.remove(output_dir +'/trimmed/'+prefix+ '_R2.chunk%d.fq' %i)
        if 'R1' in read_type:
            os.remove(output_dir +'/parse_bc/'+prefix+ '_barcode_R1.chunk%d.fq' %i)
        if 'R2' in read_type:
            os.remove(output_dir +'/parse_bc/'+prefix+ '_barcode_R2.chunk%d.fq' %i)
    '''
        
if __name__=='__main__':
    import sys
    output_dir='/Share2/home/zhangqf5/yanqiu/dscRA/output/test/RNA'
    prefix='test'
    bcset='barcodes'
    bc4='testxxx.'
    nthreads=5
    
    preprocess_fastq(output_dir, prefix, bcset, bc4, nthreads)
    #join_fq(output_dir, prefix, nthreads)
    for i in range(1,int(nthreads)+1):
        os.remove(output_dir +'/trimmed/'+prefix+ '_R1.chunk%d.fq' %i)
        os.remove(output_dir +'/trimmed/'+prefix+ '_R2.chunk%d.fq' %i)
        os.remove(output_dir +'/parse_bc/'+prefix+ '_barcode_R1.chunk%d.fq' %i)
        os.remove(output_dir +'/parse_bc/'+prefix+ '_barcode_R2.chunk%d.fq' %i)
    
