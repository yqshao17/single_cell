# make barcode dictionaries
import os
import subprocess
import sys

import pandas as pd
from collections import defaultdict
import gzip
from numpy import unique
import numpy as np
import pickle



def convert_degen_seq_to_list(seq):
    """Uses recursion to convert a degenerate sequence to a list
    For example: AGGN -> [AGGA, AGGC, AGGG, AGGT]"""
    bases=['A','G','C','T']
    seq_list = []
    N_pos = seq.find('N')
    if N_pos>=0:
        for b in bases:
            seq_list += convert_degen_seq_to_list(seq[:N_pos] + b + seq[N_pos+1:])
    else:
        seq_list.append(seq)
    return seq_list

def levenshteinDistance(s1, s2):
    if len(s1) > len(s2):
        s1, s2 = s2, s1
    distances = range(len(s1) + 1)
    for i2, c2 in enumerate(s2):
        distances_ = [i2+1]
        for i1, c1 in enumerate(s1):
            if c1 == c2:
                distances_.append(distances[i1])
            else:
                distances_.append(1 + min((distances[i1], distances[i1 + 1], distances_[-1])))
        distances = distances_
    return distances[-1]

def free_divergence_matrix(s1, s2):
    assert len(s1) == len(s2), 'Free divergence requires strings of equal length.'
    n = len(s1) + 1
    M = np.zeros((n, n))

    # First row and column
    M[0, :] = np.arange(n)
    M[:, 0] = np.arange(n)

    # Inner square
    for i in range(1, n-1):
        for j in range(1, n-1):
            if s1[i-1] == s2[j-1]:
                diag_penalty = M[i-1, j-1]
            else:
                diag_penalty = M[i-1, j-1] + 1
            M[i, j] = min(diag_penalty,
                          M[i-1, j] + 1,
                          M[i, j-1] + 1)

    # Last row
    for i in range(1, n-1):
        if s1[n-2] == s2[i-1]:
            diag_penalty = M[n-2, i-1]
        else:
            diag_penalty = M[n-2, i-1] + 1
        M[n-1, i] = min(diag_penalty,
                        M[n-2, i] + 1,
                        M[n-1, i-1])  # No penalty along last row

    # Last column
    for i in range(1, n-1):
        if s1[i-1] == s2[n-2]:
            diag_penalty = M[i-1, n-2]
        else:
            diag_penalty = M[i-1, n-2] + 1
        M[i, n-1] = min(diag_penalty,
                        M[i-1, n-1],  # No penalty along last col
                        M[i, n-2] + 1)

    # Last elt
    if s1[n-2] == s2[n-2]:
        diag_penalty = M[n-2, n-2]
    else:
        diag_penalty = M[n-2, n-2] + 1
    M[n-1, n-1] = min(diag_penalty,
                      M[n-2, n-1],  # No penalty along last col
                      M[n-1, n-2])  # No penalty along last row
    return M

def free_divergence(s1, s2):
    M = free_divergence_matrix(s1, s2)
    return M[-1, -1]

def generate_edit_dict(barcode_dict, mer, method):
    edit_dict = {0:defaultdict(list),1:defaultdict(list),2:defaultdict(list),3:defaultdict(list)}
    for bc in barcode_dict:
        for cur_mer in mer:
            if method=='Levenshtein':
                lev_dist = levenshteinDistance(bc,cur_mer)
            elif method=='FreeDivergence':
                lev_dist = free_divergence(bc,cur_mer)
            if lev_dist==0:
                edit_dict[0][cur_mer].append(bc)
            elif lev_dist==1:
                edit_dict[1][cur_mer].append(bc)
            elif lev_dist==2:
                edit_dict[2][cur_mer].append(bc)
            elif lev_dist==3:
                edit_dict[3][cur_mer].append(bc)
        print(bc,end=' ')
    return edit_dict

def generate_hamming_dict(barcode_list):
    # generate dist within 2 hamming distance from barcodes
    bc_len=len(barcode_list[0])
    edit_dict = {0:defaultdict(list),1:defaultdict(list)} #,2:defaultdict(list)
    processed=0
    for bc in barcode_list:
        edit_dict[0][bc].append(bc)
        indice = list(range(bc_len))
        for i in indice:
            bases=['A','G','C','T']
            ori_base = bc[i]
            bases.remove(ori_base)
            for base in bases:
                modi_bc = bc[:i]+base+bc[i+1:]
                #if bc not in edit_dict[1][modi_bc]:
                edit_dict[1][modi_bc].append(bc)
                '''
                if i!=indice[:-1]:
                    for j in indice[i+1:]:
                    #for j in modi_indice:
                        bases=['A','G','C','T']
                        ori_base = bc[j]
                        bases.remove(ori_base)
                        for base in bases:
                            modi_bc2 = modi_bc[:j]+base+modi_bc[j+1:]
                            if bc not in edit_dict[2][modi_bc2]:
                                edit_dict[2][modi_bc2].append(bc)
                '''
        processed+=1
        if processed%1000==0:
            print('Processed barcodes : %d'%processed)
    return edit_dict

def write_pkl(bcset_dir, dict_dir):
    bc_list = pd.read_csv(bcset_dir,names=['barcode'],index_col=False,sep='\t').barcode
    print('Total barcodes: %d'%len(bc_list))
    edit_dict=generate_hamming_dict(bc_list)
    del bc_list
    with open(dict_dir, 'wb') as f:
        pickle.dump(edit_dict, f)

if __name__ == '__main__':
    import sys
    write_pkl(sys.argv[1], sys.argv[2])
