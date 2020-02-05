trim_TSO=False
split_atac_rna=False
TSO='AAGCAGTGGTATCAACGCAGAGT'
TSO_rev='ACTCTGCGTTGATACCACTGCTT'
# ATAC
#a1='CTGTCTCTTATACACATCTCCGAGCCCACGAGAC'
#a2='CTGTCTCTTATACACATCTGACGCTGCCGACGA'

# RNA
#a1='TCACTCTGCGTTGATACCACTGCTT'
a1=''
a2='AAAAAAAAAAAAAAAAAAAA'
q5=15 # quality filter at 5'
q3=10 # quality filter at 3'
mlen1=28 # min length of read1
mlen2=20 # min length of read2

#method='FreeDivergence'

method='hamming'

data_type = 'RNA' #, ATAC, ATAC_RNA
read_type = 'R2' #R1, R1_R2

bc_type =  'R1'
bc_starts = [0]
bc_len = [16]
bc_edit_dist=[1]

umi_type = 'R1'
umi_start = 16
umi_len = 12

read_start=[0,0]

#mapping
PE=False
ATAC_PE=True
RNA_PE=False
nthreads=4

min_mapq=30
max_flen=1000
atac_mlen2=20


