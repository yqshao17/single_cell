
# ATAC
a1='CTGTCTCTTATACACATCTCCGAGCCCACGAGAC'
a2='CTGTCTCTTATACACATCTGACGCTGCCGACGA'

q5=15 # quality filter at 5'
q3=10 # quality filter at 3'
mlen1=36 # min length of read1: 16+20
mlen2=20 # min length of read2

#method='FreeDivergence'

method='hamming'

data_type = 'ATAC' #, RNA, ATAC_RNA
read_type = 'R1_R2' #R1, R2

bc_type =  'R1'
bc_starts = [0]
bc_len = [16]
bc_edit_dist=[1]

umi_type = 'R1'
umi_start = 16
umi_len = 0

read_start=[16,0]

reads_in_cells_thresh=0.92

#mapping
PE=True

min_mapq=30
max_flen=1000
atac_mlen2=20


