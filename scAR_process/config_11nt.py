trim_TSO=False
TSO='AAGCAGTGGTATCAACGCAGAGT'
TSO_rev='ACTCTGCGTTGATACCACTGCTT'
a1="CTGTCTCTTATACACATCT"
a2="CTGTCTCTTATACACATCTGACGCTGCCGACGA"
q5=15 # quality filter at 5'
q3=10 # quality filter at 3'
mlen1=25 # min length of read1
mlen2=86 # min length of read2

method='FreeDivergence'

R1A='TACTGCAGCTGAACCTC' # spacer1 
tag_seq='AGATGTGTATAAGAGACAG'
tag_edit_dist=(2,9)
tag_length=19
bc_edit_dist=(1,1,1)

umi_bc_len = [11,8,8,6]
umi_bc_starts = [0,11,36,61]
tag_start = 67

#umi_bc_len = [10,8,8,6]
#umi_bc_starts = [0,10,35,60]
#tag_start=66

reads_in_cells_thresh=0.92

#mapping
PE=False
ATAC_PE=True
RNA_PE=False
nthreads=4

min_mapq=30
max_flen=1000
atac_mlen2=20


data_config={
    'ps_cm':{'rna_path':'/Share2/home/zhangqf5/yanqiu/scAR/datasets/joint_ATAC_RNA/Paired_seq/Cell_Mix/Cell_Mix_RNA/',
            'atac_path':'/Share2/home/zhangqf5/yanqiu/scAR/datasets/joint_ATAC_RNA/Paired_seq/Cell_Mix/Cell_Mix_DNA/'},
    #'ps_ad':{'rna_path':'/Share2/home/zhangqf5/yanqiu/scAR/datasets/joint_ATAC_RNA/Paired_seq/Adult_Cerebrail_Cortex/Adult_CTX_RNA/',
    #        'atac_path':'/Share2/home/zhangqf5/yanqiu/scAR/datasets/joint_ATAC_RNA/Paired_seq/Adult_Cerebrail_Cortex/Adult_CTX_DNA/'},
    #'ps_ft':{'rna_path':'/Share2/home/zhangqf5/yanqiu/scAR/datasets/joint_ATAC_RNA/Paired_seq/Fetal_Forebrain/FB_RNA/',
    #        'atac_path':'/Share2/home/zhangqf5/yanqiu/scAR/datasets/joint_ATAC_RNA/Paired_seq/Fetal_Forebrain/FB_DNA/'},
    'ss_cm':{'rna_path':'/Share2/home/zhangqf5/yanqiu/scAR/datasets/joint_ATAC_RNA/SNARE-seq/CellLineMixture/RNA/GSE126074_CellLineMixture_SNAREseq_cDNA_counts.tsv.gz',
            'atac_path':'/Share2/home/zhangqf5/yanqiu/scAR/datasets/joint_ATAC_RNA/SNARE-seq/CellLineMixture/ATAC/GSE126074_CellLineMixture_SNAREseq_chromatin_counts.tsv.gz'},
    'ss_ad':{'rna_path':'/Share2/home/zhangqf5/yanqiu/scAR/datasets/joint_ATAC_RNA/SNARE-seq/AdBrainCortex/RNA/',
            'atac_path':'/Share2/home/zhangqf5/yanqiu/scAR/datasets/joint_ATAC_RNA/SNARE-seq/AdBrainCortex/ATAC/'},
    'ss_p0':{'rna_path':'/Share2/home/zhangqf5/yanqiu/scAR/datasets/joint_ATAC_RNA/SNARE-seq/P0_BrainCortex/RNA/',
            'atac_path':'/Share2/home/zhangqf5/yanqiu/scAR/datasets/joint_ATAC_RNA/SNARE-seq/P0_BrainCortex/ATAC/'},
    'sc_cm':{'rna_path':'/Share2/home/zhangqf5/yanqiu/scAR/datasets/joint_ATAC_RNA/sciCAR/A549/RNA',
            'atac_path':'/Share2/home/zhangqf5/yanqiu/scAR/datasets/joint_ATAC_RNA/sciCAR/A549/ATAC'},
    'sc_cm_only':{'rna_path':'/Share2/home/zhangqf5/yanqiu/scAR/datasets/joint_ATAC_RNA/sciCAR/A549/RNA_only/',
                 'atac_path':'/Share2/home/zhangqf5/yanqiu/scAR/datasets/joint_ATAC_RNA/sciCAR/A549/ATAC_only/'},
    'sc_mk':{'rna_path':'/Share2/home/zhangqf5/yanqiu/scAR/datasets/joint_ATAC_RNA/sciCAR/mouse_kidney/RNA',
            'atac_path':'/Share2/home/zhangqf5/yanqiu/scAR/datasets/joint_ATAC_RNA/sciCAR/mouse_kidney/ATAC'},
    'exp_cm':{'rna_path':'/Share2/home/zhangqf5/yanqiu/scAR/datasets/joint_ATAC_RNA/ZJS191030/ZJS191030_R/',        
            'atac_path':'/Share2/home/zhangqf5/yanqiu/scAR/datasets/joint_ATAC_RNA/ZJS191030/ZJS191030_A/'},
    'exp_cancer':{'rna_path':'/Share2/home/zhangqf5/yanqiu/scAR/datasets/joint_ATAC_RNA/WYF191116/RNA',        
            'atac_path':'/Share2/home/zhangqf5/yanqiu/scAR/datasets/joint_ATAC_RNA/WYF191116/ATAC/peak'},
    'exp_c_bin':{'atac_path':'/Share2/home/zhangqf5/yanqiu/scAR/datasets/joint_ATAC_RNA/WYF191116/ATAC/bin'}
}
