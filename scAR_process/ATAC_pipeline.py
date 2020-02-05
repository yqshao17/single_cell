

from split_sample import get_prefix

from argparse import ArgumentParser
import os 

parser = ArgumentParser(description='atac bam filtering')
parser.add_argument('mode', help='filter_bc or plot_umi_tss or  filter_umi_tss')

parser.add_argument('-l', '--inlist', help='input sample bam list for pre-filter')
parser.add_argument('-o', '--outdir', help='output directory of filtered bam file', default='./')
parser.add_argument('-t', '--threads', help='number of threads', default='4')

parser.add_argument('-i', '--inbam', help='bam file for umi_tss filtering')
parser.add_argument('-r', '--tssref', help='tss reference', default='Data/hg38_mm10_combined_epdnew.bed')
parser.add_argument('-n', '--name', help='output name prefix')
parser.add_argument('-f', '--fragref', help='fragment length ratio reference', default='Data/Fragment_length_ratio.txt')

parser.add_argument('-cr', '--cutref', help='info of UMI and TSS enrichment score')
parser.add_argument('-ct', '--cuttss', help='cutoff of tss score')
parser.add_argument('-cu', '--cutumi', help='cutoff of UMI')

args = parser.parse_args()

mode=args.mode.lower()
nthreads=int(args.threads)
outdir=args.outdir
tss_ref=args.tssref
frag_ref=args.fragref

bsub_script='bsub -q Z-ZQF -e atac.err -o atac.out -n %d ' % nthreads



if mode=='batch_plot':
    # -l --inlist
    # -o
    prefix_list=get_prefix(args.inlist)
    for prefix in prefix_list:
        os.system(bsub_script+"""python ATAC_analysis.py plot_umi_tss -i {0}/{1}_ATAC_filtered.bam -o {0} -n {1}_ATAC -r {2} -f {3}""".format(outdir, prefix, tss_ref, frag_ref))
    
elif mode=='batch_filter':
    # -l --inlist
    # -o
    # -ct, -cu
    prefix_list=get_prefix(args.inlist)
    for prefix in prefix_list:
        os.system(bsub_script+"""python ATAC_analysis.py filter_umi_tss -i {0}/{1}_ATAC_filtered.sorted.bam  -o {0} -n {1}_ATAC -cr {0}/{1}_ATAC_TssEnrichUMI.txt -ct {2} -cu {3}""".format(outdir, prefix, args.cuttss, args.cutumi))
        

elif mode=='batch_peak':
    # -l --inlist
    # -o
    prefix_list=get_prefix(args.inlist)
    for prefix in prefix_list:
        os.system(bsub_script+"""python ATAC_analysis.py peak_count -i {0}/{1}_ATAC_filtered.sorted.filtered.bam -o {0} -n {1}_ATAC""".format(outdir, prefix))
    

elif mode=='merge':
    # -n, -o
    os.system(bsub_script+"bash bash/atac_merge.sh %s %s %s %s"%(args.name, outdir, tss_ref, frag_ref))
    
