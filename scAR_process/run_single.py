import os
#from multiprocessing import Pool

import argparse
import datetime
import importlib
from parse_barcode import preprocess_fastq
from split_sample import split_atac_rna
from tools import trim, trim_RNA, run_star, sort_sam, make_combined_genome,make_gtf_annotations,generate_STAR_index
import process_rna
import process_atac
from filter_bam import filter_bam

from analysis2 import generate_single_dge_report

parser = argparse.ArgumentParser()

parser.add_argument('mode', help="""Mode: one of "all", "preproc", "postproc".
"all" runs the entire pipeline without split by bc1.
"preproc" runs trimming and parsing barcode.
"postproc" runs split atac_rna, mapping, parsing alignment and summary.
""")

# single sample mode
parser.add_argument('--fq1', help='fastq1 - mRNA or DNA reads')
parser.add_argument('--fq2', help='fastq2 - reads contain UMI and barcodes')
parser.add_argument('--bc4', help='bc4 - each sample has a unique bc4, for identification of samples')


parser.add_argument('--output_dir', help='output dir')
parser.add_argument('--prefix', help='output prefix')

parser.add_argument('--bcset', default='./barcodes/', help='directory of barcodes and barcode pkl files (bc1 is actually bc3 in experiments)')

parser.add_argument('--genome_dir', default='./', help='path containing reference genome')

parser.add_argument('--nthreads', default='4', help='number of threads to use')
# keep type
parser.add_argument('--keeptype', default='ATAC_RNA', help='kept data type')
parser.add_argument('--config', default='config', help='config file name')

args = parser.parse_args()
mode = args.mode.lower()


print(args)
print(mode)

config = importlib.import_module(args.config)

'''
if mode == 'all' or mode=='preproc':   
    ### trim ###
    print(datetime.datetime.now(), 'Trimming reads...')
    trim(args.output_dir, args.prefix, args.fq1, args.fq2, config, int(args.nthreads))
    
    ### parsing barcode ###
    print(datetime.datetime.now(), 'Parsing barcode...')
    preprocess_fastq(args.output_dir, args.prefix, args.bcset, args.bc4, config, int(args.nthreads))
'''

if mode == 'all' or mode == 'postproc':
    ### split atac and rna ###
    sample=args.prefix
    
    print(datetime.datetime.now(), 'Spliting ATAC and RNA...')
    split_atac_rna(args.output_dir, sample, config) 
    

    if 'RNA' in args.keeptype:
        print(datetime.datetime.now(),'Processing RNA...')
        print(datetime.datetime.now(),'Re-trimming RNA...')
        trim_RNA(args.output_dir, sample+'_RNA', config, int(args.nthreads))
        print(datetime.datetime.now(), 'STAR mapping...')
        run_star(args.genome_dir, args.output_dir, sample+'_RNA', int(args.nthreads), config.RNA_PE)
        print(datetime.datetime.now(), 'Sorting sam...')
        sort_sam(args.output_dir,sample+'_RNA',int(args.nthreads))
        print(datetime.datetime.now(),'Getting molecular info...')
        process_rna.molecule_info(args.genome_dir, args.output_dir, sample+'_RNA', int(args.nthreads))
        print(datetime.datetime.now(),'Getting summary statistics...')
        generate_single_dge_report(args.output_dir,args.genome_dir,'RNA',sample+'_RNA')
        print(datetime.datetime.now(),sample+' RNA is successfully done!')

    if 'ATAC' in args.keeptype:
        print(datetime.datetime.now(),'Processing ATAC...')
        ### mapping ###
        print(datetime.datetime.now(), 'STAR mapping...')
        run_star(args.genome_dir, args.output_dir, sample+'_ATAC', int(args.nthreads), config.ATAC_PE)
        print(datetime.datetime.now(), 'Sorting sam...')
        sort_sam(args.output_dir,sample+'_ATAC',int(args.nthreads))
        ### parsing alignment ###
        print(datetime.datetime.now(),'Getting molecular info...')
        process_atac.molecule_info(args.genome_dir, args.output_dir, sample+'_ATAC', config, int(args.nthreads))
        print(datetime.datetime.now(),'Getting summary statistics...')
        ### summary ###
        generate_single_dge_report(args.output_dir,args.genome_dir,'ATAC',sample+'_ATAC')
        ### filter bam file by filtered barcode ###
        if not os.path.exists(args.output_dir+'/atac_qc'):
            os.makedirs(args.output_dir+'/atac_qc')
        filter_bam(args.output_dir, sample+'_ATAC', config, int(args.nthreads)) # very slow!!
        print(datetime.datetime.now(),sample+' ATAC is successfully done!')

