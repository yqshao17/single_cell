import os
import time
import argparse
import datetime

parser = argparse.ArgumentParser()

parser.add_argument('mode', help="""Mode: one of "mkref", "resplit", "all_single", "all_batch", "all_resplit", "preproc_single", "preproc_batch", "postproc_single", "postproc_batch".
"all" runs the entire pipeline.
"preproc" runs trimming adapter and parsing barcodes.
"postproc" runs all steps from spilt RNA and ATAC, assuming that output_dir/parse_bc.
"single" input only one sample
"batch" input sample list
"resplit" resplit samples by bc1, input sample_bc_list, new sample name cannot be the same with any of the origial names.
""")

### required para
parser.add_argument('--output_dir', help='output dir')
parser.add_argument('--genome_dir', default='./', help='path containing reference genome')
parser.add_argument('--nthreads', default='4', help='number of threads to use')

### optional para
# keep type
parser.add_argument('--keeptype', default='ATAC_RNA', help='kept data type')
parser.add_argument('--config', help='config file name', default='config_RNA')


# single sample mode
parser.add_argument('--fq1', help='fastq1 - mRNA or DNA reads')
parser.add_argument('--fq2', help='fastq2 - reads contain UMI and barcodes')
parser.add_argument('--bc4', help='bc4 - each sample has a bc4')
parser.add_argument('--bcset', default='./barcodes/', help='directory of barcodes and barcode pkl files (bc1 is actually bc3 in experiments)')

parser.add_argument('--prefix', help='output prefix')
# batch mode
parser.add_argument('--raw_list', help='fastq list format: name bc4 fq1 fq2') 
parser.add_argument('--preproc_list', help='preprocessed prefix list format: prefix') # has been trimmed, barcode_headed or resplit sample list
#parser.add_argument('--resplit_list', help='preprocessed and resplit sample list format: sample')# has been trimmed, barcode_headed and resplit
parser.add_argument('--sample_bc_list', help='bc1 information format(only if mode is all_resplit): sample bc1')
# mkref mode
parser.add_argument('--genome', nargs='*', help='name(s) of genome(s)/species')
parser.add_argument('--genes', nargs='*', help='GTF file(s) with gene annotations')
parser.add_argument('--fasta', nargs='*', help='fasta file for genome(s)')
parser.add_argument('--splicing',default='True', help='whether genome has splicing')
parser.add_argument('--genomeSAindexNbases',default='14', help='set this to min(14, floor(log2(GenomeLength)/2 - 1))')


args = parser.parse_args()
mode = args.mode.lower()
output_dir=args.output_dir
genome_dir=args.genome_dir
nthreads=int(args.nthreads)
keeptype=args.keeptype
bcset=args.bcset
config_name = args.config

'''
def set_config(config_name):
    config_text = open('config.py','w')
    config_text.write('import importlib\n')
    config_text.write("config = importlib.import_module('%s')"%config_name)
    config_text.close()
    
set_config(config_name)
'''

from split_sample import resplit_bc
from tools import get_rawlist, get_prefix

outdirs=[output_dir, output_dir+'/trimmed', output_dir+'/parse_bc',
		output_dir+'/mapping',output_dir+'/molecule_info',output_dir+'/analysis']

for d in outdirs:
	if not os.path.exists(d):
		os.makedirs(d)

single_script='run_single.py'
bsub_script='bsub -q Z-ZQF -n %d' % nthreads

if mode == 'mkref':
	# required input: genome(species), genes(gtf files), fasta(reference fa)
	
	# Generate genome
	make_combined_genome(args.genome, args.fasta, args.output_dir)
	# Make a combine annotation from GTF
	make_gtf_annotations(args.genome, args.genes, args.output_dir, args.splicing)
	# Index genome with star
	generate_STAR_index(args.output_dir, args.nthreads, args.genomeSAindexNbases, args.splicing)
else:
    # write config.py
    os.system('cp %s.py %s/config.py'%(config_name, output_dir))

if mode.startswith('preproc'):
	submode=mode.split('_')[1]
	if submode=='single':
		# required input: fq1, fq2, bc4, prefix
		os.system('python %s %s --fq1 %s --fq2 %s --bc4 %s --bcset %s--output_dir %s --genome_dir %s --prefix %s --nthreads %d --config %s'
			%(single_script,'preproc',args.fq1,args.fq2,args.bc4,bcset,output_dir,genome_dir,args.prefix,nthreads,config_name))
	
	elif submode=='batch':
		# required input: raw_list
		preproc_list=[]
		raw_list=get_rawlist(args.raw_list)
		for items in raw_list:
			[prefix, bc4, fq1, fq2] = items
			preproc_list.append(prefix)
			os.system(bsub_script + ' -e preproc_batch.err -o preproc_batch.out ' +
				'python %s %s --fq1 %s --fq2 %s --bc4 %s --bcset %s --output_dir %s --genome_dir %s --prefix %s --nthreads %d --config %s'
				%(single_script,'preproc',fq1,fq2,bc4,bcset,output_dir,genome_dir,prefix,nthreads,config_name))
		f=open(output_dir+'/preproc.list','w')
		f.write(preproc_list[0])
		for prefix in preproc_list[1:]:
			f.write('\n'+prefix)
		f.close()

elif mode.startswith('postproc'):
	submode=mode.split('_')[1]	
	if submode=='single':
		# required input: prefix
		os.system('python %s %s --output_dir %s --genome_dir %s --prefix %s --nthreads %d --keeptype %s --config %s'
			%(single_script,'postproc',output_dir,genome_dir,args.prefix,nthreads,keeptype,config_name))
	
	elif submode=='batch':
		# required input: preproc_list/sample_list
		preproc_list=get_prefix(args.preproc_list)
		for prefix in preproc_list:
			os.system(bsub_script + ' -e postproc_batch.err -o postproc_batch.out ' +
				'python %s %s --output_dir %s --genome_dir %s --prefix %s --nthreads %d --keeptype %s --config %s'
				%(single_script,'postproc',output_dir,genome_dir,prefix,nthreads, keeptype,config_name))

elif mode=='resplit':
	# required input: preproc_list, sample_bc_list
	resplit_bc(output_dir, args.sample_bc_list, args.preproc_list)

elif mode.startswith('all'):
	submode=mode.split('_')[1]
	if submode=='single':
		# required input: fq1, fq2, bc4, prefix
		os.system('python %s %s --fq1 %s --fq2 %s --bc4 %s --bcset %s --output_dir %s --genome_dir %s --prefix %s --nthreads %d --keeptype %s --config %s'
			%(single_script,'all',args.fq1,args.fq2,args.bc4,bcset,output_dir,genome_dir,args.prefix,nthreads,keeptype,config_name))
	
	elif submode=='batch':
		# required input: raw_list
		raw_list=get_rawlist(args.raw_list)
		for items in raw_list:
			[prefix, bc4, fq1, fq2] = items
			os.system(bsub_script + ' -e all_batch.err -o all_batch.out ' +
				'python %s %s --fq1 %s --fq2 %s --bc4 %s --bcset %s --output_dir %s --genome_dir %s --prefix %s --nthreads %d --keeptype %s --config %s'
				%(single_script,'all',fq1,fq2,bc4,bcset,output_dir,genome_dir,prefix,nthreads,keeptype,config_name))
	
	elif submode=='resplit':
		# required input: raw_list, sample_bc_list
		# preproc
		preproc_list=[]
		raw_list=get_rawlist(args.raw_list)
		print(datetime.datetime.now(), 'Preprocessing ...')
		for items in raw_list:
			[prefix, bc4, fq1, fq2] = items
			preproc_list.append(prefix)
			#tmp_script='bsub -q Z-BNODE -n %d' % nthreads
			os.system(bsub_script + ' -e preproc_batch.err -o preproc_batch.out ' +
				'python %s %s --fq1 %s --fq2 %s --bc4 %s --bcset %s --output_dir %s --genome_dir %s --prefix %s --nthreads %d --config %s'
				%(single_script,'preproc',fq1,fq2,bc4,bcset,output_dir,genome_dir,prefix,nthreads,config_name))	
		# wait all previous samples
		done=0
		while done<len(preproc_list):
			done=0
			for prefix in preproc_list:
				if os.path.exists(output_dir+'/parse_bc/'+prefix+ '_pipeline_stats.txt'):
					done+=1
			time.sleep(15)
		f=open(output_dir+'/preproc.list','w')
		f.write(preproc_list[0])
		for prefix in preproc_list[1:]:
			f.write('\n'+prefix)
		f.close() # !!!!!!!!!
		# resplit and wait completion
		print(datetime.datetime.now(), 'Resplit samples by bc...')
		resplit_bc(output_dir, args.sample_bc_list, output_dir+'/preproc.list')
		sample_list=get_prefix(output_dir+'/sample.list')
		# postproc
		print(datetime.datetime.now(), 'Postprocessing ...')
		for prefix in sample_list:
			os.system(bsub_script + ' -e postproc_batch.err -o postproc_batch.out ' +
				'python %s %s --output_dir %s --genome_dir %s --prefix %s --nthreads %d --keeptype %s --config %s'
				%(single_script,'postproc',output_dir,genome_dir,prefix,nthreads,keeptype,config_name))

