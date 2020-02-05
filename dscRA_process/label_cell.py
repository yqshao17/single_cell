## label cell type 
import sys

def get_sample_bc(sample_bc_list):
	sample_bc_dic={}
	with open(sample_bc_list) as input:
		for line in input:
			if line.strip():
				[sample, bc]=line.strip().split()
				sample_bc_dic[bc]=sample
	return sample_bc_dic

barcode_file=sys.argv[1]
sample_bc_list=sys.argv[2]
outfile=sys.argv[3]

sample_bc=get_sample_bc(sample_bc_list)
barcodes=[el.strip() for el in open(barcode_file).readlines()]
barcodes=list(set(barcodes))

label={}
for bc in barcodes:
    try:
        label[bc]=sample_bc[bc[-6:]]
    except:
        label[bc]='Not_known'
        
with open(outfile,'w') as output:
    for bc in barcodes:
        output.write(bc+'\t'+label[bc]+'\n')