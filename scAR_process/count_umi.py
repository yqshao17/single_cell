import pysam
import sys
bamfile=sys.argv[1]
out_prefix=sys.argv[2]

samfile = pysam.AlignmentFile(bamfile,'rb')
c = 0
cell_barcodes = []
umis = []
barcodes_umi={}
NUM_READS = 0
NUM_UNIQ = 0
NUM_CHRM=0

for read in samfile:
    rname  = str(read.reference_name)
    NUM_READS += 1
    if read.is_unmapped or read.mate_is_unmapped: 
        continue
    if "chrM" in rname:
        NUM_CHRM+=1
        continue
    if (not read.is_secondary):
        NUM_UNIQ+=1
        seqname=read.qname
        barcode=':'.join(seqname.split(':')[-5:-1])
        #print(barcode)
        umi=seqname.split(':')[-1]
        barcodes_umi.setdefault(barcode,[])
        if not umi in barcodes_umi[barcode]:
            barcodes_umi[barcode].append(umi)

#print(barcodes_umi.keys())
        
with open(out_prefix+'_umi.txt','w') as output:
    for bc in barcodes_umi:
        output.write(bc+'\t'+str(len(barcodes_umi[bc]))+'\n')
    
with open(out_prefix+'_bam_stats.txt','w') as output:
    output.write("number of totol reads\t%d\n"%NUM_READS)
    output.write("number of chrM reads\t%d\n"%NUM_CHRM)
    output.write("number of uniquely mapped reads\t%d\n"%NUM_UNIQ)

    