#!/usr/bin/env python

"""
Pre-processing including following steps:
1) sorting by read names
2) convert to bed file
"""
import sys
import gzip
import pysam
import os

def is_sorted_queryname(header):
    """
    Check if bam fiel is sorted by read name.
    """
    if("HD" in header):
        if("SO" in header["HD"]):
            if(header["HD"]["SO"] == "queryname"):
                return True
    return False
    
def main():

    from argparse import ArgumentParser

    parser = ArgumentParser(description='snATAC-seq preprocessing')
    parser.add_argument('-i', '--input', help='input bam file', required=True)
    parser.add_argument('-o', '--output', help='output bed/bed.gz file', required=True)
    parser.add_argument('-m', '--mapq', help='min mappability score [30]', required=True)
    parser.add_argument('-t', '--threads', help='number of threads [3]', required=True)
    parser.add_argument('-f', '--flen', help='maximum fragment length [2000]', required=True)
    parser.add_argument('-e', '--elen', help='increase -e base pairs in each direction [75]', required=True)
    
    options = parser.parse_args()

    Num_ChrM = 0
    Num_written = 0
    
    num_threads = 1
    min_mapq = 30
    max_flen = 2000
    exlen = 75    
    # input parsing
    input_bam = options.input
    output_bed = options.output
    num_threads = int(options.threads)
    min_mapq = int(options.mapq)
    max_flen = int(options.flen)
    exlen = int(options.elen)
    
        
    if output_bed.endswith(".gz"):
        fout = gzip.open(output_bed, "wb")
    else:
        fout = open(output_bed, "w")
    
    # start reading the bam
    samfile = pysam.AlignmentFile(input_bam, "rb")
    genome_size = dict([[item["SN"], int(item["LN"])] for item in samfile.header["SQ"]])
 
    for read in samfile:
        if read.is_unmapped or read.mate_is_unmapped: continue        
        if (not read.is_secondary):
            rname  = str(read.reference_name)         
            if "chrM" in rname:
                Num_ChrM += 1
                continue          
            rstart = str(max(1, read.reference_end -5 - exlen if read.is_reverse else read.reference_start + 4 - exlen))
            rend   = str(min(genome_size[rname], read.reference_end - 5 + exlen if read.is_reverse else read.reference_start +4 + exlen))
            item=(rname, rstart, rend)
            fout.write("\t".join(list(item))+"\n")
            Num_written += 1
    samfile.close()
    
    with open(output_bed+'.stats','w') as output:
        output.write('Num_ChrM\t%d\n'%Num_ChrM)
        output.write('Num_written\t%d'%Num_written)
    
if __name__ == '__main__':
    main()
