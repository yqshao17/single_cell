#!/usr/bin/env python

"""
Pre-processing including following steps:
1) sorting by read names
2) remove duplicates
3) convert to bed file
Created by Rongxin Fang
"""
import sys
import gzip
import pysam
import os
import collections 

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

    barcode_uniq = collections.defaultdict(lambda : 0)
    barcode_total = collections.defaultdict(lambda : 0)

    from argparse import ArgumentParser
    # parameters
    
    parser = ArgumentParser(description='snATAC-seq preprocessing')
    parser.add_argument('-i', '--input', help='input bam file', required=True)
    parser.add_argument('-o', '--output', help='output bed/bed.gz file', required=True)
    parser.add_argument('-m', '--mapq', help='min mappability score [30]', required=True)
    parser.add_argument('-t', '--threads', help='number of threads [3]', required=True)
    parser.add_argument('-f', '--flen', help='maximum fragment length [2000]', required=True)
    parser.add_argument('-e', '--elen', help='increase -e base pairs in each direction [75]', required=True)
    
    options = parser.parse_args()

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

    pre_barcode = ""
    cur_list = []
    
    for read in samfile:
        cur_barcode = read.qname.split(":")[0]
        rname  = str(read.reference_name)
        rstart = str(max(1, read.reference_end -5 - exlen if read.is_reverse else read.reference_start + 4 - exlen))
        rend   = str(min(genome_size[rname], read.reference_end - 5 + exlen if read.is_reverse else read.reference_start +4 + exlen))
        #rstart = str(max(1, read.reference_start - exlen if read.is_reverse else read.reference_start + 4 - exlen))
        #rend   = str(min(genome_size[rname], read.reference_end - 5 + exlen if read.is_reverse else read.reference_end + exlen))
        if(pre_barcode == cur_barcode):
            cur_list.append((rname, rstart, rend, cur_barcode))
            
            barcode_total[cur_barcode] += 1
        else:
            for item in set(cur_list):
                        barcode_uniq[item[3]] += 1
                        fout.write("\t".join(list(item))+"\n")
            pre_barcode = cur_barcode
            cur_list = [(rname, rstart, rend, cur_barcode)]
            barcode_total[cur_barcode] += 1
                                       
    # don't forget about the last barocode
    for item in set(cur_list):
        barcode_uniq[item[3]] += 1
        fout.write("\t".join(list(item))+"\n")

    samfile.close()
    
    # write down the qc file
    with open(output_bed+".qc", "w") as fout:
        for barcode in barcode_total:
            fout.write("\t".join([barcode, str(barcode_uniq[barcode]), str(1 - float(barcode_uniq[barcode])/barcode_total[barcode])]) + "\n") 
    
if __name__ == '__main__':
    main()
