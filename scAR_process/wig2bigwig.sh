bigWigToWig ~/yanqiu/scAR/output/ZJS191030/ZJS191030_R/mapping/K562_RNA.mapped.sorted.bw ~/yanqiu/scAR/output/ZJS191030/ZJS191030_R/mapping/K562_RNA.mapped.sorted.wig
grep -v mm ~/yanqiu/scAR/output/ZJS191030/ZJS191030_R/mapping/K562_RNA.mapped.sorted.wig | grep chr | grep -v chrM > ~/yanqiu/scAR/output/ZJS191030/ZJS191030_R/mapping/K562_RNA.mapped.sorted.tmp.wig
vi
wigToBigWig ~/yanqiu/scAR/output/ZJS191030/ZJS191030_R/mapping/K562_RNA.mapped.sorted.tmp.wig Data/hg38.chrom.sizes ~/yanqiu/scAR/output/ZJS191030/ZJS191030_R/mapping/K562_RNA.mapped.sorted.bw
