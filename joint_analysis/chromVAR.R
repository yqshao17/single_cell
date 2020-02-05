#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(chromVAR)
    library(motifmatchr)
    library(SummarizedExperiment)
    library(Matrix)
})

# library(gplots) 
# library(RColorBrewer)
library(optparse)

option_list = list(
    make_option(c("-d", "--data"), type="character"),
    make_option(c("-p", "--peak"), type="character"),
    make_option(c("-o", "--outdir"), type='character', default='./'),
    make_option(c("--name"), type='character', default=''),
    make_option(c("-g", "--genome"), type='character', default='hg19'),
    make_option(c("-m", "--motif"), type='character', default='jaspar'),
    make_option(c("--process"), type='logical', default=FALSE)
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


read_count <- function(input){
    if(file_test("-f", input)){
    cellcount <- read.csv(file = input, row.names = 1, sep='\t', header=T)
    }else{
    files=list.files(path = input,full.names = TRUE)
    for (f in files){
        if (grepl("count",f)){
            cellcount <- readMM(f)
        } else if (grepl("barcode",f)){
            rna_cells = read.table(f)
        } else if (grepl("peak",f)){
            genes = read.table(f)
        }
    }
    colnames(cellcount) <- rna_cells$V1
    rownames(cellcount) <- genes$V1
    }
    return(cellcount)
}

apply_chromVAR=function(counts, peaks, motifs, genome){
    counts = SummarizedExperiment(assays=list(counts=as.matrix(counts)),rowRanges=peaks)
    counts = addGCBias(counts, genome = genome)
    counts = filterPeaks(counts, non_overlapping=T)
    print('filtered counts dim')
    print(dim(counts))
    motif_ix <- matchMotifs(motifs, counts, genome = genome)
    dev <- computeDeviations(object = counts, 
                                 annotations = motif_ix)
    return(dev)
}

process = function(counts, peaks, cutoff=800){
    index = which(width(peaks) <= cutoff)
    peaks = peaks[index]
    peaks = resize(peaks, width = 500, fix = "center")
    counts = counts[index,]
    return(list(counts, peaks))
}

if(opt$genome=='hg19'){
    library(BSgenome.Hsapiens.UCSC.hg19)
    genome = BSgenome.Hsapiens.UCSC.hg19
    species = 'Homo sapiens'
} else if (opt$genome=='mm9'){
    library(BSgenome.Mmusculus.UCSC.mm9)
    genome = BSgenome.Mmusculus.UCSC.mm9
    species = 'Mus musculus'
} else if (opt$genome=='mm10'){
    library(BSgenome.Mmusculus.UCSC.mm10)
    genome = BSgenome.Mmusculus.UCSC.mm10
    species = 'Mus musculus'
} else if (opt$genome=='hg38'){
    library(BSgenome.Hsapiens.UCSC.hg38)
    genome = BSgenome.Hsapiens.UCSC.hg38
    species = 'Homo sapiens'    
} else { print('Self defined genome')}


if(opt$motif=='jaspar'){
    motifs <- getJasparMotifs(species=species)
} else if (opt$motif=='cisBP'){
    library('chromVARmotifs')
    if(species=='Homo sapiens'){
        data("human_pwms_v1")
        motifs = human_pwms_v1
    } else if (species=='Mus musculus'){
        data("mouse_pwms_v1")
        motifs = mouse_pwms_v1
    }
}

print(paste0('Motif database: ', opt$motif))
print(paste0('Species: ', species))
    

data_file = opt$data
peak_file = opt$peak
outdir = opt$outdir
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

#counts = read.table(data_file, row.names=1, header=T)
counts = read_count(data_file)
peaks = getPeaks(peak_file, sort_peaks=T) 

peak_index = as.data.frame(peaks)
peak_index = paste0(peak_index$seqnames, ':' ,peak_index$start-1, '-', peak_index$end)

index = row.names(counts)%in%peak_index
counts = counts[index,]

dev = apply_chromVAR(counts, peaks, motifs, genome)
var = computeVariability(dev)
write.table(deviationScores(dev), paste0(outdir, '/', opt$name, 'dev.txt'), quote=F, sep='\t')
write.table(var, paste0(outdir, '/', opt$name, 'var.txt'), quote=F, sep='\t')

