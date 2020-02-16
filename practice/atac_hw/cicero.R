suppressPackageStartupMessages({
    library(BSgenome.Hsapiens.UCSC.hg19)
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    library(org.Hs.eg.db)
    library(cicero)
    library(magrittr)
    library('Matrix')
    })

###########
# import data
###########

# read in matrix data using the Matrix package
indata <- Matrix::readMM("data/Hematopoiesis-All/data/count.mtx") 
# binarize the matrix
indata@x[indata@x > 0] <- 1
# format cell info
cellinfo <- read.table("data/Hematopoiesis-All/data/barcodes.txt",comment.char = "")
row.names(cellinfo) <- cellinfo$V1
names(cellinfo) <- "cells"
# format peak info
peakinfo <- read.table("data/Hematopoiesis-All/Hematopoiesis-All_peak.bed")
names(peakinfo) <- c("chr", "bp1", "bp2")
peakinfo$site_name <- paste(peakinfo$chr, peakinfo$bp1, peakinfo$bp2, sep="_")
row.names(peakinfo) <- peakinfo$site_name

row.names(indata) <- row.names(peakinfo)
colnames(indata) <- row.names(cellinfo)

###########
# make input_cds
###########

fd <- methods::new("AnnotatedDataFrame", data = peakinfo)
pd <- methods::new("AnnotatedDataFrame", data = cellinfo)
input_cds <-  suppressWarnings(newCellDataSet(indata,
                            phenoData = pd,
                            featureData = fd,
                            expressionFamily=VGAM::binomialff(),
                            lowerDetectionLimit=0))
input_cds@expressionFamily@vfamily <- "binomialff"
input_cds <- monocle::detectGenes(input_cds)

# Ensure there are no peaks included with zero reads
input_cds <- input_cds[Matrix::rowSums(exprs(input_cds)) != 0,]
saveRDS(input_cds, "result/Hematopoiesis-All/hema_all_input.rds")

# filter peaks for lack of memory
input_cds=input_cds[Matrix::rowSums(exprs(input_cds)) >= 1000,]
saveRDS(input_cds, "result/Hematopoiesis-All/hema_all_input_f.rds")

###########
# make cicero_cds
###########
set.seed(2017)
input_cds <- detectGenes(input_cds)
input_cds <- estimateSizeFactors(input_cds)

# reduce dimension
#input_cds <- reduceDimension(input_cds, max_components = 2, num_dim=6,
#                      reduction_method = 'tSNE', norm_method = "none")
#tsne_coords <- t(reducedDimA(input_cds))
#row.names(tsne_coords) <- row.names(pData(input_cds))

# dimension reduction has already been done
cells = read.table('data/Hematopoiesis-All/GSE129785_scATAC-Hematopoiesis-All.cell_barcodes.txt', 
                   comment.char = "", head=1)
dimred <- data.frame(row.names = cells$Group_Barcode, 
                     cells$UMAP1, 
                     cells$UMAP2)
# make cicero cds
cicero_cds  <- make_cicero_cds(input_cds, k = 50, reduced_coordinates = dimred[colnames(input_cds),])
# save rds
saveRDS(cicero_rds, 'result/Hematopoiesis-All/hema_all_cicero_f.rds')


###########
#run cicero
###########

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
orgdb <- org.Hs.eg.db
bsgenome <- BSgenome.Hsapiens.UCSC.hg19
tssWindow <- 2500
flank <- 250*10^3
corCutOff <- 0.35

# choose genome region
bsgenome <- BSgenome.Hsapiens.UCSC.hg19
chromSizes <- seqlengths(bsgenome)[paste0("chr",c(1:22,"X"))]
genome <- data.frame(names(chromSizes),chromSizes)
rownames(genome) <- NULL

# run cicero, calculating connections
conns <- run_cicero(cicero_cds, genome) # Takes a few minutes to run
write.table(conns, 'result//Hematopoiesis-All/hema_all_f_conns.txt',quote = FALSE, sep='\t', row.names = FALSE)

#Annotate CDS
message("Annotating Cell Data Set...")
genes <- getTxDbGenes(txdb=txdb,orgdb=orgdb)
names(genes) <- genes$symbol
genes <- resize(genes, 1, "start") %>% resize(tssWindow * 2 + 1, "center")
geneDF <- data.frame(chromosome=seqnames(genes),start=start(genes),end=end(genes), gene=genes$symbol)
obj <- annotate_cds_by_site(input_cds, geneDF)

#Prepare for Co-Accessibility
nSites <- Matrix::colSums(assayData(obj)$exprs)
names(nSites) <- row.names(pData(obj))

#Cicero with Correlations
message("Calculating normalized gene activities...")
ciceroGA <- normalize_gene_activities(build_gene_activity_matrix(obj, conns, coaccess_cutoff = corCutOff), nSites)
saveRDS(ciceroGA,'result//Hematopoiesis-All/hema_all_f_ciceroGA.rds')

ciceroGA <- signif(ciceroGA,3)
writeMM(ciceroGA,file='result//Hematopoiesis-All/hema_all_f_ciceroGA_matrix.mtx')
write.table(row.names(ciceroGA),'result//Hematopoiesis-All/hema_all_f_ciceroGA_genes.txt',quote = FALSE, row.names=FALSE, col.names = FALSE)
write.table(colnames(ciceroGA),'result//Hematopoiesis-All/hema_all_f_ciceroGA_barcodes.txt',quote = FALSE, row.names=FALSE, col.names = FALSE)

