suppressPackageStartupMessages({
    library(BSgenome.Hsapiens.UCSC.hg19)
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    library(org.Hs.eg.db)
    library(cicero)
    library(magrittr)
    library('Matrix')
    })

getTxDbGenes <- function(txdb = NULL, orgdb = NULL, gr = NULL, ignore.strand = TRUE){
    
    if (is.null(genome)) {
        if (is.null(txdb) | is.null(orgdb)) {
            stop("If no provided genome then you need txdb and orgdb!")
        }
    }
        
    if (is.null(gr)) {
        genes <- GenomicFeatures::genes(txdb)
    }else {
        genes <- suppressWarnings(subsetByOverlaps(GenomicFeatures::genes(txdb), gr, ignore.strand = ignore.strand))
    }
    
    if (length(genes) > 1) {
        mcols(genes)$symbol <- suppressMessages(mapIds(orgdb, 
            keys = mcols(genes)$gene_id, column = "SYMBOL", keytype = "ENTREZID", 
            multiVals = "first"))
        genes <- sort(sortSeqlevels(genes), ignore.strand = TRUE)
        names(genes) <- NULL
        out <- genes
    }else {
        out <- GRanges(seqnames(gr), ranges = IRanges(0, 0), gene_id = 0, symbol = "none")[-1]
    }

    return(out)

}

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
orgdb <- org.Hs.eg.db
bsgenome <- BSgenome.Hsapiens.UCSC.hg19
tssWindow <- 2500
flank <- 250*10^3
corCutOff <- 0.35

conns = read.table( 'result//Hematopoiesis-All/hema_all_f_conns.txt', sep='\t', header = 1)
input_cds_f = readRDS('result//Hematopoiesis-All/hema_all_input_f.rds')

#Annotate CDS
message("Annotating Cell Data Set...")
genes <- getTxDbGenes(txdb=txdb,orgdb=orgdb)
names(genes) <- genes$symbol
genes <- resize(genes, 1, "start") %>% resize(tssWindow * 2 + 1, "center")
geneDF <- data.frame(chromosome=seqnames(genes),start=start(genes),end=end(genes), gene=genes$symbol)
obj <- annotate_cds_by_site(input_cds_f, geneDF)

#Prepare for Co-Accessibility
nSites <- Matrix::colSums(assayData(obj)$exprs)
names(nSites) <- row.names(pData(obj))

#Cicero with Correlations
message("Calculating normalized gene activities...")
ciceroGA <- normalize_gene_activities(build_gene_activity_matrix(obj, conns, coaccess_cutoff = corCutOff), nSites)

writeMM(ciceroGA,file='result//Hematopoiesis-All/hema_all_f_ciceroGA.mtx')
write.table(rownames(ciceroGA),'result//Hematopoiesis-All/hema_all_f_ciceroGA_genes.txt',quote = FALSE)
write.table(colnames(ciceroGA),'result//Hematopoiesis-All/hema_all_f_ciceroGA_barcodes.txt',quote = FALSE)