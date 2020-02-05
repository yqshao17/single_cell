library(Seurat)
library(Matrix)

args <- commandArgs(trailingOnly = TRUE)
input <- args[1]
outdir <- args[2]

if(file_test("-f", input)){
    cellcount <- read.csv(file = input, row.names = 1, sep='\t')
}else{
    files=list.files(path = input,full.names = TRUE)
    for (f in files){
        if (grepl("count",f)){
            cellcount <- readMM(f)
        } else if (grepl("barcode",f)){
            rna_cells = read.table(f)
        } else if (grepl("gene",f)|(grepl("peak",f))){
            genes = read.table(f)
        }
    }
colnames(cellcount) <- rna_cells$V1
rownames(cellcount) <- genes$V1
}


seuset <- CreateSeuratObject(
  counts = cellcount,
  assay = "RNA",
  min.cells = 10,
  min.features = 100
)

seuset <- NormalizeData(object = seuset,  normalization.method = "LogNormalize", scale.factor = 10000)
seuset <- FindVariableFeatures(object = seuset,selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(seuset)
seuset <- ScaleData(seuset, features = all.genes)
seuset <- RunPCA(seuset, features = VariableFeatures(object = seuset))

#Cluster the cells
seuset <- FindNeighbors(seuset, dims = 1:10)
seuset <- FindClusters(seuset, resolution = 0.5)

seuset <- RunUMAP(seuset, dims = 1:10)

pdf(paste(outdir, '/cluster_umap.pdf', sep=''))
DimPlot(seuset, reduction = "umap",label = TRUE)
dev.off()

pdf(paste(outdir, '/group_cells_umap.pdf', sep=''))
seuset <- FindClusters(seuset, resolution = 20)
DimPlot(seuset, reduction = "umap",label = TRUE)
dev.off()

as.loom(seuset, assay = 'RNA', filename = paste(outdir, '/group_labeled.loom', sep=''), 
        chunk.dims = NULL, chunk.size = NULL, overwrite = TRUE, verbose = TRUE)
