# single_cell

**scAR** and **dscRA** are two versions of processing raw sequencing data for two methods: split-seq and droplet-seq. The read structure is shown as following.

There are different modes like pre-process, repslit, post-process, etc.

Usage is
```
python pipeline.py all_resplit \
  --raw_list raw.list \
  --sample_bc_list sample_bc.list \
  --output_dir outdir \
  --genome_dir INDEX_hg \
  --bcset barcode \
  --config config_file \
  --keeptype ATAC \
  --nthreads 10

```


Retulatory is a method in trial: cells are first grouped by RNA, and then each group is considered as a pseudocell from which we can get gene-peak correlation as features. So that we convert separate RNA and ATAC features to RNA-ATAC regulatory features, which be used for downstream analysis like clustering and trajectory analysis.


Joint analysis includes clustering, cis-regulatory elements identification, transcription factor motif enrichment, disease-related snps enrichment, eqtl enrichment and trajectory analysis. 
