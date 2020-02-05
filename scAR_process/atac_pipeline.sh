indir=$1
name=$2

#python ATAC_pipeline.py batch_plot -l $indir/sample.list -o $indir/atac_qc  -r Data/hg38_mm10_combined_epdnew.bed -f Data/Fragment_length_ratio.txt

python ATAC_pipeline.py batch_filter -l $indir/sample.list  -o $indir/atac_qc  -ct 7 -cu 500
python ATAC_pipeline.py batch_peak -l $indir/sample.list  -o $indir/atac_qc
python ATAC_pipeline.py merge -n $name -o $indir/atac_qc
