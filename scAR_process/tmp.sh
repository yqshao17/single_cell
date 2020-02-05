source py2.sh

export HDF5_USE_FILE_LOCKING='FALSE'
sample=$1
indir=/Share2/home/zhangqf5/yanqiu/scAR/output/$sample/atac_qc
outdir=$indir

function snappmat(){
sample=$1
cmd="snaptools snap-add-pmat \
--snap-file=$indir/${sample}.snap \
--peak-file=$indir/${sample}_summits_extend.bed \
--verbose=True"
$cmd
}


snappmat $sample
