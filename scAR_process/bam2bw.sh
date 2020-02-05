inbam=$1
name=$2
samtools view -b -F 4 ${inbam} > ${name}.mapped.bam
#samtools view -b -F 4 ${name}.Aligned.sorted.bam > ${name}.mapped.bam
samtools sort ${name}.mapped.bam -o ${name}.mapped.sorted.bam
samtools index ${name}.mapped.sorted.bam
bamCoverage --bam ${name}.mapped.sorted.bam -o ${name}.mapped.sorted.bw
