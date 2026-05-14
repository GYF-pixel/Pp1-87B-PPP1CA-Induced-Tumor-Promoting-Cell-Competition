#!/bin/bash
#SBATCH -J Samt_KR
#SBATCH -p amd-ep2
#SBATCH -q normal
#SBATCH --mem=20480
#SBATCH -c 6
for i in KR1 KR2 KR3 KR4 KR5 KR6 KR7; do
  {
samtools sort --threads 6 -m 20G -o /storage/maxianjueLab/guoyifan/bulkKongDuRNAseq/03Mapping/hisat2/${i}.bam /storage/maxianjueLab/guoyifan/bulkKongDuRNAseq/03Mapping/hisat2/${i}.sam
samtools index /storage/maxianjueLab/guoyifan/bulkKongDuRNAseq/03Mapping/hisat2/${i}.bam
  }
done
