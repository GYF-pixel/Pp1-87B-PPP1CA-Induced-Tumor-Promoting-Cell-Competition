#!/bin/bash
#SBATCH -J Hisat2_KDR
#SBATCH -p amd-ep2
#SBATCH -q normal
#SBATCH --mem=51200
#SBATCH -c 8
for i in KR1 KR2 KR3 KR4 KR5 KR6 KR7; do
  {
hisat2 -p 8 \
-x /storage/maxianjueLab/guoyifan/RNAseq/ref/Drosophila_melanogaster.BDGP6.32.dna.toplevel.ref --rna-strandness RF \
-1 /storage/maxianjueLab/guoyifan/bulkKongDuRNAseq/01Rawdata/${i}_R1.fq.gz \
-2 /storage/maxianjueLab/guoyifan/bulkKongDuRNAseq/01Rawdata/${i}_R2.fq.gz \
-S /storage/maxianjueLab/guoyifan/bulkKongDuRNAseq/03Mapping/hisat2/${i}.sam 2>${i}.summary 
  }
done
