#!/bin/bash
#SBATCH -J SamS_KD
#SBATCH -p amd-ep2,intel-sc3,amd-ep2-short
#SBATCH -q normal
#SBATCH --mem=51200
#SBATCH -c 8
projPath="/storage/maxianjueLab/guoyifan/bulkKongDuATACseq"
cores=8
for histName in A1 A2 A3 A4 A5 KA1 KA2 KA3 KA4 KA5;do
  {
	 samtools sort -@ 8 ${projPath}/05Alignment/bam/${histName}_bowtie2_filter.bam -o ${projPath}/05Alignment/bam/${histName}_bowtie2_filter_sorted.bam
  }
done