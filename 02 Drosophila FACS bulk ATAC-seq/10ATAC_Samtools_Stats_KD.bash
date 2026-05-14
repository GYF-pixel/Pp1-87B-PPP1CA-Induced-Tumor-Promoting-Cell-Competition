#!/bin/bash
#SBATCH -J Stats_KD
#SBATCH -p amd-ep2,intel-sc3,amd-ep2-short
#SBATCH -q normal
#SBATCH --mem=51200
#SBATCH -c 8
projPath="/storage/maxianjueLab/guoyifan/bulkKongDuATACseq"
for histName in A1 A2 A3 A4 A5 KA1 KA2 KA3 KA4 KA5;do
  {
     samtools stats -@ 8 ${projPath}/05Alignment/bam/Rmdup/${histName}_bowtie2_filter_sorted_rmdup_rmmito.bam > ${projPath}/05Alignment/Final_QC/Mapping_Stats/${histName}_bowtie2_filter_sorted_rmdup.mapping_stats
  }
done
