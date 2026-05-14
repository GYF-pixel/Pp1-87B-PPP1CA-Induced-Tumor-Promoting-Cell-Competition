#!/bin/bash
#SBATCH -J SamIn_KD
#SBATCH -p amd-ep2,intel-sc3,amd-ep2-short
#SBATCH -q normal
#SBATCH --mem=60G
#SBATCH -c 12
projPath="/storage/maxianjueLab/guoyifan/bulkKongDuATACseq"
for histName in A1 A2 A3 A4 A5 KA1 KA2 KA3 KA4 KA5;do
  {
    samtools view -h -@ 12 ${projPath}/05Alignment/bam/Rmdup/${histName}_bowtie2_filter_sorted_rmdup.bam | grep -v 'mitochondrion_genome' | samtools view -b -@ 12 -o ${projPath}/05Alignment/bam/Rmdup/${histName}_bowtie2_filter_sorted_rmdup_rmmito.bam 
	samtools index -@ 12 ${projPath}/05Alignment/bam/Rmdup/${histName}_bowtie2_filter_sorted_rmdup_rmmito.bam
  }
done
