#!/bin/bash
#SBATCH -J Pic_KD
#SBATCH -p amd-ep2,intel-sc3,amd-ep2-short
#SBATCH -q normal
#SBATCH --mem=51200
#SBATCH -c 8
projPath="/storage/maxianjueLab/guoyifan/bulkKongDuATACseq"
for histName in A1 A2 A3 A4 A5 KA1 KA2 KA3 KA4 KA5;do
  {
	 picard MarkDuplicates I=${projPath}/05Alignment/bam/${histName}_bowtie2_filter_sorted.bam \
	 O=${projPath}/05Alignment/bam/Rmdup/${histName}_bowtie2_filter_sorted_rmdup.bam \
	 M=${projPath}/05Alignment/bam/Rmdup/${histName}_bowtie2_filter_sorted_rmdup.matrix \
	 ASO=coordinate REMOVE_DUPLICATES=true >${projPath}/05Alignment/bam/Rmdup/${histName}_bowtie2_filter_sorted_rmdup.log 2>&1
  }
done