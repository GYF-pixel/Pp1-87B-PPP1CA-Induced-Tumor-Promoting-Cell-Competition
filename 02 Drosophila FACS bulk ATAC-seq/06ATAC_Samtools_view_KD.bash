#!/bin/bash
#SBATCH -J SamV_KD
#SBATCH -p amd-ep2,intel-sc3,amd-ep2-short
#SBATCH -q normal
#SBATCH --mem=51200
#SBATCH -c 8
projPath="/storage/maxianjueLab/guoyifan/bulkKongDuATACseq"
cores=8
for histName in A1 A2 A3 A4 A5 KA1 KA2 KA3 KA4 KA5;do
  {
	 samtools view -h -b -@ 8 -f 3 -F 12 -F 256 -q 20 ${projPath}/05Alignment/sam/${histName}_bowtie2.sam > ${projPath}/05Alignment/bam/${histName}_bowtie2_filter.bam
  }
done