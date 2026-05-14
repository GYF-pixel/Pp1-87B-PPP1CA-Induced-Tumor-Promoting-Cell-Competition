#!/bin/bash
#SBATCH -J Samtools
#SBATCH -p amd-ep2,intel-sc3,amd-ep2-short
#SBATCH -q normal
#SBATCH --mem=120G
#SBATCH -c 12
for i in 8988t-NT1-GFP 8988t-NT1-mCherry 8988t-NT3-GFP 8988t-NT4-GFP 8988t-NT4-mCherry 8988t-NT5-mCherry 8988t-PA1-GFP 8988t-PA1-mCherry 8988t-PA2-mCherry 8988t-PA3-GFP 8988t-PA3-mCherry M1-GFP M3-GFP MF N4-2-GFP N5-2-GFP ND; do
  {
samtools sort --threads 12 -m 10G -o /storage/maxianjueLab/guoyifan/Mice_RNAseq/03Mapping/hisat2/${i}.bam /storage/maxianjueLab/guoyifan/Mice_RNAseq/03Mapping/hisat2/${i}.sam
  }
done

for i in 8988t-NT1-GFP 8988t-NT1-mCherry 8988t-NT3-GFP 8988t-NT4-GFP 8988t-NT4-mCherry 8988t-NT5-mCherry 8988t-PA1-GFP 8988t-PA1-mCherry 8988t-PA2-mCherry 8988t-PA3-GFP 8988t-PA3-mCherry M1-GFP M3-GFP MF N4-2-GFP N5-2-GFP ND; do
  {
samtools index -b -@ 12 /storage/maxianjueLab/guoyifan/Mice_RNAseq/03Mapping/hisat2/${i}.bam
  }
done
