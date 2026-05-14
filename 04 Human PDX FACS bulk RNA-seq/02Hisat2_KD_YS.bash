#!/bin/bash
#SBATCH -J Hisat2
#SBATCH -p amd-ep2,intel-sc3,amd-ep2-short
#SBATCH -q normal
#SBATCH --mem=120G
#SBATCH -c 12
for i in 8988t-NT1-GFP 8988t-NT1-mCherry 8988t-NT3-GFP 8988t-NT4-GFP 8988t-NT4-mCherry 8988t-NT5-mCherry 8988t-PA1-GFP 8988t-PA1-mCherry 8988t-PA2-mCherry 8988t-PA3-GFP 8988t-PA3-mCherry M1-GFP M3-GFP MF N4-2-GFP N5-2-GFP ND; do
  {
hisat2 -p 12 \
-x /storage/maxianjueLab/guoyifan/species_reference/Homo_sapiens/Homo_sapiens.GRCh38.110.dna.toplevel.ref --rna-strandness RF \
-1 /storage/maxianjueLab/guoyifan/Mice_RNAseq/01Rawdata/${i}_R1.fq.gz \
-2 /storage/maxianjueLab/guoyifan/Mice_RNAseq/01Rawdata/${i}_R2.fq.gz \
-S /storage/maxianjueLab/guoyifan/Mice_RNAseq/03Mapping/hisat2/${i}.sam 2>${i}.summary 
  }
done
