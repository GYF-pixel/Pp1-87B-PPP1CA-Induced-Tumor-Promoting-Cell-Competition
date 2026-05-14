#!/bin/bash
#SBATCH -J QC_RNA
#SBATCH -p amd-ep2,intel-sc3,amd-ep2-short
#SBATCH -q normal
#SBATCH --mem=20480
#SBATCH -c 8
fastqc -o /storage/maxianjueLab/guoyifan/Mice_RNAseq/02FastQC /storage/maxianjueLab/guoyifan/Mice_RNAseq/01Rawdata/*.fq.gz