#!/bin/bash
#SBATCH -J Ref
#SBATCH -p amd-ep2,intel-sc3,amd-ep2-short
#SBATCH -q normal
#SBATCH --mem-per-cpu=20480
#SBATCH -c 12
celescope rna mkref \
--thread 12 --genome_name Homo_sapiens_GFP_mCherry_insertion \
--fasta /storage/maxianjueLab/guoyifan/species_reference/Homo_sapiens_insertion_Kongdu/Homo_sapiens.GRCh38.dna.toplevel.fa \
--gtf /storage/maxianjueLab/guoyifan/species_reference/Homo_sapiens_insertion_Kongdu/Homo_sapiens.GRCh38.110.gtf