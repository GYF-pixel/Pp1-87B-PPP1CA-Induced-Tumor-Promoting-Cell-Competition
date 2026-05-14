#!/bin/bash
#SBATCH -J Hst2Ref
#SBATCH -p amd-ep2,intel-sc3,amd-ep2-short
#SBATCH -q normal
#SBATCH --mem=100G
#SBATCH -c 12
hisat2-build -p 12 \
--exon /storage/maxianjueLab/guoyifan/species_reference/Homo_sapiens/Homo_sapiens.GRCh38.110.exon \
--ss /storage/maxianjueLab/guoyifan/species_reference/Homo_sapiens/Homo_sapiens.GRCh38.110.ss \
/storage/maxianjueLab/guoyifan/species_reference/Homo_sapiens/Homo_sapiens.GRCh38.dna.toplevel.fa \
/storage/maxianjueLab/guoyifan/species_reference/Homo_sapiens/Homo_sapiens.GRCh38.110.dna.toplevel.ref > hisat2_build.log 2>&1