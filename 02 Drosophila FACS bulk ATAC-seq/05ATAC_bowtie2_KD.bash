#!/bin/bash
#SBATCH -J bwt_KD
#SBATCH -p amd-ep2,intel-sc3,amd-ep2-short
#SBATCH -q normal
#SBATCH --mem=160G
#SBATCH -c 16
projPath="/storage/maxianjueLab/guoyifan/bulkKongDuATACseq"
ref="/storage/maxianjueLab/guoyifan/species_reference/DroMelanogaster/Drosophila_melanogaster.BDGP6.32"
cores=16
for histName in A1 A2 A3 A4 A5 KA1 KA2 KA3 KA4 KA5;do
  {
		bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 0 -X 2000 --no-unal -p ${cores} -x ${ref} \
		-1 ${projPath}/03Cutadapt/${histName}_R1_cutadapt.fq.gz \
		-2 ${projPath}/03Cutadapt/${histName}_R2_cutadapt.fq.gz \
		-S ${projPath}/05Alignment/sam/${histName}_bowtie2.sam \
		> ${projPath}/05Alignment/sam/bowtie2_summary/${histName}_bowtie2_infor.log 2>&1 
  }
done
