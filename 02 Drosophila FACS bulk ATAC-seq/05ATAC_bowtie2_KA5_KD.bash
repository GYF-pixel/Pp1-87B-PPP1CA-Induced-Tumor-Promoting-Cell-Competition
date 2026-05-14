#!/bin/bash
#SBATCH -J bwt_KD
#SBATCH -p amd-ep2,intel-sc3,amd-ep2-short
#SBATCH -q normal
#SBATCH --mem=60G
#SBATCH -c 12
projPath="/storage/maxianjueLab/guoyifan/bulkKongDuATACseq"
ref="/storage/maxianjueLab/guoyifan/species_reference/DroMelanogaster/Drosophila_melanogaster.BDGP6.32"
cores=12
for histName in KA5;do
  {
		bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 0 -X 2000 --no-unal -p ${cores} -x ${ref} \
		-1 ${projPath}/03Cutadapt/${histName}_R1_cutadapt.fq.gz \
		-2 ${projPath}/03Cutadapt/${histName}_R2_cutadapt.fq.gz \
		-S ${projPath}/05Alignment/sam/${histName}_bowtie2.sam \
		> ${projPath}/05Alignment/sam/bowtie2_summary/${histName}_bowtie2_infor.log 2>&1 
  }
done
