#!/bin/bash
#SBATCH -J ctad_ATAC
#SBATCH -p amd-ep2
#SBATCH -q normal
#SBATCH --mem=30720
#SBATCH -c 8
for i in A1 A2 A3 A4 A5 KA1 KA2 KA3 KA4 KA5;do
  {
cutadapt -j 8 --times 1 -e 0.1 -O 3 --quality-cutoff 25 -m 25 \
-a CTGTCTCTTATACACATC \
-A CTGTCTCTTATACACATC \
-o /storage/maxianjueLab/guoyifan/bulkKongDuATACseq/03Cutadapt/${i}_R1_cutadapt.fq.gz \
-p /storage/maxianjueLab/guoyifan/bulkKongDuATACseq/03Cutadapt/${i}_R2_cutadapt.fq.gz \
/storage/maxianjueLab/guoyifan/bulkKongDuATACseq/01Rawdata/${i}_R1.fq.gz \
/storage/maxianjueLab/guoyifan/bulkKongDuATACseq/01Rawdata/${i}_R2.fq.gz > /storage/maxianjueLab/guoyifan/bulkKongDuATACseq/03Cutadapt/${i}_cutadapt_infor.log 2>&1
  }
done

