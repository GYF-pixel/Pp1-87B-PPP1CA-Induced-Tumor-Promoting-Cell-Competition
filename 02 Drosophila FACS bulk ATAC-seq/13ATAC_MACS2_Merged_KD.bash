#!/bin/bash
#SBATCH -J Macs_KD
#SBATCH -p amd-ep2,intel-sc3,amd-ep2-short
#SBATCH -q normal
#SBATCH --mem=80G
#SBATCH -c 8
projPath="/storage/maxianjueLab/guoyifan/bulkKongDuATACseq"
 for histName in RasPPP1CA_GFP_Pos RasPPP1CA_GFP_Neg;do
  {
macs2 callpeak -t ${projPath}/05Alignment/bam/Rmdup/${histName}_bowtie2_filter_sorted_rmdup_rmmito.bam \
      -g dm -f BAMPE --nomodel --shift -100 --extsize 200 -n macs2_${histName}_peak_q0.05 --outdir $projPath/07MACS2 -q 0.05 --keep-dup all 2>${projPath}/07MACS2/macs2Peak_${histName}_summary.txt
    }
done
