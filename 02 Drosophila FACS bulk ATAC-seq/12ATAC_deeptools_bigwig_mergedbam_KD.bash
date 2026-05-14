#!/bin/bash
#SBATCH -J merg_KD
#SBATCH -p amd-ep2,intel-sc3,amd-ep2-short
#SBATCH -q normal
#SBATCH --mem=120G
#SBATCH -c 12
projPath="/storage/maxianjueLab/guoyifan/bulkKongDuATACseq"
 for histName in RasPPP1CA_GFP_Pos RasPPP1CA_GFP_Neg;do
   {
     samtools index $projPath/05Alignment/bam/Rmdup/${histName}_bowtie2_filter_sorted_rmdup_rmmito.bam
	 bamCoverage -b $projPath/05Alignment/bam/Rmdup/${histName}_bowtie2_filter_sorted_rmdup_rmmito.bam -o ${projPath}/06BigWig/${histName}_raw.bw
	 bamCoverage -b $projPath/05Alignment/bam/Rmdup/${histName}_bowtie2_filter_sorted_rmdup_rmmito.bam -o ${projPath}/06BigWig/${histName}_normalized_BPM_20bin.bw --binSize 20 --normalizeUsing BPM	 
   }
done
