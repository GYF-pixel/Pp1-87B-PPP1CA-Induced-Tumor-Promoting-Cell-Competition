#!/bin/bash
#SBATCH -J Leng_KD
#SBATCH -p amd-ep2,intel-sc3,amd-ep2-short
#SBATCH -q normal
#SBATCH --mem=51200
#SBATCH -c 8
projPath="/storage/maxianjueLab/guoyifan/bulkKongDuATACseq"
for histName in A1 A2 A3 A4 A5 KA1 KA2 KA3 KA4 KA5;do
  {
	picard CollectInsertSizeMetrics \
	I=${projPath}/05Alignment/bam/Rmdup/${histName}_bowtie2_filter_sorted_rmdup_rmmito.bam \
	O=$projPath/05Alignment/Final_QC/Fragment_Length_Summary/${histName}_bowtie2_filter_sorted_rmdup.InsertSize.txt \
	H=$projPath/05Alignment/Final_QC/Fragment_Length_Summary/${histName}_bowtie2_filter_sorted_rmdup.InsertSize.pdf \
	METRIC_ACCUMULATION_LEVEL=ALL_READS >$projPath/05Alignment/Final_QC/Fragment_Length_Summary/${histName}_bowtie2_filter_sorted_rmdup.InsertSize.log 2>&1
  }
done
