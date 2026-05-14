#!/bin/bash
#SBATCH -J QC_KDR
#SBATCH -p amd-ep2
#SBATCH -q normal
#SBATCH --mem=51200
#SBATCH -c 8

for i in KR1 KR2 KR3 KR4 KR5 KR6 KR7; do
  {
./run_fastqc.sh -i ./raw_data -o ./qc_results -t 8
  }
done

multiqc ./qc_results -o ./qc_results/multiqc_report