#!/bin/bash
#SBATCH -J htsKDYS
#SBATCH -p amd-ep2,intel-sc3,amd-ep2-short
#SBATCH -q normal
#SBATCH --mem=50G
#SBATCH -c 12
for i in 8988t-NT1-GFP 8988t-NT1-mCherry 8988t-NT3-GFP 8988t-NT4-GFP 8988t-NT4-mCherry 8988t-NT5-mCherry 8988t-PA1-GFP 8988t-PA1-mCherry 8988t-PA2-mCherry 8988t-PA3-GFP 8988t-PA3-mCherry M1-GFP M3-GFP MF N4-2-GFP N5-2-GFP ND; do
  {
htseq-count -f bam -r pos \
--max-reads-in-buffer 50000000 \
--stranded no \
--minaqual 10 \
--type exon --idattr gene_id \
--mode union --nonunique none --secondary-alignments ignore --supplementary-alignments ignore \
--counts_output /storage/maxianjueLab/guoyifan/Mice_RNAseq/04Expression/${i}_gene.tsv \
--nprocesses 12 \
/storage/maxianjueLab/guoyifan/Mice_RNAseq/03Mapping/hisat2/${i}.bam \
/storage/maxianjueLab/guoyifan/species_reference/Homo_sapiens/Homo_sapiens.GRCh38.110.gtf > /storage/maxianjueLab/guoyifan/Mice_RNAseq/04Expression/${i}_count.HTSeq.log  2>&1
}
done
