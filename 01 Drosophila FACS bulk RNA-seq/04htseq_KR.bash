#!/bin/bash
#SBATCH -J htsKR
#SBATCH -p amd-ep2
#SBATCH -q normal
#SBATCH --mem=20480
#SBATCH -c 6
for i in KR1 KR2 KR3 KR4 KR5 KR6 KR7; do
 {
htseq-count -f bam -r pos \
--max-reads-in-buffer 50000000 \
--stranded no \
--minaqual 10 \
--type exon --idattr gene_id \
--mode union --nonunique none --secondary-alignments ignore --supplementary-alignments ignore \
--counts_output /storage/maxianjueLab/guoyifan/bulkKongDuRNAseq/04Expression/${i}_gene.tsv \
--nprocesses 6 \
/storage/maxianjueLab/guoyifan/bulkKongDuRNAseq/03Mapping/hisat2/${i}.bam \
/storage/maxianjueLab/guoyifan/RNAseq/ref/Drosophila_melanogaster.BDGP6.32.109.gtf > /storage/maxianjueLab/guoyifan/bulkKongDuRNAseq/04Expression/${i}_count.HTSeq.log  2>&1
}
done
