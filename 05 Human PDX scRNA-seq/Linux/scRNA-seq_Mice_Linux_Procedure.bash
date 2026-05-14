
####在小鼠的参考基因组里插入外源插入基因所对应的fasta信息和GTF里放入外源插入基因信息后上传到对应的文件夹里

####
module load anaconda3 
conda info -e
conda activate celescope-1.10.0

####构建参考基因组celescope rna mkref
cd /storage/maxianjueLab/guoyifan/species_reference/Mus_musculus_insertion_Kongdu
cd /storage/maxianjueLab/guoyifan/species_reference/Mus_musculus_insertion_Kongdu

#!/bin/bash
#SBATCH -J Ref
#SBATCH -p amd-ep2,intel-sc3,amd-ep2-short
#SBATCH -q normal
#SBATCH --mem-per-cpu=20480
#SBATCH -c 12
celescope rna mkref \
--thread 12 --genome_name Mus_musculus_GFP_mCherry_insertion \
--fasta /storage/maxianjueLab/guoyifan/species_reference/Mus_musculus_insertion_Kongdu/Mus_musculus.GRCm39.dna.toplevel.fa \
--gtf /storage/maxianjueLab/guoyifan/species_reference/Mus_musculus_insertion_Kongdu/Mus_musculus.GRCm39.110.gtf

####构建my.mapfile文件
CA1	/storage/maxianjueLab/guoyifan/scRNA_seq/Mice/Kongdu/Rawdata	CA1
CA2	/storage/maxianjueLab/guoyifan/scRNA_seq/Mice/Kongdu/Rawdata	CA2
CA3 /storage/maxianjueLab/guoyifan/scRNA_seq/Mice/Kongdu/Rawdata	CA3
NT1	/storage/maxianjueLab/guoyifan/scRNA_seq/Mice/Kongdu/Rawdata	NT1
NT2	/storage/maxianjueLab/guoyifan/scRNA_seq/Mice/Kongdu/Rawdata	NT2
NT3 /storage/maxianjueLab/guoyifan/scRNA_seq/Mice/Kongdu/Rawdata	NT3

####用 multi_rna 构建 celescope rna 分析的 shell 脚本
cd /storage/maxianjueLab/guoyifan/scRNA_seq/Mice/Kongdu

multi_rna \
--mapfile ./my.mapfile \
--genomeDir /storage/maxianjueLab/guoyifan/species_reference/Homo_sapiens_insertion_Kongdu \
--thread 12 \
--mod  shell \
--outdir ./

multi_rna \
--mapfile ./my.mapfile \
--genomeDir /storage/maxianjueLab/guoyifan/species_reference/Homo_sapiens_insertion_Kongdu_3X \
--thread 12 \
--mod  shell \
--outdir ./

multi_rna \
--mapfile ./my.mapfile \
--genomeDir /storage/maxianjueLab/guoyifan/species_reference/Mus_musculus_insertion_Kongdu \
--thread 12 \
--mod  shell \
--outdir ./

multi_rna \
--mapfile ./my.mapfile \
--genomeDir /storage/maxianjueLab/guoyifan/species_reference/Mus_musculus_insertion_Kongdu_3X \
--thread 12 \
--mod  shell \
--outdir ./


celescope rna sample --outdir .//NT3/00.sample --sample NT3 --thread 12 --chemistry auto  --fq1 /storage/maxianjueLab/guoyifan/scRNA_seq/Mice/Kongdu/Rawdata/NT3_R1.fastq.gz 
celescope rna barcode --outdir .//NT3/01.barcode --sample NT3 --thread 12 --chemistry auto --lowNum 2  --fq1 /storage/maxianjueLab/guoyifan/scRNA_seq/Mice/Kongdu/Rawdata/NT3_R1.fastq.gz --fq2 /storage/maxianjueLab/guoyifan/scRNA_seq/Mice/Kongdu/Rawdata/NT3_R2.fastq.gz 
celescope rna cutadapt --outdir .//NT3/02.cutadapt --sample NT3 --thread 12 --minimum_length 20 --nextseq_trim 20 --overlap 10 --insert 150  --fq .//NT3/01.barcode/NT3_2.fq 
celescope rna star --outdir .//NT3/03.star --sample NT3 --thread 12 --genomeDir /storage/maxianjueLab/guoyifan/species_reference/Homo_sapiens_insertion_Kongdu --outFilterMultimapNmax 1 --starMem 30  --fq .//NT3/02.cutadapt/NT3_clean_2.fq 
celescope rna featureCounts --outdir .//NT3/04.featureCounts --sample NT3 --thread 12 --gtf_type exon --genomeDir /storage/maxianjueLab/guoyifan/species_reference/Homo_sapiens_insertion_Kongdu  --input .//NT3/03.star/NT3_Aligned.sortedByCoord.out.bam 
celescope rna count --outdir .//NT3/05.count --sample NT3 --thread 12 --genomeDir /storage/maxianjueLab/guoyifan/species_reference/Homo_sapiens_insertion_Kongdu --expected_cell_num 3000 --cell_calling_method EmptyDrops_CR  --bam .//NT3/04.featureCounts/NT3_name_sorted.bam --force_cell_num None 
celescope rna analysis --outdir .//NT3/06.analysis --sample NT3 --thread 12 --genomeDir /storage/maxianjueLab/guoyifan/species_reference/Homo_sapiens_insertion_Kongdu  --matrix_file .//NT3/05.count/NT3_filtered_feature_bc_matrix 


####构建bash脚本
#由于上面线程是12，所以申请slot的时候可以设置成12

#!/bin/bash
#SBATCH -J NT1
#SBATCH -p amd-ep2,intel-sc3,amd-ep2-short
#SBATCH -q normal
#SBATCH --mem=300G
#SBATCH -c 12
sh /storage/maxianjueLab/guoyifan/scRNA_seq/Mice/Kongdu/shell/NT1.sh

####提交作业



