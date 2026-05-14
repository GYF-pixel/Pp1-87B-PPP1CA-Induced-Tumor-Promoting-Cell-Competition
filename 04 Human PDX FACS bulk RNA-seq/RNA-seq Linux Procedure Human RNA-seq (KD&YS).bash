#RNA-seq Linux Procedure.bash
squeue -u $(whoami)
###0.Md5 test
md5sum ./*.fq.gz

###01 FastQC
module load fastqc/0.11.9

#!/bin/bash
#SBATCH -J QC_RNA
#SBATCH -p amd-ep2,intel-sc3,amd-ep2-short
#SBATCH -q normal
#SBATCH --mem=20480
#SBATCH -c 8
fastqc -o /storage/maxianjueLab/guoyifan/Mice_RNAseq/02FastQC /storage/maxianjueLab/guoyifan/Mice_RNAseq/01Rawdata/*.fq.gz

cd fastqc
multiqc *zip #将质控结果整合

###02 hisat2
source ~/.bashrc
conda activate RNAseq

cd /storage/maxianjueLab/guoyifan/RNAseq/ref

hisat2 -h
hisat2-build -h

###构建Homo spaiens的参考基因组
cd /storage/maxianjueLab/guoyifan/species_reference/Homo_sapiens

# make exon 
hisat2_extract_exons.py Homo_sapiens.GRCh38.110.gtf> Homo_sapiens.GRCh38.110.exon &

# make splice site
hisat2_extract_splice_sites.py Homo_sapiens.GRCh38.110.gtf > Homo_sapiens.GRCh38.110.ss &

# build index

#!/bin/bash
#SBATCH -J Hst2Ref
#SBATCH -p amd-ep2,intel-sc3,amd-ep2-short
#SBATCH -q normal
#SBATCH --mem=100G
#SBATCH -c 12
hisat2-build -p 12 \
--exon /storage/maxianjueLab/guoyifan/species_reference/Homo_sapiens/Homo_sapiens.GRCh38.110.exon \
--ss /storage/maxianjueLab/guoyifan/species_reference/Homo_sapiens/Homo_sapiens.GRCh38.110.ss \
/storage/maxianjueLab/guoyifan/species_reference/Homo_sapiens/Homo_sapiens.GRCh38.dna.toplevel.fa \
/storage/maxianjueLab/guoyifan/species_reference/Homo_sapiens/Homo_sapiens.GRCh38.110.dna.toplevel.ref > hisat2_build.log 2>&1





#hisat2软件，链特异性文库设置
#--rna-strandness RF or FR
#非链特异性文库，不设置此参数

#!/bin/bash
#SBATCH -J Hisat2
#SBATCH -p amd-ep2,intel-sc3,amd-ep2-short
#SBATCH -q normal
#SBATCH --mem=120G
#SBATCH -c 12
for i in 8988t-NT1-GFP 8988t-NT1-mCherry 8988t-NT3-GFP 8988t-NT4-GFP 8988t-NT4-mCherry 8988t-NT5-mCherry 8988t-PA1-GFP 8988t-PA1-mCherry 8988t-PA2-mCherry 8988t-PA3-GFP 8988t-PA3-mCherry M1-GFP M3-GFP MF N4-2-GFP N5-2-GFP ND; do
  {
hisat2 -p 12 \
-x /storage/maxianjueLab/guoyifan/species_reference/Homo_sapiens/Homo_sapiens.GRCh38.110.dna.toplevel.ref --rna-strandness RF \
-1 /storage/maxianjueLab/guoyifan/Mice_RNAseq/01Rawdata/${i}_R1.fq.gz \
-2 /storage/maxianjueLab/guoyifan/Mice_RNAseq/01Rawdata/${i}_R2.fq.gz \
-S /storage/maxianjueLab/guoyifan/Mice_RNAseq/03Mapping/hisat2/${i}.sam 2>${i}.summary 
  }
done


#04 Samtools sort & Samtools index
module load samtools/1.11

#!/bin/bash
#SBATCH -J Samtools
#SBATCH -p amd-ep2,intel-sc3,amd-ep2-short
#SBATCH -q normal
#SBATCH --mem=120G
#SBATCH -c 12
for i in 8988t-NT1-GFP 8988t-NT1-mCherry 8988t-NT3-GFP 8988t-NT4-GFP 8988t-NT4-mCherry 8988t-NT5-mCherry 8988t-PA1-GFP 8988t-PA1-mCherry 8988t-PA2-mCherry 8988t-PA3-GFP 8988t-PA3-mCherry M1-GFP M3-GFP MF N4-2-GFP N5-2-GFP ND; do
  {
samtools sort --threads 12 -m 10G -o /storage/maxianjueLab/guoyifan/Mice_RNAseq/03Mapping/hisat2/${i}.bam /storage/maxianjueLab/guoyifan/Mice_RNAseq/03Mapping/hisat2/${i}.sam
samtools index /storage/maxianjueLab/guoyifan/Mice_RNAseq/03Mapping/hisat2/${i}.bam
  }
done


## 05 htseq

source ~/.bashrc
conda activate RNAseq

#conda install htseq
htseq-count -h

# htseq
htseq-count -f bam -r pos \
--max-reads-in-buffer 1000000 \
--stranded no \
--minaqual 10 \
--type exon --idattr gene_id \
--mode union --nonunique none --secondary-alignments ignore --supplementary-alignments ignore \
--counts_output ./count_result/test_count.tsv \
--nprocesses 1 \
./bam/test_hisat2.sort.bam  ./reference/gtf/hg38_refseq_from_ucsc.rm_XM_XR.fix_name.gtf > ./count_result/test_count.HTSeq.log  2>&1 & 

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



# feature count
featureCounts -t exon -g gene_id \
-Q 10 --primary -s 0 -p -T 1 \
-a ./reference/gtf/hg38_refseq_from_ucsc.rm_XM_XR.fix_name.gtf \
-o ./count_result/test_count.featureCounts \
./bam/test_hisat2.sort.bam \
./bam/test_hisat2.sort.2.bam > ./count_result/test_count.featureCounts.log  2>&1 & 











