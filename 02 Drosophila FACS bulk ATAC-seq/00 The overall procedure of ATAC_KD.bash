#0.Md5 test
md5sum ./*.fq.gz

projPath="/storage/maxianjueLab/guoyifan/ATAC_seq"
#1.FastQC
mkdir fastqc#创造一个文件夹存放质控结果
module load fastqc/0.11.9

fastqc -o /storage/maxianjueLab/guoyifan/bulkKongDuATACseq/02FastQC /storage/maxianjueLab/guoyifan/bulkKongDuATACseq/01Rawdata/*.fq.gz

cd fastqc
multiqc *zip #将质控结果整合

#2.Cutadapt
#cutadapt 在 CUTTag conda env里
source ~/.bashrc
conda activate CUTTag
cutadapt -h

#!/bin/bash
#SBATCH -J ctad_ATAC
#SBATCH -p amd-ep2,intel-sc3,amd-ep2-short
#SBATCH -q normal
#SBATCH --mem=20480
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
  }&
done


#3. FastQC_after_filtrate
mkdir 04FastQC_after_filtrate
module load fastqc/0.11.9

fastqc -o /storage/maxianjueLab/guoyifan/bulkKongDuATACseq/04FastQC_after_filtrate /storage/maxianjueLab/guoyifan/bulkKongDuATACseq/03Cutadapt/*_cutadapt.fq.gz

multiqc *zip  #将质控结果整合

#4. Bowtie2 Alignment
#https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
mkdir 05Alignment
mkdir 05Alignment/sam
mkdir 05Alignment/sam/bowtie2_summary
module load bowtie/2.4.2

#ATAC-seq 比对前一定要关注去除掉adaptor之后FastQC质检结果，Median Read Length，通常都是在80bp左右，最小的可以是25bp
#因此用bowtie2进行比对的时候最小的片段要为0，否则unique mapping rate会非常的低
#Bowtie2 全参数https://blog.csdn.net/herokoking/article/details/77847384

#!/bin/bash
#SBATCH -J bwt_KD
#SBATCH -p amd-ep2,intel-sc3,amd-ep2-short
#SBATCH -q normal
#SBATCH --mem=160G
#SBATCH -c 16
projPath="/storage/maxianjueLab/guoyifan/bulkKongDuATACseq"
ref="/storage/maxianjueLab/guoyifan/species_reference/DroMelanogaster/Drosophila_melanogaster.BDGP6.32"
cores=16
for histName in A1 A2 A3 A4 A5 KA1 KA2 KA3 KA4 KA5;do
  {
		bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 0 -X 2000 --no-unal -p ${cores} -x ${ref} \
		-1 ${projPath}/03Cutadapt/${histName}_R1_cutadapt.fq.gz \
		-2 ${projPath}/03Cutadapt/${histName}_R2_cutadapt.fq.gz \
		-S ${projPath}/05Alignment/sam/${histName}_bowtie2.sam \
		> ${projPath}/05Alignment/sam/bowtie2_summary/${histName}_bowtie2_infor.log 2>&1 
  }
done


#5. Samtools view
mkdir /storage/maxianjueLab/guoyifan/bulkKongDuATACseq/05Alignment/bam

module load samtools/1.11

#保留MAPQ>=20的

#!/bin/bash
#SBATCH -J SamV_KD
#SBATCH -p amd-ep2,intel-sc3,amd-ep2-short
#SBATCH -q normal
#SBATCH --mem=51200
#SBATCH -c 8
projPath="/storage/maxianjueLab/guoyifan/bulkKongDuATACseq"
cores=8
for histName in A1 A2 A3 A4 A5 KA1 KA2 KA3 KA4 KA5;do
  {
	 samtools view -h -b -@ 8 -f 3 -F 12 -F 256 -q 20 ${projPath}/05Alignment/sam/${histName}_bowtie2.sam > ${projPath}/05Alignment/bam/${histName}_bowtie2_filter.bam
  }
done


#6. Samtools sort
module load samtools/1.11

#!/bin/bash
#SBATCH -J SamS_KD
#SBATCH -p amd-ep2,intel-sc3,amd-ep2-short
#SBATCH -q normal
#SBATCH --mem=51200
#SBATCH -c 8
projPath="/storage/maxianjueLab/guoyifan/bulkKongDuATACseq"
cores=8
for histName in A1 A2 A3 A4 A5 KA1 KA2 KA3 KA4 KA5;do
  {
	 samtools sort -@ 8 ${projPath}/05Alignment/bam/${histName}_bowtie2_filter.bam -o ${projPath}/05Alignment/bam/${histName}_bowtie2_filter_sorted.bam
  }
done


#7. Picard remove Duplicates
module load R/4.2.1   #绘制Insetsize
module load picard/2.25.1
picard -h                   #查看可用工具
picard MarkDuplicates -h    #查看MarkDuplicates使用方法

mkdir 05Alignment/bam/Rmdup

#!/bin/bash
#SBATCH -J Pic_KD
#SBATCH -p amd-ep2,intel-sc3,amd-ep2-short
#SBATCH -q normal
#SBATCH --mem=51200
#SBATCH -c 8
projPath="/storage/maxianjueLab/guoyifan/bulkKongDuATACseq"
for histName in A1 A2 A3 A4 A5 KA1 KA2 KA3 KA4 KA5;do
  {
	 picard MarkDuplicates I=${projPath}/05Alignment/bam/${histName}_bowtie2_filter_sorted.bam \
	 O=${projPath}/05Alignment/bam/Rmdup/${histName}_bowtie2_filter_sorted_rmdup.bam \
	 M=${projPath}/05Alignment/bam/Rmdup/${histName}_bowtie2_filter_sorted_rmdup.matrix \
	 ASO=coordinate REMOVE_DUPLICATES=true >${projPath}/05Alignment/bam/Rmdup/${histName}_bowtie2_filter_sorted_rmdup.log 2>&1
  }
done

#8. Samtools Remove Mitocondrial Reads & Bam File Index
#过滤线粒体reads并添加索引

module load samtools/1.11
#!/bin/bash
#SBATCH -J SamIn_KD
#SBATCH -p amd-ep2,intel-sc3,amd-ep2-short
#SBATCH -q normal
#SBATCH --mem=120G
#SBATCH -c 12
projPath="/storage/maxianjueLab/guoyifan/bulkKongDuATACseq"
for histName in A1 A2 A3 A4 A5 KA1 KA2 KA3 KA4 KA5;do
  {
    samtools view -h -@ 12 ${projPath}/05Alignment/bam/Rmdup/${histName}_bowtie2_filter_sorted_rmdup.bam | grep -v 'mitochondrion_genome' | samtools view -b -@ 12 -o ${projPath}/05Alignment/bam/Rmdup/${histName}_bowtie2_filter_sorted_rmdup_rmmito.bam 
	samtools index -@ 12 ${projPath}/05Alignment/bam/Rmdup/${histName}_bowtie2_filter_sorted_rmdup_rmmito.bam
  }
done


#9. Final-align stats 
mkdir 05Alignment/Final_QC/Mapping_Stats

#!/bin/bash
#SBATCH -J Stats_KD
#SBATCH -p amd-ep2,intel-sc3,amd-ep2-short
#SBATCH -q normal
#SBATCH --mem=51200
#SBATCH -c 8
projPath="/storage/maxianjueLab/guoyifan/bulkKongDuATACseq"
for histName in A1 A2 A3 A4 A5 KA1 KA2 KA3 KA4 KA5;do
  {
     samtools stats -@ 8 ${projPath}/05Alignment/bam/Rmdup/${histName}_bowtie2_filter_sorted_rmdup_rmmito.bam > ${projPath}/05Alignment/Final_QC/Mapping_Stats/${histName}_bowtie2_filter_sorted_rmdup.mapping_stats
  }
done

#10. Picard InsertSize/Fragment_Length_Summary
mkdir 05Alignment/Final_QC
mkdir 05Alignment/Final_QC/Fragment_Length_Summary


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


#11.  bam to bw

source ~/.bashrc
conda activate deeptools

mkdir 06BigWig

#!/bin/bash
#SBATCH -J bCov_KD
#SBATCH -p amd-ep2,intel-sc3,amd-ep2-short
#SBATCH -q normal
#SBATCH --mem=120G
#SBATCH -c 12
projPath="/storage/maxianjueLab/guoyifan/bulkKongDuATACseq"
for histName in A1 A2 A3 A4 A5 KA1 KA2 KA3 KA4 KA5;do
  {
	 bamCoverage -b ${projPath}/05Alignment/bam/Rmdup/${histName}_bowtie2_filter_sorted_rmdup_rmmito.bam -o ${projPath}/06BigWig/${histName}_raw.bw
	 bamCoverage -b ${projPath}/05Alignment/bam/Rmdup/${histName}_bowtie2_filter_sorted_rmdup_rmmito.bam -o ${projPath}/06BigWig/${histName}_normalized_BPM_20bin.bw --binSize 20 --normalizeUsing BPM
     bamCoverage -b ${projPath}/05Alignment/bam/Rmdup/${histName}_bowtie2_filter_sorted_rmdup_rmmito.bam -o ${projPath}/06BigWig/${histName}_normalized_RPKM_20bin.bw --binSize 20 --normalizeUsing RPKM	 
  }
done


#12. Merge bam files across biological replicates
source ~/.bashrc
conda activate deeptools

cd /storage/maxianjueLab/guoyifan/bulkKongDuATACseq/05Alignment/bam/Rmdup

#Group information 
#RasPPP1CA_GFP_Pos A4 A5 KA4 KA5
#RasPPP1CA_GFP_Neg A1 A2 A3 KA1 KA2 KA3
samtools merge RasPPP1CA_GFP_Neg_bowtie2_filter_sorted_rmdup_rmmito.bam A1_bowtie2_filter_sorted_rmdup_rmmito.bam A2_bowtie2_filter_sorted_rmdup_rmmito.bam A3_bowtie2_filter_sorted_rmdup_rmmito.bam KA1_bowtie2_filter_sorted_rmdup_rmmito.bam KA2_bowtie2_filter_sorted_rmdup_rmmito.bam KA3_bowtie2_filter_sorted_rmdup_rmmito.bam &
samtools merge RasPPP1CA_GFP_Pos_bowtie2_filter_sorted_rmdup_rmmito.bam A4_bowtie2_filter_sorted_rmdup_rmmito.bam A5_bowtie2_filter_sorted_rmdup_rmmito.bam  KA4_bowtie2_filter_sorted_rmdup_rmmito.bam KA5_bowtie2_filter_sorted_rmdup_rmmito.bam &

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
     bamCoverage -b ${projPath}/05Alignment/bam/Rmdup/${histName}_bowtie2_filter_sorted_rmdup_rmmito.bam -o ${projPath}/06BigWig/${histName}_normalized_RPKM_20bin.bw --binSize 20 --normalizeUsing RPKM	 
   }
done


#13. MACS2 Call Peak
mkdir 07MACS2

source ~/.bashrc
conda activate macs2

# No control
ATAC-seq关心的是在哪切断，断点才是peak的中心，所以使用shift模型，–shift -75或-100
对人细胞系ATAC-seq 数据call peak的参数设置如下：
macs2 callpeak -t H1hesc.final.bam -n sample --shift -100 --extsize 200 --nomodel -B --SPMR -g hs --outdir Macs2_out 2> sample.macs2.log

#!/bin/bash
#SBATCH -J Macs_KD
#SBATCH -p amd-ep2,intel-sc3,amd-ep2-short
#SBATCH -q normal
#SBATCH --mem=80G
#SBATCH -c 8
projPath="/storage/maxianjueLab/guoyifan/bulkKongDuATACseq"
for histName in A1 A2 A3 A4 A5 KA1 KA2 KA3 KA4 KA5;do
  {
macs2 callpeak -t ${projPath}/05Alignment/bam/Rmdup/${histName}_bowtie2_filter_sorted_rmdup_rmmito.bam \
      -g dm -f BAMPE --nomodel --shift -100 --extsize 200 -n macs2_${histName}_peak_q0.05 --outdir $projPath/07MACS2 -q 0.05 --keep-dup all 2>${projPath}/07MACS2/macs2Peak_${histName}_summary.txt
    }
done

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



#14. Deeptools visulization

#14.1 Deeptools_BPM
cd /storage/maxianjueLab/guoyifan/bulkKongDuATACseq
mkdir 08Deeptools_BPM

module load R/4.2.1
source ~/.bashrc
conda activate deeptools


####################TSS
#!/bin/bash
#SBATCH -J depTSS_KD
#SBATCH -p amd-ep2,intel-sc3,amd-ep2-short
#SBATCH -q normal
#SBATCH --mem=120G
#SBATCH -c 12
projPath="/storage/maxianjueLab/guoyifan/bulkKongDuATACseq"
ref_gtf="/storage/maxianjueLab/guoyifan/species_reference/Drosophila_melanogaster_BDGP6_32/Drosophila_melanogaster.BDGP6.32.109.gtf"
for histName in RasPPP1CA_GFP_Pos RasPPP1CA_GFP_Neg A1 A2 A3 A4 A5 KA1 KA2 KA3 KA4 KA5;do
{
     computeMatrix reference-point -S ${projPath}/06BigWig/${histName}_normalized_BPM_20bin.bw \
							   -p 12 \
							  --referencePoint TSS \
							  --afterRegionStartLength 3000 \
							  --beforeRegionStartLength 3000 \
							  -R $ref_gtf \
							  --skipZeros  --missingDataAsZero \
							  -o $projPath/08Deeptools_BPM/${histName}_refPoint_TSS_data.gz
     plotHeatmap -m $projPath/08Deeptools_BPM/${histName}_refPoint_TSS_data.gz \
            --missingDataColor 1 \
            --colorList 'white,#925E9F' \
            --heatmapHeight 12 \
			--sortUsing sum --startLabel "TSS" \
            --endLabel "TES" --xAxisLabel "" \
            --regionsLabel "Peaks" \
            --samplesLabel "${histName}" \
            -o $projPath/08Deeptools_BPM/${histName}_refPoint_TSS_heatmap.pdf
    }
done


####################TSS_TES
#!/bin/bash
#SBATCH -J depGB_KD
#SBATCH -p amd-ep2,intel-sc3,amd-ep2-short
#SBATCH -q normal
#SBATCH --mem=120G
#SBATCH -c 12
projPath="/storage/maxianjueLab/guoyifan/bulkKongDuATACseq"
ref_gtf="/storage/maxianjueLab/guoyifan/species_reference/Drosophila_melanogaster_BDGP6_32/Drosophila_melanogaster.BDGP6.32.109.gtf"
for histName in RasPPP1CA_GFP_Pos RasPPP1CA_GFP_Neg A1 A2 A3 A4 A5 KA1 KA2 KA3 KA4 KA5;do
{
     computeMatrix reference-point -S ${projPath}/06BigWig/${histName}_normalized_BPM_20bin.bw \
							   -p 8 \
							  --binSize 20 \
							  --regionBodyLength 5000 \
							  --referencePoint TSS \
							  --afterRegionStartLength 3000 \
							  --beforeRegionStartLength 3000 \
							  -R $ref_gtf \
							  --skipZeros  --missingDataAsZero \
							  -o $projPath/08Deeptools_BPM/${histName}_refPoint_TSS_TES_data.gz
     plotHeatmap -m $projPath/08Deeptools_BPM/${histName}_refPoint_TSS_TES_data.gz \
            --missingDataColor 1 \
            --colorList 'white,#339933' \
            --heatmapHeight 12 \
			--sortUsing sum --startLabel "TSS" \
            --endLabel "TES" --xAxisLabel "" \
            --regionsLabel "Peaks" \
            --samplesLabel "${histName}" \
            -o $projPath/08Deeptools_BPM/${histName}_refPoint_TSS_TES_heatmap.pdf
    }
done


####################TSS Summary (Whole genome input)
cd /storage/maxianjueLab/guoyifan/bulkKongDuATACseq/08Deeptools_BPM
mkdir TFs_of_interested

#####GFP negative
#!/bin/bash
#SBATCH -J depWG_KD
#SBATCH -p amd-ep2,intel-sc3,amd-ep2-short
#SBATCH -q normal
#SBATCH --mem=80G
#SBATCH -c 8
projPath="/storage/maxianjueLab/guoyifan/bulkKongDuATACseq"
ref_gtf="/storage/maxianjueLab/guoyifan/species_reference/Drosophila_melanogaster_BDGP6_32/Drosophila_melanogaster.BDGP6.32.109.gtf"
TFs_of_interested="/storage/maxianjueLab/guoyifan/bulkKongDuATACseq/08Deeptools_BPM/TFs_of_interested"
for histName in WT_Ras_RasPPP1CA_GFP_neg_Summary;do
{
computeMatrix reference-point -S /storage/maxianjueLab/guoyifan/FTH/02Analysis/02ATAC-seq/06BigWig/E_normalized_BPM_20bin.bw \
	                                                     /storage/maxianjueLab/guoyifan/FTH/02Analysis/02ATAC-seq/06BigWig/A_normalized_BPM_20bin.bw \
														 ${projPath}/06BigWig/RasPPP1CA_GFP_Neg_normalized_BPM_20bin.bw \
							   -p 8 \
							  --referencePoint TSS \
							  --afterRegionStartLength 3000 \
							  --beforeRegionStartLength 3000 \
							  -R $ref_gtf \
							  --skipZeros  --missingDataAsZero \
							  -o $projPath/08Deeptools_BPM/${histName}_TSS_whole_genome_data.gz

plotHeatmap -m $projPath/08Deeptools_BPM/${histName}_TSS_whole_genome_data.gz \
            --missingDataColor 1 \
            --heatmapHeight 12 \
			--sortUsing sum --startLabel "TSS" \
            --endLabel "TES" --xAxisLabel "" \
            --regionsLabel "Whole genome" \
            --samplesLabel "WT Peaks"  "Ras Peaks" "RasPPP1CA Peaks" \
            -o $projPath/08Deeptools_BPM/${histName}_TSS_whole_genome_heatmap.pdf
	}
done



#####GFP positive
#!/bin/bash
#SBATCH -J depWG_KD
#SBATCH -p amd-ep2,intel-sc3,amd-ep2-short
#SBATCH -q normal
#SBATCH --mem=80G
#SBATCH -c 8
projPath="/storage/maxianjueLab/guoyifan/bulkKongDuATACseq"
ref_gtf="/storage/maxianjueLab/guoyifan/species_reference/Drosophila_melanogaster_BDGP6_32/Drosophila_melanogaster.BDGP6.32.109.gtf"
TFs_of_interested="/storage/maxianjueLab/guoyifan/bulkKongDuATACseq/08Deeptools_BPM/TFs_of_interested"
for histName in WT_Ras_RasPPP1CA_GFP_pos_Summary;do
{
computeMatrix reference-point -S /storage/maxianjueLab/guoyifan/FTH/02Analysis/02ATAC-seq/06BigWig/F_normalized_BPM_20bin.bw \
	                                                     /storage/maxianjueLab/guoyifan/FTH/02Analysis/02ATAC-seq/06BigWig/B_normalized_BPM_20bin.bw \
														 ${projPath}/06BigWig/RasPPP1CA_GFP_Pos_normalized_BPM_20bin.bw \
							   -p 8 \
							  --referencePoint TSS \
							  --afterRegionStartLength 3000 \
							  --beforeRegionStartLength 3000 \
							  -R $ref_gtf \
							  --skipZeros  --missingDataAsZero \
							  -o $projPath/08Deeptools_BPM/${histName}_TSS_whole_genome_data.gz

plotHeatmap -m $projPath/08Deeptools_BPM/${histName}_TSS_whole_genome_data.gz \
            --missingDataColor 1 \
            --heatmapHeight 12 \
			--sortUsing sum --startLabel "TSS" \
            --endLabel "TES" --xAxisLabel "" \
            --regionsLabel "Whole genome" \
            --samplesLabel "WT Peaks"  "Ras Peaks" "RasPPP1CA Peaks" \
            -o $projPath/08Deeptools_BPM/${histName}_TSS_whole_genome_heatmap.pdf
	}
done


####################TSS Summary (MACS2 bed input)
##########MACS2 Summit.bed

#####GFP positive
#!/bin/bash
#SBATCH -J depPK_KD
#SBATCH -p amd-ep2,intel-sc3,amd-ep2-short
#SBATCH -q normal
#SBATCH --mem=80G
#SBATCH -c 8
projPath="/storage/maxianjueLab/guoyifan/bulkKongDuATACseq"
ref_gtf="/storage/maxianjueLab/guoyifan/species_reference/Drosophila_melanogaster_BDGP6_32/Drosophila_melanogaster.BDGP6.32.109.gtf"
TFs_of_interested="/storage/maxianjueLab/guoyifan/bulkKongDuATACseq/08Deeptools_BPM/TFs_of_interested"
for histName in WT_Ras_RasPPP1CA_GFP_pos_Summary;do
{
computeMatrix reference-point -S /storage/maxianjueLab/guoyifan/FTH/02Analysis/02ATAC-seq/06BigWig/F_normalized_BPM_20bin.bw \
	                                                     /storage/maxianjueLab/guoyifan/FTH/02Analysis/02ATAC-seq/06BigWig/B_normalized_BPM_20bin.bw \
														 ${projPath}/06BigWig/RasPPP1CA_GFP_Pos_normalized_BPM_20bin.bw \
							   -p 8 \
							  --referencePoint TSS \
							  --afterRegionStartLength 3000 \
							  --beforeRegionStartLength 3000 \
							  -R $TFs_of_interested/macs2_WT_GFP_Pos_peak_q0.05_summits.bed \
							  $TFs_of_interested/macs2_Ras_GFP_Pos_peak_q0.05_summits.bed \
							  $TFs_of_interested/macs2_RasPPP1CA_GFP_Pos_peak_q0.05_summits.bed \
							  --skipZeros  --missingDataAsZero \
							  -o $projPath/08Deeptools_BPM/${histName}_TSS_Macs2_Peak_Summit_data.gz

plotHeatmap -m $projPath/08Deeptools_BPM/${histName}_TSS_Macs2_Peak_Summit_data.gz \
            --heatmapHeight 12 \
			--sortUsing sum --startLabel "TSS" \
            --endLabel "TES" --xAxisLabel "" \
            --regionsLabel "WT open-region" "Ras open-region" "RasPPP1CA open-region"\
            --samplesLabel "WT Peaks"  "Ras Peaks" "RasPPP1CA Peaks" \
            -o $projPath/08Deeptools_BPM/${histName}_TSS_Macs2_Peak_Summit_heatmap.pdf
	}
done



#####GFP negative
#!/bin/bash
#SBATCH -J depPK_KD
#SBATCH -p amd-ep2,intel-sc3,amd-ep2-short
#SBATCH -q normal
#SBATCH --mem=80G
#SBATCH -c 8
projPath="/storage/maxianjueLab/guoyifan/bulkKongDuATACseq"
ref_gtf="/storage/maxianjueLab/guoyifan/species_reference/Drosophila_melanogaster_BDGP6_32/Drosophila_melanogaster.BDGP6.32.109.gtf"
TFs_of_interested="/storage/maxianjueLab/guoyifan/bulkKongDuATACseq/08Deeptools_BPM/TFs_of_interested"
for histName in WT_Ras_RasPPP1CA_GFP_neg_Summary;do
{
computeMatrix reference-point -S /storage/maxianjueLab/guoyifan/FTH/02Analysis/02ATAC-seq/06BigWig/E_normalized_BPM_20bin.bw \
	                                                     /storage/maxianjueLab/guoyifan/FTH/02Analysis/02ATAC-seq/06BigWig/A_normalized_BPM_20bin.bw \
														 ${projPath}/06BigWig/RasPPP1CA_GFP_Neg_normalized_BPM_20bin.bw \
							   -p 8 \
							  --referencePoint TSS \
							  --afterRegionStartLength 3000 \
							  --beforeRegionStartLength 3000 \
							  -R $TFs_of_interested/macs2_WT_GFP_Neg_peak_q0.05_summits.bed \
							  $TFs_of_interested/macs2_Ras_GFP_Neg_peak_q0.05_summits.bed \
							  $TFs_of_interested/macs2_RasPPP1CA_GFP_Neg_peak_q0.05_summits.bed \
							  --skipZeros  --missingDataAsZero \
							  -o $projPath/08Deeptools_BPM/${histName}_TSS_Macs2_Peak_Summit_data.gz

plotHeatmap -m $projPath/08Deeptools_BPM/${histName}_TSS_Macs2_Peak_Summit_data.gz \
            --heatmapHeight 12 \
			--sortUsing sum --startLabel "TSS" \
            --endLabel "TES" --xAxisLabel "" \
            --regionsLabel "WT open-region" "Ras open-region" "RasPPP1CA open-region"\
            --samplesLabel "WT Peaks"  "Ras Peaks" "RasPPP1CA Peaks" \
            -o $projPath/08Deeptools_BPM/${histName}_TSS_Macs2_Peak_Summit_heatmap.pdf
	}
done


##########MACS2 narrowPeak to bed

#####GFP positive
#!/bin/bash
#SBATCH -J depnarPK_KD
#SBATCH -p amd-ep2,intel-sc3,amd-ep2-short
#SBATCH -q normal
#SBATCH --mem=80G
#SBATCH -c 8
projPath="/storage/maxianjueLab/guoyifan/bulkKongDuATACseq"
ref_gtf="/storage/maxianjueLab/guoyifan/species_reference/Drosophila_melanogaster_BDGP6_32/Drosophila_melanogaster.BDGP6.32.109.gtf"
TFs_of_interested="/storage/maxianjueLab/guoyifan/bulkKongDuATACseq/08Deeptools_BPM/TFs_of_interested"
for histName in WT_Ras_RasPPP1CA_GFP_pos_Summary;do
{
computeMatrix reference-point -S /storage/maxianjueLab/guoyifan/FTH/02Analysis/02ATAC-seq/06BigWig/F_normalized_BPM_20bin.bw \
	                                                     /storage/maxianjueLab/guoyifan/FTH/02Analysis/02ATAC-seq/06BigWig/B_normalized_BPM_20bin.bw \
														 ${projPath}/06BigWig/RasPPP1CA_GFP_Pos_normalized_BPM_20bin.bw \
							   -p 8 \
							  --referencePoint TSS \
							  --afterRegionStartLength 3000 \
							  --beforeRegionStartLength 3000 \
							  -R $TFs_of_interested/macs2_WT_GFP_Pos_peak_q0.05_peaks_narrowPeak.bed \
							  $TFs_of_interested/macs2_Ras_GFP_Pos_peak_q0.05_peaks_narrowPeak.bed \
							  $TFs_of_interested/macs2_RasPPP1CA_GFP_Pos_peak_q0.05_peaks_narrowPeak.bed \
							  --skipZeros  --missingDataAsZero \
							  -o $projPath/08Deeptools_BPM/${histName}_TSS_Macs2_narrowPeak_data.gz

plotHeatmap -m $projPath/08Deeptools_BPM/${histName}_TSS_Macs2_narrowPeak_data.gz \
            --heatmapHeight 12 \
			--sortUsing sum --startLabel "TSS" \
            --endLabel "TES" --xAxisLabel "" \
            --regionsLabel "WT open-region" "Ras open-region" "RasPPP1CA open-region"\
            --samplesLabel "WT Peaks"  "Ras Peaks" "RasPPP1CA Peaks" \
            -o $projPath/08Deeptools_BPM/${histName}_TSS_Macs2_narrowPeak_heatmap.pdf
	}
done



#####GFP negative
#!/bin/bash
#SBATCH -J depnarPK_KD
#SBATCH -p amd-ep2,intel-sc3,amd-ep2-short
#SBATCH -q normal
#SBATCH --mem=80G
#SBATCH -c 8
projPath="/storage/maxianjueLab/guoyifan/bulkKongDuATACseq"
ref_gtf="/storage/maxianjueLab/guoyifan/species_reference/Drosophila_melanogaster_BDGP6_32/Drosophila_melanogaster.BDGP6.32.109.gtf"
TFs_of_interested="/storage/maxianjueLab/guoyifan/bulkKongDuATACseq/08Deeptools_BPM/TFs_of_interested"
for histName in WT_Ras_RasPPP1CA_GFP_neg_Summary;do
{
computeMatrix reference-point -S /storage/maxianjueLab/guoyifan/FTH/02Analysis/02ATAC-seq/06BigWig/E_normalized_BPM_20bin.bw \
	                                                     /storage/maxianjueLab/guoyifan/FTH/02Analysis/02ATAC-seq/06BigWig/A_normalized_BPM_20bin.bw \
														 ${projPath}/06BigWig/RasPPP1CA_GFP_Neg_normalized_BPM_20bin.bw \
							   -p 8 \
							  --referencePoint TSS \
							  --afterRegionStartLength 3000 \
							  --beforeRegionStartLength 3000 \
							  -R $TFs_of_interested/macs2_WT_GFP_Neg_peak_q0.05_peaks_narrowPeak.bed \
							  $TFs_of_interested/macs2_Ras_GFP_Neg_peak_q0.05_peaks_narrowPeak.bed \
							  $TFs_of_interested/macs2_RasPPP1CA_GFP_Neg_peak_q0.05_peaks_narrowPeak.bed \
							  --skipZeros  --missingDataAsZero \
							  -o $projPath/08Deeptools_BPM/${histName}_TSS_Macs2_narrowPeak_data.gz

plotHeatmap -m $projPath/08Deeptools_BPM/${histName}_TSS_Macs2_narrowPeak_data.gz \
            --heatmapHeight 12 \
			--sortUsing sum --startLabel "TSS" \
            --endLabel "TES" --xAxisLabel "" \
            --regionsLabel "WT open-region" "Ras open-region" "RasPPP1CA open-region"\
            --samplesLabel "WT Peaks"  "Ras Peaks" "RasPPP1CA Peaks" \
            -o $projPath/08Deeptools_BPM/${histName}_TSS_Macs2_narrowPeak_heatmap.pdf
	}
done



###############################################################
#14.2 Deeptools_RPKM
cd /storage/maxianjueLab/guoyifan/bulkKongDuATACseq
mkdir 08Deeptools_RPKM

module load R/4.2.1
source ~/.bashrc
conda activate deeptools

####################TSS
#!/bin/bash
#SBATCH -J depTSS_KD
#SBATCH -p amd-ep2,intel-sc3,amd-ep2-short
#SBATCH -q normal
#SBATCH --mem=120G
#SBATCH -c 12
projPath="/storage/maxianjueLab/guoyifan/bulkKongDuATACseq"
ref_gtf="/storage/maxianjueLab/guoyifan/species_reference/Drosophila_melanogaster_BDGP6_32/Drosophila_melanogaster.BDGP6.32.109.gtf"
for histName in RasPPP1CA_GFP_Pos RasPPP1CA_GFP_Neg A1 A2 A3 A4 A5 KA1 KA2 KA3 KA4 KA5;do
{
     computeMatrix reference-point -S ${projPath}/06BigWig/${histName}_normalized_RPKM_20bin.bw \
							   -p 12 \
							  --referencePoint TSS \
							  --afterRegionStartLength 3000 \
							  --beforeRegionStartLength 3000 \
							  -R $ref_gtf \
							  --skipZeros  --missingDataAsZero \
							  -o $projPath/08Deeptools_RPKM/${histName}_refPoint_TSS_data.gz
     plotHeatmap -m $projPath/08Deeptools_RPKM/${histName}_refPoint_TSS_data.gz \
            --missingDataColor 1 \
            --colorList 'white,#925E9F' \
            --heatmapHeight 12 \
			--sortUsing sum --startLabel "TSS" \
            --endLabel "TES" --xAxisLabel "" \
            --regionsLabel "Peaks" \
            --samplesLabel "${histName}" \
            -o $projPath/08Deeptools_RPKM/${histName}_refPoint_TSS_heatmap.pdf
    }
done


####################TSS_TES
#!/bin/bash
#SBATCH -J depGB_KD
#SBATCH -p amd-ep2,intel-sc3,amd-ep2-short
#SBATCH -q normal
#SBATCH --mem=120G
#SBATCH -c 12
projPath="/storage/maxianjueLab/guoyifan/bulkKongDuATACseq"
ref_gtf="/storage/maxianjueLab/guoyifan/species_reference/Drosophila_melanogaster_BDGP6_32/Drosophila_melanogaster.BDGP6.32.109.gtf"
for histName in RasPPP1CA_GFP_Pos RasPPP1CA_GFP_Neg A1 A2 A3 A4 A5 KA1 KA2 KA3 KA4 KA5;do
{
     computeMatrix reference-point -S ${projPath}/06BigWig/${histName}_normalized_RPKM_20bin.bw \
							   -p 8 \
							  --binSize 20 \
							  --regionBodyLength 5000 \
							  --referencePoint TSS \
							  --afterRegionStartLength 3000 \
							  --beforeRegionStartLength 3000 \
							  -R $ref_gtf \
							  --skipZeros  --missingDataAsZero \
							  -o $projPath/08Deeptools_RPKM/${histName}_refPoint_TSS_TES_data.gz
     plotHeatmap -m $projPath/08Deeptools_RPKM/${histName}_refPoint_TSS_TES_data.gz \
            --missingDataColor 1 \
            --colorList 'white,#339933' \
            --heatmapHeight 12 \
			--sortUsing sum --startLabel "TSS" \
            --endLabel "TES" --xAxisLabel "" \
            --regionsLabel "Peaks" \
            --samplesLabel "${histName}" \
            -o $projPath/08Deeptools_RPKM/${histName}_refPoint_TSS_TES_heatmap.pdf
    }
done


####################TSS Summary (Whole genome input)
cd /storage/maxianjueLab/guoyifan/bulkKongDuATACseq/08Deeptools_RPKM
mkdir TFs_of_interested

#####GFP negative
#!/bin/bash
#SBATCH -J depWG_KD
#SBATCH -p amd-ep2,intel-sc3,amd-ep2-short
#SBATCH -q normal
#SBATCH --mem=80G
#SBATCH -c 8
projPath="/storage/maxianjueLab/guoyifan/bulkKongDuATACseq"
ref_gtf="/storage/maxianjueLab/guoyifan/species_reference/Drosophila_melanogaster_BDGP6_32/Drosophila_melanogaster.BDGP6.32.109.gtf"
TFs_of_interested="/storage/maxianjueLab/guoyifan/bulkKongDuATACseq/08Deeptools_RPKM/TFs_of_interested"
for histName in WT_Ras_RasPPP1CA_GFP_neg_Summary;do
{
computeMatrix reference-point -S /storage/maxianjueLab/guoyifan/FTH/02Analysis/02ATAC-seq/06BigWig/E_normalized_RPKM_20bin.bw \
	                                                     /storage/maxianjueLab/guoyifan/FTH/02Analysis/02ATAC-seq/06BigWig/A_normalized_RPKM_20bin.bw \
														 ${projPath}/06BigWig/RasPPP1CA_GFP_Neg_normalized_RPKM_20bin.bw \
							   -p 8 \
							  --referencePoint TSS \
							  --afterRegionStartLength 3000 \
							  --beforeRegionStartLength 3000 \
							  -R $ref_gtf \
							  --skipZeros  --missingDataAsZero \
							  -o $projPath/08Deeptools_RPKM/${histName}_TSS_whole_genome_data.gz

plotHeatmap -m $projPath/08Deeptools_RPKM/${histName}_TSS_whole_genome_data.gz \
            --missingDataColor 1 \
            --heatmapHeight 12 \
			--sortUsing sum --startLabel "TSS" \
            --endLabel "TES" --xAxisLabel "" \
            --regionsLabel "Whole genome" \
            --samplesLabel "WT Peaks"  "Ras Peaks" "RasPPP1CA Peaks" \
            -o $projPath/08Deeptools_RPKM/${histName}_TSS_whole_genome_heatmap.pdf
	}
done



#####GFP positive
#!/bin/bash
#SBATCH -J depWG_KD
#SBATCH -p amd-ep2,intel-sc3,amd-ep2-short
#SBATCH -q normal
#SBATCH --mem=80G
#SBATCH -c 8
projPath="/storage/maxianjueLab/guoyifan/bulkKongDuATACseq"
ref_gtf="/storage/maxianjueLab/guoyifan/species_reference/Drosophila_melanogaster_BDGP6_32/Drosophila_melanogaster.BDGP6.32.109.gtf"
TFs_of_interested="/storage/maxianjueLab/guoyifan/bulkKongDuATACseq/08Deeptools_RPKM/TFs_of_interested"
for histName in WT_Ras_RasPPP1CA_GFP_pos_Summary;do
{
computeMatrix reference-point -S /storage/maxianjueLab/guoyifan/FTH/02Analysis/02ATAC-seq/06BigWig/F_normalized_RPKM_20bin.bw \
	                                                     /storage/maxianjueLab/guoyifan/FTH/02Analysis/02ATAC-seq/06BigWig/B_normalized_RPKM_20bin.bw \
														 ${projPath}/06BigWig/RasPPP1CA_GFP_Pos_normalized_RPKM_20bin.bw \
							   -p 8 \
							  --referencePoint TSS \
							  --afterRegionStartLength 3000 \
							  --beforeRegionStartLength 3000 \
							  -R $ref_gtf \
							  --skipZeros  --missingDataAsZero \
							  -o $projPath/08Deeptools_RPKM/${histName}_TSS_whole_genome_data.gz

plotHeatmap -m $projPath/08Deeptools_RPKM/${histName}_TSS_whole_genome_data.gz \
            --missingDataColor 1 \
            --heatmapHeight 12 \
			--sortUsing sum --startLabel "TSS" \
            --endLabel "TES" --xAxisLabel "" \
            --regionsLabel "Whole genome" \
            --samplesLabel "WT Peaks"  "Ras Peaks" "RasPPP1CA Peaks" \
            -o $projPath/08Deeptools_RPKM/${histName}_TSS_whole_genome_heatmap.pdf
	}
done


####################TSS Summary (MACS2 bed input)
##########MACS2 Summit.bed

#####GFP positive
#!/bin/bash
#SBATCH -J depPK_KD
#SBATCH -p amd-ep2,intel-sc3,amd-ep2-short
#SBATCH -q normal
#SBATCH --mem=80G
#SBATCH -c 8
projPath="/storage/maxianjueLab/guoyifan/bulkKongDuATACseq"
ref_gtf="/storage/maxianjueLab/guoyifan/species_reference/Drosophila_melanogaster_BDGP6_32/Drosophila_melanogaster.BDGP6.32.109.gtf"
TFs_of_interested="/storage/maxianjueLab/guoyifan/bulkKongDuATACseq/08Deeptools_RPKM/TFs_of_interested"
for histName in WT_Ras_RasPPP1CA_GFP_pos_Summary;do
{
computeMatrix reference-point -S /storage/maxianjueLab/guoyifan/FTH/02Analysis/02ATAC-seq/06BigWig/F_normalized_RPKM_20bin.bw \
	                                                     /storage/maxianjueLab/guoyifan/FTH/02Analysis/02ATAC-seq/06BigWig/B_normalized_RPKM_20bin.bw \
														 ${projPath}/06BigWig/RasPPP1CA_GFP_Pos_normalized_RPKM_20bin.bw \
							   -p 8 \
							  --referencePoint TSS \
							  --afterRegionStartLength 3000 \
							  --beforeRegionStartLength 3000 \
							  -R $TFs_of_interested/macs2_WT_GFP_Pos_peak_q0.05_summits.bed \
							  $TFs_of_interested/macs2_Ras_GFP_Pos_peak_q0.05_summits.bed \
							  $TFs_of_interested/macs2_RasPPP1CA_GFP_Pos_peak_q0.05_summits.bed \
							  --skipZeros  --missingDataAsZero \
							  -o $projPath/08Deeptools_RPKM/${histName}_TSS_Macs2_Peak_Summit_data.gz

plotHeatmap -m $projPath/08Deeptools_RPKM/${histName}_TSS_Macs2_Peak_Summit_data.gz \
            --heatmapHeight 12 \
			--sortUsing sum --startLabel "TSS" \
            --endLabel "TES" --xAxisLabel "" \
            --regionsLabel "WT open-region" "Ras open-region" "RasPPP1CA open-region"\
            --samplesLabel "WT Peaks"  "Ras Peaks" "RasPPP1CA Peaks" \
            -o $projPath/08Deeptools_RPKM/${histName}_TSS_Macs2_Peak_Summit_heatmap.pdf
	}
done



#####GFP negative
#!/bin/bash
#SBATCH -J depPK_KD
#SBATCH -p amd-ep2,intel-sc3,amd-ep2-short
#SBATCH -q normal
#SBATCH --mem=80G
#SBATCH -c 8
projPath="/storage/maxianjueLab/guoyifan/bulkKongDuATACseq"
ref_gtf="/storage/maxianjueLab/guoyifan/species_reference/Drosophila_melanogaster_BDGP6_32/Drosophila_melanogaster.BDGP6.32.109.gtf"
TFs_of_interested="/storage/maxianjueLab/guoyifan/bulkKongDuATACseq/08Deeptools_RPKM/TFs_of_interested"
for histName in WT_Ras_RasPPP1CA_GFP_neg_Summary;do
{
computeMatrix reference-point -S /storage/maxianjueLab/guoyifan/FTH/02Analysis/02ATAC-seq/06BigWig/E_normalized_RPKM_20bin.bw \
	                                                     /storage/maxianjueLab/guoyifan/FTH/02Analysis/02ATAC-seq/06BigWig/A_normalized_RPKM_20bin.bw \
														 ${projPath}/06BigWig/RasPPP1CA_GFP_Neg_normalized_RPKM_20bin.bw \
							   -p 8 \
							  --referencePoint TSS \
							  --afterRegionStartLength 3000 \
							  --beforeRegionStartLength 3000 \
							  -R $TFs_of_interested/macs2_WT_GFP_Neg_peak_q0.05_summits.bed \
							  $TFs_of_interested/macs2_Ras_GFP_Neg_peak_q0.05_summits.bed \
							  $TFs_of_interested/macs2_RasPPP1CA_GFP_Neg_peak_q0.05_summits.bed \
							  --skipZeros  --missingDataAsZero \
							  -o $projPath/08Deeptools_RPKM/${histName}_TSS_Macs2_Peak_Summit_data.gz

plotHeatmap -m $projPath/08Deeptools_RPKM/${histName}_TSS_Macs2_Peak_Summit_data.gz \
            --heatmapHeight 12 \
			--sortUsing sum --startLabel "TSS" \
            --endLabel "TES" --xAxisLabel "" \
            --regionsLabel "WT open-region" "Ras open-region" "RasPPP1CA open-region"\
            --samplesLabel "WT Peaks"  "Ras Peaks" "RasPPP1CA Peaks" \
            -o $projPath/08Deeptools_RPKM/${histName}_TSS_Macs2_Peak_Summit_heatmap.pdf
	}
done


##########MACS2 narrowPeak to bed

#####GFP positive
#!/bin/bash
#SBATCH -J depnarPK_KD
#SBATCH -p amd-ep2,intel-sc3,amd-ep2-short
#SBATCH -q normal
#SBATCH --mem=80G
#SBATCH -c 8
projPath="/storage/maxianjueLab/guoyifan/bulkKongDuATACseq"
ref_gtf="/storage/maxianjueLab/guoyifan/species_reference/Drosophila_melanogaster_BDGP6_32/Drosophila_melanogaster.BDGP6.32.109.gtf"
TFs_of_interested="/storage/maxianjueLab/guoyifan/bulkKongDuATACseq/08Deeptools_RPKM/TFs_of_interested"
for histName in WT_Ras_RasPPP1CA_GFP_pos_Summary;do
{
computeMatrix reference-point -S /storage/maxianjueLab/guoyifan/FTH/02Analysis/02ATAC-seq/06BigWig/F_normalized_RPKM_20bin.bw \
	                                                     /storage/maxianjueLab/guoyifan/FTH/02Analysis/02ATAC-seq/06BigWig/B_normalized_RPKM_20bin.bw \
														 ${projPath}/06BigWig/RasPPP1CA_GFP_Pos_normalized_RPKM_20bin.bw \
							   -p 8 \
							  --referencePoint TSS \
							  --afterRegionStartLength 3000 \
							  --beforeRegionStartLength 3000 \
							  -R $TFs_of_interested/macs2_WT_GFP_Pos_peak_q0.05_peaks_narrowPeak.bed \
							  $TFs_of_interested/macs2_Ras_GFP_Pos_peak_q0.05_peaks_narrowPeak.bed \
							  $TFs_of_interested/macs2_RasPPP1CA_GFP_Pos_peak_q0.05_peaks_narrowPeak.bed \
							  --skipZeros  --missingDataAsZero \
							  -o $projPath/08Deeptools_RPKM/${histName}_TSS_Macs2_narrowPeak_data.gz

plotHeatmap -m $projPath/08Deeptools_RPKM/${histName}_TSS_Macs2_narrowPeak_data.gz \
            --heatmapHeight 12 \
			--sortUsing sum --startLabel "TSS" \
            --endLabel "TES" --xAxisLabel "" \
            --regionsLabel "WT open-region" "Ras open-region" "RasPPP1CA open-region"\
            --samplesLabel "WT Peaks"  "Ras Peaks" "RasPPP1CA Peaks" \
            -o $projPath/08Deeptools_RPKM/${histName}_TSS_Macs2_narrowPeak_heatmap.pdf
	}
done



#####GFP negative
#!/bin/bash
#SBATCH -J depnarPK_KD
#SBATCH -p amd-ep2,intel-sc3,amd-ep2-short
#SBATCH -q normal
#SBATCH --mem=80G
#SBATCH -c 8
projPath="/storage/maxianjueLab/guoyifan/bulkKongDuATACseq"
ref_gtf="/storage/maxianjueLab/guoyifan/species_reference/Drosophila_melanogaster_BDGP6_32/Drosophila_melanogaster.BDGP6.32.109.gtf"
TFs_of_interested="/storage/maxianjueLab/guoyifan/bulkKongDuATACseq/08Deeptools_RPKM/TFs_of_interested"
for histName in WT_Ras_RasPPP1CA_GFP_neg_Summary;do
{
computeMatrix reference-point -S /storage/maxianjueLab/guoyifan/FTH/02Analysis/02ATAC-seq/06BigWig/E_normalized_RPKM_20bin.bw \
	                                                     /storage/maxianjueLab/guoyifan/FTH/02Analysis/02ATAC-seq/06BigWig/A_normalized_RPKM_20bin.bw \
														 ${projPath}/06BigWig/RasPPP1CA_GFP_Neg_normalized_RPKM_20bin.bw \
							   -p 8 \
							  --referencePoint TSS \
							  --afterRegionStartLength 3000 \
							  --beforeRegionStartLength 3000 \
							  -R $TFs_of_interested/macs2_WT_GFP_Neg_peak_q0.05_peaks_narrowPeak.bed \
							  $TFs_of_interested/macs2_Ras_GFP_Neg_peak_q0.05_peaks_narrowPeak.bed \
							  $TFs_of_interested/macs2_RasPPP1CA_GFP_Neg_peak_q0.05_peaks_narrowPeak.bed \
							  --skipZeros  --missingDataAsZero \
							  -o $projPath/08Deeptools_RPKM/${histName}_TSS_Macs2_narrowPeak_data.gz

plotHeatmap -m $projPath/08Deeptools_RPKM/${histName}_TSS_Macs2_narrowPeak_data.gz \
            --heatmapHeight 12 \
			--sortUsing sum --startLabel "TSS" \
            --endLabel "TES" --xAxisLabel "" \
            --regionsLabel "WT open-region" "Ras open-region" "RasPPP1CA open-region"\
            --samplesLabel "WT Peaks"  "Ras Peaks" "RasPPP1CA Peaks" \
            -o $projPath/08Deeptools_RPKM/${histName}_TSS_Macs2_narrowPeak_heatmap.pdf
	}
done


#15. 重复样本的处理
#15.1 Overlapping peaks using bedtools
module load bedtools/2.30.0

bedtools intersect \
-a macs2/Nanog-rep1_peaks.narrowPeak \
-b macs2/Nanog-rep2_peaks.narrowPeak \
-wo > bedtools/Nanog-overlaps.bed

#-a: 参数后加重复样本1（A）
#-b：参数后加重复样本2（B），也可以加多个样本
#如果只有-a和-b参数，返回的是相对于A的overlaps。加上参数-wo返回A和B原始的记录加上overlap的记录。参数-wa返回每个overlap中A的原始记录。
#-wo：Write the original A and B entries plus the number of base pairs of overlap between the two features.
#-wa：Write the original entry in A for each overlap.
#-v：Only report those entries in A that have no overlaps with 

#15.2 Irreproducibility Discovery Rate (IDR)
#https://github.com/nboley/idr
#https://www.jianshu.com/p/d8a7056b4294

#conda create -n idr
#conda install -c bioconda idr
cd /storage/maxianjueLab/guoyifan/FTH/02Analysis/02ATAC-seq

source ~/.bashrc
conda activate IDR

mkdir 09IDR
#Sort peak by -log10(p-value)
cd /storage/maxianjueLab/guoyifan/CUT_Tag/peakCalling/MACS2

#复制narrowPeak文件到09IDR

ref_gtf="/storage/maxianjueLab/guoyifan/species_reference/Drosophila_melanogaster_BDGP6_32/Drosophila_melanogaster.BDGP6.32.109.gtf"
projPath="/storage/maxianjueLab/guoyifan/bulkKongDuATACseq"
for histName in RasPPP1CA_GFP_Pos RasPPP1CA_GFP_Neg A1 A2 A3 A4 A5 KA1 KA2 KA3 KA4 KA5;do
{
 cp $projPath/07MACS2/macs2_${histName}_peak_q0.05_peaks.narrowPeak $projPath/09IDR/macs2_${histName}_peak_q0.05_peaks.narrowPeak 
 sort -k8,8nr $projPath/09IDR/macs2_${histName}_peak_q0.05_peaks.narrowPeak > $projPath/09IDR/macs2_${histName}_IDR_peaks.narrowPeak 
	}&
done	

#IDR寻找重复Peak
cd /storage/maxianjueLab/guoyifan/bulkKongDuATACseq/09IDR
mkdir Results
mkdir Second_IDR
#https://cloud.tencent.com/developer/article/1624533
idr --samples ../idr/test/data/peak1 ../idr/test/data/peak2 --peak-list ../idr/test/data/merged_peaks

#RasPPP1CA_GFP_Pos A4 A5 KA4 KA5

#A4 A5
idr --samples macs2_KA4_IDR_peaks.narrowPeak  macs2_KA5_IDR_peaks.narrowPeak  --peak-list macs2_RasPPP1CA_GFP_Pos_IDR_peaks.narrowPeak \
--input-file-type narrowPeak \
--rank p.value \
--output-file ./Results/RasPPP1CA_GFP_Pos_KA4_KA5 \
--plot \
--log-output-file ./Results/RasPPP1CA_GFP_Pos_KA4_KA5.log

idr --samples macs2_KA4_IDR_peaks.narrowPeak  macs2_KA5_IDR_peaks.narrowPeak  \
--input-file-type narrowPeak \
--rank p.value \
--output-file ./Results/KA4_KA5 \
--plot \
--log-output-file ./Results/KA4_KA5.log

#KA4 KA5
idr --samples macs2_A4_IDR_peaks.narrowPeak  macs2_A5_IDR_peaks.narrowPeak  --peak-list macs2_RasPPP1CA_GFP_Pos_IDR_peaks.narrowPeak \
--input-file-type narrowPeak \
--rank p.value \
--output-file ./Results/RasPPP1CA_GFP_Pos_A4_A5 \
--plot \
--log-output-file ./Results/RasPPP1CA_GFP_Pos_A4_A5.log


idr --samples macs2_A4_IDR_peaks.narrowPeak  macs2_A5_IDR_peaks.narrowPeak  \
--input-file-type narrowPeak \
--rank p.value \
--output-file ./Results/A4_A5 \
--plot \
--log-output-file ./Results/A4_A5.log

cd /storage/maxianjueLab/guoyifan/bulkKongDuATACseq/09IDR/Second_IDR

idr --samples A4_A5.narrowPeak  KA4_KA5.narrowPeak \
--input-file-type narrowPeak \
--rank p.value \
--output-file ./A4_A5_KA4_KA5 \
--plot \
--log-output-file ./A4_A5_KA4_KA5.log

idr --samples A4_A5.narrowPeak  KA4_KA5.narrowPeak --peak-list macs2_RasPPP1CA_GFP_Pos_IDR_peaks.narrowPeak \
--input-file-type narrowPeak \
--rank p.value \
--output-file ./RasPPP1CA_GFP_Pos_A4_A5_KA4_KA5 \
--plot \
--log-output-file ./RasPPP1CA_GFP_Pos_A4_A5_KA4_KA5.log

wc -l *-idr 计算下common peaks的个数，接着可再计算下与总peaks的比率。
如果想看IDR<0.05的，可以通过第5列信息过滤：
awk '{if($5 >= 540) print $0}' sample-idr | wc -l

#RasPPP1CA_GFP_Neg A1 A2 A3 KA1 KA2 KA3

#KA1 A1
idr --samples macs2_KA1_IDR_peaks.narrowPeak  macs2_A1_IDR_peaks.narrowPeak  --peak-list macs2_RasPPP1CA_GFP_Neg_IDR_peaks.narrowPeak \
--input-file-type narrowPeak \
--rank p.value \
--output-file ./Results/RasPPP1CA_GFP_Neg_KA1_A1 \
--plot \
--log-output-file ./Results/RasPPP1CA_GFP_Neg_KA1_A1.log &

idr --samples macs2_KA1_IDR_peaks.narrowPeak  macs2_A1_IDR_peaks.narrowPeak  \
--input-file-type narrowPeak \
--rank p.value \
--output-file ./Results/KA1_A1 \
--plot \
--log-output-file ./Results/KA1_A1.log &

#KA2 A2
idr --samples macs2_KA2_IDR_peaks.narrowPeak  macs2_A2_IDR_peaks.narrowPeak  --peak-list macs2_RasPPP1CA_GFP_Neg_IDR_peaks.narrowPeak \
--input-file-type narrowPeak \
--rank p.value \
--output-file ./Results/RasPPP1CA_GFP_Neg_KA2_A2 \
--plot \
--log-output-file ./Results/RasPPP1CA_GFP_Neg_KA2_A2.log &

idr --samples macs2_KA2_IDR_peaks.narrowPeak  macs2_A2_IDR_peaks.narrowPeak  \
--input-file-type narrowPeak \
--rank p.value \
--output-file ./Results/KA2_A2 \
--plot \
--log-output-file ./Results/KA2_A2.log &

#KA3 A3
idr --samples macs2_KA3_IDR_peaks.narrowPeak  macs2_A3_IDR_peaks.narrowPeak  --peak-list macs2_RasPPP1CA_GFP_Neg_IDR_peaks.narrowPeak \
--input-file-type narrowPeak \
--rank p.value \
--output-file ./Results/RasPPP1CA_GFP_Neg_KA3_A3 \
--plot \
--log-output-file ./Results/RasPPP1CA_GFP_Neg_KA3_A3.log &

idr --samples macs2_KA3_IDR_peaks.narrowPeak  macs2_A3_IDR_peaks.narrowPeak  \
--input-file-type narrowPeak \
--rank p.value \
--output-file ./Results/KA3_A3 \
--plot \
--log-output-file ./Results/KA3_A3.log &

cd /storage/maxianjueLab/guoyifan/bulkKongDuATACseq/09IDR/Second_IDR

idr --samples KA1_A1.narrowPeak  KA2_A2.narrowPeak \
--input-file-type narrowPeak \
--rank p.value \
--output-file ./KA1_A1_KA2_A2 \
--plot \
--log-output-file ./KA1_A1_KA2_A2.log

idr --samples KA1_A1_KA2_A2.narrowPeak  KA3_A3.narrowPeak \
--input-file-type narrowPeak \
--rank p.value \
--output-file ./KA1_A1_KA2_A2_KA3_A3 \
--plot \
--log-output-file ./KA1_A1_KA2_A2_KA3_A3.log

idr --samples KA1_A1_KA2_A2.narrowPeak  KA3_A3.narrowPeak --peak-list macs2_RasPPP1CA_GFP_Neg_IDR_peaks.narrowPeak \
--input-file-type narrowPeak \
--rank p.value \
--output-file ./RasPPP1CA_GFP_Neg_KA1_A1_KA2_A2_KA3_A3 \
--plot \
--log-output-file ./RasPPP1CA_GFP_Neg_KA1_A1_KA2_A2_KA3_A3.log

#16 ChIPSeeker in R
mkdir 10ChIPSeeker

library(dplyr)
library(dbplyr)
library(stringr)
library(ggplot2)
library(viridis)
library(clusterProfiler)
library(ChIPseeker)
library(GenomicFeatures)
library(GenomicRanges)
library(chromVAR) ## For FRiP analysis and differential analysis
library(DESeq2) ## For differential analysis section
library(ggpubr) ## For customizing figures
library(corrplot) ## For correlation plot
library(biomaRt)
library(curl)
library(org.Dm.eg.db)


txdb <- makeTxDbFromGFF("Drosophila_melanogaster.BDGP6.32.109.gtf",
                        format="gtf")    #可以使用gtf和gff3
						
						
keytypes(txdb)    #感兴趣的话，可以用以下方法探索txdb都包含了什么内容
keys(txdb)

list.files()
macs2_Ras_GFP_Neg_peak_q0.05_summits.bed
macs2_Ras_GFP_Pos_peak_q0.05_summits.bed
macs2_RasPPP1CA_GFP_Neg_peak_q0.05_summits.bed
macs2_RasPPP1CA_GFP_Pos_peak_q0.05_summits.bed
macs2_WT_GFP_Neg_peak_q0.05_summits.bed
macs2_WT_GFP_Pos_peak_q0.05_summits.bed

#读入单个summits文件
peaks <- readPeakFile("macs2_Ras_GFP_Neg_peak_q0.05_summits.bed")
peaks <- readPeakFile("macs2_Ras_GFP_Pos_peak_q0.05_summits.bed")
peaks <- readPeakFile("macs2_RasPPP1CA_GFP_Neg_peak_q0.05_summits.bed")
peaks <- readPeakFile("macs2_RasPPP1CA_GFP_Pos_peak_q0.05_summits.bed")
peaks <- readPeakFile("macs2_WT_GFP_Neg_peak_q0.05_summits.bed")
peaks <- readPeakFile("macs2_WT_GFP_Pos_peak_q0.05_summits.bed")

peaks <- readPeakFile("02IDR_RasPPP1CA_GFP_Neg_KA1_A1_KA2_A2_KA3_A3.bed")
peaks <- readPeakFile("02IDR_RasPPP1CA_GFP_Pos_A4_A5_KA4_KA5.bed")

#结构注释
peakAnno <- annotatePeak(peaks, TxDb=txdb, tssRegion=c(-2000, 2000))

#最后将我们的注释结果转为数据框，便于查看
df <- as.data.frame(peakAnno)
head(df)
#将注释到的基因提取出来（MACS2，第14列），用于后续功能分析
#将注释到的基因提取出来（IDR.bed，第29列），用于后续功能分析
#gene <- df[,14]
gene <- df[,29]
###对基因进行注释-获取gene_symbol
gene.df <- bitr(gene, fromType = "ENSEMBL", 
                toType = c("SYMBOL","ENTREZID"),
                OrgDb = org.Dm.eg.db) 
# geneID    需要转换的ID
# fromType  当前ID类型
# toType    转换成什么ID，使用keytypes()查看有哪些类型
# OrgDb     注释数据库
# 注意相同的FBgn号会注释出不同的基因，所以以下代码会报错
#row.names(gene.df) <- gene.df$ENSEMBL重复 FBgn0004828对应 His3.3B和His3.3A
#row.names(gene.df) <- gene.df$SYMBOL重复 His2B:CG33872对应 FBgn0053910和FBgn0053872
write.csv(gene.df,file="peakAnno_Ras_GFP_Neg.csv",sep="/t")
write.csv(gene.df,file="peakAnno_Ras_GFP_Pos.csv",sep="/t")
write.csv(gene.df,file="peakAnno_RasPPP1CA_GFP_Neg.csv",sep="/t")
write.csv(gene.df,file="peakAnno_RasPPP1CA_GFP_Pos.csv",sep="/t")
write.csv(gene.df,file="peakAnno_WT_GFP_Neg.csv",sep="/t")
write.csv(gene.df,file="peakAnno_WT_GFP_Pos.csv",sep="/t")

write.csv(gene.df,file="02IDR_peakAnno_RasPPP1CA_GFP_Neg_KA1_A1_KA2_A2_KA3_A3.csv",sep="/t")
write.csv(gene.df,file="02IDR_peakAnno_RasPPP1CA_GFP_Pos_A4_A5_KA4_KA5.csv",sep="/t")


a <- read.table(file = "04Venn_Only_RasPPP1CA_GFP_Pos_568.txt",sep = "\t",header = T)
gene.df2 <- bitr(a[,1], fromType = "ENSEMBL", 
                toType = c("SYMBOL","ENTREZID"),
                OrgDb = org.Dm.eg.db) 
write.csv(gene.df2,file="04Venn_Anno_Only_RasPPP1CA_GFP_Pos_568.csv",sep="/t")

a <- read.table(file = "04Venn_Only_RasPPP1CA_GFP_Neg_737.txt",sep = "\t",header = T)
gene.df2 <- bitr(a[,1], fromType = "ENSEMBL", 
                toType = c("SYMBOL","ENTREZID"),
                OrgDb = org.Dm.eg.db) 
write.csv(gene.df2,file="04Venn_Anno_Only_RasPPP1CA_GFP_Neg_737.csv",sep="/t")



#一次也可以读入多个summits文件，使用list存储，然后使用lapply注释
files = list(E75_Flag = ("macs2_E75_Flag_peak_q0.05_summits.bed"), 
				 NY_N = ("macs2_NY_N_peak_q0.05_summits.bed"), 
				 NY_Y = ("macs2_NY_Y_peak_q0.05_summits.bed"),
				 ENY_N = ("macs2_ENY_N_peak_q0.05_summits.bed"), 
				 ENY_Y = ("macs2_ENY_Y_peak_q0.05_summits.bed"))
peakAnnoList <- lapply(files, 
                       annotatePeak,
                       TxDb=txdb,
                       tssRegion=c(-2000, 2000))
plotAnnoBar(peakAnnoList)
plotDistToTSS(peakAnnoList)


#注释完，进行可视化，多种图可供选择
plotAnnoBar(peakAnno)
plotDistToTSS(peakAnno)
vennpie(peakAnno)
plotAnnoPie(peakAnno)
#install.packages("ggupset")
library(ggupset)
upsetplot(peakAnno)
#install.packages("ggimage")
library(ggimage)
upsetplot(peakAnno, vennpie=TRUE)
















