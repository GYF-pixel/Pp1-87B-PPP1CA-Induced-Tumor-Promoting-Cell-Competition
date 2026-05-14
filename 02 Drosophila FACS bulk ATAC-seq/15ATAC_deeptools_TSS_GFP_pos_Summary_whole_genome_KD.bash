#!/bin/bash
#SBATCH -J depWG_KD
#SBATCH -p amd-ep2,intel-sc3,amd-ep2-short
#SBATCH -q normal
#SBATCH --mem=80G
#SBATCH -c 8
projPath="/storage/maxianjueLab/guoyifan/bulkKongDuATACseq"
ref_gtf="/storage/maxianjueLab/guoyifan/species_reference/Drosophila_melanogaster_BDGP6_32/Drosophila_melanogaster.BDGP6.32.109.gtf"
TFs_of_interested="/storage/maxianjueLab/guoyifan/bulkKongDuATACseq/08Deeptools/TFs_of_interested"
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
							  -o $projPath/08Deeptools/${histName}_TSS_whole_genome_data.gz

plotHeatmap -m $projPath/08Deeptools/${histName}_TSS_whole_genome_data.gz \
            --missingDataColor 1 \
            --heatmapHeight 12 \
			--sortUsing sum --startLabel "TSS" \
            --endLabel "TES" --xAxisLabel "" \
            --regionsLabel "Whole genome" \
            --samplesLabel "WT Peaks"  "Ras Peaks" "RasPPP1CA Peaks" \
            -o $projPath/08Deeptools/${histName}_TSS_whole_genome_heatmap.pdf
	}
done
