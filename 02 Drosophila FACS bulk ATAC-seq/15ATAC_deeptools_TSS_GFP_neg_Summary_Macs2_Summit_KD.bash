#!/bin/bash
#SBATCH -J depPK_KD
#SBATCH -p amd-ep2,intel-sc3,amd-ep2-short
#SBATCH -q normal
#SBATCH --mem=80G
#SBATCH -c 8
projPath="/storage/maxianjueLab/guoyifan/bulkKongDuATACseq"
ref_gtf="/storage/maxianjueLab/guoyifan/species_reference/Drosophila_melanogaster_BDGP6_32/Drosophila_melanogaster.BDGP6.32.109.gtf"
TFs_of_interested="/storage/maxianjueLab/guoyifan/bulkKongDuATACseq/08Deeptools/TFs_of_interested"
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
							  -o $projPath/08Deeptools/${histName}_TSS_Macs2_Peak_Summit_data.gz

plotHeatmap -m $projPath/08Deeptools/${histName}_TSS_Macs2_Peak_Summit_data.gz \
            --heatmapHeight 12 \
			--sortUsing sum --startLabel "TSS" \
            --endLabel "TES" --xAxisLabel "" \
            --regionsLabel "WT open-region" "Ras open-region" "RasPPP1CA open-region"\
            --samplesLabel "WT Peaks"  "Ras Peaks" "RasPPP1CA Peaks" \
            -o $projPath/08Deeptools/${histName}_TSS_Macs2_Peak_Summit_heatmap.pdf
	}
done
