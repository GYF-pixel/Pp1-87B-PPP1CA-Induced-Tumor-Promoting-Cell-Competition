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
							  -o $projPath/08Deeptools/${histName}_refPoint_TSS_data.gz
     plotHeatmap -m $projPath/08Deeptools/${histName}_refPoint_TSS_data.gz \
            --missingDataColor 1 \
            --colorList 'white,#925E9F' \
            --heatmapHeight 12 \
			--sortUsing sum --startLabel "TSS" \
            --endLabel "TES" --xAxisLabel "" \
            --regionsLabel "Peaks" \
            --samplesLabel "${histName}" \
            -o $projPath/08Deeptools/${histName}_refPoint_TSS_heatmap.pdf
    }
done
