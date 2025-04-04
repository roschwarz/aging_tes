#!/usr/bin/env bash

####### Brain ####

computeMatrix scale-regions -S ../../results/cage_RS/brain_downsampled_gclipped/coverage/brain_downsampled.cage_peaks.one.merged.bw \
                    -R ../../data/processed/s03_brain_multi_TE_island.bed \
                    -p 20 \
                    -o ../../data/processed/s03_brain_multi_TE_island_matrix.gz \
            --beforeRegionStartLength 1000 \
            --afterRegionStartLength 1000 \
            --regionBodyLength 2000 \
            --missingDataAsZero 

plotHeatmap -m  ../../data/processed/s03_brain_multi_TE_island_matrix.gz\
                    --dpi 300 \
                    -out ../../results/figures/s03_brain_multi_TE_island_density.svg \
                    --colorList 'white, #373F43' \
                    --startLabel "TE island start" \
                    --endLabel "TE island end" \
                    --linesAtTickMarks \
                    --heatmapWidth 10 \
                    --heatmapHeight 34 \
                    --regionsLabel "TE island"


####### skin ####

computeMatrix scale-regions -S ../../results/cage_RS/skinII_gclipped/coverage/skinII.cage_peaks.one.merged.bw \
                    -R ../../data/processed/s03_skin_multi_TE_island.bed \
                    -p 20 \
                    -o ../../data/processed/s03_skin_multi_TE_island_matrix.gz \
            --beforeRegionStartLength 1000 \
            --afterRegionStartLength 1000 \
            --regionBodyLength 2000 \
            --missingDataAsZero 

plotHeatmap -m  ../../data/processed/s03_skin_multi_TE_island_matrix.gz\
                    --dpi 300 \
                    -out ../../results/figures/s03_skin_multi_TE_island_density.svg \
                    --colorList 'white, #373F43' \
                    --startLabel "TE island start" \
                    --endLabel "TE island end" \
                    --linesAtTickMarks \
                    --heatmapWidth 10 \
                    --heatmapHeight 34 \
                    --regionsLabel "TE island"

####### blood ####

computeMatrix scale-regions -S ../../results/cage_RS/blood_gclipped/coverage/blood.cage_peaks.one.merged.bw \
                    -R ../../data/processed/s03_blood_multi_TE_island.bed \
                    -p 20 \
                    -o ../../data/processed/s03_blood_multi_TE_island_matrix.gz \
            --beforeRegionStartLength 1000 \
            --afterRegionStartLength 1000 \
            --regionBodyLength 2000 \
            --missingDataAsZero 

plotHeatmap -m  ../../data/processed/s03_blood_multi_TE_island_matrix.gz\
                    --dpi 300 \
                    -out ../../results/figures/s03_blood_multi_TE_island_density.svg \
                    --colorList 'white, #373F43' \
                    --startLabel "TE island start" \
                    --endLabel "TE island end" \
                    --linesAtTickMarks \
                    --heatmapWidth 10 \
                    --heatmapHeight 34 \
                    --regionsLabel "TE island"

