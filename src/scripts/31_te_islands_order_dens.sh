
echo 'Compute Matrix'

computeMatrix scale-regions \
    -S ../../results/te_islands/SINE_ignore_strand.bw ../../results/te_islands/LINE_ignore_strand.bw ../../results/te_islands/LTR_ignore_strand.bw ../../results/te_islands/DNA_ignore_strand.bw  \
    -R ../../data/shared/brain_can_TE_island.bed ../../data/shared/skin_can_TE_island.bed ../../data/shared/blood_can_TE_island.bed -b 500 -a 500 -p 56 --missingDataAsZero \
    -o ../../results/te_islands/deeptools_matrix.gz

echo 'Make The Plot'
plotHeatmap -m ../../results/te_islands/deeptools_matrix.gz \
    --dpi 300 \
    -out ../../results/figures/31_te_island_dens.svg \
    --colorList 'white, #373F43' \
    --startLabel 'TE island start' \
    --endLabel 'TE island end' \
    --labelRotation 45 \
    --linesAtTickMarks \
    --heatmapWidth 1.5 \
    --heatmapHeight 10

