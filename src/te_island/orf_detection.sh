#!/usr/bin/bash
#

TE_ISLAND_PATH="../../results/te_island/"

for tissue in blood brain_downsampled skin; do

    bedtools getfasta -fi /misc/paras/data/genomes/GRCm38_v102/GRCm38.fa \
        -fo "$TE_ISLAND_PATH"/"$tissue"_indie_te_island.fa \
        -bed "$TE_ISLAND_PATH"/"$tissue"_indie_te_island.bed \
        -nameOnly

    orfipy "$TE_ISLAND_PATH"/"$tissue"_indie_te_island.fa \
        --bed "$tissue"_te_islands_orfs.bed \
        --min 100 \
        --max 10000 \
        --procs 20 \
        --table 1 \
        --outdir "$TE_ISLAND_PATH"/orfs_"$tissue" \
        --longest
done
