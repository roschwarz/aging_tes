#!/usr/bin/bash

# merges TE instances in close proximity to TE regions and extends the created 
# TE regions 500 bp towards the 5' end.

TE_annotation=$1
ChromSizes=$2
output=$3

bedtools merge -i "$TE_annotation" -s -d 500 -c 6 -o distinct | \
    awk '{print $1"\t"$2"\t"$3"\tTE_Cluster_" NR "\t42\t"$4}' > "$output"mm10_TE_region.bed 

bedtools slop -i <(grep -v "^chr[1-9a-zA-Z]*_" "$output"mm10_TE_region.bed) \
    -g "$ChromSizes" -l 500 -r 0 -s > "$output"mm10_TE_region_5prime_extended.bed 

bedtools slop -i <(grep -v "^chr[1-9a-zA-Z]*_" "$output"mm10_TE_region.bed) \
    -g "$ChromSizes" -l 0 -r 500 -s > "$output"mm10_TE_region_3prime_extended.bed 

bedtools slop -i <(grep -v "^chr[1-9a-zA-Z]*_" "$output"mm10_TE_region.bed) \
    -g "$ChromSizes" -l 500 -r 500 -s > "$output"mm10_TE_region_full_extended.bed 
