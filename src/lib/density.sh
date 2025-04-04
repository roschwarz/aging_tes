#!/usr/bin/bash

calc_matrix(){

    local regionOfInterest=$1
    local matrixName=$2
    shift 2 # shift all arguments to the left
    local bwFiles=("$@") # rebuild the array with all arguments that are left

    matrix_command=(computeMatrix scale-regions -S "${bwFiles[@]}" -p 20)
        matrix_command+=(-R "$regionOfInterest")
        matrix_command+=(--beforeRegionStartLength 1000)
        matrix_command+=(--afterRegionStartLength 1000)
        matrix_command+=(--regionBodyLength 2000)
        matrix_command+=(--missingDataAsZero)
        matrix_command+=(-o "$matrixName")


    echo "${matrix_command[@]}"
    "${matrix_command[@]}"

}

calc_matrix_center(){

    local regionOfInterest=$1
    local matrixName=$2
    shift 2 # shift all arguments to the left
    local bwFiles=("$@") # rebuild the array with all arguments that are left

    matrix_command=(computeMatrix reference-point -S "${bwFiles[@]}" -p 40)
        matrix_command+=(-R "$regionOfInterest")
        matrix_command+=(-b 1000)
        matrix_command+=(-a 1000)
        matrix_command+=(--referencePoint center)
        matrix_command+=(--missingDataAsZero)
        matrix_command+=(-o "$matrixName")


    echo "${matrix_command[@]}"
    "${matrix_command[@]}"

}

draw_plot(){

    local matrix=$1
    local fig_name=$2

    plot_command=(plotHeatmap -m "$matrix")
    plot_command+=(-out "$fig_name")
    plot_command+=(--sortUsing sum)

    echo "${plot_command[@]}"
    "${plot_command[@]}"

}

prep_density_plot(){

    narrowPath=$1
    chromSize=$2
    destiny=$3

    bedGraphPath=$destiny/

    mkdir -p "$bedGraphPath"

    for file in "$narrowPath"/*.narrowPeak; do

        sra=$(basename "$file" .unique.sorted.rmdup.narrowPeak).bedgraph

        awk 'BEGIN {OFS="\t"}{print $1, $2, $3, $8}' "$file" > $bedGraphPath/$sra

        bwFile=$bedGraphPath/$(basename "$sra" .bedgraph).bw

        bedGraphToBigWig "$bedGraphPath"/"$sra" "$chromSize" "$bwFile"

        bwFiles+=($(basename "$sra" .bedgraph).bw)

    done

    echo "${bwFiles[@]}"
}
