#!/usr/bin/env bash

source ../lib/cageFileHandler.sh
ChromSizes=../../data/shared/GRCm38.p6.chrom.sizes

OUTPUTROOT=../../results/cage_RS/ 
READS=clipped

create_bigWigs(){

    local OUTPUTROOT=$1
    local read_type=$2
    local ChromSizes=$3

    # merge peaks if they overlap and create bigWig files as bw don't have
    # strand specificity?
    for data_set in $(getDirectory "$OUTPUTROOT" "$read_type"); do

        tissue=$(basename "$data_set" _gclipped)
        
        if [[ "$tissue" != "brain" && "$tissue" != "skinI" ]]; then

            bedGraph_file="$data_set"/coverage/${tissue}.cage_peaks.one.bedgraph
            bedGraph_M_file="$data_set"/coverage/${tissue}.cage_peaks.one.merged.bedgraph

            bedtools merge -i "$bedGraph_file" | \
                     awk -v OFS='\t' '{print $1,$2,$3,1}' | \
                     sort -k1,1 -k2,2n > "$bedGraph_M_file"

            bedGraphToBigWig "$bedGraph_M_file" "$ChromSizes" "$data_set"/coverage/${tissue}.cage_peaks.one.merged.bw


        fi

    done

}


calc_density_matrix(){

    local OUTPUTROOT=$1
    local read_type=$2

    for data_set in $(getDirectory "$OUTPUTROOT" "$read_type"); do

        tissue=$(basename "$data_set" _gclipped)
        
        if [[ "$tissue" != "brain" && "$tissue" != "skinI" ]]; then


            bed_file=../../data/shared/${tissue}_independent_TE_regions.bed
            bw_file="$data_set"/coverage/${tissue}.cage_peaks.one.merged.bw
            matrix=../../data/shared/${tissue}_independent_TE_regions.dens_matrix.gz

            # Compute the matrix that is needed for the figure
            echo computeMatrix scale-regions -S $bw_file -R ${bed_file} -p 20 \
            --beforeRegionStartLength 1000 --afterRegionStartLength 1000 --regionBodyLength 2000 --missingDataAsZero -o $matrix


            computeMatrix scale-regions -S $bw_file \
                    -R ${bed_file} \
                    -p 20 \
                    -o $matrix \
            --beforeRegionStartLength 1000 \
            --afterRegionStartLength 1000 \
            --regionBodyLength 2000 \
            --missingDataAsZero 
        fi
    done


}

create_density_plots(){

    local OUTPUTROOT=$1
    local read_type=$2

    for data_set in $(getDirectory "$OUTPUTROOT" "$read_type"); do

        tissue=$(basename "$data_set" _gclipped)
        
        if [[ "$tissue" != "brain" && "$tissue" != "skinI" ]]; then

            matrix=../../data/shared/${tissue}_independent_TE_regions.dens_matrix.gz
            plotName=../../results/figures/$tissue.density_plot.svg


           echo plotHeatmap -m $matrix \
            --dpi 300 \
            -out $plotName \
            --colorList 'white, #373F43' \
            --startLabel "TE-region start" \
            --endLabel "TE-region end" \
            --labelRotation 90 \
            --linesAtTickMarks \
            --heatmapWidth 10 \
            --heatmapHeight 34 \
            --regionsLabel "TE-regions"

            plotHeatmap -m $matrix \
                    --dpi 300 \
                    -out $plotName \
                    --colorList 'white, #373F43' \
                    --startLabel "TE-region start" \
                    --endLabel "TE-region end" \
                    --labelRotation 90 \
                    --linesAtTickMarks \
                    --heatmapWidth 10 \
                    --heatmapHeight 34 \
                    --regionsLabel "TE-regions"



        fi


    done

}

# create_bigWigs $OUTPUTROOT $READS $ChromSizes
# calc_density_matrix $OUTPUTROOT $READS
# create_density_plots $OUTPUTROOT $READS
#
#
create_bigWigs_quant(){

    local OUTPUTROOT=$1
    local ChromSizes=$2

    # merge peaks if they overlap and create bigWig files as bw don't have
    # strand specificity?
    mapfile -t dirs < <(find "$OUTPUTROOT" -mindepth 1 -maxdepth 1 -type d \( -name "brain" -o -name "blood" -o -name "skinI" -o -name "skinII" \))

    for data_set in "${dirs[@]}"; do

        tissue=$(basename "$data_set")

        if [[ "$tissue" != "skinI" ]]; then

            bedGraph_file="$data_set"/coverage/${tissue}.quant_peaks.one.bedgraph
            bedGraph_M_file="$data_set"/coverage/${tissue}.quant_peaks.one.merged.bedgraph

            echo $bedGraph_file
            echo $bedGraph_M_file

            bedtools merge -i $bedGraph_file | \
                     awk -v OFS='\t' '{print $1,$2,$3,1}' | \
                     sort -k1,1 -k2,2n > $bedGraph_M_file

            bedGraphToBigWig $bedGraph_M_file $ChromSizes "$data_set"/coverage/${tissue}.quant_peaks.one.merged.bw


        fi

    done

}

calc_density_matrix_quant(){

    local OUTPUTROOT=$1

    mapfile -t dirs < <(find "$OUTPUTROOT" -mindepth 1 -maxdepth 1 -type d \( -name "brain" -o -name "blood" -o -name "skinI" -o -name "skinII" \))

    for data_set in "${dirs[@]}"; do

        tissue=$(basename "$data_set")
        
        if [[ "$tissue" != "skinI" ]]; then

            if [[ "$tissue" == "brain" ]]; then
                
                bed_file=../../data/shared/${tissue}_downsampled_independent_TE_regions.bed

            else
                bed_file=../../data/shared/${tissue}_independent_TE_regions.bed
            fi

            bw_file="$data_set"/coverage/${tissue}.quant_peaks.one.merged.bw
            matrix=../../data/shared/${tissue}_independent_TE_regions.quant.dens_matrix.gz

            # Compute the matrix that is needed for the figure
            echo computeMatrix scale-regions -S $bw_file -R ${bed_file} -p 20 \
            --beforeRegionStartLength 1000 --afterRegionStartLength 1000 --regionBodyLength 2000 --missingDataAsZero -o $matrix


            computeMatrix scale-regions -S $bw_file \
                    -R ${bed_file} \
                    -p 20 \
                    -o $matrix \
            --beforeRegionStartLength 1000 \
            --afterRegionStartLength 1000 \
            --regionBodyLength 2000 \
            --missingDataAsZero 
        fi
    done


}


create_density_plots_quant(){

    local OUTPUTROOT=$1

    mapfile -t dirs < <(find "$OUTPUTROOT" -mindepth 1 -maxdepth 1 -type d \( -name "brain" -o -name "blood" -o -name "skinI" -o -name "skinII" \))

    for data_set in "${dirs[@]}"; do

        tissue=$(basename "$data_set")
        
        if [[ "$tissue" != "skinI" ]]; then

            matrix=../../data/shared/${tissue}_independent_TE_regions.quant.dens_matrix.gz
            plotName=../../results/figures/$tissue.quant.density_plot.svg


           echo plotHeatmap -m $matrix \
            --dpi 300 \
            -out $plotName \
            --colorList 'white, #373F43' \
            --startLabel "TE-region start" \
            --endLabel "TE-region end" \
            --labelRotation 90 \
            --linesAtTickMarks \
            --heatmapWidth 2 \
            --heatmapHeight 7 \
            --regionsLabel "TE-regions"

            plotHeatmap -m $matrix \
                    --dpi 300 \
                    -out $plotName \
                    --colorList 'white, #373F43' \
                    --startLabel "TE-region start" \
                    --endLabel "TE-region end" \
                    --labelRotation 90 \
                    --linesAtTickMarks \
                    --heatmapWidth 2 \
                    --heatmapHeight 7 \
                    --regionsLabel "TE-regions"
        fi
    done

}
OUTPUTROOT=../../results/quant_RS/
#create_bigWigs_quant $OUTPUTROOT $ChromSizes
#calc_density_matrix_quant $OUTPUTROOT
create_density_plots_quant $OUTPUTROOT

############# Quant Density for annotation based TE regions ###################

#calc_density_matrix_quant(){
#
#    local OUTPUTROOT=$1
#
#    mapfile -t dirs < <(find "$OUTPUTROOT" -mindepth 1 -maxdepth 1 -type d \( -name "brain" -o -name "blood" -o -name "skinI" -o -name "skinII" \))
#
#    for data_set in "${dirs[@]}"; do
#
#        tissue=$(basename "$data_set")
#        
#        if [[ "$tissue" != "skinI" ]]; then
#
#            echo Calc density matrix of $tissue
#
#            bed_file=../../data/shared/mm10_TE_region_3prime_extended.bed
#
#            #if [[ "$tissue" == "brain" ]]; then
#            #    
#            #    bed_file=../../data/shared/${tissue}_downsampled_independent_TE_regions.bed
#
#            #else
#            #    bed_file=../../data/shared/${tissue}_independent_TE_regions.bed
#            #fi
#
#            bw_file="$data_set"/coverage/${tissue}.quant_peaks.one.merged.bw
#            matrix=../../data/shared/${tissue}_TE_regions.quant.dens_matrix.gz
#
#            # Compute the matrix that is needed for the figure
#            echo computeMatrix scale-regions -S $bw_file -R ${bed_file} -p 20 \
#            --beforeRegionStartLength 1000 --afterRegionStartLength 1000 --regionBodyLength 2000 --missingDataAsZero -o $matrix
#
#
#            computeMatrix scale-regions -S $bw_file \
#                    -R ${bed_file} \
#                    -p 20 \
#                    -o $matrix \
#            --beforeRegionStartLength 1000 \
#            --afterRegionStartLength 1000 \
#            --regionBodyLength 2000 \
#            --missingDataAsZero 
#        fi
#    done
#
#
#}
#
#
#create_density_plots_quant(){
#
#    local OUTPUTROOT=$1
#
#    mapfile -t dirs < <(find "$OUTPUTROOT" -mindepth 1 -maxdepth 1 -type d \( -name "brain" -o -name "blood" -o -name "skinI" -o -name "skinII" \))
#
#    for data_set in "${dirs[@]}"; do
#
#        tissue=$(basename "$data_set")
#        
#        if [[ "$tissue" != "skinI" ]]; then
#
#            echo Plot density plot of $tissue
#
#            matrix=../../data/shared/${tissue}_TE_regions.quant.dens_matrix.gz
#            plotName=../../results/figures/$tissue.TE_region_3prime_extended.quant.density_plot.svg
#
#
#           echo plotHeatmap -m $matrix \
#            --dpi 300 \
#            -out $plotName \
#            --colorList 'white, #373F43' \
#            --startLabel "TE-region start" \
#            --endLabel "TE-region end" \
#            --labelRotation 90 \
#            --linesAtTickMarks \
#            --heatmapWidth 2\
#            --heatmapHeight 9 \
#            --regionsLabel "TE-regions"
#
#            plotHeatmap -m $matrix \
#                    --dpi 300 \
#                    -out $plotName \
#                    --colorList 'white, #373F43' \
#                    --startLabel "TE-region start" \
#                    --endLabel "TE-region end" \
#                    --labelRotation 90 \
#                    --linesAtTickMarks \
#                    --heatmapWidth 2 \
#                    --heatmapHeight 7 \
#                    --regionsLabel "TE-regions"
#        fi
#    done
#
#}
#OUTPUTROOT=../../results/quant_RS/
#create_bigWigs_quant $OUTPUTROOT $ChromSizes
#calc_density_matrix_quant $OUTPUTROOT
#create_density_plots_quant $OUTPUTROOT
