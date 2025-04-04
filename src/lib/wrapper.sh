
########################
# Peak specific things #
########################

callPeaks(){

    # Call peakachu to call peaks in stranded .bam-files, named as following
    # <name>[reverse|forward].bam.
    #
    #
    # @param
    #
    # <- directory with .bam files
    # <- directory name of file destiny
    #
    # @return
    #
    # -> stores peak files in the given peak result directory

    echo Call peaks with peakachu using bashbone.
    local bam_dir=$1
    local peak_res=$2

    for strand in forward reverse; do
        mapper=(segemehl);
        bamFiles=$(find "$bam_dir"/*"$strand".bam -type f )
        segemehl=(${bamFiles[@]})
        echo peaks::peakachu -t 56 -f 150 -r mapper -o "$peak_res" -z true
        peaks::peakachu -t 56 -f 150 -r mapper -o "$peak_res" -z true
    done

}

callPeaks_tissue(){


    local OUTPUTROOT=$1
    local read_type=$2
    
    for data_set in $(getDirectory "$OUTPUTROOT" "$read_type"); do

        bam_dir="$data_set"/bams_stranded_segemehl/
        peak_res=$data_set/raw_peaks/
        
        mkdir -p "$peak_res"

        callPeaks "$bam_dir" "$peak_res"
    done
}

merge_narrowPeaks(){

    # merge peaks that are within a certain proximity (size).
    # The input are merged narrowPeak files, so that they are transformed before
    # merging the peaks. Be careful, this works only when the peaks are called in
    # strand specific manner. The bedtools call does not use the -s flag as the 
    # peak files are already strand specific. Names of the merged peaks are 
    # stored in the 4. column as a comma separated list.
    #
    # @param (optional)
    #
    # <- directory with peaks
    # <- window size for merging
    # (<-) - output directory, if not set it will be set to the peak directory
    #
    # @return
    #
    # -> .bed file that contains coordinates of merged peaks
    #

	local peaks_res=$1    
	local size=$2
    local output=$3

    if [[ ! "$output" ]]; then
        echo Output is not given and is set to "$peaks_res"
        output="$peaks_res"
    fi

    echo merge peaks

    for peaks in "$peaks_res"/*.narrowPeak; do

        [[ $(getStrand "$peaks") == 'forward' ]] && strand="+" || strand="-"

        merged_peaks=$(basename "$peaks" .narrowPeak).bed

        bedtools merge -i "$peaks" -d "$size" -c 4,5,6 -o distinct,distinct,distinct | \
        sed "s/\./$strand/" > "$output"/"${merged_peaks}"

    done

}
