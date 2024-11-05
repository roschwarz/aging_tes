
peaksBedtoGTF(){

	# Takes two bed files and merges them and crate a gtf file
	# To-Do: 
    #   -make it more generic, to read more than 2 bed files or only one.
    #   - provide opportunity to set the attribute in collumn 3 (tts, tss)
	# If multiple peaks were merged up-stream of that function only the first 
	# peak name is used.
    # In the current pipeline the quant-peak bed-files are annotated in an 
    # anti-sense manner. Now, the strand will be switched in the GTF file.

	bedFile1=$1
	bedFile2=$2
	resultName=$3
    type=$4


    if [[ $type == 'quant' ]]; then

        cat "$bedFile1" "$bedFile2" | \
            sort -k1,1 -k2,2n | \
            sed 's/,peak_[0-9]*//g' | \
        awk 'BEGIN{OFS="\t";}{if($6=="+") print $1,"quant_peak","tts",$2,$3,$5,"-",".","gene_id \""$4"_rev\";" ; else print $1,"quant_peak","tts",$2,$3,$5,"+",".","gene_id \""$4"_fwd\";"}' > $resultName


    else

        cat "$bedFile1" "$bedFile2" | \
            sort -k1,1 -k2,2n | \
            sed 's/,peak_[0-9]*//g' | \
        awk 'BEGIN{OFS="\t";}{if($6=="+") print $1,"cage_peak","tss",$2,$3,$5,$6,".","gene_id \""$4"_fwd\";" ; else print $1,"cage_peak","tss",$2,$3,$5,$6,".","gene_id \""$4"_rev\";"}' > $resultName

    fi

}


