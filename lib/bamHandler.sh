#!/usr/bin/env bash

#source /misc/paras/data/rschwarz/coding/biohacker/bash/bin/boneCalls.sh
source /misc/paras/data/rschwarz/projects/mCQuaRna/src/lib/cageFileHandler.sh

###########################
# Separation of bam files #
###########################

splitBams(){

	cmds=()

	bamDir=$1
	result=$2
	OUTPUTROOT=$3

	for f in "$bamDir"/*.unique.sorted.bam ; do

        f_name=$(basename "$f" .unique.sorted.bam)
        mkdir -p "$result"
        
        commander::makecmd -a cmds -s '|' -c {COMMANDER[0]}<<- CMD
                samtools view -F 16 $f -@ 56 -h -b > $result/$f_name.forward.bam;
                samtools view -f 16 $f -@ 56 -h -b > $result/$f_name.reverse.bam;
                samtools index $result/$f_name.forward.bam  $result/$f_name.forward.bam.bai;
                samtools index $result/$f_name.reverse.bam  $result/$f_name.reverse.bam.bai;
CMD

	done

	commander::printcmd -a cmds
	commander::qsubcmd -n splitBams -o $OUTPUTROOT/sge/splitBams -r -w -p threads -t 56 -a cmds 
	
}

merge_bams_universal(){

    # merges provided bam files and creates an index of the merged bam file
    #
    # @paramters
    # <- string of bam files that will be merged
    # <- name of the result file including path
    #
    # @result
    # -> stores a merged bam file
    # -> stores an index of the merged bam file
    #
    # !!! Don't put $bam_files in quotes, it will work !!!
    bam_files=$1
    output=$2

    echo merge following files "$bam_files".
    samtools merge -@ 56 $output $bam_files

    echo create index for "$output"
    samtools index -@ 56 "$output"

}


########################
# Merging of bam files #
########################

merge_bams_age(){

    dataSet=$1 
    result_dir=$2 
    tissue=$3

    for file in "$dataSet"/*.bam; do
        age+=($(getAge "$file"))
        strand+=($(getStrand "$file"))
    done


    ages=$(getUniqArrEntries "${age[@]}")
    strands=$(getUniqArrEntries "${strand[@]}")

    for age in $ages; do

        for strand in $strands; do

            bam_files=$(giveGroupMembers "$dataSet" "$tissue" "$age" "$strand")
            result_file="$result_dir"/"${tissue}"."${age}"."${strand}".bam

            merge_bams_universal "$bam_files" "$result_file"

        done
    done

}

tissue_dict(){

    tissue=$1

    if [ "$tissue" == "brain" ]; then
        echo "[bB]rain"
    elif [ "$tissue" == "skinI" ]; then
        echo "_Skin_"
    elif [ "$tissue" == "brain_downSampled" ]; then
        echo "brain_downsampled"
    elif [ "$tissue" == "skinII" ]; then
        echo "[0-9]Skin_"
    else
        echo "blubb"
    fi


}

merge_bams_tissue(){

    dataSet=$1 
    result_dir=$2 
    tissue=$3

    for strand in forward reverse; do

        bam_files=$(giveGroupMembers_II "$dataSet" "" "" "$strand")
        result_file="$result_dir"/"${tissue}"."${strand}".bam
        merge_bams_universal "$bam_files" "$result_file"
        

    done

    return 0

}


merge_bams(){
    #
    # Merges .bam-files of same tissue, age, and strand.
    # Merges .bam_files of same tissue and strand where the age is not 
    # considered.
    # 

    rippchen_res=$1 
    read_type=$2

    dirs=$(getDirectory "$rippchen_res" "$read_type")
    
    for dataSet in $dirs; do

        tissue=$(basename "$dataSet" _gclipped)
        result_dir="$dataSet"/bams_merged_segemehl
        dataSet="$dataSet"/bams_stranded_segemehl

        mkdir -p "$result_dir"         

        #merge_bams_age "$dataSet" "$result_dir" "$tissue"

        echo ==== Merge bamfiles of $tissue ====
        merge_bams_tissue "$dataSet" "$result_dir" "$tissue"

    done
}

