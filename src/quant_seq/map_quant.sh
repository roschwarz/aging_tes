#!/usr/bin/env bash

source /misc/paras/data/programs/bashbone/latest/rippchen/activate.sh -l true -c false 
source ../lib/cageFileHandler.sh
source ../lib/bamHandler.sh

mapping=false # map the data with rippchen (used mappers: segemehl and star)
stats=false
splitbams=false
mergebams=true

echo Map quant data

FASTQDIR=../../data/raw/quant_seq/
OUTPUTROOT=../../results/quant_seq/ 

mkdir -p $OUTPUTROOT 

map_data(){

    echo Prepare and map that using rippchen.
    
    local fastq_dir=$1
    local OUTPUTROOT=$2
    cmds=()

    dirs=( $(find "$fastq_dir" ! -path "$fastq_dir" -type d) )

    mkdir -p $OUTPUTROOT

    for tissue_dir in "${dirs[@]}"; do

        tissue=$(basename $tissue_dir)

        result_dir=$OUTPUTROOT/$tissue 

        mkdir -p $result_dir

        for file in "$tissue_dir"/*.fastq.gz; do

            tissue=$(getTissue "$file")
            file_name=$(basename "$file")
            log=$(basename "$file" .fastq.gz).log

            commander::makecmd -a cmds -s '|' -c {COMMANDER[0]}<<- CMD
                rippchen.sh
                -v 2
                -o $result_dir
                -l $result_dir/logs/prep/$log
                -tmp /data/tmp/rschwarz
                -t 56
                -xmem 200000
                -g /misc/paras/data/genomes/GRCm38_v102/GRCm38.fa
                -gtf /misc/paras/data/genomes/GRCm38_v102/GRCm38.fa.gtf
                -a1 CTGTCTCTTATACACATCT,AGATCGGAAGAGC
                -no-rrm
                -no-stats
                -no-split
                -no-bwa
                -no-star
                -no-quant
                -1 $file
                -skip md5
		CMD

        done
    done

    commander::printcmd -a cmds
    commander::qsubcmd -n quant_map -o "$result_dir"/sge/mapping -r -p threads -t 56 -a cmds
    

}

stats(){

    echo Doing stats for the mapped data.

    
    local fastq_dir=$1
    local OUTPUTROOT=$2
    cmds=()

    dirs=( $(find "$fastq_dir" ! -path "$fastq_dir" -type d) )

    for tissue_dir in "${dirs[@]}"; do

        tissue=$(basename $tissue_dir)

        result_dir=$OUTPUTROOT/$tissue 

	    ls -v "$tissue_dir"/*.fastq.gz > "$result_dir"/R1.list

        commander::makecmd -a cmds -s '|' -c {COMMANDER[0]}<<- CMD
                rippchen.sh
                -v 2
                -o $result_dir
                -l $result_dir/logs/stats/stat.log
                -tmp /ssd/tmp/rschwarz -r
                -t 56
                -xmem 200000
                -g /misc/paras/data/genomes/GRCm38_v102/GRCm38.fa
                -gtf /misc/paras/data/genomes/GRCm38_v102/GRCm38.fa.gtf
                -a1 CTGTCTCTTATACACATCT,AGATCGGAAGAGC
                -no-rrm
                -no-split
                -no-bwa
                -no-star
                -no-quant
                -1 $result_dir/R1.list
                -redo stats 
                -skip md5
CMD

    done

	commander::printcmd -a cmds
	commander::qsubcmd -n stats -o "$OUTPUTROOT"/sge/stats -r -p threads -t 56 -i 5 -a cmds


}

split_bams(){
    # separates bam files into forward and reverse strand    
    
    local OUTPUTROOT=$1
    
    dirs=( $(find "$OUTPUTROOT" -mindepth 1 -maxdepth 1 -type d \( -name "brain" -o -name "blood" -o -name "skinI" -o -name "skinII" \)))

    for rippchen_dir in "${dirs[@]}"; do

        echo "$rippchen_dir"

        result_dir="$rippchen_dir"/bams_stranded_segemehl
        bam_dir="$rippchen_dir"/mapped/segemehl/

        splitBams "$bam_dir" "$result_dir" "$OUTPUTROOT"

    done

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
    OUTPUTROOT=$1 
    #read_type=$2

    #dirs=$(getDirectory "$rippchen_res" "$read_type")
    
    dirs=( $(find "$OUTPUTROOT" -mindepth 1 -maxdepth 1 -type d \( -name "brain" -o -name "blood" -o -name "skinI" -o -name "skinII" \)))
    
    for dataSet in ${dirs[@]}; do

        tissue=$(basename "$dataSet" _gclipped)
        result_dir="$dataSet"/bams_merged_segemehl
        dataSet="$dataSet"/bams_stranded_segemehl

        mkdir -p "$result_dir"         

        echo ==== Merge bamfiles of $tissue ====
        merge_bams_tissue "$dataSet" "$result_dir" "$tissue"

        echo ==== Merge bamfiles of $tissue AND age ====
        echo merge_bams_age "$dataSet" "$result_dir" "$tissue"


    done
}



#====== MAPPING =====
#
# Data is mapped with the rippchen pipeline and subsequently some stats are
# calculated with the rippchen pipeline.
# 
if $mapping; then

    echo Map reads with segemehl using rippchen.
    map_data "$FASTQDIR" "$OUTPUTROOT"

fi


if $stats; then

    echo Doing stats using rippchen.
    stats "$FASTQDIR" "$OUTPUTROOT"
    
fi

if $splitbams; then

    echo Split bams into forward and reverse bam files.
    split_bams "$OUTPUTROOT"

fi

if $mergebams; then

    echo Merge bams of same tissue, age, and strand.
    merge_bams "$OUTPUTROOT" "$READS"

fi
