#!/usr/bin/env bash

source /misc/paras/data/rschwarz/coding/biohacker/bash/bin/boneCalls.sh
source ../lib/cageFileHandler.sh
source ../lib/bamHandler.sh


raw=false
gclipped=true
mapping=true # map the data with rippchen (used mappers: segemehl and star)
stats=false
splitbams=false
mergebams=false

SIZE=50

callPeaks=false # call peaks with peakachu

if $raw; then

	echo Run the pipeline for raw data

	FASTQDIR=../../data/cage_RS/
	OUTPUTROOT=../../results/cage_RS/ 
    READS=raw

fi

if $gclipped; then

	echo Run the pipeline for G-clipped data

	FASTQDIR=../../data/processed/cage_seq/
	OUTPUTROOT=../../results/cage_seq_gclipped
    READS=clipped

    mkdir -p $OUTPUTROOT

fi

GENEANNOTATION_GTF=/misc/paras/data/genomes/GRCm38_v102/GRCm38.fa.aggregated.gtf

map_data(){

    echo Prepare and map that using rippchen.
    
    fastq_dir=$1
    read_type=$2

    cmds=()

    if [[ $read_type == "clipped" ]]; then
        dirs=( $(find "$fastq_dir" -maxdepth 1 -type d -regex '.*_gclipped$') )
    else
        dirs=( $(find "$fastq_dir" -maxdepth 1 -type d ! -regex '.*_gclipped$') )
    fi

    for tissue_dir in "${dirs[@]}"; do

        result_dir=$(echo "$tissue_dir" | sed "{s/data/results/}")

        for fi in "$tissue_dir"/*.fastq.gz; do

            tissue=$(getTissue "$fi")
            file_name=$(basename "$fi")
            log=$(basename "$fi" .fastq.gz).log

            commander::makecmd -a cmds -s '|' -c {COMMANDER[0]}<<- CMD
                rippchen.sh
                -v 2
                -o $result_dir
                -l $result_dir/logs/prep/$log
                -tmp /data/tmp/rschwarz -r
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
                -1 $fi
                -skip md5
		CMD

        done
    done

    result_dir=$(echo "$fastq_dir" | sed "{s/data/results/}")

    commander::printcmd -a cmds
    commander::qsubcmd -n cage_prep -o "$result_dir"/sge/prep -r -p threads -t 56 -a cmds
    

}

stats(){

    echo Doing stats for the mapped data.

	cmds=()

	fastq_dir=$1
    read_type=$2
    
    if [[ $read_type == "clipped" ]]; then
        dirs=( $(find "$fastq_dir" -maxdepth 1 -type d -regex '.*_gclipped$') )
    else
        dirs=( $(find "$fastq_dir" -maxdepth 1 -type d ! -regex '.*_gclipped$') )
    fi

    for tissue_dir in "${dirs[@]}"; do

        result_dir=$(echo "$tissue_dir" | sed "{s/data/results/}")

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

    result_dir=$(echo "$fastq_dir" | sed "{s/data/results/}")

	commander::printcmd -a cmds
	commander::qsubcmd -n stats -o $result_dir/sge/stats -r -p threads -t 56 -i 5 -a cmds

}


#====== MAPPING =====
#
# Data is mapped with the rippchen pipeline and subsequently some stats are
# calculated with the rippchen pipeline.
# 
if $mapping; then

    echo Map reads with segemehl using rippchen.
    
    map_data "$FASTQDIR" "$READS"

fi


if $stats; then

    echo Doing stats using rippchen.
	
    stats "$FASTQDIR" "$READS"
    
fi


if $splitbams; then

    echo Split bams into forward and reverse bam files.
    split_bams "$OUTPUTROOT" "$READS"

fi

if $mergebams; then

    echo Merge bams of same tissue, age, and strand.
    merge_bams "$OUTPUTROOT" "$READS"

fi

callPeaks(){

    echo Call peaks with peakachu using bashbone functionality.
    OUTPUTROOT=$1
    read_type=$2
    
    
    #mkdir -p $result
    #
    for data_set in $(getDirectory "$OUTPUTROOT" "$read_type"); do

        peak_res=$data_set/raw_peaks/

        for strand in forward reverse; do
            mapper=(segemehl);
            mkdir -p "$peak_res"
            tissue=$(basename "$data_set" _gclipped)
            bamFiles=$(find "$data_set"/bams_stranded_segemehl/*"$strand".bam -type f )
            segemehl=(${bamFiles[@]})
            echo peaks::peakachu -t 56 -f 150 -r mapper -o "$peak_res" -z true
            peaks::peakachu -t 56 -f 150 -r mapper -o "$peak_res" -z true
        done

    done
   
   # for tissue in _Skin_ [0-9]Skin_ [bB]rain Blood; do
   # for strand in forward reverse; do
   #     for age in [oO] [yY]; do
   # 	bamFiles=$(collect_bams_grouped $SPLITBAMS_OUTPUT $tissue $age $strand)
   # 		segemehl=($bamFiles)
   # 		echo peaks::peakachu -t 56 -f 150 -r mapper -o $result -z true
   # 		peaks::peakachu -t 20 -f 150 -r mapper -o $result -z true
   #     done
   # done
   # done


}

if $callPeaks; then
	
	echo Call peaks with peakachu and merge called peaks that are within "$SIZE"

	#callPeaks "$OUTPUTROOT" "$READS"


    for data_set in $(getDirectory "$OUTPUTROOT" "$read_type"); do

        raw_peaks=$data_set/raw_peaks/segemehl/peakachu
        echo $data_set
	    #merge_narrowPeaks $raw_peaks $SIZE

    done
    
fi
