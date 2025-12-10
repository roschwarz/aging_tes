#!/usr/bin/env bash

source ../lib/cageFileHandler.sh
source ../lib/functions.sh

SIZE=50

echo Call peaks for quant data

OUTPUTROOT=../../results/quant_seq/


callPeaks(){

    # Call peaks using peakachu that is a function within bashbone.
    
    echo Call peaks with peakachu using bashbone functionality.

    source /misc/paras/data/programs/bashbone/latest/bashbone/activate.sh -l true -c false

    OUTPUTROOT=$1
    
    dirs=( $(find "$OUTPUTROOT" -mindepth 1 -maxdepth 1 -type d \( -name "brain" -o -name "blood" -o -name "skinI" -o -name "skinII" \)))

    for data_set in "${dirs[@]}"; do

        echo Call peaks for $data_set
        peak_res=$data_set/raw_peaks/

        echo $peak_res

        for strand in forward reverse; do
            mapper=(segemehl);
            mkdir -p "$peak_res"
            tissue=$(basename "$data_set")
            bamFiles=$(find "$data_set"/bams_stranded_segemehl/*"$strand".bam -type f )
            segemehl=(${bamFiles[@]})
            bamIndices=(${bamIndex[@]})
            echo peaks::peakachu -t 56 -f 150 -r mapper -o "$peak_res" -z true
            peaks::peakachu -t 56 -f 150 -r mapper -m 2000000 -g /misc/paras/data/genomes/GRCm38_v102/GRCm38.fa -o "$peak_res" -z true
        done

    done
   
}


merge_narrowPeaks(){
# merge peaks that are within a certain proximity (size).
# The input are narrowPeak files created with peakachu. Be careful, this works 
# only when the peaks are called in strand specific manner. 

	local OUTPUTROOT=$1    
	local size=$2

    echo Merge peaks located within "$size"
    mapfile -t dirs < <(find "$OUTPUTROOT" -mindepth 1 -maxdepth 1 -type d \( -name "brain" -o -name "blood" -o -name "skinI" -o -name "skinII" \))

    for data_set in "${dirs[@]}"; do

        peak_res=$data_set/raw_peaks/segemehl/peakachu/

        for peaks in $peak_res/*.narrowPeak; do

            strand=$(getStrand $peaks)
            tissue=$(getTissue $peaks)
            output="$tissue"."$strand".bed

            [[ $(getStrand $peaks) == 'forward' ]] && strand="+" || strand="-"

            bedtools merge -i $peaks -d $size -c 4,5,6 -o distinct,distinct,distinct | \
            sed "s/\./$strand/" > $peak_res/${output}

        done


    done

}

create_peak_gtf(){

    # IMPORTANT: The peaks are called on the opposite strand due to the usage
    # of function that I usually used to call the cage peaks. However, the 
    # quant data was sequenced in the opposite direction. Therefore, I 
    # implemented I switch the strand in the gtf file, but for the moment
    # they are not switched in the bed files!!!
	local OUTPUTROOT=$1    

    mapfile -t dirs < <(find "$OUTPUTROOT" -mindepth 1 -maxdepth 1 -type d \( -name "brain" -o -name "blood" -o -name "skinI" -o -name "skinII" \))
    for data_set in "${dirs[@]}"; do

        
        tissue=$(basename "$data_set")
        peak_res=$data_set/raw_peaks/segemehl/peakachu/

        echo Create a gtf file for QUANT-peaks of $tissue.

        bedFile1=$peak_res/$tissue.forward.bed
        bedFile2=$peak_res/$tissue.reverse.bed
        result=$peak_res/$tissue.quant_peaks.gtf       

        peaksBedtoGTF $bedFile1 $bedFile2 $result quant

    done
}


count_peak_reads(){

	local OUTPUTROOT=$1    

    local threads=10

    source /misc/paras/data/programs/bashbone/latest/bashbone/activate.sh

    mapfile -t dirs < <(find "$OUTPUTROOT" -mindepth 1 -maxdepth 1 -type d \( -name "brain" -o -name "blood" -o -name "skinI" -o -name "skinII" \))

    for data_set in "${dirs[@]}"; do
        
        count_res=$data_set/counts
        bam_dir=$data_set/mapped/segemehl/
        tissue=$(basename "$data_set")

        gtfFile=$data_set/raw_peaks/segemehl/peakachu/$tissue.quant_peaks.gtf

        mkdir -p $count_res
        
        mapfile -t bamFiles < <(find $bam_dir -type f -name "*.unique.sorted.bam")
        
        commander::makecmd -a cmds -s '|' -c {COMMANDER[0]}<<- CMD
            featureCounts \
            -a $gtfFile \
            -s 2 \
            -o $count_res/${tissue}_peak_counts.csv \
            -T $threads \
            -t tts \
            ${bamFiles[@]}
		CMD

    done

    commander::printcmd -a cmds
    commander::qsubcmd -n count_reads -o "$OUTPUTROOT"/sge/counting -r -p threads -t $threads -a cmds

}

echo Call peaks with peakachu and merge called peaks that are within "$SIZE"
#callPeaks "$OUTPUTROOT"
#merge_narrowPeaks "$OUTPUTROOT" "$SIZE"
#create_peak_gtf "$OUTPUTROOT"
count_peak_reads $OUTPUTROOT
