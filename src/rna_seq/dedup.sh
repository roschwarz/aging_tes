#! /usr/bin/env bash
export RIPPCHEN=/misc/paras/data/programs/bashbone/latest/rippchen
source $RIPPCHEN/activate.sh

tissue="skinII"

dedup(){

    declare -a cmds

    for i in $(ls $tissue/*UMI.fastq.gz); do
	
	R1=$(echo $i | sed s/UMI/R1/)
	R2=$(echo $i | sed s/UMI/R2/)
	UMI=$i
	
	cmds+=("fq1=($R1); fq2=($R2); fq3=($UMI); preprocess::dedup -1 fq1 -2 fq2 -3 fq3 -o $tissue/deduped -t 20 -p /data/tmp/rschwarz/")

    done

    commander::qsubcmd -r -p threads -t 20 -i 10 -n dedub -o $tissue/deduped -a cmds -w
}
#dedup

counter(){
    
    echo -e "file\toriginal\tdeduped\tper" 
    echo -e "file\toriginal\tdeduped\tper" > skinII/counts.tsv

    for i in $(ls skinII/*R1.fastq.gz); do
	
	de=$(echo $i | sed 's/skinII/skinII\/deduped/')	
	original=$(echo $(zcat $i | wc -l)/4 | bc)
	deduped=$(echo $(zcat $de | wc -l)/4 | bc)
	per=$(python -c "print($deduped.0/$original.0)")

	echo -e "$i\t$original\t$deduped\t$per" 
	echo -e "$i\t$original\t$deduped\t$per" >> skinII/counts.tsv

    done
}
#counter

