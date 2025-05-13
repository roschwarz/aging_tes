#! /usr/bin/env bash
export RIPPCHEN=/misc/paras/data/programs/bashbone/latest/rippchen
source $RIPPCHEN/activate.sh

###########
# Naming  #
###########

RAW_DATA=../../data/raw/rna_seq/female/20250226_929_RS/
RENAMED_DATA=../../data/raw/rna_seq/female/
PROCESSED_DATA=../../data/processed/rna_seq/deduplicated/female/

quality_check(){

    mkdir -p $RAW_DATA/fastqc/
    fastqc -t 40 -o $RAW_DATA/fastqc/ $(ls -f $RAW_DATA/*.fastq.gz)
    multiqc $RAW_DATA/fastqc/ -o $RAW_DATA

}
#quality_check

dedup(){

    declare -a cmds

    mkdir -p $PROCESSED_DATA

    for i in $(ls $RENAMED_DATA/*UMI.fastq.gz); do
        
        R1=$(echo $i | sed s/UMI/R1/)
        R2=$(echo $i | sed s/UMI/R2/)
        UMI=$i
        
        cmds+=("fq1=($R1); fq2=($R2); fq3=($UMI); preprocess::dedup -1 fq1 -2 fq2 -3 fq3 -o $PROCESSED_DATA -t 20")

        done

    commander::printcmd -a cmds
    commander::qsubcmd -r -p threads -t 20 -i 10 -n dedub -o $PROCESSED_DATA -a cmds
}
#dedup

quality_check_dedup(){

    mkdir -p $PROCESSED_DATA/fastqc/
    fastqc -t 40 -o $PROCESSED_DATA/fastqc/ $(ls -f $PROCESSED_DATA/*.fastq.gz)
    multiqc $PROCESSED_DATA/fastqc/ -o $PROCESSED_DATA

}
#quality_check_dedup
#

detector_old(){

    DETEctorEnv=/home/lakatos/rschwarz/.cache/pypoetry/virtualenvs/detector-w0DMliGL-py3.9/bin/activate
    results=../../results/rna_seq/female/
    threads=56

    mkdir -p $results
    
    cmds=("source $DETEctorEnv;\
          DETEctor.py quant \
          -t $threads \
          -f $PROCESSED_DATA \
          -r mm10.TE.transcriptome.v102 \
          -o $results")

    commander::printcmd -a cmds
    commander::qsubcmd -r -p threads -t "$threads" -i 1 -n sal_fem -o "$results" -a cmds

}
#detector_old
#
#

salmon_map(){

    results_dir="../../results/rna_seq/female/"

    mkdir -p $results_dir

    export salmon=/misc/paras/data/rschwarz/programs/salmon.v0.8.2/salmon
    export index=/misc/paras/data/rschwarz/common_data/salmonTE_indices/references/mm10.TE.transcriptome.v102
    export results="$results_dir"

    find $PROCESSED_DATA -type f -name "*.fastq.gz" | \
        sort | paste -d ' ' - - | \
        parallel --jobs 10 --colsep ' ' --env salmon --env index --env results \
        --joblog "$results_dir/parallel_joblog.txt" --bar \
        'output_dir=$(basename {1} _R1.fastq.gz); $salmon quant -l A -1 {1} -2 {2} -p 5 -i $index -o $results/${output_dir/}'

}
salmon_map
