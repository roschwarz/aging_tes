#! /usr/bin/env bash
export RIPPCHEN=/misc/paras/data/programs/bashbone/latest/rippchen
source $RIPPCHEN/activate.sh -l true -c false -r true -a "$@"

###########
# Naming  #
###########

RAW_DATA=/misc/paras/data/rschwarz/projects/aging_tes/data/raw/rna_seq/female/20250728_965_RS/
RENAMED_DATA=/misc/paras/data/rschwarz/projects/aging_tes/data/raw/rna_seq/female/
PROCESSED_DATA=/misc/paras/data/rschwarz/projects/aging_tes/data/processed/rna_seq/deduplicated/female/

quality_check(){

    mkdir -p $RAW_DATA/fastqc/
    fastqc -t 40 -o $RAW_DATA/fastqc/ $(ls -f $RAW_DATA/*.fastq.gz)
    multiqc $RAW_DATA/fastqc/ -o $RAW_DATA

}
#quality_check

# Deduplicate based on UMI using the preprocess tool from the rippchen pipeline
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
    commander::qsubcmd -c bashbone -r -p threads -t 20 -i 10 -n dedub -o $PROCESSED_DATA -a cmds
}
#dedup

quality_check_dedup(){

    mkdir -p $PROCESSED_DATA/fastqc/
    fastqc -t 40 -o $PROCESSED_DATA/fastqc/ $(ls -f $PROCESSED_DATA/*.fastq.gz)
    multiqc $PROCESSED_DATA/fastqc/ -o $PROCESSED_DATA

}
#quality_check_dedup

detector_old(){

    #DETEctorEnv=/home/lakatos/rschwarz/.cache/pypoetry/virtualenvs/detector-w0DMliGL-py3.9/bin/activate
    DETEctorEnv=/misc/paras/data/rschwarz/coding/DETEctor/.venv/bin/activate
    DETEctor=/misc/paras/data/rschwarz/coding/DETEctor/DETEctor.py
    results=/misc/misc/data/rschwarz/projects/aging_tes/results/rna_seq/female/detector/
    threads=56

    mkdir -p $results
    
    cmds=("source $DETEctorEnv;\
          python3 $DETEctor map-salmon \
          -t $threads \
          -f $PROCESSED_DATA \
          -r mm10.TE.transcriptome.v102 \
          -l paired \
          -o $results")

    commander::printcmd -a cmds
    commander::qsubcmd -r -p threads -t "$threads" -i 1 -n sal_fem -o "$results" -a cmds

}
#detector_old


# Map the raw data with the rippchen pipeline from Konstantin
standard(){

    results_dir="/misc/paras/data/rschwarz/projects/aging_tes/results/rna_seq/female/rippchen/"

    mkdir -p $results_dir
    threads=56

    declare -a cmds

    for fi in $(find $RAW_DATA -type f -name "*R1.fastq.gz" | rev | cut -c 13- | rev | uniq); do
        fastqR1=${fi}_R1.fastq.gz
        fastqR2=${fi}_R2.fastq.gz
        fastqR3=${fi}_UMI.fastq.gz

        log=$(basename $fi .fastq).log

        commander::makecmd -a cmds -s '|' -c {COMMANDER[0]}<<- CMD
            $RIPPCHEN/scripts/rippchen.sh
                -v 2
                -t $threads
                -1 $fastqR1
                -2 $fastqR2
                -3 $fastqR3
                -g /misc/paras/data/genomes/GRCm38_v102/GRCm38.fa
                -gtf /misc/paras/data/genomes/GRCh38_v102/GRCh38.fa.gtf
                -a1 AGATCGGAAGAGC -a2 AGATCGGAAGAGC
                -o $results_dir
                -l $results_dir/logs/prep/resume_$log
                -tmp /data/tmp/rschwarz
                -no-bwa
                -no-star
                -no-stats
                -no-cor
                -skip md5
CMD
        done
        
        commander::printcmd -a cmds
        commander::qsubcmd -r -p threads -t $threads -i 4 -n female -o $results_dir -a cmds

}
standard
