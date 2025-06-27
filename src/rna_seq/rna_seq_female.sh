#! /usr/bin/env bash
export RIPPCHEN=/misc/paras/data/programs/bashbone/latest/rippchen
source $RIPPCHEN/activate.sh -l true -c false -r true -a "$@"

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

salmon_map(){

    results_dir="../../results/rna_seq/female/skin_dup"
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
#salmon_map


merge_tables(){
    
    for tissue in skin; do #, brain; do

        awk '
        FNR==1 {
            file_num++
            if (file_num==1) header="TE"
            header = header "\t" FILENAME
            next
        }
        {
            key=$1
            value=$NF
            data[key,file_num]=value
            keys[key]=1
        }
        END {
            print header
            for (k in keys) {
                line = k
                for (i=1; i<=file_num; i++) {
                    line = line "\t" ( (k,i) in data ? data[k,i] : "NA")
                }
                print line
            }
        }
        ' ../../results/rna_seq/female/skin_dup/*/quant.sf > ../../results/rna_seq/female/skin_dup/EXPR.csv # ../../results/rna_seq/female/*${tissue}*/quant.sf > ../../results/rna_seq/female/EXPR_${tissue}.csv


    done
}
#merge_tables


salmonTE_map(){


    results_dir="../../results/rna_seq/female/salmonTE/"

    mkdir -p $results_dir
    declare -a fastqs

	for files in $(ls "$PROCESSED_DATA"*R*.fastq.gz | rev | cut -c 12- | rev | uniq); do
    		
	    fastqs+="${files}R1.fastq.gz ${files}R2.fastq.gz "
	    
	done

	    cmds=("cd $resultDir/; source ~/programs/conda/bin/activate salmonTE; \
		/misc/paras/data/rschwarz/programs/SalmonTE/SalmonTE.py quant \
		    --reference=mm10.TE.transcriptome.v102 \
		    --outpath=$resultDir \
		    --num_threads=56 \
            --exprtype=count \
    		    $fastqs")
	
    commander::printcmd -a cmds
	commander::qsubcmd -r -p threads -t 56 -i 1 -n sal_"${tissue}" -o "$results_dir" -a cmds
}
#salmonTE_map


# Map the raw data with the rippchen pipeline from Konstantin
standard(){

    results_dir="../../results/rna_seq/female/rippchen"

    mkdir -p $results_dir
    threads=56

    declare -a cmds

    for fi in $(find $RAW_DATA -type f -name "*skin*R1.fastq.gz" | rev | cut -c 13- | rev | uniq); do
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
