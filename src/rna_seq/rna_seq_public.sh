#!/usr/bin/env bash

# Update path to ../../data/public/PMCID_PMC11529319/male/ when you want to run the script for the male data
RAW_DATA=../../data/public/PMCID_PMC11529319/female/

quality_check_dedup(){

    mkdir -p $RAW_DATA/fastqc/
    fastqc -t 40 -o $RAW_DATA/fastqc/ $(ls -f $RAW_DATA/*.fastq.gz)
    multiqc $RAW_DATA/fastqc/ -o $RAW_DATA

}
#quality_check_dedup

salmon_map(){

    results_dir="../../results/rna_seq/public/PMCID_PMC11529319/"

    mkdir -p $results_dir

    export salmon=/misc/paras/data/rschwarz/programs/salmon.v0.8.2/salmon
    export index=/misc/paras/data/rschwarz/common_data/salmonTE_indices/references/mm10.TE.transcriptome.v102
    export results="$results_dir"

    find $RAW_DATA -type f -name "*.fastq.gz" | \
        sort | paste -d ' ' - - | \
        parallel --jobs 5 --colsep ' ' --env salmon --env index --env results \
        --joblog "$results_dir/parallel_joblog.txt" --bar \
        'output_dir=$(basename {1} _R1.fastq.gz); $salmon quant -l A -1 {1} -2 {2} -p 10 -i $index -o $results/${output_dir/}'

}
#salmon_map
#
#
merge_tables(){
    
    echo "merge count tables"
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
    ' ../../results/rna_seq/public/PMCID_PMC11529319/female/*/quant.sf > ../../results/rna_seq/public/PMCID_PMC11529319/female/EXPR.csv

}
merge_tables
