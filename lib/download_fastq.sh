#!/usr/bin/bash

source /misc/paras/data/programs/bashbone/latest/bashbone/activate.sh -r true -a "$@"

THREADS=56

sra_list=$1
fastq_dir=$2

mkdir -p "$fastq_dir"


# Get information if data is paired- and single-end sequenced 
# Source bashbone as our common installed efetch produces errors.
# bashbone version of efetch=19.0
# common version of efetch=16.2
collect_seq_protocol(){
    
    bashbone -c
    local sra_id=$1

    seq_protocol="$(efetch -db sra -id "$sra_id" -format runinfo -mode xml | \
        grep '<LibraryLayout>' | sed -nE 's/^\s*<([^>]+>)(.+)<\/\1/\2/p')"

    bashbone -c

    echo "$seq_protocol"

}


# zip data with pigz; advantaged parallization is possible
compression(){

    local sra_id=$1
    local seq_protocol=$2

    if [[ $seq_protocol == "PAIRED" ]]; then
        sra_ids=("${sra_id}"_1.fastq "${sra_id}"_2.fastq)
    else
        sra_ids=("${sra_id}".fastq)
    fi

    zip=("pigz" "-p" "$THREADS" "${sra_ids[@]}")

    echo "${zip[@]}"
    "${zip[@]}"

}

# Download data with aws because it is much faster than fasterq-dump
# After download the data is converted from sra to fastq with fasterq-dump
# and zipped with pigz
download_via_aws(){

    local sra_id=$1
    local fastq_dir=$2
    local seq_protocol=$3
    
    echo "Download $sra_id"
    aws s3 sync s3://sra-pub-run-odp/sra/"$sra_id" "$fastq_dir" --no-sign-request
    
    fasterq=/misc/paras/data/programs/sra-toolkit/latest/bin/fasterq-dump

    toFastq=("$fasterq" "./$sra_id" "--progres" "-e" "$THREADS")

    if [[ $seq_protocol == "PAIRED" ]]; then

        toFastq+=("--split-files")
    fi
    
    echo "${toFastq[@]}"

    cd "$fastq_dir"

    "${toFastq[@]}"
    
    compression "$sra_id" "$seq_protocol"

    rm "$sra_id"

    cd ../
    
}

# get info if paired or single end data
first_id=$(head -n 1 "$sra_list")
seq_protocol=$(collect_seq_protocol "$first_id")

# Download the fastq files
sra_id_list=()
while IFS= read -r line; do
    download_via_aws "$line" "$fastq_dir" "$seq_protocol"
    sra_id_list+=("$line")
done < "$sra_list"
