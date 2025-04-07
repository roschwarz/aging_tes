#!/usr/bin/env bash
source /misc/paras/data/programs/bashbone/latest/rippchen/activate.sh

FASTQ_PATH="../../data/processed/rna_seq/deduplicated/"
RESULTS_DIR="../../results/rna_seq/"

runSalmonTE(){
 
    tissue=$1
	fastqDir=$2/${tissue}/

    resultDir="$RESULTS_DIR"/"${tissue}"/alignment_SalmonTE
	mkdir -p "$resultDir"
    declare -a fastqs

	for files in $(ls "$fastqDir"/*R*.fastq.gz | rev | cut -c 12- | rev | uniq); do
    		
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
	commander::qsubcmd -r -p threads -t 56 -i 1 -n sal_"${tissue}" -o "$resultDir" -a cmds

}


for tissue in brain blood skinII; do
    runSalmonTE $tissue $FASTQ_PATH
done
