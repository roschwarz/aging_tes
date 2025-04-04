#!/usr/bin/env bash

#=========== Get information about files and directories =======================
#
# The majority of the contained functions are specific for the mCQuaRna project.
#
# Provides helper function to easily get files and the respective group 
# assignment due to their file names in the cage data set. 

getDirectory(){

    root_dir=$1
    read_type=$2

    if [[ $read_type == "clipped" ]]; then
        dirs=( $(find "$root_dir" -maxdepth 1 -type d -regex '.*_gclipped$') )
    else
        dirs=( $(find "$root_dir" -maxdepth 1 -type d ! -regex '.*_gclipped$') )
    fi

    echo "${dirs[@]}"
}

getTissue(){

    file=$1

    file=$(basename $file)

    if echo $file | grep -q [Bb]lood; then
	    tissue=blood
    elif echo $file | grep -q '_Skin_\|skinI\.'; then
	    tissue=skinI
    elif echo $file | grep -q '[0-9]Skin_\|skinII\.'; then
	    tissue=skinII
    elif echo $file | grep -q [bB]rain; then
        if echo $file | grep -q down[sS]ampled; then
	        tissue=brain_downsampled
        else
            tissue=brain
        fi
    else
	    tissue="NA"
    fi

    echo $tissue

}

getAge(){


    file=$1

    file=$(basename "$file")


   if echo $file | grep -Eq "_[oO]_|_MO|_Mo_" ; then
       age='old'
   elif echo $file | grep -Eq "_[yY]_|_MY|_My_" ; then  
       age='young'
   else
       age='NA'
   fi       

   echo $age
}


getStrand(){
    
    file=$1

    file=$(basename $file)

    if echo $file | grep -Eq "forward|fwd"; then

	    strand="forward"

    elif echo $file | grep -Eq "reverse|rev"; then

	strand="reverse"
    else
	strand="NA"
    fi

    echo $strand


}

giveGroupMembers(){

    # takes a file path and the group definitions

    filesPath=$1
    tissue=$2
    age=$3
    strand=$4

    
    declare -a members

    if [[ ! ($tissue || $age || $strand) ]]; then
	    echo at least one condition is missing
    fi

    for file in "$filesPath"/*.bam; do

	    if [[ $(getTissue "$file") == "$tissue" && \
                $(getAge "$file") == "$age" && \
                $(getStrand "$file") == "$strand" ]]; then

	        members+=( "$file")

	    fi 
    done


    echo "${members[@]}"

}


giveGroupMembers_II(){

    # takes a file path and the group definitions

    filesPath=$1
    tissue=$2
    age=$3
    strand=$4

    #echo Search for files of $tissue, $age and $strand
    
    declare -a members

    if [ ! $tissue ] && [ ! $age ] && [ ! $strand ]; then
	    echo No condition provided, but you need at least one condition!
        return 1
    fi

    members=$(find "$filesPath" -maxdepth 1 -type f -name "*$age*$tissue*$strand.bam" -print0 | xargs -0 -I {} echo -n "{} ")

    echo "$members"

}

getConditions(){
    # Returns information about tissue, age and strand of files of that 
    # project.
    file=$1

    file=$(basename "$file")

    tissue=$(getTissue "$file")
    age=$(getAge "$file")
    strand=$(getStrand "$file")

    echo "$tissue"\;"$age"\;"$strand"
}

getUniqArrEntries(){
    #takes an array and returns unique entries of that array
	
    arr=$*

    uniques=$( echo "${arr[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' ')

    echo "$uniques"

}

