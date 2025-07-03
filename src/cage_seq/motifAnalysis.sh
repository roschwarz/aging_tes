#!/usr/bin/env bash

source ../lib/cageFileHandler.sh

OUTPUTROOT=../../results/cage_seq_gclipped/
READS=clipped

ChromSizes=../../data/shared/GRCm38.p6.chrom.sizes
GenomeFA=/misc/paras/data/genomes/GRCm38_v102/GRCm38.fa
path=./data/motif/



getPromoter(){

    local TSS=$1 
    local ChromSizes=$2
    local basePairs=$3
    local promoter=$4
    

    bedtools slop \
            -i  $TSS \
            -g $ChromSizes \
            -s \
            -l $basePairs \
            -r 0 > $promoter


}

mergeTss(){

    local path=$1
    local tissue=$2
    local result=$3
        
    cat <(sed 's/,peak_[0-9]*//g' $path/$tissue.forward.tss.bed | awk -v OFS="\t" '{$4 = $4"_fwd"; print $0}') \
                <(sed 's/,peak_[0-9]*//g' $path/$tissue.reverse.tss.bed | awk  -v OFS="\t" '{$4 = $4"_rev"; print $0}') \
                > $result



}

createPromoter(){

    # Define the promoter region of iTEs (TEs that intersect with a CAGE peak)
	local OUTPUTROOT=$1    
    local read_type=$2
    local ChromSizes=$3
    local basePairs=$4

    for data_set in $(getDirectory "$OUTPUTROOT" "$read_type"); do

        tissue=$(basename "$data_set" _gclipped)

        if [[ $tissue != "brain" && $tissue != "skinI" ]]; then

            tssDir=$data_set/tss
            tssFile=$tssDir/$tissue.tss.bed

            promoterFile=$(echo "$tssFile" | sed 's/tss.bed/promoter.bed/')
            independentTERegion=../../data/shared/${tissue}_independent_TE_regions.bed  
            tssIndependentTERegion=../../data/shared/${tissue}_tss_independent_TE_regions.bed  
            promoterIndependentTERegion=../../data/shared/${tissue}_promoter_independent_TE_regions.bed  
            echo XXXXX "$independentTERegion"
            echo Merge forward and reverse ctss into one bed \("$tssFile"\)
            mergeTss "$tssDir" "$tissue" "$tssFile"

            echo Extract coordinates of promoter for all called tss \("$promoterFile"\)
            getPromoter "$tssFile" "$ChromSizes" "$basePairs" "$promoterFile"


            echo Get tss of independently expressed TE regions \("$tssIndependentTERegion"\)
            bedtools intersect -a "$independentTERegion" \
                -b <(sort -k1,1 -k2,2n "$tssFile") -s > "$tssIndependentTERegion"


            echo Extract promtoer coordinates of independently expressed TE regions \("$tssIndependentTERegion"\)
            promoterFile=$(echo "$tssFile" | sed 's/tss.bed/promoter.bed/')
            getPromoter "$tssIndependentTERegion" "$ChromSizes" "$basePairs" "$promoterIndependentTERegion"
        fi

    done



}

findMotifs(){

    # apply findMotifsGenome of the HOMER suite to identifiy potential TFBS 
    # motifs in independently expressed TE regions
	local OUTPUTROOT=$1    
    local read_type=$2
    local GenomeFA=$3


    echo Search for motifs

    for data_set in $(getDirectory "$OUTPUTROOT" "$read_type"); do

        tissue=$(basename "$data_set" _gclipped)

        if [[ $tissue != "brain" && $tissue != "skinI" ]]; then

            promoterIndependentTERegion=../../data/shared/${tissue}_promoter_independent_TE_regions.bed  
            result=$data_set/homer
            mkdir -p "$result"

            echo Run Homer for "$tissue"
            findMotifsGenome.pl "$promoterIndependentTERegion" \
                    "$GenomeFA" \
                    "$result" \
                    -p 56
        fi
    done

}

collectMotifCoordinates(){

    # apply scanMotifGenomeWide to collect the coordinates of the called motifs
	local OUTPUTROOT=$1    
    local read_type=$2
    local GenomeFA=$3
    local numberOfMotifs=$4


    echo Search for motifs

    for data_set in $(getDirectory "$OUTPUTROOT" "$read_type"); do

        
        tissue=$(basename "$data_set" _gclipped)

        if [[ $tissue != "brain" && $tissue != "skinI" ]]; then


            for i in $(seq 1 "$numberOfMotifs"); do 

                motifFile=$data_set/homer/homerResults/motif$i.motif 
                motifName=$(grep '^>'  "$motifFile" | awk -F":" '{print $2}' | awk -F"/" '{print $1}')
                resultFile=$data_set/homer/${motifName}.coordinates.bed

                scanMotifGenomeWide.pl \
                    "$motifFile" "$GenomeFA" -p 30 | \
                    awk -v OFS='\t' '{print $2,$3,$4,$7,$6,$5}' | \
                    sort -k1,1 -k2,2n > "$resultFile"
            done
        fi
    done
     


}

motifTEIntersection_brain(){
    
    # Intersecting TSS and Sox5 motif with iTEs and their flanking regions
    echo intersect sox5 motif with independently expressed TEs
    local ChromSizes=$1
    results=../../results/cage_RS/brain_downsampled_gclipped/annotations

    mkdir -p "$results"

    tss=../../results/cage_RS/brain_downsampled_gclipped/tss/brain_downsampled.tss.bed
    sox5Coordinates=../../results/cage_RS/brain_downsampled_gclipped/homer/Sox5.coordinates.bed

    independentTERegion=../../data/shared/brain_downsampled_independent_TE_regions.bed
    up_flank="$results"/brain_downsampled_TE_region_up_stream.bed
    down_flank="$results"/brain_downsampled_TE_region_down_stream.bed
    regions=("$independentTERegion" "$up_flank" "$down_flank")


    bedtools flank -i "$independentTERegion" -g "$ChromSizes" -l 500 -r 0 -s > "$up_flank"
    bedtools flank -i "$independentTERegion" -g "$ChromSizes" -r 500 -l 0 -s > "$down_flank"

    # Intersect TE.region and flanks with TSS
    for file in "${regions[@]}"; do
            res=$(basename "$file" bed)tss.intersection.bed
            bedtools intersect -a "$file" -b "$tss" -s -wa -wb |
                    awk -v OFS="," 'BEGIN{print "chromosome","region.start","region.end","region.id","strand","feature.start","feature.end","feature.id","region.width"}; {print $1,$2,$3,$4,$6,$8,$9,$10,$3-$2}'\
                    > "$results"/"$res"
    done

    # Intersect TE.region with sox5 motif
    for file in "${regions[@]}"; do
            res=$(basename "$file" bed)sox.intersection.bed
            bedtools intersect -a "$file" -b "$sox5Coordinates" -s -wa -wb | \
            awk -v OFS='\t' '{if($12=="-") $8=($9-1); else $9=($8+1); print $0}' |
            awk -v OFS="," 'BEGIN{print "chromosome","region.start","region.end","region.id","strand","feature.start","feature.end","feature.id","region.width"}; {print $1,$2,$3,$4,$6,$8,$9,$10,$3-$2}'\
            > "$results"/"$res"
    done


}

motifTEIntersection_skin(){
    
    # Intersecting TSS and Sox5 motif with iTEs and their flanking regions
    echo intersect sox10 motif with independently expressed TEs
    local ChromSizes=$1
    results=../../results/cage_RS/skinII_gclipped//annotations

    mkdir -p "$results"

    tss=../../results/cage_RS/skinII_gclipped/tss/skinII.tss.bed
    sox10Coordinates=../../results/cage_RS/skinII_gclipped/homer/SOX10.coordinates.bed

    independentTERegion=../../data/shared/skinII_independent_TE_regions.bed
    up_flank="$results"/skin_TE_region_up_stream.bed
    down_flank="$results"/skin_TE_region_down_stream.bed
    regions=("$independentTERegion" "$up_flank" "$down_flank")


    bedtools flank -i "$independentTERegion" -g "$ChromSizes" -l 500 -r 0 -s > "$up_flank"
    bedtools flank -i "$independentTERegion" -g "$ChromSizes" -r 500 -l 0 -s > "$down_flank"

    # Intersect TE.region and flanks with TSS
    for file in "${regions[@]}"; do
            res=$(basename "$file" bed)tss.intersection.bed
            bedtools intersect -a "$file" -b "$tss" -s -wa -wb |
                    awk -v OFS="," 'BEGIN{print "chromosome","region.start","region.end","region.id","strand","feature.start","feature.end","feature.id","region.width"}; {print $1,$2,$3,$4,$6,$8,$9,$10,$3-$2}'\
                    > "$results"/"$res"
    done

    # Intersect TE.region with sox10 motif
    for file in "${regions[@]}"; do
            res=$(basename "$file" bed)sox.intersection.bed
            bedtools intersect -a "$file" -b "$sox10Coordinates" -s -wa -wb | \
            awk -v OFS='\t' '{if($12=="-") $8=($9-1); else $9=($8+1); print $0}' |
            awk -v OFS="," 'BEGIN{print "chromosome","region.start","region.end","region.id","strand","feature.start","feature.end","feature.id","region.width"}; {print $1,$2,$3,$4,$6,$8,$9,$10,$3-$2}'\
            > "$results"/"$res"
    done


}

motifTEIntersection_blood(){
    
    # Intersecting TSS and Sox5 motif with iTEs and their flanking regions
    echo intersect sox4 motif with independently expressed TEs
    local ChromSizes=$1
    results=../../results/cage_RS/blood_gclipped/annotations

    mkdir -p "$results"

    tss=../../results/cage_RS/blood_gclipped/tss/blood.tss.bed
    sox4Coordinates=../../results/cage_RS/blood_gclipped/homer/RBP1.coordinates.bed

    independentTERegion=../../data/shared/blood_independent_TE_regions.bed
    up_flank="$results"/blood_TE_region_up_stream.bed
    down_flank="$results"/blood_TE_region_down_stream.bed
    regions=("$independentTERegion" "$up_flank" "$down_flank")


    bedtools flank -i "$independentTERegion" -g "$ChromSizes" -l 500 -r 0 -s > "$up_flank"
    bedtools flank -i "$independentTERegion" -g "$ChromSizes" -r 500 -l 0 -s > "$down_flank"

    # Intersect TE.region and flanks with TSS
    for file in "${regions[@]}"; do
            res=$(basename "$file" bed)tss.intersection.bed
            bedtools intersect -a "$file" -b "$tss" -s -wa -wb |
                    awk -v OFS="," 'BEGIN{print "chromosome","region.start","region.end","region.id","strand","feature.start","feature.end","feature.id","region.width"}; {print $1,$2,$3,$4,$6,$8,$9,$10,$3-$2}'\
                    > "$results"/"$res"
    done

    # Intersect TE.region with sox4 motif
    for file in "${regions[@]}"; do
            res=$(basename "$file" bed)sox.intersection.bed
            bedtools intersect -a "$file" -b "$sox4Coordinates" -s -wa -wb | \
            awk -v OFS='\t' '{if($12=="-") $8=($9-1); else $9=($8+1); print $0}' |
            awk -v OFS="," 'BEGIN{print "chromosome","region.start","region.end","region.id","strand","feature.start","feature.end","feature.id","region.width"}; {print $1,$2,$3,$4,$6,$8,$9,$10,$3-$2}'\
            > "$results"/"$res"
    done


}

createPromoter $OUTPUTROOT $READS $ChromSizes 500
#findMotifs $OUTPUTROOT $READS $GenomeFA
#collectMotifCoordinates "$OUTPUTROOT" "$READS" "$GenomeFA" 5
#motifTEIntersection_brain "$ChromSizes"
#motifTEIntersection_skin "$ChromSizes"
#motifTEIntersection_blood "$ChromSizes"
