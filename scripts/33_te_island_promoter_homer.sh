

canonicalTEIslands(){
    bed_path=../../data/shared/
    bed_suffix=_promoter_can_TE_island.bed
    genome=/misc/paras/data/genomes/GRCm38_v102/GRCm38.fa


    for tissue in brain skin blood; do

        bed_file=$bed_path$tissue$bed_suffix
        output=../../results/te_islands/"$tissue"_homer/

        echo findMotifsGenome.pl $bed_file $genome $output -p 56 -size 200     
        findMotifsGenome.pl $bed_file $genome $output -p 56 -size 200
    done
}

down_regulated_Island(){

    bed_path=../../results/TEItx/homer/
    bed_suffix=_promoter_top500_down.bed
    genome=/misc/paras/data/genomes/GRCm38_v102/GRCm38.fa

    for tissue in brain skin blood; do

        bed_file=$bed_path$tissue$bed_suffix
        output=../../results/TEItx/homer/"$tissue"_homer/

        echo findMotifsGenome.pl $bed_file $genome $output -p 56 -size 200     
        findMotifsGenome.pl $bed_file $genome $output -p 56 -size 200
    done


}
down_regulated_Island
