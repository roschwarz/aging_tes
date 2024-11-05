#!/env/usr/bin bash
source /misc/paras/data/rschwarz/workspace/phd/projects/mCQuaRna/bin/boneCalls.sh


filterBackground(){

    local target=$1
    local regionFile=$2
    local region=$3 
    
    local file=$(basename $target) 

    if [ ! -f "./storage/motif_analysis/regions_background/.info.txt" ]; then

	echo "The files contain all instances detected as expressed subtracted by \n
    	the target instances (detected DETEs) as recommended by HOMER." > ./storage/motif_analysis/regions_background/.info.txt


    fi
    
    grep -vFw -f <(cut -f 4 $target) $regionFile > ./storage/motif_analysis/regions_background/$region.background.$file

    echo ./storage/motif_analysis/regions_background/$region.background.$file

}

filterForDETEs(){

    local target=$1
    local regionFile=$2
    local region=$3 
    
    local file=$(basename $target) 

    if [ ! -f "./storage/motif_analysis/dete_regions/.info.txt" ]; then

	echo "These files contain the respective regions of my target TE set. \n
    	The file name is created as following region.typeOfExpressedTEs.expressionDirection.tissue.bed .\n
	The region is either promoter - up-stream of CTSSs, down - down-stream of CTSSs or center - window that \n
	has the calle CTSS in  the middle." > ./storage/motif_analysis/regions_background/.info.txt


    fi
    

    grep -Fw -f <(cut -f 4 $target) $regionFile > ./storage/motif_analysis/dete_regions/$region.$file


    echo ./storage/motif_analysis/dete_regions/$region.$file


}


deteRegionFiles=./storage/motif_analysis/dete_regions
regionsBackground=./storage/motif_analysis/regions_background
regionFiles=../cage/rippchen/results_Gclipped/peaks/intersection
genome=/misc/paras/data/genomes/GRCm38_v102/GRCm38.fa
resultsDir=./storage/motif_analysis/Homer_runs_1k/

mkdir -p $deteRegionFiles
mkdir -p $regionsBackground
mkdir -p $resultsDir

cmds=()

echo Run HOMER for 1k windows instead of 500 bp. Run for 3 types of window \
    promoter, center and down. The TE background took expressed TEs as background. \
    For the third homer background run the lists of TEsi, which are considered, were updated \
from dete_lists to te_lists. Now, all expressed TEs are also considered. > $resultsDir/.info.txt

for file in ./storage/motif_analysis/te_lists/*bed; do

    fileName=$(basename $file)
    tissue=$(getTissue $file)
    output=homer_background_third

    for region in center promoter down; do
    
	target=$(filterForDETEs $file $regionFiles/$tissue.TE.1k.$region.bed $region)
	#background=$(filterBackground $file $regionFiles/$tissue.TE.1k.$region.bed $region)
   	results=$resultsDir/${output}/$(basename $target .bed)
	
	echo $results

    	#findMotifsGenome.pl $target $genome $resultsDir -p 20 -bg $background
    
	commander::makecmd -a cmds -s '|' -c {COMMANDER[0]}<<- CMD
    		findMotifsGenome.pl $target $genome $results -p 20
	CMD
    done 
    
done
	
commander::printcmd -a cmds
commander::qsubcmd -n homer -o $output/sge/stats.log -r -p threads -t 20 -i 10 -a cmds
