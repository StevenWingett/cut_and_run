##== linux command ==##

# Setup
## depending on how you load picard and your server environment, the picardCMD can be different. Adjust accordingly.
picardCMD="java -jar /usr/local/bin/picard.jar"
projPath=$PWD
declare -a sampleList=("K27me3_rep1" "K27me3_rep2" "K4me3_rep1" "K4me3_rep2" "IgG_rep1" "IgG_rep2")

# Processing
echo Project path set to $projPath
mkdir -p $projPath/alignment/removeDuplicate/picard_summary

for sample in ${sampleList[@]}
do
    echo "##################### Processing $sample ##############################"

    ## Sort by coordinate
    $picardCMD SortSam I=$projPath/alignment/sam/${sample}_bowtie2.sam O=$projPath/alignment/sam/${sample}_bowtie2.sorted.sam SORT_ORDER=coordinate

    ## mark duplicates
    $picardCMD MarkDuplicates I=$projPath/alignment/sam/${sample}_bowtie2.sorted.sam O=$projPath/alignment/removeDuplicate/${sample}_bowtie2.sorted.dupMarked.sam METRICS_FILE=$projPath/alignment/removeDuplicate/picard_summary/${sample}_picard.dupMark.txt

    ## remove duplicates
    $picardCMD MarkDuplicates I=$projPath/alignment/sam/${sample}_bowtie2.sorted.sam O=$projPath/alignment/removeDuplicate/${sample}_bowtie2.sorted.rmDup.sam REMOVE_DUPLICATES=true METRICS_FILE=$projPath/alignment/removeDuplicate/picard_summary/${sample}_picard.rmDup.txt
done

echo Done