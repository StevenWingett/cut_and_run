##== linux command ==##

# Setup
projPath=$PWD
declare -a sampleList=("K27me3_rep1" "K27me3_rep2" "K4me3_rep1" "K4me3_rep2" "IgG_rep1" "IgG_rep2")


# Processing
minQualityScore=2
for sample in ${sampleList[@]}
do
    echo Processing $sample
    samtools view -q $minQualityScore ${projPath}/alignment/sam/${sample}_bowtie2.sam >${projPath}/alignment/sam/${sample}_bowtie2.qualityScore$minQualityScore.sam
done

echo Done