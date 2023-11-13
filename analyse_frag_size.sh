##== linux command ==##

# Setup
projPath=$PWD
declare -a sampleList=("K27me3_rep1" "K27me3_rep2" "K4me3_rep1" "K4me3_rep2" "IgG_rep1" "IgG_rep2")

# Processing
mkdir -p $projPath/alignment/sam/fragmentLen

for sample in ${sampleList[@]}
do
    echo Processing $sample
    ## Extract the 9th column from the alignment sam file which is the fragment length
    samtools view -F 0x04 $projPath/alignment/sam/${sample}_bowtie2.sam | awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($9)}' | sort | uniq -c | awk -v OFS="\t" '{print $2, $1/2}' >$projPath/alignment/sam/fragmentLen/${sample}_fragmentLen.txt
done

echo Done