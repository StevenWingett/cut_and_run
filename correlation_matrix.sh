##== linux command ==##

# Setup
projPath=$PWD
declare -a sampleList=("K27me3_rep1" "K27me3_rep2" "K4me3_rep1" "K4me3_rep2" "IgG_rep1" "IgG_rep2")

# Processing
for sample in ${sampleList[@]}
do
    echo Processing $sample
    ## We use the mid point of each fragment to infer which 500bp bins does this fragment belong to.
    binLen=500
    awk -v w=$binLen '{print $1, int(($2 + $3)/(2*w))*w + w/2}' $projPath/alignment/bed/${sample}_bowtie2.fragments.bed | sort -k1,1V -k2,2n | uniq -c | awk -v OFS="\t" '{print $2, $3, $1}' |  sort -k1,1V -k2,2n  >$projPath/alignment/bed/${sample}_bowtie2.fragmentsCount.bin$binLen.bed

done

echo Done
