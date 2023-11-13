##== linux command ==##

# Setup
projPath=$PWD
declare -a sampleList=("K27me3_rep1" "K27me3_rep2" "K4me3_rep1" "K4me3_rep2" "IgG_rep1" "IgG_rep2")
chromSize="/home/swingett/data1_wingett/JP/cut_and_run/test_dataset/tutorial/mapping/hg38.chrom.size"

# Processing
for sample in ${sampleList[@]}
do
    echo Processing $sample

    # Get seqDepth
    seqDepth=$(cat ${projPath}/alignment/sam/bowtie2_summary/${sample}_bowtie2_spikeIn.seqDepth)
    echo -e "\tSequence depth: $seqDepth"

    if [[ "$seqDepth" -gt "1" ]]; then
        
        mkdir -p $projPath/alignment/bedgraph

        scale_factor=`echo "10000 / $seqDepth" | bc -l`
        echo -e "\tScaling factor for $sample is: $scale_factor!"
        bedtools genomecov -bg -scale $scale_factor -i $projPath/alignment/bed/${sample}_bowtie2.fragments.bed -g $chromSize > $projPath/alignment/bedgraph/${sample}_bowtie2.fragments.normalized.bedgraph
        
    fi
done

echo Done