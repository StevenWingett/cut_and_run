##== linux command ==##

# Setup
projPath=$PWD
declare -a sampleList=("K27me3_rep1" "K27me3_rep2" "K4me3_rep1" "K4me3_rep2" "IgG_rep1" "IgG_rep2")

# Processing

for sample in ${sampleList[@]}
do
    echo Processing $sample
    ## Filter and keep the mapped read pairs
    samtools view -bS -F 0x04 $projPath/alignment/sam/${sample}_bowtie2.sam >$projPath/alignment/bam/${sample}_bowtie2.mapped.bam

    ## Convert into bed file format
    bedtools bamtobed -i $projPath/alignment/bam/${sample}_bowtie2.mapped.bam -bedpe >$projPath/alignment/bed/${sample}_bowtie2.bed

    ## Keep the read pairs that are on the same chromosome and fragment length less than 1000bp.
    awk '$1==$4 && $6-$2 < 1000 {print $0}' $projPath/alignment/bed/${sample}_bowtie2.bed >$projPath/alignment/bed/${sample}_bowtie2.clean.bed

    ## Only extract the fragment related columns
    cut -f 1,2,6 $projPath/alignment/bed/${sample}_bowtie2.clean.bed | sort -k1,1 -k2,2n -k3,3n  >$projPath/alignment/bed/${sample}_bowtie2.fragments.bed
done

echo Done