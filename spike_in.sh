##== linux command ==##

# Setup
projPath=$PWD
declare -a sampleList=("IgG_rep1" "IgG_rep2" "K27me3_rep1" "K27me3_rep2" "K4me3_rep1" "K4me3_rep2")
cores=8
spikeInRef="/home/swingett/data1_wingett/JP/cut_and_run/test_dataset/tutorial/mapping/ecoli_U00096.3_genome/e_coli_U00096.3"
chromSize="/home/swingett/data1_wingett/JP/cut_and_run/test_dataset/tutorial/mapping/hg38.chrom.size"

bowtie_option="--end-to-end"
#bowtie_option="--local"    # For reads > 25bp length - really?


# Processing
for sample in ${sampleList[@]}
do
    echo Processing $sample
    bowtie2 ${bowtie_option} --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -p ${cores} -x ${spikeInRef} -1 ${projPath}/fastq/${sample}_R1.fastq.gz -2 ${projPath}/fastq/${sample}_R2.fastq.gz -S $projPath/alignment/sam/${sample}_bowtie2_spikeIn.sam &> $projPath/alignment/sam/bowtie2_summary/${sample}_bowtie2_spikeIn.txt

    seqDepthDouble=`samtools view -F 0x04 $projPath/alignment/sam/${sample}_bowtie2_spikeIn.sam | wc -l`
    seqDepth=$((seqDepthDouble/2))
    echo $seqDepth >$projPath/alignment/sam/bowtie2_summary/${sample}_bowtie2_spikeIn.seqDepth
done

echo Done