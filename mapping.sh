##== linux command ==##

# Setup
projPath=$PWD
declare -a sampleList=("IgG_rep1" "IgG_rep2" "K27me3_rep1" "K27me3_rep2" "K4me3_rep1" "K4me3_rep2")

cores=8
ref="/data1/Genome_References/Ensembl/Homo_sapiens/GRCh38/Bowtie2/Homo_sapiens.GRCh38.dna.primary_assembly"

bowtie_option="--end-to-end"
#bowtie_option="--local"    # For reads > 25bp length


#Processing
mkdir -p ${projPath}/alignment/sam/bowtie2_summary
mkdir -p ${projPath}/alignment/bam
mkdir -p ${projPath}/alignment/bed
mkdir -p ${projPath}/alignment/bedgraph

# Processing
for sample in ${sampleList[@]}
do
    echo Processing $sample
    bowtie2 ${bowtie_option} --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -p ${cores} -x ${ref} -1 ${projPath}/fastq/${sample}_R1.fastq.gz -2 ${projPath}/fastq/${sample}_R2.fastq.gz -S ${projPath}/alignment/sam/${sample}_bowtie2.sam &> ${projPath}/alignment/sam/bowtie2_summary/${sample}_bowtie2.txt

done

echo Done