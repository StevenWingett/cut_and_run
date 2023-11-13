##== linux command ==##

# Setup
projPath=$PWD
declare -a sampleList=("K27me3_rep1" "K27me3_rep2" "K4me3_rep1" "K4me3_rep2" "IgG_rep1.1" "IgG_rep1.2" "IgG_rep2")
declare -a nonControlSampleList=("K27me3_rep1" "K27me3_rep2" "K4me3_rep1" "K4me3_rep2")
genome_annotation_gtf_file=$projPath/Homo_sapiens.GRCh38.102.gtf
cores=8

# Processing
# Heatmap visualization on specific regions (requires https://deeptools.readthedocs.io/en/develop/)
mkdir -p $projPath/alignment/bigwig      

echo Creating heatmaps on specific regions
for sample in ${sampleList[@]}
do
    echo Processing $sample
    samtools sort -o $projPath/alignment/bam/${sample}.sorted.bam $projPath/alignment/bam/${sample}_bowtie2.mapped.bam                                                     
    samtools index $projPath/alignment/bam/${sample}.sorted.bam                                                                                                              
    bamCoverage -b $projPath/alignment/bam/${sample}.sorted.bam -o $projPath/alignment/bigwig/${sample}_raw.bw   
done

# Heatmap over transcription units
echo Creating heatmaps on transcription units
nonControlsamplesString=''

for nonControlsample in ${nonControlSampleList[@]}
do
    echo Processing $nonControlsample
    nonControlsamplesString="$nonControlsamplesString $projPath/alignment/bigwig/${nonControlsample}_raw.bw"
done

# echo $subcommand

mkdir -p $projPath/data/matrix

computeMatrix scale-regions -S $nonControlsamplesString \
                            -R $genome_annotation_gtf_file \
                            --beforeRegionStartLength 3000 \
                            --regionBodyLength 5000 \
                            --afterRegionStartLength 3000 \
                            --skipZeros -o $projPath/data/matrix/matrix_gene.mat.gz -p $cores

plotHeatmap -m $projPath/data/matrix/matrix_gene.mat.gz -out $projPath/data/matrix/Histone_gene.png --sortUsing sum

# Heatmap on CUT&Tag peaks 
for nonControlsample in ${nonControlSampleList[@]}
do
    echo Processing $nonControlsample
    
    awk '{split($6, summit, ":"); split(summit[2], region, "-"); print summit[1]"\t"region[1]"\t"region[2]}' $projPath/peakCalling/SEACR/${nonControlsample}_seacr_control.peaks.stringent.bed >$projPath/peakCalling/SEACR/${nonControlsample}_seacr_control.peaks.summitRegion.bed

    computeMatrix reference-point -S $projPath/alignment/bigwig/${nonControlsample}_raw.bw \
                 -R $projPath/peakCalling/SEACR/${nonControlsample}_seacr_control.peaks.summitRegion.bed \
                 --skipZeros -o $projPath/peakCalling/SEACR/${nonControlsample}_SEACR.mat.gz -p $cores -a 3000 -b 3000 --referencePoint center

    plotHeatmap -m $projPath/peakCalling/SEACR/${nonControlsample}_SEACR.mat.gz -out $projPath/peakCalling/SEACR/${nonControlsample}_SEACR.heatmap.png --sortUsing sum --startLabel "Peak Start" --endLabel "Peak End" --xAxisLabel "" --regionsLabel "Peaks" --samplesLabel "${nonControlsample}"
done

echo Done