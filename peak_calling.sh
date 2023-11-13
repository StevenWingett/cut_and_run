##== linux command ==##

# See https://github.com/yezhengSTAT/CUTTag_tutorial/issues/2


# Setup
projPath=$PWD
declare -a sampleList=("K27me3_rep1" "K27me3_rep2" "K4me3_rep1" "K4me3_rep2")
declare -a matchingControls=("IgG_rep1" "IgG_rep2" "IgG_rep1" "IgG_rep2")  # Control to match with each sample

#seacr="${projPath}/SEACR-1.3/SEACR_1.3.sh"
#histControl=$2

#Processing 
# Get the location of THIS script
scriptDir=$(cd "$(dirname "${BASH_SOURCE[0]}")" &> /dev/null && (pwd -W 2> /dev/null || pwd))
seacr="${scriptDir}/SEACR-1.3/SEACR_1.3.sh"   # Folder in same directory as this script

mkdir -p $projPath/peakCalling/SEACR

arrayLength=${#sampleList[@]}
for (( i=0; i<$arrayLength; i++ ))
do 
    #echo "${matchingControls[$i]}"

    sampleName=${sampleList[$i]}
    controlName=${matchingControls[$i]}

    echo Processing sample ${sampleName} with matching control ${controlName}

    bash $seacr $projPath/alignment/bedgraph/${sampleName}_bowtie2.fragments.normalized.bedgraph \
    $projPath/alignment/bedgraph/${controlName}_bowtie2.fragments.normalized.bedgraph \
    non stringent $projPath/peakCalling/SEACR/${sampleName}_seacr_control.peaks

    bash $seacr $projPath/alignment/bedgraph/${sampleName}_bowtie2.fragments.normalized.bedgraph 0.01 non stringent $projPath/peakCalling/SEACR/${sampleName}_seacr_top0.01.peaks
done
