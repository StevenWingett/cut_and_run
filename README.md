# Overview
**These are collected scripts that constitute a pipeline to analyse Cut and Run / Cut and Tag datasets.  The scripts have been taken from the [tutorial](https://yezhengstat.github.io/CUTTag_tutorial/) put together by Ye Zheng, Kami Ahmad and Steven Henikoff.  The scripts have been modified to facilitate processing - see below for further details.**

# Pre-pipeline steps

## QC
Use FASTQC and FastQ Screen.


## Trimming
Use trim galore if reads >25bp


# Pipeline

Pipeline details:
https://yezhengstat.github.io/CUTTag_tutorial/

## Prep for pipeline

### Filenames
Input fastq files need to be of the format when passing to the pipeline

[sample_name]_rep[1,2,3,...n].fastq.gz

No underscores are allowed in the sample_name.

### Create a chromosome sizes file for the genome (not the spike-in genome)
`wget https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/faSize`

`./faSize -detailed -tab file.fasta`


## Mapping
Place FASTQ files in a folder named 'fastq'

`nohup bash ./mapping.sh > log.mapping.sh.log &`

 
## Spike in

`nohup bash ./spike_in.sh > log.spike_in.sh.log &`


## Summarise mapping

`Rscript alignment_summary.R`


## Identify / remove duplicates
`nohup bash analyse_duplicates.sh > log.analyse_duplicates.sh.log &`

`Rscript summarise_duplicates.R`


## Assess mapped fragment size distribution
`nohup bash analyse_frag_size.sh > log.analyse_frag_size.out &`

`Rscript summarise_frag_length.R`


## Filter by read quality ##
`nohup bash filter_quality.sh > log_filter_quality.sh.log &`


## Convert to different file formats
`nohup bash file_format_conversion.sh > log.file_format_conversion.sh.out &`


## Assess replicate reproducibility
`nohup bash correlation_matrix.sh > log.correlation_matrix.sh.out &`

`Rscript correlation_matrix.R`


## Spike-in calibration
`nohup bash spike_in_calibration.sh > log.spike_in_calibration.sh.out &`

`Rscript spike_in_calibration.R`


## Peak calling
`nohup bash peak_calling.sh > log.peak_calling.sh.out &`

`Rscript summarise_peaks.R`


## Visualisation
`nohup bash visualisation.sh > log.visualisation.sh.out &`