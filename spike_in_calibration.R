##=== R command ===## 
library(tidyverse)
library(viridis)
library(ggpubr)


# Setup
projPath = getwd()
sampleList = c("K27me3_rep1", "K27me3_rep2", "K4me3_rep1", "K4me3_rep2", "IgG_rep1", "IgG_rep2")
histList = c("K27me3", "K4me3", "IgG")
multiplier = 10000

# Processing
# Make an output directory
outdir = paste0(projPath, '/summary_results')
if(dir.exists(outdir) == FALSE) {
    dir.create(outdir, showWarnings = TRUE, recursive = FALSE, mode = "0777")
}

# Read in the duplication summary file
infile = paste0(projPath, "/summary_results/duplication_summary.tsv")
print(paste("Reading in", infile))
alignDupSummary = read_delim(infile, delim = "\t",  col_types = cols(Replicate = col_factor()))


# Read in the alignment summary file
infile = paste0(projPath, "/summary_results/alignment_summary.tsv")
print(paste("Reading in", infile))
alignResult = read_delim(infile, delim = "\t",  col_types = cols(Replicate = col_factor()))

# Edit so it is the same as the alignResult tible in a previous script
# We haven't removed the % from the final column as this does not appear necessary
alignResult$Histone <- as.factor(alignResult$Histone)
alignResult <- alignResult[1:(length(alignResult)-2)]    # Drop last 2 columns


scaleFactor = c()

for(sample in sampleList){
  spikeDepth = read.table(paste0(projPath, "/alignment/sam/bowtie2_summary/", sample, "_bowtie2_spikeIn.seqDepth"), header = FALSE, fill = TRUE)$V1[1]
  
  histInfo = strsplit(sample, "_")[[1]]
  scaleFactor = data.frame(scaleFactor = multiplier/spikeDepth, Histone = histInfo[1], Replicate = histInfo[2])  %>% rbind(scaleFactor, .)
}
scaleFactor$Histone = factor(scaleFactor$Histone, levels = histList)
calibrationSummary <- left_join(alignDupSummary, scaleFactor, by = c("Histone", "Replicate"))

outfile = paste0(outdir, "/spike_in_calibration_summary.tsv")
write_delim(calibrationSummary, outfile, delim = "\t")

## Generate sequencing depth boxplot
fig6A = scaleFactor %>% ggplot(aes(x = Histone, y = scaleFactor, fill = Histone)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 20) +
    ylab("Spike-in Scalling Factor") +
    xlab("")

normDepth = inner_join(scaleFactor, alignResult, by = c("Histone", "Replicate")) %>% mutate(normDepth = MappedFragNum_hg38 * scaleFactor)

fig6B = normDepth %>% ggplot(aes(x = Histone, y = normDepth, fill = Histone)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 20) +
    ylab("Normalization Fragment Count") +
    xlab("") + 
    coord_cartesian(ylim = c(1000000, 130000000))


outfile = paste0(outdir, "/spike_in_calibration_summary_plots.png")
png(outfile, width = 2000, height = 1000, res = 120)
ggarrange(fig6A, fig6B, ncol = 2, common.legend = TRUE, legend="bottom")
dev.off()


print('Done')