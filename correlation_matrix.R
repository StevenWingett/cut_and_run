##=== R command ===##
## Path to the project and histone list
library(tidyverse)
library(corrplot)

# Setup
projPath = getwd()
sampleList = c("K27me3_rep1", "K27me3_rep2", "K4me3_rep1", "K4me3_rep2", "IgG_rep1", "IgG_rep2")
histList = c("K27me3", "K4me3", "IgG")


# Processing
# Make an output directory
outdir = paste0(projPath, '/summary_results')
if(dir.exists(outdir) == FALSE) {
    dir.create(outdir, showWarnings = TRUE, recursive = FALSE, mode = "0777")
}

reprod = c()
fragCount = NULL
for(sample in sampleList){
  
  if(is.null(fragCount)){
    
    fragCount = read.table(paste0(projPath, "/alignment/bed/", sample, "_bowtie2.fragmentsCount.bin500.bed"), header = FALSE) 
    colnames(fragCount) = c("chrom", "bin", sample)
  
  }else{
    
    fragCountTmp = read.table(paste0(projPath, "/alignment/bed/", sample, "_bowtie2.fragmentsCount.bin500.bed"), header = FALSE)
    colnames(fragCountTmp) = c("chrom", "bin", sample)
    fragCount = full_join(fragCount, fragCountTmp, by = c("chrom", "bin"))
    
  }
}

M = cor(fragCount %>% select(-c("chrom", "bin")) %>% log2(), use = "complete.obs") 


outfile = paste0(outdir, "/correlation_matrix.png")
png(outfile, width = 2000, height = 1000, res = 120)
corrplot(M, method = "color", outline = T, addgrid.col = "darkgray", order="hclust", addrect = 3, rect.col = "black", rect.lwd = 3,cl.pos = "b", tl.col = "indianred4", tl.cex = 1, cl.cex = 1, addCoef.col = "black", number.digits = 2, number.cex = 1, col = colorRampPalette(c("midnightblue","white","darkred"))(100))
dev.off()

outfile = paste0(outdir, "/correlation_matrix.pdf")
pdf(outfile)
corrplot(M, method = "color", outline = T, addgrid.col = "darkgray", order="hclust", addrect = 3, rect.col = "black", rect.lwd = 3,cl.pos = "b", tl.col = "indianred4", tl.cex = 1, cl.cex = 1, addCoef.col = "black", number.digits = 2, number.cex = 1, col = colorRampPalette(c("midnightblue","white","darkred"))(100))
dev.off()

print("Done")