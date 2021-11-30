args <- commandArgs(TRUE)

config_file_run <- args[1] #file with path to different variables (working directory, outputdirectory, samplesheet information)

library(data.table)
library(methods)
library(Rtsne)
library(magrittr)
library(multipanelfigure)
library(ggplot2)

setwd(workdir)

#------Data Loading----#

cat("Loading data...")

count_data = read.table(file=list.files(pattern = "expr.groups.*cpm_log2.*.txt$"), header = TRUE, sep = "\t")
row.names(count_data) = count_data$ID
count_data$ID <- NULL
conds = fread(samplesheet, select = c(1:2)) 
conds = conds[order(conds$Sample_Project), ]

#-----Performing tSNE-----#

cat("Performing tSNE...")

max_perplexity = round(ncol(count_data)/3) #perplexity parameter should not be bigger than 3*perplexity < nrow(x)

set.seed(42)
pdf(file=paste0(outputdir, "/tSNE_result.pdf"), width=6, height=9)
par(mfrow=c(3,3))
for(i in seq(from=5, to=max_perplexity)){
  x <- t(count_data)
  res <- Rtsne(x, perplexity = i, max_iter = 10000)
  plot(res$Y, col = c("#3C5488B2", "#DC0000B2")[as.factor(conds$Group)],
             pch = 19, xlab="tSNE 1", ylab="tSNE 2", cex.lab=0.5)
  mtext(paste0("Perplexity:", i))
} 

dev.off()

