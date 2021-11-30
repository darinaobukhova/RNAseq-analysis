#!/usr/bin/env Rscript

args <- commandArgs(TRUE)

config_file_run <- args[1] #file with path to different variables

source(as.character(config_file_run))

#---Check and install missing libraries---#
list.of.packages <- c("tidyverse", "methods", "Matrix",
                      "DESeq2", "reshape2", "RColorBrewer", "readr", "dplyr", "AnnotationDbi",
                      "org.Hs.eg.db", "data.table", "ggsci", "magrittr", "multipanelfigure",
                      "Rtsne", "ggplot2", "ellipse", "scales", "sva", "DESeq2")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")

#--- Load libraries---#

lapply(list.of.packages, FUN = function (X) {
  do.call("library", list(X))
})


#---------Data upload--------#

cat("Loading data...")

setwd(wd)

count_data = read.table(file=list.files(pattern="^.*with_counts*.txt$"), header = TRUE, sep = "\t")
row.names(count_data) = count_data$ID
colnames(count_data) <- gsub("^.*\\.", "", colnames(count_data)) #For uniformity of naming
count_data$ID <- NULL
conds = fread(samplesheet) 
conds = conds[order(conds$Group), ]
idx <- match(conds$SampleID, colnames(count_data))
count_data <- count_data[, idx] #Ordering according to the group

#---Batch correction---#

cat("Performing Batch Correction...")

batch <- conds$Sample_Project
group <- conds$Group

corrected_data <- ComBat_seq(as.matrix(count_data), batch=batch, group=group)

save(count_data, corrected_data, conds, file=paste0(outputdir, "/count_data.RData"))

#Applying a variance stabilizing transformation

count_data = vst(as.matrix(count_data))
corrected_data = vst(corrected_data)

#Performing PCA on corrected and uncorrected data (for comparison)

pcDat_raw = prcomp(t(as.matrix(count_data)), scale = F, center = T) 

var_explained_raw <- pcDat_raw$sdev^2/sum(pcDat_raw$sdev^2)

write.csv2(as.data.frame(round(pcDat_raw$x, 3)), 
           file=paste(outputdir, "RawPCA_coords.csv", sep="/"), quote=F)

pcDat_corrected = prcomp(t(as.matrix(corrected_data)), scale = F, center = T) 

var_explained_corrected <- pcDat_corrected$sdev^2/sum(pcDat_corrected$sdev^2)

write.csv2(as.data.frame(round(pcDat_corrected$x, 3)), 
           file=paste(outputdir, "CorrectedPCA_coords.csv", sep="/"), quote=F)

#Scale and center columns
d=svd(scale(as.matrix(corrected_data))) #apply SVD
assays=t(d$u) %*% scale(as.matrix(corrected_data)) #projection on eigenassays, representing the sample space

#--------t-SNE-------#

set.seed(42)

tsne_res <- Rtsne(t(corrected_data), perplexity = 5)

write.csv2(as.data.frame(tsne_res$Y), file=paste(outputdir, "tSNE_coords.csv", sep="/"), quote=F)

#--------Generating plots-------#

cat("Plotting...")

theme_set(theme_bw())

#---PCA for raw vs batch-corrected data---#

p1 <- pcDat_raw$x %>% 
  as.data.frame() %>% 
  ggplot(aes(x = PC1, y = PC2, color = factor(conds$Group), shape=factor(conds$Sample_Project)))+
  geom_point(size = 3.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.text = element_text(size = 6)) + 
  labs(x = paste0("PC1: ", round(var_explained_raw[1]*100, 1), "%"),
       y = paste0("PC2: ", round(var_explained_raw[2]*100, 1), "%"),
       title = "Uncorrected data") +
  guides(color = guide_legend(title="Group"), shape = guide_legend(title="Sequencing Run"))+
  stat_ellipse(linetype=2)+
  scale_color_npg()

p2 <- pcDat_corrected$x %>% 
  as.data.frame() %>% 
  ggplot(aes(x = PC1, y = PC2, color = factor(conds$Group), shape=factor(conds$Sample_Project))) +
  geom_point(size = 3.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.text = element_blank(), legend.position = "none") + 
  labs(x = paste0("PC1: ", round(var_explained_corrected[1]*100, 1), "%"),
       y = paste0("PC2: ", round(var_explained_corrected[2]*100, 1), "%"), 
       title = "Batch corrected data") +
  stat_ellipse(linetype=2) +
  scale_color_npg()

figure <- multi_panel_figure(columns = 1, rows = 2, panel_label_type = "none")

figure %<>%
  fill_panel(p1, column = 1, row = 1) %<>%
  fill_panel(p2, column = 1, row = 2) 

save_multi_panel_figure(figure, filename=paste0(outputdir, 
                                                "/Uncorrected-BatchCorrected-PCA.pdf"))


#---PCA plots with sample names--#

pdf(file=paste0(outputdir, "/PCA_samples.labeled.pdf"), height=7, width=11)

pcDat_corrected$x %>% 
  as.data.frame() %>% 
  ggplot(aes(x = PC1, y = PC2, color = factor(conds$Group), shape=factor(conds$Sample_Project))) + 
  ggrepel::geom_text_repel(aes(label = rownames(pcDat_corrected$x)), max.overlaps = Inf) +
  geom_point(size = 3.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.text = element_text(size = 6)) + 
  labs(x = paste0("PC1: ", round(var_explained_corrected[1]*100, 1), "%"),
      y = paste0("PC2: ", round(var_explained_corrected[2]*100, 1), "%")) +
  scale_color_npg() +
  guides(color = guide_legend(title="Group"), shape = guide_legend(title="Sequencing Run"))+
  stat_ellipse(type="norm", linetype=2)

dev.off()

#---tSNE with sample names---#

pdf(file=paste0(outputdir, "/tSNE_samples.labeled.pdf"), height=7, width=9)

tsne_res$Y %>% 
  as.data.frame() %>% 
  ggplot(aes(x = V1, y = V2, color = factor(conds$Group), shape=factor(conds$Sample_Project))) + geom_point(size = 4) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "tSNE 1", y = "tSNE 2") +
  ggrepel::geom_text_repel(aes(label = conds$SampleID)) +
  scale_color_npg()

dev.off()


#--PCA plot without sample names for a combined figure---# 

f1 <- pcDat_corrected$x %>% 
  as.data.frame() %>% 
  ggplot(aes(x = PC1, y = PC2, color = factor(conds$Group))) + geom_point(size = 3) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "right",
        legend.text = element_text(size = 6), legend.title = element_blank()) + 
  labs(x = paste0("PC1: ", round(var_explained_corrected[1]*100, 1), "%"),
       y = paste0("PC2: ", round(var_explained_corrected[2]*100, 1), "%")) +
  scale_color_npg() +
  guides(shape = guide_legend(override.aes = list(size = 0.5)))

#---Projection of the samples over the top 2 eigenarrays---#

f2 <- assays %>% 
  as.data.frame() %>% 
  ggplot(aes(x = assays[1,], y = assays[2,], colour = factor(conds$Group))) + geom_point(size = 3) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "Eigenarray 1", y = "Eigenarray 2") +
  scale_color_npg()


#---t-SNE plot---#

f3 <- tsne_res$Y %>% 
  as.data.frame() %>% 
  ggplot(aes(x = V1, y = V2, color = factor(conds$Group))) + geom_point(size = 4) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "tSNE 1", y = "tSNE 2") +
  scale_color_npg()
 
#---Combining all the plots into one figure---#

figure1 <- multi_panel_figure(columns = 2, rows = 2, panel_label_type = "none")

figure1 %<>%
  fill_panel(f1, column = 1, row = 1) %<>%
  fill_panel(f2, column = 2, row = 1) %<>%
  fill_panel(f3, column = 2, row = 2)

save_multi_panel_figure(figure1, filename=paste0(outputdir, "/2Dplot.pdf"))

