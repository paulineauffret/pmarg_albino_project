#!/usr/bin/Rscript
#Parameter file for DESeq2 analysis

###################################################################################################
#
#     Set working directories
#
###################################################################################################
input_dir <-"~/juvenile"              #Input data directory
output_dir <- "~/juvenile_de_output"     #Output directory


###################################################################################################
#
#     Set metadata parameters
#
###################################################################################################
sampleSheet <- paste(input_dir,"/sampleSheet.info",sep="")                              #SampleSheet. /!\ Need an Id column
colonnes_interet <- c(1,2,3)  #columns to keep
grep_motif <- ".txt"   #input count files motifs
sample_info_name <- "Sample"
#Color parameters
tb <- "grey20"
an <- "antiquewhite3"
colors <- c(tb,an,an,tb,an,tb,an,tb,an,tb)

###################################################################################################
#
#     Set annotations parameters
#
###################################################################################################
annot41k <- read.csv("~/sequence_annotation.csv",h=T, sep="\t") ; head(annot41k) ; dim(annot41k)
annotFuc <- read.csv("~/blastx_transcriptome_combined_e-4_unique.txt",h=T, sep="\t") ; head(annotFuc)
annotFuc <- annotFuc[,c(1,2)] ; head(annotFuc)
colnames(annotFuc) <- c("transcript_name","fucata_hit") ; head(annotFuc) ; dim(annotFuc)


###################################################################################################
#
#     Set filtering parameters
#
###################################################################################################
threshold <- 10 #min expression value accross samples
threshold_nb_samples <- 2
outlier <- ""

###################################################################################################
#
#     Set DESeq2 parameters
#
###################################################################################################
design_formula <- as.formula("~Phenotype")
pvalueAdjustMethod = "BH"
alpha=0.05
FC=1.5
nb_comp=1
contrast_lst <- c("Phenotype","Albino","Black wild-type")
fitParams <- "parametric"
title_graphs <- "Albino vs Black wild-type"

###################################################################################################
#
#     Set PCA parameters
#
###################################################################################################
pca_intergroup <- c("Phenotype")
pca_colors <- c(an, tb)
pca_shape <- c(16,17)
point_size=4
aes_pca_color <- "Phenotype"
aes_pca_shape <- ""

###################################################################################################
#
#     Set Clustering parameters
#
###################################################################################################
nb_max_genes_heatmap <- 75
heatmap_color_variable <- aes_pca_color
dist1 <- "euclidean"
clust <- "average"
ann_colors = list(Phenotype = c("Albino" = an, "Black wild-type" = tb))


