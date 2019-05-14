#!/usr/bin/Rscript
#ADifferential analysis, P. margaritifera albino phenotype versus black wild-type phenotype 
#INPUT = HTseq-count output files (txt files, one per sample, raw) + sample sheet with samples metadata (txt file, rows = samples ; col = attributes)
#see https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.pdf

###################################################################################################
#
#     Load libraries
#
###################################################################################################
#source("http://bioconductor.org/biocLite.R")
library(contrast)
library(RpngDevice)
library(DESeq2)
library(tximport)
library(gplots)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(VennDiagram)
library(readODS)
library(assertthat)
library(scales)
library(WGCNA)
library(genefilter)
library(ggrepel)
library(rafalib)
library(egg)

rm(list=ls()) 
sessionInfo <- sessionInfo()

###################################################################################################
#
#     I. Set parameters - Load parameters file
#
###################################################################################################
#Choose the dataset : 0 = juvenile ; 1 = mantle ; 2 = pearl sac
dataset_id <- 0
if (dataset_id == 0) {
	source("params_juvenile.r")
} else if (dataset_id == 1) {
	source("params_mantle.r")
} else {
	source("params_pearlsac.r")
}

source("de_analysis_functions.r")

system(command=paste("mkdir -p ",output_dir,sep=""))

setwd(output_dir) ; getwd()

###################################################################################################
#
#     II. DESeq2 analysis from HTSeq-count data
#
###################################################################################################
#Imports samples infos
samplesInfo <- read.csv(sampleSheet, h=T,sep="\t") ; head(samplesInfo) ; dim(samplesInfo)
target_file <- samplesInfo[,colonnes_interet] ; head(target_file) ; dim(target_file)

sampleFiles <- grep(grep_motif,list.files(input_dir),value=TRUE)
sampleCondition <- samplesInfo
sampleTable <- data.frame(sampleName = sampleFiles, fileName = sampleFiles)
sampleTable <- merge(sampleTable,samplesInfo,by.x="sampleName",by.y=sample_info_name)
sampleTable <- sampleTable[which(sampleTable$fileName != outlier),]

#Creates DEseq object
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = input_dir, design = design_formula) ; ddsHTSeq

#Trims too low transcripts : only keep rows with sum > threshold
ddsHTSeq_f <- ddsHTSeq[ rowSums(counts(ddsHTSeq) >= threshold) >= threshold_nb_samples, ] ; ddsHTSeq_f
sum(colSums(counts(ddsHTSeq_f)))
nb_genes <- dim(ddsHTSeq_f)[1] ; nb_genes
bf <- rownames(ddsHTSeq) ; length(bf)
af <- rownames(ddsHTSeq_f) ; length(af)
trash <- setdiff(bf, af) ; length(trash) ; test <- length(trash)+length(af) ; test

#Look at count distribution
print(colData(ddsHTSeq_f)$Id)
boxplot_bf <- boxplot(counts(ddsHTSeq),outline=F,names=colData(ddsHTSeq)$Id, las=1, col=colors, cex.names=0.7, cex=0.7, main="Counts distribution before filtering")
boxplot_af <- boxplot(counts(ddsHTSeq_f),outline=F,names=colData(ddsHTSeq_f)$Id, las=1, col=colors, cex.names=0.7,main="Counts distribution after filtering")

#Estimate size factors
dds <- estimateSizeFactors(ddsHTSeq_f)
sizeFactors(dds)

#Run DE analysis with DESeq2
dds <- DESeq(dds, fitType=fitParams)
resultsNames(dds)

#Plot dispersion
plotDispEsts(dds)

###################################################################################################
#
#     III. Graphic exploration
#
###################################################################################################
###################################################################################################
#     III.1. Principal Components Analysis
###################################################################################################
#PCA axis 1 & 2
rld <- rlog(ddsHTSeq_f, blind=FALSE)
data <- plotPCA(rld, intgroup=pca_intergroup, returnData=TRUE)
pointLabel <- colData(ddsHTSeq_f)$Id
percentVar <- round(100*attr(data, "percentVar"),2)
ggplot(data, aes(PC1, PC2)) +
  geom_point(size=point_size,aes_string(color=aes_pca_color)) +#, shape=aes_pca_shape)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  #ggtitle(paste("Principal Component Analysis plot ",nb_genes," genes", sep="")) +
  geom_hline(aes(yintercept=0), col="grey") +
  geom_vline(aes(xintercept=0), col="grey") +
  #scale_shape_manual(values=pca_shape) +
  scale_color_manual(values=pca_colors) +
  theme_light() +
  ggtitle("Adult Stage, Pearl Sac") +
  theme(plot.title = element_text(face="bold", size=12)) +
  #theme(legend.position="left")
  theme(legend.position='none')
  #geom_text(aes(label=pointLabel),hjust=0.5, vjust=2, size=2.7)
dev.off()


##PCA axis 1 & 3
data_2 <- plotPCA.san(rld, intgroup=pca_intergroup, returnData=TRUE)
pointLabel <- colData(ddsHTSeq_f)$Id
percentVar <- round(100*attr(data_2, "percentVar"))
ggplot(data_2, aes(PC1, PC3)) +
  geom_point(size=point_size,aes_string(color=aes_pca_color)) +#, shape=aes_pca_shape)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC3: ",percentVar[3],"% variance")) +
  ggtitle(paste("Principal Component Analysis plot ",nb_genes," genes", sep="")) +
  geom_hline(aes(yintercept=0), col="grey") +
  geom_vline(aes(xintercept=0), col="grey") +
  #scale_shape_manual(values=pca_shape) +
  scale_color_manual(values=pca_colors) +
  geom_text(aes(label=pointLabel),hjust=0.5, vjust=2, size=3) 


##PCA axis 2 & 3
data_3 <- plotPCA.san2(rld, intgroup=pca_intergroup, returnData=TRUE)
pointLabel <- colData(ddsHTSeq_f)$Id
percentVar <- round(100*attr(data_3, "percentVar"))
ggplot(data_3, aes(PC2, PC3)) +
  geom_point(size=point_size,aes_string(color=aes_pca_color)) + #, shape=aes_pca_shape)) +
  xlab(paste0("PC2: ",percentVar[2],"% variance")) +
  ylab(paste0("PC3: ",percentVar[3],"% variance")) +
  ggtitle(paste("Principal Component Analysis plot ",nb_genes," genes", sep="")) +
  geom_hline(aes(yintercept=0), col="grey") +
  geom_vline(aes(xintercept=0), col="grey") +
  #scale_shape_manual(values=pca_shape) +
  scale_color_manual(values=pca_colors) +
  geom_text(aes(label=colData(ddsHTSeq_f)$Id),hjust=0.5, vjust=2, size=3)


###################################################################################################
#     III.2. Heatmap & Clustering
###################################################################################################
#Selects only most abundant transcripts
nCounts <- counts(dds, normalized=TRUE)
select <- order(rowMeans(nCounts),decreasing=TRUE)[1:nb_max_genes_heatmap]

#Selects corresponding norm counts
nt <- normTransform(dds)
log2.norm.counts <- assay(nt)
log2.norm.counts <- assay(nt)[select,]

#Gets the metadata
colnames(log2.norm.counts) <- colData(dds)$Id
df <- as.data.frame(colData(dds)[,heatmap_color_variable])
rownames(df) <- colData(dds)$Id
colnames(df) <- heatmap_color_variable

#Draw heatmap
pheatmap(log2.norm.counts,
         clustering_distance_cols = dist1, 
         clustering_method = clust, 
         cluster_rows=FALSE,
         cluster_cols=TRUE, 
         annotation_col=df,
         show_rownames=TRUE,
         fontsize_row=5,
         fontsize_col=10,
         fontsize=8,
         annotation_colors = ann_colors, 
         main=paste("Clustered heatmap of ",nb_max_genes_heatmap," most abundant genes\n",dist1," distance with ",clust, " clustering method",sep=""),
         las=1)

#Heatmap sample to sample
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds), Id)
hc <- hclust(distsRL)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(mat, Rowv=as.dendrogram(hc),
          symm=TRUE, trace='none',
          col=rev(hmcol),margin=c(10, 10), colRow=colors, colCol=colors)


#Clustering
distance_matrix <- dist(t(log2.norm.counts), method = dist1)
hh <- hclust(distance_matrix, method = clust)
myplclust(hh, labels=rownames(df), lab.col=colors)


###################################################################################################
#
#     IV. DE genes results
#
###################################################################################################
###################################################################################################
#     IV.1. Get DE genes
###################################################################################################
res <- results(dds,contrast=contrast_lst, pAdjustMethod = pvalueAdjustMethod, alpha = alpha) ; dim(res)

#Selects results with p-value <= 0.05 et |log2FC| >= 2
res_f <- res[which(res$padj<=alpha),] ; res_f<-res_f[which(abs(res_f$log2FoldChange) >= FC),] ; dim(res_f)

#Merge DE genes and annotations
res_f <- as.data.frame(res_f)
res_f$tr_name <- rownames(res_f)
res_f_annot <- merge(res_f,annot41k,by.x="tr_name",by.y="Name", all.x=TRUE, all.y=FALSE) ; head(res_f_annot)
res_f_annot_fuc <- merge(res_f_annot,annotFuc,by.x="tr_name",by.y="transcript_name", all.x=TRUE, all.y=FALSE) ; head(res_f_annot_fuc)

###################################################################################################
#     IV.2. MA plots
###################################################################################################
plotMA.DESeqResults(res, main=title_graphs, alpha=alpha)


write.table(set1_set2_set3_annot_fuc,"juv_mant_sac", row.names = FALSE, quote=F, sep="\t")




#Uniq1 - Uniq2
inter12 <- intersect(uniq1, uniq2) ; length(inter12)
inter12 <- intersect(tmp_set, uniq1) ; length(inter12)
draw.pairwise.venn(length(uniq1),length(uniq2),length(inter12),category=c(uniq1_name,uniq2_name),
                   fill=c(col_uniq1,col_uniq2), cat.col=c(col_uniq1,col_uniq2), cat.cex=0.8, cat.pos=c(3,6), cat.dist=c(0.02,0.02),
                   cat.fontfamily=c("sans","sans"))
dev.off()

all <- merge(uniq1_tab, uniq2_tab, by.x="tr_name", by.y="tr_name",all.x=T, all.y=T)

head(all)
###################################################################################################
#
#     Print tables and figures
#
###################################################################################################
getwd()
png("boxplot_count_data_bf_filtering.png")
boxplot(counts(ddsHTSeq),outline=F,names=colData(ddsHTSeq)$Id, las=1, col=colors, cex.names=0.7)
dev.off()

png("boxplot_count_data_after_filtering.png")
boxplot(counts(ddsHTSeq_f),outline=F,names=colData(ddsHTSeq_f)$Id, las=1, col=colors, cex.names=0.7)
dev.off()

png("plotDispersion.png")
plotDispEsts(dds)
dev.off()

png("PCA_PC1_PC2.png")
data <- plotPCA(rld, intgroup=pca_intergroup, returnData=TRUE)
pointLabel <- colData(ddsHTSeq_f)$Id
percentVar <- round(100*attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2)) +
  geom_point(size=point_size,aes_string(color=aes_pca_color))+#, shape=aes_pca_shape)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle(paste("Principal Component Analysis plot ",nb_genes," genes", sep="")) +
  geom_hline(aes(yintercept=0), col="grey") +
  geom_vline(aes(xintercept=0), col="grey") +
  #scale_shape_manual(values=pca_shape) +
  scale_color_manual(values=pca_colors) +
  geom_text(aes(label=pointLabel),hjust=0.5, vjust=2, size=3)
dev.off()

png("PCA_PC1_PC3.png")
data_2 <- plotPCA.san(rld, intgroup=pca_intergroup, returnData=TRUE)
pointLabel <- colData(ddsHTSeq_f)$Id
percentVar <- round(100*attr(data_2, "percentVar"))
ggplot(data_2, aes(PC1, PC3)) +
  geom_point(size=point_size,aes_string(color=aes_pca_color))+#, shape=aes_pca_shape)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC3: ",percentVar[3],"% variance")) +
  ggtitle(paste("Principal Component Analysis plot ",nb_genes," genes", sep="")) +
  geom_hline(aes(yintercept=0), col="grey") +
  geom_vline(aes(xintercept=0), col="grey") +
  #scale_shape_manual(values=pca_shape) +
  scale_color_manual(values=pca_colors) +
  geom_text(aes(label=pointLabel),hjust=0.5, vjust=2, size=3)
dev.off()

png("PCA_PC2_PC3.png")
data_3 <- plotPCA.san2(rld, intgroup=pca_intergroup, returnData=TRUE)
pointLabel <- colData(ddsHTSeq_f)$Id
percentVar <- round(100*attr(data_3, "percentVar"))
ggplot(data_3, aes(PC2, PC3)) +
  geom_point(size=point_size,aes_string(color=aes_pca_color, shape=aes_pca_shape)) +
  xlab(paste0("PC2: ",percentVar[2],"% variance")) +
  ylab(paste0("PC3: ",percentVar[3],"% variance")) +
  ggtitle(paste("Principal Component Analysis plot ",nb_genes," genes", sep="")) +
  geom_hline(aes(yintercept=0), col="grey") +
  geom_vline(aes(xintercept=0), col="grey") +
  scale_shape_manual(values=pca_shape) +
  scale_color_manual(values=pca_colors) +
  geom_text(aes(label=colData(ddsHTSeq_f)$Id),hjust=0.5, vjust=2, size=3)
dev.off()


png("heatmap.png")
pheatmap(log2.norm.counts,
         clustering_distance_cols = dist1, 
         clustering_method = clust, 
         cluster_rows=FALSE,
         cluster_cols=TRUE, 
         annotation_col=df,
         show_rownames=TRUE,
         fontsize_row=5,
         fontsize_col=10,
         fontsize=8,
         annotation_colors = ann_colors, 
         main=paste("Clustered heatmap of ",nb_max_genes_heatmap," most abundant genes\n",dist1," distance with ",clust, " clustering method",sep=""),
         las=1)
dev.off()

png("heatmap_samplebysample.png")
heatmap.2(mat, Rowv=as.dendrogram(hc),
          symm=TRUE, trace='none',
          col=rev(hmcol),margin=c(10, 10), colRow=colors, colCol=colors)
dev.off()

png("clustering.png")
myplclust(hh, labels=rownames(df), lab.col=colors)
dev.off()

png("MAplot_res.png")
plotMA.DESeqResults(res, main=title_graphs, alpha=alpha)
dev.off()

write.table(trash,"excluded_transcripts_from_filtering.txt",quote=F,row.names=F)
write.table(res,"de_analysis_complete.txt", quote=F, sep="\t")
write.table(res_f,paste("de_analysis_pv_",alpha,"_FC_",FC,".txt", sep=""), quote=F, sep="\t")
write.table(mcols(res,use.names=TRUE),"res_description.txt",quote=F,row.names=F)
write.table(res_f_annot, paste("de_analysis_pv_",alpha,"_FC_",FC,"_annot.txt", sep=""),row.names = FALSE, quote=F, sep="\t")
write.table(res_f_annot_fuc, paste("de_analysis_pv_",alpha,"_FC_",FC,"_annot_fuc.txt", sep="") ,row.names = FALSE, quote=F, sep="\t")

