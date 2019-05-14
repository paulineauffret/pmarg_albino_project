rm(list=ls()) 
#Load libraries
library(VennDiagram)
library(RSvgDevice)
library(PerformanceAnalytics)
library(dplyr)
library(magrittr)
library(vcfR)
library(adegenet)
library(genepop)
library(PopGenome)
library(stackr)
library(Rcpp)
require("devtools")
library(netview)
#library(pegas)
library(miscTools)
library(Matrix)
library(ggplot2)
library(networkD3)
library(scales)
library(igraph)
library(radiator)
library(ape)
library(parallel)
#install.packages("pcadapt")
library(pcadapt)
#source("https://bioconductor.org/biocLite.R")
#biocLite("qvalue")
library(qvalue)


#################################################################################
##    STEP 1 : Load data
#################################################################################
wd <- "~"
setwd(wd)
sessionInfo()
citation()

#Input files
juv_man_DP20 <-"juv_mantle_merge_freebayes_nAlleles_4_minMapQ_20_minCov_5_filter_maxmiss1_minmaf0.1_minDP20_rmvindels-alb_blck.vcf"

#Set current file to analyse
path_to_file <- juv_man_DP20
filename="juvenile_mantle_DP20"
result_dir=paste(wd,"JUV_MANTLE",sep="/")
path_to_file <- juv_DP20
filename="juv_DP20"
result_dir=paste(wd,"MANTLE",sep="/")

#Import VCF
vcf <- read.vcfR(path_to_file) ; head(vcf)

#Only keeps SNPS, not indels
vcf2 <- extract.indels(vcf)

#Creating genlight object
data <- vcfR2genlight(vcf2)

#Sample names
data$ind.names <- unlist(lapply(as.vector(data$ind.names), function(x) unlist(strsplit(   unlist(strsplit(x, ".", fixed=TRUE))[5], "_", fixed=TRUE))[1]  )) ; data$ind.names

#Set populations
data$pop <-as.factor(c("B","A","A","A","B","B","B","A","B","A","A","A","B","B","A","B","A","B"))

#################################################################################
##    STEP 2 : Dataset description
#################################################################################
#glplot
glPlot(data)

#ACP
pca1 <- glPca(data)
svg(paste(paste(result_dir,filename,sep="/"),"pca.svg",sep="_"))
scatter(pca1, posi="topright", clabel=1.3,col=c("antiquewhite3", "grey20"))
dev.off()

#DAPC 
dapc1 <- dapc(data, n.pca=npc, n.da=nda)
svg(paste(paste(result_dir,filename,sep="/"),"dapc.svg",sep="_"))
scatter(dapc1,scree.da=FALSE, bg="white", posi.pca="topright", legend=TRUE, col=c("antiquewhite3", "grey20", "darkred"))
dev.off()

#Compoplot
svg(paste(paste(result_dir,filename,sep="/"),"compoplot.svg",sep="_"))
compoplot(dapc1, col=c("antiquewhite3","grey20"),txt.leg=levels(data$pop), ncol=1, posi="topright", cleg=1, lab=data$ind.names, show.lab=TRUE, main=paste("Compoplot\nnumber of pca axes=",npc," ; number of da axes=",nda,sep=""))  
dev.off()

#Loadings plot
s <- 0.00001
loadp <- loadingplot(dapc1$var.contr, thres=s) ; loadp
abline(s,0,col="red")


###################################################################
#  STEP 3 : PCAdapt
###################################################################
path_to_file <- juv_man_DP20
pcafile <- read.pcadapt(path_to_file, type = "vcf")
plot(x, option = "screeplot")
plot(x, K=2, option = "pca")
svg(paste(paste(result_dir,filename,sep="/"),"pcadapt.svg",sep="_"))
plot(x, option = "scores", pop = data$pop, col=c("antiquewhite3","grey20"))
plot(x, option = "scores", i = 2, j = 3, pop = data$pop, col=c("antiquewhite3","grey20"))
dev.off()

summary(x)
plot(x , option = "manhattan", K=5)
plot(x, option = "qqplot", K=5)
hist(x$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
plot(x, option = "stat.distribution", K=5)
qval <- qvalue(x$pvalues)$qvalues
alpha <- 0.01
outliers <- which(qval < alpha) ; length(outliers)
padj <- p.adjust(x$pvalues,method="BH")
outliers <- which(padj < alpha) ; length(outliers)
length(x$pvalues)

#place your significant ID in a vector, or specify wi-hich column of a file to look for:
outliers_v <- as.vector(outliers)
#Here, my significant values are in the table res, column outliers  (res$outliers)

#Get VCF:
vcf<-vcf2
#Extract names of markers
names<-cbind(c(1:nrow(vcf@fix)),vcf@fix[,1],vcf@fix[,2],vcf@fix[,4], vcf@fix[,5]) ; head(names)
sig<- names[names[,1] %in% outliers,]
colnames(sig)<-c("index","CHROM","pos", "ref","alt")
dim(sig) ; length(sig) ; head(sig)
sig[1,]
# Merge with annotations
annot <- read.csv("~/sequence_annotation.csv",h=T, sep="\t") ; head(annot) ; dim(annot)
sig_annot <- merge(sig, annot, by.x="CHROM", by.y="Name")

write.table(sig_annot, "juv_mant_pcadapt_0.01_K3", quote=F, row.names = F, sep="\t")

