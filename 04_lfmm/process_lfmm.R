library(readr)
library(dplyr, quietly = TRUE, warn.conflicts = FALSE)
library(lfmm)
library(methods)

## Define command line parser ----

parser = argparse::ArgumentParser()

parser$add_argument('--genotype', dest = 'genotype_file', type="character", help='vcf file, freebayes format')
parser$add_argument('--phenotype', dest = 'phenotype_file', type="character", help='')
parser$add_argument('-o', '--output', dest='output_file', default="processed/processed_vcf.txt", help='output file name')

args = parser$parse_args()

## Define Y, the genotype matrix ----

vcf.lfmm <- read_tsv(args$genotype_file)

tY <- vcf.lfmm %>% select(-CHROM, -POS, -ID, -REF, -ALT, -QUAL, -FILTER, -INFO, -FORMAT)
Y <- t(tY %>% select(sort(colnames(tY))))

## Define X, the phenotype vector ----

phenotype <- read_tsv(args$phenotype_file) %>% 
  mutate(Phenotype = as.numeric(factor(Phenotype))) %>% 
  arrange(Sample_name)

X <- phenotype$Phenotype

## Define K, the latent variables count ----

#pc <- prcomp(Y)
#plot(pc$sdev[1:20]^2, xlab = 'PC', ylab = "Variance explained")

K <- 1

## Compute lfmm model ----

mod.lfmm <- lfmm_ridge(Y = Y, X = X, K = K, algorithm = "alternated")

pv <- lfmm_test(Y = Y, 
                X = X, 
                lfmm = mod.lfmm, 
                calibrate = "gif")

pvalues <- pv$calibrated.pvalue 

vcf.lfmm %>% 
  mutate(pvalues = pvalues) %>%
  select(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, pvalues) %>%
  write_tsv(args$output_file)
