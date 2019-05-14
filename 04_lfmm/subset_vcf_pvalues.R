library(readr)
library(dplyr)

library(argparse)

## Define command line parser ----

parser = argparse::ArgumentParser()

parser$add_argument('--vcf', dest = 'vcf_file', type = "character", help = 'vcf file, freebayes format')
parser$add_argument('--pvalues', dest = 'pvalues_file', type = "character", help = 'pseudo vcf file, with a \"pvalues\" column')
parser$add_argument('--threshold', dest = 'threshold', type = "double", default = 0.01, help = 'pvalue threshold')
parser$add_argument('-o', '--output', dest='output_file', default="subset_vcf.txt", help='output file name')

#args = parser$parse_args(c("--vcf", "manteaux_poches_maxmiss0.9_minmaf0.1-alb_blck_DP20.recode.vcf", "--pvalues", "processed_vcf.txt"))
args = parser$parse_args()

vcf.pval <- read_tsv(args$pvalues_file)
vcf.raw <- read_tsv(args$vcf_file, comment = "#", col_names = colnames(vcf.pval))

left_join(vcf.pval %>% 
            filter(pvalues < args$threshold) %>% select(CHROM, POS), 
          vcf.raw) %>%
  write_tsv(args$output_file)