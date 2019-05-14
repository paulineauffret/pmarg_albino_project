# pmarg_albino_project
Tracing key genes associated with Pinctada margaritifera albino phenotype from juvenile to cultured pearl harvest stages by multiple whole transcriptome sequencing

#Pipeline
Many thanks to Leo Milhade and Jeremy Le luyer for help with lfmm pipeline and alternative splicing analysis, respectively.

01_from_raw_to_count_matrix_and_vcf
	01_1_juvenile
		01_fastqc.sh --> 01_fastqc.qsub
		02_trimmomatic.sh --> 02_trimmomatic.qsub
		03_mapping_PE_bwa.sh
		04_filtering_bamfiles.sh
		05_counting_tag_htseq.sh
		06_filtering_bam_snp.sh --> 06_filtering_bam_snp.qsub
		07_snp_calling_freebayes.qsub
		08_vcftools.qsub
	01_2_mantle
		01_trimming_PE_cutadapt.qsub
		02_mapping_PE_bwa.sh
		03_filtering_bamfiles.sh
		04_counting_tag_htseq.sh
		05_1_markdup.sh --> 05_2_markdup.qsub
		06_freebayes.qsub
		07_vcftools.qsub
	01_3_pearlsac
		01_1_fastqc_raw_reads.sh --> 01_2_fastqc_raw_reads.qsub
		02_1_trimming_trimmomatic_PE.sh --> 02_2_trimming_trimmomatic_PE.qsub
		03_mapping_PE_bwa.sh
		04_filtering_bamfiles.sh
		05_counting_tag_htseq.sh
	01_4_merge_vcf.qsub
02_DE_analysis
03_goatools
04_lfmm
05_PCAdapt.r
06_annotate_snp
07_splicing

#Material and Methods
##Differential expression analysis
Read quality was assessed with fastqc v0.11.5 (“Babraham Bioinformatics - FastQC A Quality Control Tool for High Throughput Sequence Data” n.d.) and multiQC (Ewels et al. 2016). Raw reads were filtered for sequencing adapters removal and quality trimming (Q=28) using i. Cutadapt v1.13 (Martin 2011) for Mantle (M) M dataset and ii. Trimmomatic v0.36 (Bolger, Lohse, and Usadel 2014) for Pearl Sac (PS) and Juvenile (J) datasets. For each dataset, only surviving paired-end reads were retained. Filtered reads were mapped on the multi-tissue reference transcriptome of Pinctada margaritifera (41,075 contigs) (Le Luyer et al. 2019) using bwa v0.7.15 (Li and Durbin 2009) with standard parameters. Low mapping quality (Q≥5), mispaired and multi-mapping were removed using Samtools v1.4.1 (Li et al. 2009). A matrix of raw counts was built using HTSeq-count v0.6.1 (Anders, Pyl, and Huber 2015). To minimize the false-positive rate, the count matrices were filtered for low expressed transcripts. For each dataset, all transcripts with less than 10 counts in at least 2 samples were discarded. We identified differentially expressed genes (DEGs) between albino and black wild-type individuals using DESeq2 v1.16.2 (Love, Huber, and Anders 2014) with R v3.4.0 (https://www.R-project.org/) following the standard workflow. The DESeq2 method internally corrects for library size and uses negative binomial generalized linear models to test for differential expression. In this study the statistical models were built using 'counts ~ Phenotype' design formula, where Phenotype qualitative variable indicates oyster phenotype (albino/black wild-type). All features with absolute log2 fold change greater than 1.5 and adjusted p-value smaller than 0.05 (Benjamini-Hochberg method) were reported as differentially expressed (DEGs). Overlapping DEGs between the three datasets were vizualised using VennDiagram R package v1.6.20 (https://github.com/cran/VennDiagram). To test the overrepresentation of gene ontology (GO) terms in resulting DEGs, we used Goatools v0.6.10 (Klopfenstein et al. 2018) throught ‘go_enrichment’ pipeline (https://github.com/enormandeau/go_enrichment) with go-basic.obo database (release 2017-04-14). The resulting lists of significant GO enriched terms were filtered for Fisher’s Test p-value < 0.05 and used for semantic-based clustering using REVIGO with allowed similarity=0.5 (Supek et al. 2011). DEGs lists were sumitted to KASS server ((Moriya et al. 2007), Last updated: April 3, 2015) to visualize related KEGG pathways. 
##Population genetics analysis
In J and M datasets, we investigated variable sites between white albino and black wild-type populations. We did not include PS dataset because of possible contamination with recipient RNA during sampling of pearl sacs. Single nucleotide polymorphisms (SNPs) were called from preprocessed aligned reads using Freebayes v1.1.0-3-g961e5f3 (Garrison and Marth 2012) with a required minimum mapping quality of 20. Preprocessing of aligned reads included marking and removing duplicates, correcting N cigar reads, sorting and indexing bam files using gatk v4.0.2.1 (McKenna et al. 2010) and Picard tools suite v1.119 (“Picard Tools - By Broad Institute” n.d.). Resulting variant calling files (VCF) were filtered for missing data and indels (none authorized), allele frequency (≥0.1) and depth (≥20) using vcftools v0.1.14 (Danecek et al. 2011). Filtered VCF for J and M were merged in order to minimize family effect and focus on phenotype-associated events. To investigate the structure of albino and black wild-type individuals within the two populations, we performed population genetics analysis on the filtered VCF files using the following R packages : vcfR v1.8.0 (Knaus and Grünwald 2017), adegenet v2.1.1 (Jombart and Ahmed 2011) and genepop v1.0.5 (Rousset 2008) with R v3.4.0. We then looked for outlier SNPs (adjusted pvalue<0.01, Benjamini-Hochberg procedure) with PCAdapt v4.0.3 (Luu, Bazin, and Blum 2017) PCAdapt tests for outliers using correlations between SNPs and first principal components (K) of PCA analysis. We used K=2 (J and M separately) and K=3 (M and J merged) and adjusted p-value<0.01according to Benjamini-Hochberg procedure. We also performed phenotype-associated SNPs analysis using lfmm v0.0 (https://bcm-uga.github.io/lfmm/index.html) to highlight potential genetic markers of albinism. Lfmm software constructs latent factor mixed models (LFMMs), which are statistical regression models to test associations between a multidimensional set of response variables (here, genotypes) and a variable of interest (here, phenotype). LFMMs include unobserved variables, called latent factors, that correct the model for confounding effects due to population structure and other hidden causes. We selected phenotype-associated SNPs according to adjusted pvalue<0.01. Phenotype-associated SNPs were annotated using using ‘LongOrfs’ function implemented in Transdecoder v3.0.1 (Haas et al. 2013) and personal python script.
##Alternative gene splicing and exon usage
To detect differential splicing events in the three RNA-seq datasets, the filtered reads were mapped on the draft genome of P. margaritifera using GSNAP aligner v2017-03-17 (Wu et al. 2016) allowing five mismatches, splicing and using the ‘splitting-output’ function to retain only concordant and unique mapped paired-end reads, as described previously (Le Luyer et al. 2019). We used the QORTs (Hartley and Mullikin 2015) and JunctionSeq R packages (Hartley and Mullikin 2016) to detect significant differences in exon usage. Only exons and junctions with a minimal coverage of six were used for the analysis and only differences with FDR < 0.01 were considered significant.
