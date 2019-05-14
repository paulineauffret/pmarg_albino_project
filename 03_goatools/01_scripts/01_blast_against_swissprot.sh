#!/bin/bash

# Path to software directory
PATH_DIR=~/bin/go_enrichment-master

# Global variables
SEQUENCE_FILE=${PATH_DIR}/03_sequences/Trinity.transcriptome.300118.fa
SWISSPROT_RESULT=${PATH_DIR}/04_blast_results/Trinity.transcriptome.300118.fa_swissprot_results.txt
SWISSPROT_HITS=${PATH_DIR}/04_blast_results/Trinity.transcriptome.300118.fa_swissprot_hits.txt
SWISSPROT_DB=/db/blast/all/uniprot_swissprot	

# Blast all sequences against swissprot (must be installed locally)
# WARNING use `-j N` if you need to limit the number of CPUs used to N <integer>
cat $SEQUENCE_FILE | parallel --gnu -k --block 1k -j 16 --recstart '>' --pipe 'blastx -db '$SWISSPROT_DB' -query - -evalue 1e-3 -outfmt 6 -max_target_seqs 1' > $SWISSPROT_RESULT

#blastx -db $SWISSPROT_DB -query $SEQUENCE_FILE -evalue 1e-3 -outfmt 6 -max_target_seqs 1 > $SWISSPROT_RESULT

# TODO filter blasts on similarity, evalue, length...

# Extract analyzed_genes.hits
awk '{print $1,$2}' $SWISSPROT_RESULT | uniq > $SWISSPROT_HITS

