#!/bin/bash


PATH_DIR=~/bin/go_enrichment-master

# Global variables
SWISSPROT_HITS=${PATH_DIR}/04_blast_results/Trinity.transcriptome.300118.fa_swissprot_hits.txt
ANNOTATION_FOLDER=${PATH_DIR}/05_annotations
FISHER_FOLDER=${PATH_DIR}/06_fisher_tests

# Get info from uniprot for each hit in parallel
cat $SWISSPROT_HITS | while read i; do feature=$(echo $i | cut -d " " -f 1) ; hit=$(echo $i | cut -d "|" -f 2 | cut -d "." -f 1); echo "wget -q -O - http://www.uniprot.org/uniprot/${hit}.txt > $ANNOTATION_FOLDER/${feature}.info"; done > temp_wget_commands.txt



#cat temp_wget_commands.txt | parallel "eval {}"
bash temp_wget_commands.txt
rm temp_wget_commands.txt

# Extract wanted info
for i in $ANNOTATION_FOLDER/*.info; do echo -e "$(basename $i | perl -pe 's/\.info//')\t$(cat $i | grep -E '^DR\s+GO' | awk '{print $3}' | perl -pe 's/\n//')"; done > $FISHER_FOLDER/all_go_annotations.csv

