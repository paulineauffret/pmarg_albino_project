#!/bin/bash

PATH_DIR=~/bin/go_enrichment-master

# Global variables
GOATOOLS=find_enrichment.py
FISHER_FOLDER=${PATH_DIR}/06_fisher_tests
GO_DATABASE=${PATH_DIR}/02_go_database/go-basic.obo
#SIGNIFICANT_IDS=${PATH_DIR}/05_annotations/significant_ids_MT_oihana.txt
SIGNIFICANT_IDS=${PATH_DIR}/05_annotations/significant_ids_SP_oihana.txt
ALL_IDS=${PATH_DIR}/05_annotations/all_ids_oihana.txt
ANNOTATIONS=${PATH_DIR}/06_fisher_tests/all_go_annotations_oihana.csv
#ENRICHMENT=${PATH_DIR}/06_fisher_tests/go_enrichment_oihana_MT.csv
ENRICHMENT=${PATH_DIR}/06_fisher_tests/go_enrichment_oihana_SP.csv

# Running goa tools
echo "Running enrichment analysis..."

source $CONDA2/activate goatools-0.6.10

$GOATOOLS --pval=0.05 --indent --obo $GO_DATABASE $SIGNIFICANT_IDS $ALL_IDS $ANNOTATIONS > $ENRICHMENT

source deactivate

echo "  --- Please find your results in '$FISHER_FOLDER/go_enrichment.csv' ---"

