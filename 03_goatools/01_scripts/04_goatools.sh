#!/bin/bash

PATH_DIR=/~/bin/go_enrichment-master

# Global variables
GOATOOLS=find_enrichment.py
FISHER_FOLDER=${PATH_DIR}/06_fisher_tests
GO_DATABASE=${PATH_DIR}/02_go_database/go-basic.obo

ALL_IDS=${PATH_DIR}/05_annotations/all_ids.txt
ANNOTATIONS=${PATH_DIR}/06_fisher_tests/all_go_annotations.csv

SIGNIFICANT_ID_1=~/bin/go_enrichment-master/05_annotations/to_goatools/significant_ids_Jdown.txt
SIGNIFICANT_ID_2=~/bin/go_enrichment-master/05_annotations/to_goatools/significant_ids_JPSdown.txt
SIGNIFICANT_ID_3=~/bin/go_enrichment-master/05_annotations/to_goatools/significant_ids_JPSup.txt
SIGNIFICANT_ID_4=~/bin/go_enrichment-master/05_annotations/to_goatools/significant_ids_Jup.txt
SIGNIFICANT_ID_5=~/bin/go_enrichment-master/05_annotations/to_goatools/significant_ids_Mdown.txt
SIGNIFICANT_ID_6=~/bin/go_enrichment-master/05_annotations/to_goatools/significant_ids_MJdown.txt
SIGNIFICANT_ID_7=~/bin/go_enrichment-master/05_annotations/to_goatools/significant_ids_MJPSdown.txt
SIGNIFICANT_ID_8=~/bin/go_enrichment-master/05_annotations/to_goatools/significant_ids_MJPSup.txt
SIGNIFICANT_ID_9=~/bin/go_enrichment-master/05_annotations/to_goatools/significant_ids_MJup.txt
SIGNIFICANT_ID_10=~/bin/go_enrichment-master/05_annotations/to_goatools/significant_ids_MPSdown.txt
SIGNIFICANT_ID_11=~/bin/go_enrichment-master/05_annotations/to_goatools/significant_ids_MPSup.txt
SIGNIFICANT_ID_12=~/bin/go_enrichment-master/05_annotations/to_goatools/significant_ids_Mup.txt
SIGNIFICANT_ID_13=~/bin/go_enrichment-master/05_annotations/to_goatools/significant_ids_PSdown.txt
SIGNIFICANT_ID_14=~/bin/go_enrichment-master/05_annotations/to_goatools/significant_ids_PSup.txt
ENRICHMENT_1=~/bin/go_enrichment-master/06_fisher_tests/go_enrichment_significant_ids_Jdown.txt
ENRICHMENT_2=~/bin/go_enrichment-master/06_fisher_tests/go_enrichment_significant_ids_JPSdown.txt
ENRICHMENT_3=~/bin/go_enrichment-master/06_fisher_tests/go_enrichment_significant_ids_JPSup.txt
ENRICHMENT_4=~/bin/go_enrichment-master/06_fisher_tests/go_enrichment_significant_ids_Jup.txt
ENRICHMENT_5=~/bin/go_enrichment-master/06_fisher_tests/go_enrichment_significant_ids_Mdown.txt
ENRICHMENT_6=~/bin/go_enrichment-master/06_fisher_tests/go_enrichment_significant_ids_MJdown.txt
ENRICHMENT_7=~/bin/go_enrichment-master/06_fisher_tests/go_enrichment_significant_ids_MJPSdown.txt
ENRICHMENT_8=~/bin/go_enrichment-master/06_fisher_tests/go_enrichment_significant_ids_MJPSup.txt
ENRICHMENT_9=~/bin/go_enrichment-master/06_fisher_tests/go_enrichment_significant_ids_MJup.txt
ENRICHMENT_10=~/bin/go_enrichment-master/06_fisher_tests/go_enrichment_significant_ids_MPSdown.txt
ENRICHMENT_11=~/bin/go_enrichment-master/06_fisher_tests/go_enrichment_significant_ids_MPSup.txt
ENRICHMENT_12=~/bin/go_enrichment-master/06_fisher_tests/go_enrichment_significant_ids_Mup.txt
ENRICHMENT_13=~/bin/go_enrichment-master/06_fisher_tests/go_enrichment_significant_ids_PSdown.txt
ENRICHMENT_14=~/bin/go_enrichment-master/06_fisher_tests/go_enrichment_significant_ids_PSup.txt


# Running goa tools
echo "Running enrichment analysis..."

source $CONDA2/activate goatools-0.6.10

for i in {1..14}
do
	SI=SIGNIFICANT_ID_$i ;
	EN=ENRICHMENT_$i ;
	$GOATOOLS --pval=0.05 --method "bonferroni,sidak,holm,fdr" --indent --obo $GO_DATABASE ${!SI} $ALL_IDS $ANNOTATIONS > ${!EN} ;
done

#$GOATOOLS --pval=0.05 --indent --obo $GO_DATABASE $SIGNIFICANT_IDS $ALL_IDS $ANNOTATIONS > $ENRICHMENT

source deactivate

echo "  --- Please find your results in '$FISHER_FOLDER/go_enrichment.csv' ---"

