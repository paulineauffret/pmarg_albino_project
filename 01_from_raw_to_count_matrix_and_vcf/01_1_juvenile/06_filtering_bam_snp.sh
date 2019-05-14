#!/bin/bash

data=~/juvenile/results/04_filtering_bamfiles						#path to input bam files (resulting from step 04)
script=~/juvenile/scripts/06_filtering_bam_snp_scripts/06_filtering_bam_snp.qsub	#path to 06_filtering_bam_snp.qsub script

for file in $(ls ${data}/*.bam)
do
	base=${file##*/}
	name=${base%%_F4*}
        toEval="cat ${script} | sed 's/__BASE__/$base/g'"; eval $toEval > ${script%.*}_$name.qsub 
	qsub ${script%.*}_$name.qsub
done
