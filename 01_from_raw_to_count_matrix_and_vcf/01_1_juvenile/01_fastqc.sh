#!/bin/bash

DATA_DIRECTORY=~/juvenile/raw				#path to raw data directory (raw fastq.gz files, each pair R1 - R2 placed into one subdirectory)
QC_SCRIPT=~/juvenile/scripts/01_fastqc.qsub		#path to 01_2_fastqc_raw_reads.qsub script
SCRIPT=~/juvenile/scripts				#path to scripts directory
HEADER=${SCRIPT}/input_files/header.txt			#path to header.txt file


#Loop into each subdirectory of DATA_DIRECTORY containing pairs of raw fastq files and run QC script
cd $DATA_DIRECTORY
for dir in $(ls) ;
do
	cd $dir ;
	FILE_R1=$(ls *_R1.fastq.gz) ;
	FILE_R2=$(ls *_R2.fastq.gz) ;
	echo "cd ${DATA_DIRECTORY}/${dir}" >> ${SCRIPT}/tmp ;
	echo "FILE_R1=${DATA_DIRECTORY}/${dir}/${FILE_R1}" >> ${SCRIPT}/tmp ;
	echo "FILE_R2=${DATA_DIRECTORY}/${dir}/${FILE_R2}" >> ${SCRIPT}/tmp ;
	cat $HEADER ${SCRIPT}/tmp $QC_SCRIPT > $SCRIPT/${QC_SCRIPT##*/}_${dir}.qsub
	qsub ${SCRIPT}/${QC_SCRIPT##*/}_${dir}.qsub ;
	rm ${SCRIPT}/tmp ;
	cd .. ;
done
