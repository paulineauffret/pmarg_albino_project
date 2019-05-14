#!/bin/bash

#Global variables
DATA=~/juvenile/results/03_mapping_bwa/mapping_BWA_Trinity.transcriptome.300118.fa	#path to working directory containing bam files (resulting from step 03)
WORKING_DIRECTORY=~/juvenile/results/04_filtering_bamfiles
SCRIPTBASE=~/juveniles/scripts								#path to scripts directory
SCRIPT=${SCRIPTBASE}/04_filtering_bamfiles
HEADER=${SCRIPTBASE}/input_files/header_samtools.txt        				#path to header_samtools.txt file (containing PBS parameters)
SAMTOOLS_ENV=". ~/samtools/latest/env.sh"						#path to samtools version 1.4.1 conda environment if needed
SAMTOOLS="samtools"									#path to samtools version 1.4.1
SAMTOOLS_FILTERING_PARAMS="-F 4 -F 256 -q 5 -f2"					#samtools filtering parameters
EXT="F4_F256_q5_f2"									#filtered bam file extension
LOGDIR=~/juvenile/log

mkdir -p ${WORKING_DIRECTORY}
mkdir -p ${SCRIPT}

cd ${DATA}
for file in $(ls *.bam)
do
	cp ${HEADER} ${SCRIPT}/filter_${file}.qsub ;
	echo "${SAMTOOLS_ENV}" >> ${SCRIPT}/filter_${file}.qsub ;
	echo "${SAMTOOLS} view -b ${SAMTOOLS_FILTERING_PARAMS} ${DATA}/${file} >& ${WORKING_DIRECTORY}/${file%.*}_${EXT}.bam 2> ${LOGDIR}/filter_${file}.log" >> ${SCRIPT}/filter_${file}.qsub ;
	echo "${SAMTOOLS} flagstat ${DATA}/${file%.*}_${EXT}.bam >& ${WORKING_DIRECTORY}/${file%.*}_${EXT}.flagstat 2> ${LOGDIR}/filter_${file}.log" >> ${SCRIPT}/filter_${file}.qsub ;
	echo "${SAMTOOLS} flagstat ${DATA}/${file} >& ${DATA}/${file}.flagstat 2> ${LOGDIR}/filter_${file}.log" >> ${SCRIPT}/filter_${file}.qsub ;
	qsub ${SCRIPT}/filter_${file}.qsub ;
done

