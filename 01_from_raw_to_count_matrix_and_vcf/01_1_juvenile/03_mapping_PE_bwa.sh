#!/bin/bash

#Global variables
ASSEMBLY=~/Trinity.transcriptome.300118.fa				#path to reference transcriptome assembly
INPUT=~/juvenile/scripts/input_files/lst_trimmed_fastq_files.txt	#path to trimmed fastq list file containing absolute path to the 50 trimmed fastq files (resulting from step 02)
WORKING_DIRECTORY=~/juvenile/results/03_mapping_bwa             	#path to working/output directory 
SCRIPT=~/juvenile/scripts                         			#path to scripts directory
MAP_SCRIPT=${SCRIPT}/03_mapping_PE_bwa_scripts				#name of bwa scripts directory
HEADER=${SCRIPT}/input_files/header_bwa.txt                         	#path to header.txt file (containing PBS parameters)
BWA="bwa"                       					#path to bwa version 0.7.15
BWA_ENV="~/bwa_env.sh"                        				#path to bwa version 0.7.15 conda environment if needed
INDEX=1      	                  					#boolean : 1 if transcriptome index exists, else 0       
NB_CPU=16								#number of cpus
NB_SAMPLES=3                      					#number of samples
TAG="mapping_BWA_${ASSEMBLY##*/}"					#output and log tag
LOGDIR=~/juvenile/log							#log directory

#Source trimmed fastq list file
source ${INPUT}

#Creatind index if not existing
if [[ $INDEX == 0 ]]
then
	cp ${HEADER} ${SCRIPT}/create_bwa_index_${ASSEMBLY##*/}.qsub ;
	echo "${BWA_ENV}" >> ${SCRIPT}/create_bwa_index_${ASSEMBLY##*/}.qsub ;
	echo "time ${BWA} index $ASSEMBLY" >> ${SCRIPT}/create_bwa_index_${ASSEMBLY##*/}.qsub ;
	qsub ${SCRIPT}/create_bwa_index_${ASSEMBLY##*/}.qsub ;
else

#Running bwa
mkdir -p $WORKING_DIRECTORY
mkdir -p $WORKING_DIRECTORY/$TAG
mkdir -p $MAP_SCRIPT

cd $WORKING_DIRECTORY/$TAG

for i in `seq 1 $NB_SAMPLES`
do
	R2=LEFT_$i ;
	R1=RIGHT_$i ;
	temp=${!R1##*/} ;
	prefix=${temp%%_R1*} ;
	cp ${HEADER} ${MAP_SCRIPT}/mapping_BWA_${ASSEMBLY##*/}_${prefix}.qsub ;
	echo "${BWA_ENV}" >> ${MAP_SCRIPT}/mapping_BWA_${ASSEMBLY##*/}_${prefix}.qsub ;
	echo "time ${BWA} mem -t ${NB_CPU} -M ${ASSEMBLY} ${!R1} ${!R2} >& ${WORKING_DIRECTORY}/${TAG}/${prefix}.bam 2> ${LOGDIR}/${TAG}_${prefix}.log" >> ${MAP_SCRIPT}/mapping_BWA_${ASSEMBLY##*/}_${prefix}.qsub ;
	qsub ${MAP_SCRIPT}/mapping_BWA_${ASSEMBLY##*/}_${prefix}.qsub ;
done

fi

 
