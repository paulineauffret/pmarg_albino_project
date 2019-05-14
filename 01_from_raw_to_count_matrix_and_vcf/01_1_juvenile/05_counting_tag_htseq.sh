#!/bin/bash

WD=~/juvenile/results/03_mapping_bwa/mapping_BWA_Trinity.transcriptome.300118.fa			#path to filtered bam files directory (resulting from step 03)
OUT=~/juvenile/results/03_mapping_bwa/mapping_BWA_Trinity.transcriptome.300118.fa/counting_htseq	#path to ouptut directory
SCRIPT=~/juvenile/scripts                 								#path to scripts directory
HEADER=~/juvenile/scripts/input_files/header.txt                 					#path to header.txt file (containing PBS parameters)
GFF=~/juvenile/scripts/input_files/Trinity.100aaorf.minexpr0.5.gff3					#path to gff file
EXT="F4_F256_q5_f2.bam"											#filtered bam extension
SAMTOOLS_ENV=". ~/samtools/latest/env.sh"   								#path to samtools version 1.4.1 conda environment if needed
SAMTOOLS="samtools"       										#path to samtools version 1.4.1
HTSEQ_ENV=". ~/htseq-count/latest/env.sh"								#path to htseq version 0.6.1 conda environment if needed
HTSEQ_COUNT="htseq-count"										#path to htseq-count version 0.6.1

cd ${WD}

mkdir -p ${OUT}

for file in $(ls *${EXT})
do
	cp ${HEADER} ${SCRIPT}/counting_htseq_${file}.qsub ;
	echo "cd ${WD}" >> ${SCRIPT}/counting_htseq_${file}.qsub ;
	echo "${SAMTOOLS_ENV}" >> ${SCRIPT}/counting_htseq_${file}.qsub ;
	echo "${SAMTOOLS} sort ${file} > ${file%.*}_sorted.bam" >> ${SCRIPT}/counting_htseq_${file}.qsub ;
	echo "${SAMTOOLS} index ${file%.*}_sorted.bam > ${file%.*}_sorted.bam.bai" >> ${SCRIPT}/counting_htseq_${file}.qsub ;
	echo "${HTSEQ_ENV}" >> ${SCRIPT}/counting_htseq_${file}.qsub ;
	echo "${HTSEQ_COUNT} -f 'bam' -s 'no' -r 'pos' -t 'CDS' -i 'Name' ${file%.*}_sorted.bam ${GFF} >> ${OUT}/${file}_htseq-count.txt" >> ${SCRIPT}/counting_htseq_${file}.qsub ;
	qsub ${SCRIPT}/counting_htseq_${file}.qsub ;
done

