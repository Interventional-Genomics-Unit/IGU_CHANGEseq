#!/bin/bash
####
#  QC BAM alignment
# input; analysis folder path, sample name
#
####

SAMPLE = $1
DIR = $2 #analysis dir
echo "${SAMPLE}"
echo "${DIR}"

TXT_IN = ${DIR}/identified/${SAMPLE}_identified_matched.txt
BAM_IN = ${DIR}/aligned/${SAMPLE}_sorted.bam
BAM_OUT = ${DIR}/identified/${SAMPLE}_identified_matched_sorted.bam
HIST = ${DIR}/identified/${SAMPLE}_identified_matched_coverage.txt

cut -f1,2,3,4,5 $TEXT_IN > ${DIR}/identified/${SAMPLE}_identified_matched.bed
samtools view -b -h -L ${DIR}/identified/${SAMPLE}_identified_matched.bed $BAM_IN | samtools sort > $BAM_OUT
samtools coverage -A --histogram $BAM_OUT

