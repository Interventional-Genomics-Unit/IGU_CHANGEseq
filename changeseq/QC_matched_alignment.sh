#!/bin/bash
####
#  QC BAM alignment
# input; analysis folder path, sample name
#
####


echo "SAMPLE: $1"
echo "CONTROL: $2"
echo "DIR: $3"


#input
TXT_IN=${3}/raw_results/tables/${1}_identified_matched.txt
STBAM_IN=${3}/aligned/${1}_sorted.bam
BAM_IN=${3}/aligned/${1}.bam
CTL_BAM_IN=${3}/aligned/${2}.bam

#output
BAM_OUT=${3}/qc/${1}_identified_matched_sorted.bam
HIST=${3}/qc/${1}_identified_matched_coverage.txt
STATS=${3}/qc/${1}_aligned_stats.txt
CTL_STATS=${3}/qc/${2}_aligned_stats.txt

cut -f1,2,3,4,5 $TXT_IN > ${3}/raw_results/tables/${1}_identified_matched.bed
samtools view -b -h -L ${3}/raw_results/tables/${1}_identified_matched.bed $STBAM_IN | samtools sort > $BAM_OUT
samtools index $BAM_OUT
#samtools coverage -A --histogram $BAM_OUT > $HIST
samtools stats $BAM_IN | grep ^SN | cut -f 2- > $STATS
samtools stats $CTL_BAM_IN | grep ^SN | cut -f 2- > $CTL_STATS
