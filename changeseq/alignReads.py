"""
alignReads
"""

from __future__ import print_function
import sys
import subprocess
import os
import logging

logger = logging.getLogger('root')
logger.propagate = False

def alignReads(sample_name,BWA_path, HG19_path, read1, read2, outfile):

    #outputfile
    output_folder = os.path.dirname(outfile)
    bam_filename = f"{output_folder}/{sample_name}.bam"
    sorted_bam_file = f"{output_folder}/{sample_name}_sorted.bam"

    # Check if genome is already indexed by bwa
    index_files_extensions = ['.pac', '.amb', '.ann', '.bwt', '.sa']

    genome_indexed = True
    for extension in index_files_extensions:
        if not os.path.isfile(HG19_path + extension):
            genome_indexed = False
            break

    # If the genome is not already indexed, index it
    if not genome_indexed:
        logger.info('Genome index files not detected. Running BWA to generate indices.')
        bwa_index_command = '{0} index {1}'.format(BWA_path, HG19_path)
        logger.info('Running bwa command: %s', bwa_index_command)
        subprocess.call(bwa_index_command.split())
        logger.info('BWA genome index generated')
    else:
        logger.info('BWA genome index found.')

    # Run paired end alignment against the genome
    logger.info('Running paired end mapping for {0}'.format(sample_name))
    bwa_alignment_command = f'{BWA_path} mem -t 48 {HG19_path} {read1} {read2} | samtools view -bS - > {bam_filename}'
    subprocess.call(bwa_alignment_command, shell=True, stdout=sys.stdout, stderr=sys.stderr)


    # samtools sort
    samtools_sort_command = f"samtools view -@ 48 -h -F 4 {bam_filename} | awk '$3 != \"chrM\" || $1 ~ /^@/' | samtools sort -@ 48 -o {sorted_bam_file}"
    samtools_index_command = f'samtools index {sorted_bam_file}'


    # sort bam and index bam
    logger.info('samtools sort bam file')
    logger.info(samtools_sort_command)
    subprocess.call(samtools_sort_command, shell=True,stdout=sys.stdout,stderr=sys.stderr)

    logger.info('samtools index sorted bam file')
    logger.info(samtools_index_command)
    subprocess.call(samtools_index_command, shell=True, stdout=sys.stdout, stderr=sys.stderr)