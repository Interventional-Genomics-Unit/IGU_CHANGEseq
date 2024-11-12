from __future__ import print_function
import subprocess
import os
import logging
import sys

logger = logging.getLogger('root')
logger.propagate = False


def fastqQC(read1, read2,logfile):
    sample_name = os.path.basename(logfile).split('.')[0]
    # Run QC
    logger.info('Running QC for {0}'.format(sample_name))
    fastp_command = 'fastp -i {0} -I {1} --overrepresentation_analysis -h {2}'.format(
        read1,
        read2,
        logfile
    )
    logger.info(fastp_command)
    subprocess.check_call(fastp_command, shell=True)
    logger.info('QC for {0} completed.'.format(sample_name))

def summarize_cutadapt(logfile_base):
    pass
    report = {logfile_base}.cutadapt_info_file.tsv
    with open(report) as f:
        for line in f:
            line = line.split("\t")




def trim_tn5(read1, read2, read1_out, read2_out, Tn5, cutadapt,logfile_base):
    sample_name = os.path.basename(logfile_base)

    # Run QC
    trimming_input = f"-a {Tn5} -A {Tn5}"
    trimmed_reads_PE_output = f"-o {read1_out} -p {read2_out}"
    # too_short_reads_PE_output = f"--too-short-output {output_dir}/{label}.too_short.R1.fastq.gz --too-short-paired-output {output_dir}/{label}.too_short.R2.fastq.gz"
    # untrimmed_reads_PE_output = f"--untrimmed-output {output_dir}/{label}.untrimmed.R1.fastq.gz --untrimmed-paired-output {output_dir}/{label}.untrimmed.R2.fastq.gz"
    # other_options = f"-j {njobs} -Z -m 10 --overlap 20 -e 0.1 --info-file {output_dir}/{label}.cutadapt_info_file.tsv --pair-filter both"
    other_options = f" -Z --overlap 40 -e 0.15 --info-file {logfile_base}.cutadapt_info_file.tsv"

    # cutadapt_command = f"{cutadapt} {trimming_input} {trimmed_reads_PE_output} {too_short_reads_PE_output} {untrimmed_reads_PE_output} {other_options} {R1} {R2}"
    cutadapt_command = f"{cutadapt} {trimming_input} {trimmed_reads_PE_output} {other_options} {read1} {read2} > {logfile_base}.txt"
    logger.info(cutadapt_command)
    subprocess.call(cutadapt_command, shell=True, stdout=sys.stdout, stderr=sys.stderr)
    logger.info('Running QC for {0}'.format(sample_name))

    logger.info('QC for {0} completed.'.format(sample_name))


#def unmerged_fastqQC(read1, read2,read1_out,read2_out,adapter_list,logfile):
#    sample_name = os.path.basename(logfile).split('.')[0]
    # Run QC
    #logger.info('Running QC for {0}'.format(sample_name))
    #fastp_command = 'fastp --adapter_fasta {0} -i {1} -I {2} -o {3} -O {4} --overrepresentation_analysis -h {5}'.format(
    #    adapter_list,
    #    read1,
    #    read2,
    #    read1_out,
    #    read2_out,
    #    logfile
    #)
    #logger.info(fastp_command)
    #subprocess.check_call(fastp_command, shell=True)
    #logger.info('QC for {0} completed.'.format(sample_name))