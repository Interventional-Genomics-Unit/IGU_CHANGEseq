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
    logger.info('Running Fast_QC for {0}'.format(sample_name))
    fastp_command = 'fastp -i {0} -I {1} --overrepresentation_analysis -h {2}'.format(
        read1,
        read2,
        logfile
    )
    logger.info(fastp_command)
    subprocess.check_call(fastp_command, shell=True)
    logger.info('QC for {0} completed.'.format(sample_name))

def trim_tn5(read1, read2, read1_out, read2_out, Tn5, cutadapt,cutadapt_logfile):

    # Run QC
    trimming_input = f"-a {Tn5} -A {Tn5}"
    trimmed_reads_PE_output = f"-o {read1_out} -p {read2_out}"
    other_options = f" -Z --overlap 40 -e 0.15" # --info-file {logfile_base}.cutadapt_info_file.tsv"

    cutadapt_command = f"{cutadapt} {trimming_input} {trimmed_reads_PE_output} {other_options} {read1} {read2} > {cutadapt_logfile}"
    logger.info(cutadapt_command)
    subprocess.call(cutadapt_command, shell=True, stdout=sys.stdout, stderr=sys.stderr)



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