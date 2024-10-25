from __future__ import print_function
import subprocess
import os
import logging

logger = logging.getLogger('root')
logger.propagate = False


def merged_fastqQC(read1, read2,logfile):
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



def unmerged_fastqQC(read1, read2,read1_out,read2_out,adapter_list,logfile):
    sample_name = os.path.basename(logfile).split('.')[0]
    # Run QC
    logger.info('Running QC for {0}'.format(sample_name))
    fastp_command = 'fastp --adapter_fasta {0} -i {1} -I {2} -o {3} -O {4} --overrepresentation_analysis -h {5}'.format(
        adapter_list,
        read1,
        read2,
        read1_out,
        read2_out,
        logfile
    )
    logger.info(fastp_command)
    subprocess.check_call(fastp_command, shell=True)
    logger.info('QC for {0} completed.'.format(sample_name))