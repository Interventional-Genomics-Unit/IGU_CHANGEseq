from __future__ import print_function
import subprocess
import os
import logging

logger = logging.getLogger('root')
logger.propagate = False


def fastqQC(read1, read2,outfile):
    sample_name = os.path.basename(outfile).split('.')[0]
    # Run QC
    logger.info('Running QC for {0}'.format(sample_name))
    fastp_command = 'fastp -i {0} -I {1} --overrepresentation_analysis -h {2}'.format(
        read1,
        read2,
        outfile
    )
    logger.info(fastp_command)
    subprocess.check_call(fastp_command, shell=True)
    logger.info('QC for {0} completed.'.format(sample_name))