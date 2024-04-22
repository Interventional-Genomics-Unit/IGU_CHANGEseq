from __future__ import print_function
import subprocess
import os
import logging

logger = logging.getLogger('root')
logger.propagate = False


def deduplicateReads(bc_pattern,read1, read2, processed_read1,processed_read2):

    sample_name = os.path.basename(processed_read1).split('_preprocessed.')[0]
    sample_name = sample_name[:-3]
    output_folder = os.path.dirname(processed_read1)
    base_name = os.path.join(output_folder, sample_name)

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Run umi extraction
    logger.info('Running umi extraction for {0}'.format(sample_name))
    extract_command ='umi_tools extract -I {0} --bc-pattern={1} --read2-in={2} --stdout={3} --read2-out={4} --log={5}'.format(
        read1,
        bc_pattern,
        read2,
        processed_read1,
        processed_read2,
        base_name + ".log"
    )

    logger.info(extract_command)
    subprocess.check_call(extract_command, shell=True)
    logger.info('UMI extracting for {0} completed.'.format(sample_name))

