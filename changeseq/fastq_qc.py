from __future__ import print_function
import subprocess
import os
import logging
import sys
import regex as re

logger = logging.getLogger('root')
logger.propagate = False

def parse_cutadapt_log(sample_name,logfile_path):
    """
    Parse an existing Cutadapt logfile to extract read statistics
    to log how many were trimmed, filtered etc.
    """

    try:
        with open(logfile_path, "r") as f:
            log_text = f.read()
    except Exception as e:
        logger.error(f"Could not read cutadapt logfile: {e}")
        return None

    # Regex patterns used by Cutadapt’s summary section
    reads_input = re.search(r"Total reads processed for:\s+([\d,]+)", log_text)
    reads_written = re.search(r"Reads written \(passing filters\):\s+([\d,]+)", log_text)

    if not reads_input or not reads_written:
        logger.error("Could not parse read count information from the cutadapt logfile.")
        return None

    input_reads = int(reads_input.group(1).replace(",", ""))
    output_reads = int(reads_written.group(1).replace(",", ""))

    removed = input_reads - output_reads
    percent_removed = removed / input_reads if input_reads > 0 else 0.0

    # Logging (your logger will handle formatting)
    logger.info(f"Reads before trimming {sample_name}: {input_reads}")
    logger.info(f"Reads after trimming {sample_name}: {output_reads}")
    logger.info(f"Reads removed {sample_name}: {removed} ({percent_removed:.2%})")

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

def trim_tn5(read1, read2, read1_out, read2_out, Tn5, cutadapt,cutadapt_logfile, 
             trim_options = "-Z --overlap 35 -e 0.15 -m 30 -j 48"):

    # Run QC
    # -Z : Enables "zero-cap" mode for trimming with 3' adapters (no re-trimming).
    # --overlap 35 : Minimum 35 bp overlap with adapter for trimming to occur. (adapter are 42bp)
    # -e 0.15 : Maximum allowed error rate (15%).
    # -m 25 : discard reads shorter than 25bp after trimming
    sample_name = os.path.basename(cutadapt_logfile).split('.')[0]
    trimming_input = f"-a {Tn5} -A {Tn5}"
    trimmed_reads_PE_output = f"-o {read1_out} -p {read2_out} "

    cutadapt_command = f"{cutadapt} {trimming_input} {trimmed_reads_PE_output} {trim_options} {read1} {read2} > {cutadapt_logfile}"
    logger.info(cutadapt_command)
    subprocess.call(cutadapt_command, shell=True, stdout=sys.stdout, stderr=sys.stderr)
    parse_cutadapt_log(sample_name, cutadapt_logfile)





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