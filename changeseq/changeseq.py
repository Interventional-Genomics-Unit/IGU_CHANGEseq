#!/usr/bin/env python
#-*- coding: utf-8 -*-

"""
circleseq.py as the wrapper for CIRCLE-seq analysis
"""

from alignReads import alignReads
from visualize import visualizeOfftargets
from mergeReads import mergeReads
import argparse
import os
import sys
import subprocess
import traceback
import log
from utility import get_parameters
from copy import deepcopy as dp
import findCleavageSites
import callVariants
from annotate import annotate
from report_qc import write_qc
from multi_sample_combiner import process_results
from fastq_qc import fastqQC,trim_tn5

logger = log.createCustomLogger('root')
global p_dir
p_dir = os.path.dirname(os.path.realpath(__file__))

class CircleSeq:
    def __init__(self):
        self.parameters = {}
        self.output_dir = {}
        self.findCleavageSites_input_bam = {}
        self.vis_input_tsv = {}
        self.annotation_file = {}

    def parseManifest(self,analysis_folder,fq_dir,manifest,settings, sample='all'):
        logger.info('Loading manifest...')

        try:
            parameters = get_parameters(analysis_folder,fq_dir,manifest,settings)
            self.parameters = dp(parameters)

            if sample != 'all':
                self.parameters['samples'] = {}
                self.parameters['samples'][sample] = parameters['samples'][sample]
            # print (self.parameters)
            # Make folders for output
            for folder in ['preprocessed','aligned', 'raw_output', 'fastq', 'variants','qc','processed_output']:
                self.output_dir[folder] = os.path.join(self.parameters["analysis_folder"], folder)
                if not os.path.exists(self.output_dir[folder]):
                    os.makedirs(self.output_dir[folder])

            # Just to initialize some default input file names for the identify and visualization steps
            for sample in self.parameters['samples']:
                self.findCleavageSites_input_bam[sample] = [f"{self.output_dir['aligned']}/{sample}_sorted.bam",
                                                            f"{self.output_dir['aligned']}/control_{sample}_sorted.bam"]

            # initialize some default input file names for annotation files
            for sample in self.parameters['samples']:
                self.annotation_file[sample] = f"{self.output_dir[ 'raw_output']}/{sample}_annotated_results.csv"

        except Exception as e:
            logger.error('Incorrect or malformed manifest file. Please ensure your manifest contains all required fields.')
            sys.exit()

    def processReads(self):

        for sample in self.parameters['samples']:

            # Trim tn5 adapter
            if str(self.parameters['merged_analysis'][0]) != 'T':
                sample_read1_outfile = os.path.join(self.parameters["analysis_folder"], 'preprocessed', sample + '_R1_processesed.fastq.gz')
                sample_read2_outfile = os.path.join(self.parameters["analysis_folder"], 'preprocessed',sample + '_R2_processesed.fastq.gz')
                cutadapt_logfile = os.path.join(self.parameters["analysis_folder"], 'preprocessed', sample + '_trim_log.txt')

                trim_tn5(self.parameters['samples'][sample]['read1'],
                         self.parameters['samples'][sample]['read2'],
                         sample_read1_outfile,
                         sample_read2_outfile,
                         Tn5=self.parameters['changeseq_adapter'],
                         cutadapt=self.parameters['cutadapt'],
                         cutadapt_logfile = cutadapt_logfile)

                self.parameters['samples'][sample]['read1'] = sample_read1_outfile
                self.parameters['samples'][sample]['read2'] = sample_read2_outfile
                cutadapt_logfile = os.path.join(self.parameters["analysis_folder"], 'preprocessed', sample + '_control_trim_log.txt')
                control_read1_outfile = os.path.join(self.parameters["analysis_folder"], 'preprocessed',
                                                    sample + '_R1_processesed_control.fastq.gz')
                control_read2_outfile = os.path.join(self.parameters["analysis_folder"], 'preprocessed',
                                                    sample + '_R2_processesed_control.fastq.gz')
                trim_tn5(self.parameters['samples'][sample]['controlread1'],
                         self.parameters['samples'][sample]['controlread2'],
                         control_read1_outfile,
                         control_read2_outfile,
                         Tn5 = self.parameters['changeseq_adapter'],
                         cutadapt = self.parameters['cutadapt'],
                         cutadapt_logfile = cutadapt_logfile)
                self.parameters['samples'][sample]['controlread1'] = control_read1_outfile
                self.parameters['samples'][sample]['controlread2'] = control_read2_outfile

            if self.parameters['dedup_umi']:

                from deduplicate import deduplicateReads
                logger.info('Deduplicating UMIs...')

                for sample in self.parameters['samples']:
                    sample_processed_read1_path = os.path.join(self.parameters["analysis_folder"], 'preprocessed', sample + '_R1_preprocessed.fastq.gz')
                    sample_processed_read2_path = os.path.join(self.parameters["analysis_folder"], 'preprocessed',sample + '_R2_preprocessed.fastq.gz')
                    control_processed_read1_path = os.path.join(self.parameters["analysis_folder"], 'preprocessed', 'control_' + sample +
                                                               '_R1_preprocessed.fastq.gz')
                    control_processed_read2_path = os.path.join(self.parameters["analysis_folder"], 'preprocessed','control_' +
                                                               sample + '_R2_preprocessed.fastq.gz')
                    deduplicateReads(self.parameters['bc_pattern'],
                                     self.parameters['samples'][sample]['read1'],
                                     self.parameters['samples'][sample]['read2'],
                                     sample_processed_read1_path,
                                     sample_processed_read2_path)
                    self.parameters['samples'][sample]['read1'] = sample_processed_read1_path
                    self.parameters['samples'][sample]['read2'] = sample_processed_read2_path


                    deduplicateReads(self.parameters['bc_pattern'],
                                     self.parameters['samples'][sample]['controlread1'],
                                     self.parameters['samples'][sample]['controlread2'],
                                     control_processed_read1_path,
                                     control_processed_read2_path)
                    self.parameters['samples'][sample]['controlread1'] = control_processed_read1_path
                    self.parameters['samples'][sample]['controlread2'] = control_processed_read2_path


    def alignReads(self):
        if self.parameters['merged_analysis']:
            logger.info('Merging reads...')
            try:
                self.merged = {}
                for sample in self.parameters['samples']:
                    sample_merge_path = os.path.join(self.parameters["analysis_folder"], 'fastq', sample + '_merged.fastq.gz')
                    control_sample_merge_path = os.path.join(self.parameters["analysis_folder"], 'fastq', 'control_' + sample + '_merged.fastq.gz')

                    mergeReads(self.parameters['samples'][sample]['read1'],
                                self.parameters['samples'][sample]['read2'],
                               sample_merge_path)
                    mergeReads(self.parameters['samples'][sample]['controlread1'],
                                 self.parameters['samples'][sample]['controlread2'],
                               control_sample_merge_path)

                    sample_alignment_path = os.path.join(self.parameters["analysis_folder"], 'aligned', sample + '.sam')
                    control_sample_alignment_path = os.path.join(self.parameters["analysis_folder"], 'aligned', 'control_' + sample + '.sam')

                    alignReads(self.parameters['bwa'],
                               self.parameters['reference_genome'],
                               sample_merge_path,
                               '',
                               sample_alignment_path)

                    alignReads(self.parameters['bwa'],
                               self.parameters['reference_genome'],
                               control_sample_merge_path,
                               '',
                               control_sample_alignment_path)

                    self.merged[sample] = sample_alignment_path
                    logger.info('Finished merging and aligning reads.')

            except Exception as e:
                logger.error('Error aligning')
                logger.error(traceback.format_exc())
                quit()
        else:
            try:
                self.aligned = {}
                self.aligned_sorted = {}
                for sample in self.parameters['samples']:
                    sample_alignment_path = os.path.join(self.parameters["analysis_folder"], 'aligned', sample + '.sam')
                    control_sample_alignment_path = os.path.join(self.parameters["analysis_folder"], 'aligned', 'control_' + sample + '.sam')
                    alignReads(self.parameters['bwa'],
                               self.parameters['reference_genome'],
                               self.parameters['samples'][sample]["read1"],
                               self.parameters['samples'][sample]["read2"],
                               sample_alignment_path)
                    alignReads(self.parameters['bwa'],
                               self.parameters['reference_genome'],
                               self.parameters['samples'][sample]['controlread1'],
                               self.parameters['samples'][sample]['controlread2'],
                               control_sample_alignment_path)
                    self.aligned[sample] = sample_alignment_path
                    self.aligned_sorted[sample] = os.path.join(self.parameters["analysis_folder"], 'aligned', sample + '_sorted.bam')
                    logger.info('Finished aligning reads to genome.')
                    self.findCleavageSites_input_bam[sample] = [sample_alignment_path.replace(".sam","_sorted.bam"),
                                                                control_sample_alignment_path.replace(".sam","_sorted.bam")]

            except Exception as e:
                logger.error('Error aligning')
                logger.error(traceback.format_exc())
                quit()

    def findCleavageSites(self):
        logger.info('Identifying off-target cleavage sites.')

        try:
            for sample in self.parameters['samples']:
                bam, control_bam = self.findCleavageSites_input_bam[sample]

                identified_sites_file = os.path.join(self.parameters["analysis_folder"],  'raw_output', sample)

                findCleavageSites.compare(self.parameters['reference_genome'], bam, control_bam, self.parameters['samples'][sample]['target'],
                                          self.parameters['search_radius'], self.parameters['window_size'], self.parameters['mapq_threshold'], self.parameters['gap_threshold'],
                                          self.parameters['start_threshold'], self.parameters['mismatch_threshold'],
                                          sample, self.parameters['samples'][sample]['description'],
                                          identified_sites_file, False,
                                          merged=self.parameters['merged_analysis'],read_count_cutoff=self.parameters['read_threshold'],read_length=151)

        except Exception as e:
            logger.error('Error identifying off-target cleavage site.')
            logger.error(traceback.format_exc())
            quit()

    def addAnnotations(self):
        for sample in self.parameters['samples']:
            try:
                matched_file = f"{self.output_dir[ 'raw_output']}/{sample}_identified_matched.txt"
                print(f"Annotating {matched_file}")
                annotate(matched_file, self.parameters['annotate_path'])
            except Exception as e:
                print('Error Annotating for sample %s.' % (sample))

    def visualize(self):
        logger.info('Visualizing off-target sites')
        for sample in self.parameters['samples']:
            try:
                infile = os.path.join(self.parameters["analysis_folder"],  'raw_output',sample + '_identified_matched_annotated.csv')
                outfile = os.path.join(self.parameters["analysis_folder"],  'raw_output', sample + '_offtargets')
                visualizeOfftargets(infile, outfile, title=sample,PAM=self.parameters["PAM"])
            except Exception as e:
                logger.error('Error visualizing off-target sites.')
                logger.error(traceback.format_exc())

    def analyze(self):
        logger.info('Normalizing and Joining reads')
        for rep_group_name,replicates in self.parameters['replicates'].items():
            infiles = []
            qcfiles = []
            for sample in self.parameters['replicates'][rep_group_name]['sample_name']:
                infiles.append(os.path.join(self.parameters["analysis_folder"], 'raw_output',
                                            sample + '_identified_matched_annotated.csv'))
                qcfiles.append(os.path.join(self.parameters["analysis_folder"], 'qc', sample + '_qc_report.txt'))
            logger.info('Normalizing {rep_group_name}')

            outfolder = os.path.join(self.parameters["analysis_folder"],'processed_output')
            process_results(rep_group_name, replicates, infiles, qcfiles,
                            outfolder=outfolder,
                            normalization_method=self.parameters['normalize'],
                            read_threshold = int(self.parameters['read_threshold']), PAM=self.parameters["PAM"])


    def QC(self):
        logger.info('Starting QC')
        try:
            for sample in self.parameters['samples']:

                logger.info('Running fastq quality and adapter analysis for {0}'.format(sample))
                fastqc_logfile = os.path.join(self.parameters["analysis_folder"], 'preprocessed', sample + '.html')

                fastqQC(self.parameters['samples'][sample]['read1'],
                        self.parameters['samples'][sample]['read2'],
                        fastqc_logfile)

        except Exception as e:
            logger.error('Error with Fastqc')

        try:
            for sample in self.parameters['samples']:
                #if self.merged_analysis:
                script_path = p_dir + "/QC_matched_alignment.sh"
                logger.info('Running alignment coverage QC for {0}'.format(sample))
                coverage_command = 'sh {0} {1} {2}'.format(script_path,
                     sample, self.parameters["analysis_folder"])
                logger.info(coverage_command)
                subprocess.check_call(coverage_command, shell=True)
                logger.info('OT Coverage for {0} completed.'.format(sample))

        except Exception as e:
            logger.error('Error with alignment coverage QC ')

        try:
            for sample in self.parameters['samples']:
                logger.info('Running  report  for {0}'.format(sample))
                preprocessed_logfile =  os.path.join(self.parameters["analysis_folder"], 'preprocessed', sample + '_trim_log.txt')
                coverage_stat_file  = os.path.join(self.parameters["analysis_folder"], 'qc', sample + '_aligned_stats.txt')
                qc_file = os.path.join(self.parameters["analysis_folder"], 'qc', sample + '_qc_report.txt')

                write_qc(qc_file, preprocessed_logfile, coverage_stat_file)

        except Exception as e:
            logger.error('Error with QC report')


    def callVariants(self):

        try:
            if self.parameters['variant_analysis']:
                logger.info('Identifying genomic variants')

                for sample in self.parameters['samples']:
                    sorted_bam_file = os.path.join(self.parameters['analysis_folder'], 'coverage', sample + '_identified_matched_sorted.bam')
                    identified_sites_file = os.path.join(self.parameters['analysis_folder'], 'identified', sample + '_identified_matched.txt')
                    variants_basename = os.path.join(self.parameters['analysis_folder'], 'variants', sample)

                    callVariants.getVariants(identified_sites_file, self.parameters['reference_genome'], sorted_bam_file,
                                             variants_basename, self.parameters['search_radius'], self.parameters['mismatch_threshold'])

                    try:
                        annotated_file = variants_basename +  '_Variants.txt'
                        annotate(annotated_file, annotate_path=self.parameters['annotate_file'])
                    except Exception as e:
                        logger.error('Error annotating genomic variants.')


                logger.info('Finished identifying genomic variants')

        except Exception as e:
            logger.error('Error identifying genomic variants.')
            logger.error(traceback.format_exc())
            quit()

    def parallel(self, manifest_path, lsf, run='all'):
        logger.info('Submitting parallel jobs')
        current_script = __file__

        try:
            for sample in self.parameters['samples']:
                cmd = 'python {0} {1} --manifest {2} --sample {3}'.format(current_script, run, manifest_path, sample)
                logger.info(cmd)
                subprocess.call(lsf.split() + [cmd])
            logger.info('Finished job submission')

        except Exception as e:
            logger.error('Error submitting jobs.')
            logger.error(traceback.format_exc())

    def referenceFree(self):
        pass

def parse_args():
    parser = argparse.ArgumentParser()

    subparsers = parser.add_subparsers(description='Individual Step Commands',
                                       help='Use this to run individual steps of the pipeline',
                                       dest='command')

    all_parser = subparsers.add_parser('all', help='Run all steps of the pipeline')
    all_parser.add_argument('--manifest', '-m', help='The file name of the sample manifest(.csv)', required=True)
    all_parser.add_argument('--raw_fastq_folder', '-fq', help='Directory where raw fastq file are stored (abs path)', required=True)
    all_parser.add_argument('--analysis_folder', '-dir', help='Directory in which all pipeline outputs will be saved (aabs path)', required=True)
    all_parser.add_argument('--settings', '-stg', help='The file name of the parameters and settings are (.csv)', default = 'default')
    all_parser.add_argument('--sample', '-s', help='Specify sample to process (default is all)', default='all')

    parallel_parser = subparsers.add_parser('parallel', help='Run all steps of the pipeline in parallel')
    parallel_parser.add_argument('--manifest', '-m', help='The file name of the sample manifest(.csv)', required=True)
    parallel_parser.add_argument('--raw_fastq_folder', '-fq', help='Directory where raw fastq file are stored (abs path)', required=True)
    parallel_parser.add_argument('--analysis_folder', '-dir', help='Directory in which all pipeline outputs will be saved (aabs path)', required=True)
    parallel_parser.add_argument('--settings', '-stg', help='The file name of the parameters and settings are (.csv)', default = 'default')
    parallel_parser.add_argument('--lsf', '-l', help='Specify LSF CMD', default='bsub -R rusage[mem=32000] -P Genomics -q standard')
    parallel_parser.add_argument('--run', '-r', help='Specify which steps of pipepline to run (all, align, identify, visualize, variants)', default='all')

    align_parser = subparsers.add_parser('align', help='Run alignment only')
    align_parser.add_argument('--manifest', '-m', help='The file name of the sample manifest(.csv)', required=True)
    align_parser.add_argument('--raw_fastq_folder', '-fq', help='Directory where raw fastq file are stored (abs path)', required=True)
    align_parser.add_argument('--analysis_folder', '-dir', help='Directory in which all pipeline outputs will be saved (aabs path)', required=True)
    align_parser.add_argument('--settings', '-stg', help='The file name of the parameters and settings are (.csv)', default = 'default')
    align_parser.add_argument('--sample', '-s', help='Specify sample to process (default is all)', default='all')

    data_parser = subparsers.add_parser('makefiles', help='Combine and normalize replicates, produces vizualations')
    data_parser.add_argument('--outdir', '-o',help='Specify the final annnotation table location. Default is changeseq directory',default=os.path.dirname(os.path.realpath(__file__)) + "/data/")
    data_parser.add_argument('--ftp_path', '-f', help='RefSeq FTP Path. Must be in refseq.txt format', default="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ncbiRefSeq.txt.gz")
    data_parser.add_argument('--reset_output', '-r', help='Change path changeseq uses for stored annotation file',default=False)

    fa_parser = subparsers.add_parser('set_fasta', help='Sets the default genome fasta path so it is not needed in sample manifest')
    fa_parser.add_argument('-fasta', '-fa',
                             help='Specify full fasta path. make sure this is unzipped')

    merge_parser = subparsers.add_parser('merge', help='Merge paired end reads')
    merge_parser.add_argument('--manifest', '-m', help='The file name of the sample manifest(.csv)', required=True)
    merge_parser.add_argument('--raw_fastq_folder', '-fq', help='Directory where raw fastq file are stored (abs path)', required=True)
    merge_parser.add_argument('--analysis_dir', '-dir', help='Directory in which all pipeline outputs will be saved (aabs path)', required=True)
    merge_parser.add_argument('--settings', '-stg', help='The file name of the parameters and settings are (.csv)', default = 'default')
    merge_parser.add_argument('--sample', '-s', help='Specify sample to process (default is all)', default='all')


    identify_parser = subparsers.add_parser("identify", help='Run identification only')
    identify_parser.add_argument('--manifest', '-m', help='The file name of the sample manifest(.csv)', required=True)
    identify_parser.add_argument('--raw_fastq_folder', '-fq', help='Directory where raw fastq file are stored (abs path)', required=True)
    identify_parser.add_argument('--analysis_folder', '-dir', help='Directory in which all pipeline outputs will be saved (aabs path)', required=True)
    identify_parser.add_argument('--settings', '-stg', help='The file name of the parameters and settings are (.csv)', default = 'default')
    identify_parser.add_argument('--sample', '-s', help='Specify sample to process (default is all)', default='all')

    visualize_parser = subparsers.add_parser('visualize', help='Run visualization only')
    visualize_parser.add_argument('--manifest', '-m', help='The file name of the sample manifest(.csv)', required=True)
    visualize_parser.add_argument('--raw_fastq_folder', '-fq', help='Directory where raw fastq file are stored (abs path)', required=True)
    visualize_parser.add_argument('--analysis_folder', '-dir', help='Directory in which all pipeline outputs will be saved (abs path)', required=True)
    visualize_parser.add_argument('--settings', '-stg', help='The file name of the parameters and settings are (.csv)', default = 'default')
    visualize_parser.add_argument('--sample', '-s', help='Specify sample to process (default is all)', default='all')

    variants_parser = subparsers.add_parser('variants', help='Run variants analysis only')
    variants_parser.add_argument('--manifest', '-m', help='The file name of the sample manifest(.csv)', required=True)
    variants_parser.add_argument('--raw_fastq_folder', '-fq', help='Directory where raw fastq file are stored (abs path)', required=True)
    variants_parser.add_argument('--analysis_folder', '-dir', help='Directory in which all pipeline outputs will be saved (abs path)', required=True)
    variants_parser.add_argument('--settings', '-stg', help='The file name of the parameters and settings are (.csv)', default = 'default')
    variants_parser.add_argument('--sample', '-s', help='Specify sample to process (default is all)', default='all')

    analyzer_parser = subparsers.add_parser('analyze', help='Run analysis  only')
    analyzer_parser.add_argument('--manifest', '-m', help='The file name of the sample manifest(.csv)', required=True)
    analyzer_parser.add_argument('--raw_fastq_folder', '-fq', help='Directory where raw fastq file are stored (abs path)',required=True)
    analyzer_parser.add_argument('--analysis_folder', '-dir',
                            help='Directory in which all pipeline outputs will be saved (aabs path)', required=True)
    analyzer_parser.add_argument('--settings', '-stg', help='The file name of the parameters and settings are (.csv)',default='default')
    analyzer_parser.add_argument('--sample', '-s', help='Specify sample to process (default is all)', default='all')

    coverage_parser = subparsers.add_parser('qc', help='Run QC analysis of matched sites')
    coverage_parser.add_argument('--manifest', '-m', help='The file name of the sample manifest(.csv)', required=True)
    coverage_parser.add_argument('--raw_fastq_folder', '-fq', help='Directory where raw fastq file are stored (abs path)', required=True)
    coverage_parser.add_argument('--analysis_folder', '-dir', help='Directory in which all pipeline outputs will be saved (aabs path)', required=True)
    coverage_parser.add_argument('--settings', '-stg', help='The file name of the parameters and settings are (.csv)', default = 'default')
    coverage_parser.add_argument('--sample', '-s', help='Specify sample to process (default is all)', default='all')

    reference_free_parser = subparsers.add_parser('reference-free', help='Run reference-free discovery only')
    reference_free_parser.add_argument('--manifest', '-m', help='The file name of the sample manifest(.csv)', required=True)
    reference_free_parser.add_argument('--raw_fastq_folder', '-fq', help='Directory where raw fastq file are stored (abs path)', required=True)
    reference_free_parser.add_argument('--analysis_dir', '-dir', help='Directory in which all pipeline outputs will be saved (aabs path)', required=True)
    reference_free_parser.add_argument('--settings', '-stg', help='The file name of the parameters and settings are (.csv)', default = 'default')
    reference_free_parser.add_argument('--sample', '-s', help='Specify sample to process (default is all)', default='all')

    return parser.parse_args()

def main():
    args = parse_args()

    if args.command == 'all':
        c = CircleSeq()
        c.parseManifest(args.analysis_folder,args.raw_fastq_folder,args.manifest,args.settings, args.sample)
        c.processReads()
        c.alignReads()
        c.findCleavageSites()
        c.addAnnotations()
        c.visualize()
        c.analyze()
        c.QC()
        c.callVariants()
    elif args.command == 'parallel':
        c = CircleSeq()
        c.parseManifest(args.analysis_folder,args.raw_fastq_folder,args.manifest,args.settings, args.sample)
        c.parallel(args.manifest, args.lsf, args.run)
    elif args.command == 'merge':
        c = CircleSeq()
        c.parseManifest(args.analysis_folder,args.raw_fastq_folder,args.manifest,args.settings, args.sample)
        c.alignReads()
    elif args.command == 'align':
        c = CircleSeq()
        c.parseManifest(args.analysis_folder,args.raw_fastq_folder,args.manifest,args.settings, args.sample)
        c.alignReads()
    elif args.command == 'identify':
        c = CircleSeq()
        c.parseManifest(args.analysis_folder,args.raw_fastq_folder,args.manifest,args.settings, args.sample)
        c.findCleavageSites()
        c.addAnnotations()
        c.visualize()
    elif args.command == 'visualize':
        c = CircleSeq()
        c.parseManifest(args.analysis_folder,args.raw_fastq_folder,args.manifest,args.settings, args.sample)
        c.addAnnotations()
        c.visualize()
    elif args.command == 'analyze':
        c = CircleSeq()
        c.parseManifest(args.analysis_folder,args.raw_fastq_folder,args.manifest,args.settings, args.sample)
        c.analyze()
    elif args.command == 'variants':
        c = CircleSeq()
        c.parseManifest(args.analysis_folder,args.raw_fastq_folder,args.manifest,args.settings, args.sample)
        c.callVariants()
    elif args.command == 'qc':
        c = CircleSeq()
        c.parseManifest(args.analysis_folder,args.raw_fastq_folder,args.manifest,args.settings, args.sample)
        c.QC()
    elif args.command == 'makefiles':
        from make_annotation_table import makefiles
        makefiles(args.ftp_path,args.outdir,args.reset_output,p_dir)

if __name__ == '__main__':
    main()
