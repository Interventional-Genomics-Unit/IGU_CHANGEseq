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
import yaml
import validation
import findCleavageSites
import callVariants
from annotate import annotate
from deduplicate import deduplicateReads
from fastq_qc import merged_fastqQC,unmerged_fastqQC

logger = log.createCustomLogger('root')
global p_dir
p_dir = os.path.dirname(os.path.realpath(__file__))

class CircleSeq:

    def __init__(self):
        self.search_radius = 20
        self.window_size = 3
        self.mapq_threshold = 50
        self.start_threshold = 1
        self.gap_threshold = 3
        self.mismatch_threshold = 6
        self.read_threshold = 6
        self.merged_analysis = True
        self.all_chromosomes = False
        self.variant_analysis = False
        self.dedup_umi = False
        self.genome = None
        self.refseq_names = None


    def parseManifest(self, manifest_path, sample='all'):
        logger.info('Loading manifest...')

        with open(manifest_path, 'r') as f:
            manifest_data = yaml.load(f,Loader=yaml.FullLoader)

        try:
            # Validate manifest data
            validation.validateManifest(manifest_data)

            self.BWA_path  = manifest_data['bwa']
            self.reference_genome = manifest_data['reference_genome']
            self.analysis_folder = manifest_data['analysis_folder']
            self.raw_fastq_folder = manifest_data['raw_fastq_folder']
            # Allow the user to specify read threshold, window_size and search_radius if they'd like
            if os.path.isfile(p_dir + "/changeseq/data/paths.txt"):
                self.annotate_file = open(p_dir + "/changeseq/data/paths.txt", "r").readlines()[0]

            elif os.path.isfile(p_dir + "/data/paths.txt"):
                self.annotate_file = open(p_dir + "/data/paths.txt", "r").readlines()[0]
            else:
                print("No Annotation file path found in /changeseq/data/paths.txt")
                self.annotate_file = None

            if os.path.isfile(p_dir + "/changeseq/data/adapter.fa"):
                self.adapter_list = p_dir + "/changeseq/data/adapter.fa"
            elif os.path.isfile(p_dir + "/data/adapter.fa"):
                self.adapter_list = p_dir + "/data/adapter.fa"
            else:
                print("No adapter.fa for in  in /changeseq/data/")
                quit()

            if 'search_radius' in manifest_data:
                self.search_radius = manifest_data['search_radius']
            if 'window_size' in manifest_data:
                self.window_size = manifest_data['window_size']
            if 'mapq_threshold' in manifest_data:
                self.mapq_threshold = manifest_data['mapq_threshold']
            if 'start_threshold' in manifest_data:
                self.start_threshold = manifest_data['start_threshold']
            if 'gap_threshold' in manifest_data:
                self.gap_threshold = manifest_data['gap_threshold']
            if 'mismatch_threshold' in manifest_data:
                self.mismatch_threshold = manifest_data['mismatch_threshold']
            if 'read_threshold' in manifest_data:
                self.read_threshold = manifest_data['read_threshold']
            if 'merged_analysis' in manifest_data:
                self.merged_analysis = manifest_data['merged_analysis']
            if 'all_chromosomes' in manifest_data:
                self.all_chromosomes = manifest_data['all_chromosomes']
            if 'variant_analysis' in manifest_data:
                self.variant_analysis = manifest_data['variant_analysis']
            if 'read_length' in manifest_data:
                self.read_length = manifest_data['read_length']
            else:
                self.read_length = 151
            if 'dedup_umi' in manifest_data:
                self.dedup_umi = manifest_data['dedup_umi']
                if manifest_data['dedup_umi'] == True:
                    self.bc_pattern = manifest_data['bc_pattern']
                    self.read_length -= len(self.bc_pattern)

            if 'genome' in manifest_data:
                self.genome = manifest_data['genome']
                if self.genome in ['hg38','hg19']:
                    self.refseq_names = p_dir+"/refseq_gene_name.py"
            # Allow the user to specify PAM seq. Yichao 4/29/2020
            if 'PAM' in manifest_data:
                self.PAM = manifest_data['PAM']
            else:
                self.PAM = "NGG"

            # Allow the user to specify Read Count cutoff. Yichao 4/29/2020
            if 'read_count_cutoff' in manifest_data:
                self.read_count_cutoff = manifest_data['read_count_cutoff']
            else:
                self.read_count_cutoff = 6

            # Do not allow to run variant_analysis with merged_analysis
            if self.merged_analysis and self.variant_analysis:
                logger.error('merged_analysis is not compatible with variant_analysis. Please remove one option.')
                sys.exit()

            if sample == 'all':
                self.samples = manifest_data['samples']
            else:
                self.samples = {}

                try:
                    self.samples[sample] = manifest_data['samples'][sample]
                except KeyError as e:
                    logger.error(sample + " is not a sample in the manifest file")
                    logger.error(traceback.format_exc())
                    quit()
            # Make folders for output
            for folder in ['preprocessed','aligned', 'identified', 'fastq', 'visualization', 'variants','coverage']:
                output_folder = os.path.join(self.analysis_folder, folder)
                if not os.path.exists(output_folder):
                    os.makedirs(output_folder)

        except Exception as e:
            logger.error('Incorrect or malformed manifest file. Please ensure your manifest contains all required fields.')
        #    sys.exit()

    def processReads(self):

        for sample in self.samples:
            logfile = os.path.join(self.analysis_folder, 'preprocessed', sample + '.html')
            if self.merged_analysis:

                merged_fastqQC(self.samples[sample]['read1'],
                        self.samples[sample]['read2'],
                        logfile)
            else:
                sample_read1_outfile = os.path.join(self.analysis_folder, 'preprocessed', sample + '_R1_processesed.fastq.gz')
                sample_read2_outfile = os.path.join(self.analysis_folder, 'preprocessed',sample + '_R2_processesed.fastq.gz')
                unmerged_fastqQC(self.samples[sample]['read1'],
                        self.samples[sample]['read2'],
                                 sample_read1_outfile,
                                 sample_read2_outfile,
                                 self.adapter_list,
                                 logfile)
                self.samples[sample]['read1'] = sample_read1_outfile
                self.samples[sample]['read2'] = sample_read2_outfile


                control_read1_outfile = os.path.join(self.analysis_folder, 'preprocessed',
                                                    sample + '_R1_processesed_control.fastq.gz')
                control_read2_outfile = os.path.join(self.analysis_folder, 'preprocessed',
                                                    sample + '_R2_processesed_control.fastq.gz')
                unmerged_fastqQC(self.samples[sample]['controlread1'],
                                 self.samples[sample]['controlread2'],
                                 control_read1_outfile,
                                 control_read2_outfile,
                                 self.adapter_list,
                                 logfile)
                self.samples[sample]['controlread1'] = control_read1_outfile
                self.samples[sample]['controlread2'] = control_read2_outfile

        if self.dedup_umi:
            logger.info('Deduplicating UMIs...')

            for sample in self.samples:
                sample_processed_read1_path = os.path.join(self.analysis_folder, 'preprocessed', sample + '_R1_preprocessed.fastq.gz')
                sample_processed_read2_path = os.path.join(self.analysis_folder, 'preprocessed',sample + '_R2_preprocessed.fastq.gz')
                control_processed_read1_path = os.path.join(self.analysis_folder, 'preprocessed', 'control_' + sample +
                                                           '_R1_preprocessed.fastq.gz')
                control_processed_read2_path = os.path.join(self.analysis_folder, 'preprocessed','control_' +
                                                           sample + '_R2_preprocessed.fastq.gz')
                deduplicateReads(self.bc_pattern,
                                 self.samples[sample]['read1'],
                                 self.samples[sample]['read2'],
                                 sample_processed_read1_path,
                                 sample_processed_read2_path)


                deduplicateReads(self.bc_pattern,
                                 self.samples[sample]['controlread1'],
                                 self.samples[sample]['controlread2'],
                                 control_processed_read1_path,
                                 control_processed_read2_path)


    def alignReads(self):
        if self.merged_analysis:
            logger.info('Merging reads...')
            try:
                self.merged = {}
                for sample in self.samples:
                    sample_merge_path = os.path.join(self.analysis_folder, 'fastq', sample + '_merged.fastq.gz')
                    control_sample_merge_path = os.path.join(self.analysis_folder, 'fastq', 'control_' + sample + '_merged.fastq.gz')
                    read1 = self.samples[sample]['read1']
                    read2 = self.samples[sample]['read2']
                    control1 = self.samples[sample]['controlread1']
                    control2 = self.samples[sample]['controlread2']
                    if self.dedup_umi:
                        read1 = os.path.join(self.analysis_folder, 'preprocessed',
                                                                  sample + '_R1_preprocessed.fastq.gz')
                        read2 = os.path.join(self.analysis_folder, 'preprocessed',
                                                                   sample + '_R2_preprocessed.fastq.gz')
                        control1 = os.path.join(self.analysis_folder, 'preprocessed',
                                                                    'control_' + sample +
                                                                    '_R1_preprocessed.fastq.gz')
                        control2 = os.path.join(self.analysis_folder, 'preprocessed', 'control_' +
                                                                    sample + '_R2_preprocessed.fastq.gz')

                    mergeReads(read1,
                               read2,
                               sample_merge_path)
                    mergeReads(control1,
                               control2,
                               control_sample_merge_path)

                    sample_alignment_path = os.path.join(self.analysis_folder, 'aligned', sample + '.sam')
                    control_sample_alignment_path = os.path.join(self.analysis_folder, 'aligned', 'control_' + sample + '.sam')

                    alignReads(self.BWA_path,
                               self.reference_genome,
                               sample_merge_path,
                               '',
                               sample_alignment_path)

                    alignReads(self.BWA_path,
                               self.reference_genome,
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
            logger.info('Aligning reads...')
            try:
                self.aligned = {}
                self.aligned_sorted = {}
                for sample in self.samples:
                    sample_alignment_path = os.path.join(self.analysis_folder, 'aligned', sample + '.sam')
                    control_sample_alignment_path = os.path.join(self.analysis_folder, 'aligned', 'control_' + sample + '.sam')
                    alignReads(self.BWA_path,
                               self.reference_genome,
                               self.samples[sample]['read1'],
                               self.samples[sample]['read2'],
                               sample_alignment_path)
                    alignReads(self.BWA_path,
                               self.reference_genome,
                               self.samples[sample]['controlread1'],
                               self.samples[sample]['controlread2'],
                               control_sample_alignment_path)
                    self.aligned[sample] = sample_alignment_path
                    self.aligned_sorted[sample] = os.path.join(self.analysis_folder, 'aligned', sample + '_sorted.bam')
                    logger.info('Finished aligning reads to genome.')

            except Exception as e:
                logger.error('Error aligning')
                logger.error(traceback.format_exc())
                quit()

    def findCleavageSites(self):
        logger.info('Identifying off-target cleavage sites.')

        try:
            for sample in self.samples:
                if self.merged_analysis:
                    sorted_bam_file = os.path.join(self.analysis_folder, 'aligned', sample + '.bam')
                    control_sorted_bam_file = os.path.join(self.analysis_folder, 'aligned', 'control_' + sample + '.bam')
                else:
                    sorted_bam_file = os.path.join(self.analysis_folder, 'aligned', sample + '_sorted.bam')
                    control_sorted_bam_file = os.path.join(self.analysis_folder, 'aligned', 'control_' + sample + '_sorted.bam')
                identified_sites_file = os.path.join(self.analysis_folder, 'identified', sample)
                logger.info('Window: {0}, MAPQ: {1}, Gap: {2}, Start {3}, Mismatches {4}, Search_Radius {5}'.format(self.window_size, self.mapq_threshold, self.gap_threshold, self.start_threshold, self.mismatch_threshold, self.search_radius))
                findCleavageSites.compare(self.reference_genome, sorted_bam_file, control_sorted_bam_file, self.samples[sample]['target'],
                                          self.search_radius, self.window_size, self.mapq_threshold, self.gap_threshold,
                                          self.start_threshold, self.mismatch_threshold, sample, self.samples[sample]['description'],
                                          identified_sites_file, self.all_chromosomes,merged=self.merged_analysis,read_count_cutoff=self.read_threshold,read_length=self.read_length)

        except Exception as e:
            logger.error('Error identifying off-target cleavage site.')
            logger.error(traceback.format_exc())
            quit()

    def visualize(self):
        logger.info('Visualizing off-target sites')
        try:
            for sample in self.samples:
                 if sample != 'control':
                     if self.annotate_file:
                         infile = os.path.join(self.analysis_folder, 'identified', sample + '_identified_matched.txt')
                         annotate(infile, annotate_path=self.annotate_file)
                         print("Annotating Sample;", sample)

                         outfile = os.path.join(self.analysis_folder, 'visualization', sample + '_offtargets')
                         viz_infile = os.path.join(self.analysis_folder, 'identified', sample + '_identified_matched_annotated.csv')
                         visualizeOfftargets(infile=viz_infile, outfile=outfile, title=sample, PAM=self.PAM, annotation=True)
            logger.info('Finished visualizing off-target sites')

        except Exception as e:
            logger.error('Error visualizing off-target sites.')
            logger.error(traceback.format_exc())

    def OTCoverage(self):
        logger.info('QC Coverage')
        try:
            for sample in self.samples:
                #if self.merged_analysis:
                script_path = p_dir + "/QC_matched_alignment.sh"
                logger.info('Running OT Coverage for {0}'.format(sample))
                coverage_command = 'sh {0} {1} {2}'.format(script_path,
                     sample, self.analysis_folder)
                logger.info(coverage_command)
                subprocess.check_call(coverage_command, shell=True)
                logger.info('OT Coverage for {0} completed.'.format(sample))

        except Exception as e:
            logger.error('Error with OT coverage')


    def callVariants(self):

        try:
            if self.variant_analysis:
                logger.info('Identifying genomic variants')

                for sample in self.samples:
                    sorted_bam_file = os.path.join(self.analysis_folder, 'coverage', sample + '_identified_matched_sorted.bam')
                    identified_sites_file = os.path.join(self.analysis_folder, 'identified', sample + '_identified_matched.txt')
                    variants_basename = os.path.join(self.analysis_folder, 'variants', sample)
                    logger.info('Mismatches {0}, Search_Radius {1}'.format(self.mismatch_threshold, self.search_radius))
                    callVariants.getVariants(identified_sites_file, self.reference_genome, sorted_bam_file, variants_basename, self.search_radius, self.mismatch_threshold)

                    try:
                        annotated_file = variants_basename +  '_Variants.txt'
                        annotate(annotated_file, annotate_path=self.annotate_file)
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
            for sample in self.samples:
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
    all_parser.add_argument('--manifest', '-m', help='Specify the manifest Path', required=True)
    all_parser.add_argument('--sample', '-s', help='Specify sample to process (default is all)', default='all')

    parallel_parser = subparsers.add_parser('parallel', help='Run all steps of the pipeline in parallel')
    parallel_parser.add_argument('--manifest', '-m', help='Specify the manifest Path', required=True)
    parallel_parser.add_argument('--lsf', '-l', help='Specify LSF CMD', default='bsub -R rusage[mem=32000] -P Genomics -q standard')
    parallel_parser.add_argument('--run', '-r', help='Specify which steps of pipepline to run (all, align, identify, visualize, variants)', default='all')

    align_parser = subparsers.add_parser('align', help='Run alignment only')
    align_parser.add_argument('--manifest', '-m', help='Specify the manifest Path', required=True)
    align_parser.add_argument('--sample', '-s', help='Specify sample to process (default is all)', default='all')

    merge_parser = subparsers.add_parser('merge', help='Merge paired end reads')
    merge_parser.add_argument('--manifest', '-m', help='Specify the manifest Path', required=True)
    merge_parser.add_argument('--sample', '-s', help='Specify sample to process (default is all)', default='all')

    identify_parser = subparsers.add_parser('identify', help='Run identification only')
    identify_parser.add_argument('--manifest', '-m', help='Specify the manifest Path', required=True)
    identify_parser.add_argument('--sample', '-s', help='Specify sample to process (default is all)', default='all')

    visualize_parser = subparsers.add_parser('visualize', help='Run visualization only')
    visualize_parser.add_argument('--manifest', '-m', help='Specify the manifest Path', required=True)
    visualize_parser.add_argument('--sample', '-s', help='Specify sample to process (default is all)', default='all')

    variants_parser = subparsers.add_parser('variants', help='Run variants analysis only')
    variants_parser.add_argument('--manifest', '-m', help='Specify the manifest Path', required=True)
    variants_parser.add_argument('--sample', '-s', help='Specify sample to process (default is all)', default='all')

    variants_parser = subparsers.add_parser('coverage', help='Run coverage analysis of matched sites')
    variants_parser.add_argument('--manifest', '-m', help='Specify the manifest Path', required=True)
    variants_parser.add_argument('--sample', '-s', help='Specify sample to process (default is all)', default='all')

    data_parser = subparsers.add_parser('makefiles', help='Combine and normalize replicates, produces vizualations')
    data_parser.add_argument('--outdir', '-o', help='Specify the final annnotation table location. Default is changeseq directory', default=p_dir +"/data/")
    data_parser.add_argument('--ftp_path', '-f', help='RefSeq FTP Path. Must be in refseq.txt format', default ="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ncbiRefSeq.txt.gz")
    data_parser.add_argument('--reset_output', '-r', help='Change path changeseq uses for stored annotation file',default = False)


    reference_free_parser = subparsers.add_parser('reference-free', help='Run reference-free discovery only')
    reference_free_parser.add_argument('--manifest', '-m', help='Specify the manifest Path', required=True)
    reference_free_parser.add_argument('--sample', '-s', help='Specify sample to process (default is all)', default='all')

    return parser.parse_args()

def main():
    args = parse_args()

    if args.command == 'all':
        c = CircleSeq()
        print(args.sample)
        c.parseManifest(args.manifest, args.sample)
        c.processReads()
        c.alignReads()
        c.findCleavageSites()
        c.visualize()
        c.OTCoverage()
        c.callVariants()
    elif args.command == 'parallel':
        c = CircleSeq()
        c.parseManifest(args.manifest)
        c.parallel(args.manifest, args.lsf, args.run)
    elif args.command == 'align':
        c = CircleSeq()
        c.parseManifest(args.manifest, args.sample)
        c.alignReads()
    elif args.command == 'identify':
        c = CircleSeq()
        c.parseManifest(args.manifest, args.sample)
        c.findCleavageSites()
    elif args.command == 'merge':
        c = CircleSeq()
        c.parseManifest(args.manifest, args.sample)
        c.alignReads()
    elif args.command == 'visualize':
        c = CircleSeq()
        c.parseManifest(args.manifest, args.sample)
        c.visualize()
    elif args.command == 'variants':
        c = CircleSeq()
        c.parseManifest(args.manifest, args.sample)
        c.callVariants()
    elif args.command == 'coverage':
        c = CircleSeq()
        c.parseManifest(args.manifest, args.sample)
        c.OTCoverage()
    elif args.command == 'makefiles':
        from make_annotation_table import makefiles
        makefiles(args.ftp_path,args.outdir,args.reset_output,p_dir)

if __name__ == '__main__':
    main()
