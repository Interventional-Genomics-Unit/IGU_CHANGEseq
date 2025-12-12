import string
import re
import gzip
import validation
import yaml
import os
import logging
import pandas as pd
import subprocess
import sys
from collections import OrderedDict
"""
FASTQ generator function from umi package
"""
logger = logging.getLogger('root')
logger.propagate = False

def write_default_yaml(param_name,param):
    default_yaml = os.path.dirname(os.path.realpath(validation.__file__)) + "/default.yaml"
    with open(default_yaml, 'r') as f:
        default = yaml.load(f, Loader=yaml.FullLoader)
    default[param_name] = param
    print(f"Storing {param_name} as {param}")
    with open(default_yaml, "w") as f:
        yaml.dump(default,f,default_flow_style=False)

def make_folders(analysis_folder):
    output_dir ={}
    for folder in ['preprocessed', 'aligned', 'raw_results', 'fastq', 'variants', 'qc',
                   'post-process_results', 'post-process_results/visualization', 'post-process_results/tables',
                   'raw_results/visualizations', 'raw_results/tables']:
        output_dir[folder] = os.path.join(analysis_folder, folder)
        if not os.path.exists(output_dir[folder]):
            os.makedirs(output_dir[folder])
    return output_dir

def fq(file):
    if re.search('.gz$', file):
        fastq = gzip.open(file, 'rb')
    else:
        fastq = open(file, 'r')
    with fastq as f:
        while True:
            l1 = f.readline().rstrip('\n')
            if not l1:
                break
            l2 = f.readline().rstrip('\n')
            l3 = f.readline().rstrip('\n')
            l4 = f.readline().rstrip('\n')
            yield [l1, l2, l3, l4]

def reverseComplement(sequence):
    transtab = string.maketrans("ACGT","TGCA")
    return sequence.translate(transtab)[::-1]


def get_parameters(analysis_folder,fq_dir,sample_manifest,settings='default'):
    sample_manifest = f'{analysis_folder}/{os.path.basename(sample_manifest)}'
    yaml_fname = sample_manifest.replace(".csv", ".yaml")
    default_yaml = os.path.dirname(os.path.realpath(validation.__file__)) + "/default.yaml"

    with open(default_yaml, 'r') as f:
        default = yaml.load(f, Loader=yaml.FullLoader)

    validation.exists(filepath=sample_manifest)


    manifest_df = pd.read_csv(sample_manifest)
    if settings != 'default':
        settings = f'{analysis_folder}/{os.path.basename(settings)}'
        logger.info(f"custom settings loading from file {settings}")
        settings_dict = pd.read_csv(settings, names=['parameter', 'setting']).set_index('parameter').to_dict()['setting']
    else:
        logger.info(f"using default settings")
        settings_dict = {}

    return_dict = {}
    for param,default_s in default.items():

        if param in settings_dict.keys():
            try:
                return_dict[param] = int(settings_dict[param]) if type(settings_dict[param]) != bool else settings_dict[param]
            except Exception:
                return_dict[param] = settings_dict[param]
        else:
            try:

                return_dict[param] = int(default[param]) if type(default[param]) != bool else default[param]
            except Exception:
                return_dict[param] = default[param]
        if type(return_dict[param]) != dict:
            logger.info(f'{param} is set to {str(return_dict[param])}')



    for x in ['merged_analysis','variant_analysis']:
        return_dict[x] =False if str(return_dict[x]) == 'False' or str(return_dict[x]) == 'FALSE' else True
        logger.info(f'{x} is set to {str(return_dict[x])}')


    return_dict['analysis_folder'] = analysis_folder
    return_dict['raw_fastq_folder'] = fq_dir

    validation.validateManifest(return_dict)

    return_dict['samples'] = {}
    return_dict['replicates'] = {}

    num_samples = len(manifest_df['sequencing_sample_name'])
    if num_samples == 0:
        logger.error("ERROR; NO SAMPLES DETECTED IN MANIFEST")
        exit()

    fq_files = [x for x in os.listdir(fq_dir) if "fastq" in x]

    columns = manifest_df.columns
    if 'replicate_group_name' in columns:
        for name in set(manifest_df['replicate_group_name']):
            replicate_samples = manifest_df.loc[manifest_df['replicate_group_name'] == name, ['sample_name', 'target']].to_dict('list')

            return_dict['replicates'][name] = replicate_samples

    for i in range(num_samples):
        #sample_basename = [x[:x.find('_001.f') - 2] for x in fq_files if manifest_df.iloc[i]['sequencing_sample_name'] in x][0]
        sample_basename = [x for x in fq_files if manifest_df.iloc[i]['sequencing_sample_name'] in x if '_R1_0' in x][0]
        control_basename = [x for x in fq_files if manifest_df.iloc[i]['control_sequencing_sample_name'] in x if '_R1_0' in x][0]
        if len(sample_basename)==0:
            logger.error('fastq file for {0} is not detected'.format(manifest_df.iloc[i]['sequencing_sample_name']))
            sys.exit()
        if len(control_basename)==0:
            logger.error('fastq file for {0} is not detected'.format(manifest_df.iloc[i]['control_sequencing_sample_name']))
            sys.exit()

        return_dict['samples'][ manifest_df.iloc[i]['sample_name']] = {}
        return_dict['samples'][manifest_df.iloc[i]['sample_name']]['target'] = manifest_df.iloc[i]['target']
        return_dict['samples'][manifest_df.iloc[i]['sample_name']]['read1'] = fq_dir + sample_basename
        return_dict['samples'][manifest_df.iloc[i]['sample_name']]['read2'] = fq_dir + sample_basename.replace("_R1_0","_R2_0")
        return_dict['samples'][manifest_df.iloc[i]['sample_name']]['controlread1'] = fq_dir + control_basename
        return_dict['samples'][manifest_df.iloc[i]['sample_name']]['controlread2'] = fq_dir + control_basename.replace("_R1_0","_R2_0")
        return_dict['samples'][manifest_df.iloc[i]['sample_name']]['description'] = str(manifest_df.iloc[i]['description'])

    with open(yaml_fname, 'w') as yf:
        yaml.dump( return_dict, yf, default_flow_style=False)
    logger.info("yaml manifest file created")
    #logger.info(return_dict)


    return return_dict

def copy_control_samples(samples):
    # changeseq is built to run every control with every sample except, usually theres only one control for all samples
    # ro by pass this just copy one control for all
    representative_controls = {}
    controlreads = {}
    for sample, info in samples.items():
        if samples[sample]['controlread1'] not in controlreads.keys():
            controlreads[samples[sample]['controlread1']] = []
        controlreads[samples[sample]['controlread1']].append(sample)
    for controlread, samples in controlreads.items():
        representative = samples[0]
        for s in samples:
            representative_controls[s] = representative
    return representative_controls


def check_control_exists(sample,representative_control,control_outfile):
    run_control_flag = True
    representative_control_file = control_outfile.replace(sample,representative_control)
    if sample != representative_control:
        if validation.exists(representative_control_file):
            run_control_flag = False
            make_control_copy(representative_control_file, control_outfile)
    else:
        logger.info(f'representative control does not exsists, running control for {sample}')
    return run_control_flag

def make_control_copy(representative_control_file,control_outfile):
    # copies control files that the sample
    logger.info(f'skipping re-run of {control_outfile}')
    cmd = f'cp {representative_control_file} {control_outfile}'
    logger.info(cmd)
    subprocess.check_call(cmd, shell=True)
    logger.info('done.')




