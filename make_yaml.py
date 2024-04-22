import pandas as pd
import os
import sys

parameters_csv = sys.argv[1]
samples_csv = sys.argv[2]


parameters_csv =  "/groups/clinical/projects/Assay_Dev/IGU_CHANGEseq/CS_05/12878/cs_medit_parameters.csv"
samples_csv = "/groups/clinical/projects/Assay_Dev/IGU_CHANGEseq/CS_05/12878/12878_medit.csv"
yaml_fname = samples_csv.replace(".csv",".yaml")
print("yaml file out: ",yaml_fname)

#manifest = pd.read_csv(manifest_csv,names = ['parameter','setting'])
#i = manifest[manifest['parameter'] == 'name'].index[0]
#settings = manifest.iloc[0:i].set_index('parameter').to_dict()['setting']
#samples = manifest.iloc[i:]
settings = pd.read_csv(parameters_csv,names = ['parameter','setting']).set_index('parameter').to_dict()['setting']
manifest = pd.read_csv(samples_csv)


with open(yaml_fname,'w') as yf:
    yf.write('reference_genome: ' +settings['reference_genome'] + '\n')
    yf.write('analysis_folder: ' +settings['analysis_folder'] +'\n')
    yf.write('annotate_file: ' + settings['annotate_file'] + '\n')
    yf.write('raw_fastq_folder: ' + settings['raw_fastq_folder'] + '\n')
    yf.write('\n')

    yf.write('bwa: ' + settings['bwa'] + '\n')
    yf.write('samtools: ' + settings['samtools'] + '\n\n')

    yf.write('read_threshold: %s' %settings['read_threshold'] +'\n')
    yf.write('window_size: %s' %settings['window_size'] + '\n')
    yf.write('mapq_threshold: %s' %settings['mapq_threshold'] + '\n')
    yf.write('start_threshold: %s' %settings['start_threshold'] + '\n')
    yf.write('gap_threshold: %s' %settings['gap_threshold'] + '\n')
    yf.write('mismatch_threshold: %s' %settings['mismatch_threshold'] +'\n')
    yf.write('search_radius : ' + str(30) + '\n')
    yf.write('merged_analysis: ' + settings['merged_analysis'] + '\n')
    yf.write('variant_analysis: ' + settings['variant_analysis'] + '\n')
    yf.write('dedup_umi: ' + settings['dedup_umi'] + '\n')
    yf.write('bc_pattern: ' + settings['bc_pattern'] + '\n')

    yf.write('samples:\n')
    #for samplename in manifest.loc[manifest['parameter'] == 'name', 'setting']:
    num_samples = len(manifest['sequencing_sample_name'])
    fq_dir = settings['raw_fastq_folder']
    fq_files = [x for x in os.listdir(fq_dir) if "fastq" in x]
    for i in range(num_samples):
        sample_basename = [x[:x.find('_001.f')-2] for x in fq_files if manifest.iloc[i]['sequencing_sample_name'] in x][0]
        control_basename = [x[:x.find('_001.f')-2] for x in fq_files if manifest.iloc[i]['control_sequencing_sample_name'] in x][0]

        yf.write('  ' + manifest.iloc[i]['sample_name'] + ':\n')
        yf.write('    target: ' + manifest.iloc[i]['target'] + '\n')
        yf.write('    read1: ' +fq_dir + sample_basename + 'R1_001.fastq.gz'+'\n')
        yf.write('    read2: ' + fq_dir + sample_basename + 'R2_001.fastq.gz'+'\n')
        yf.write('    controlread1: ' + fq_dir + control_basename + 'R1_001.fastq.gz'+ '\n')
        yf.write('    controlread2: ' + fq_dir + control_basename + 'R2_001.fastq.gz'+ '\n')
        yf.write('    description: '+ str(manifest.iloc[i]['description']) + '\n')
yf.close()


