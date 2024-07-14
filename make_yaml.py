import pandas as pd
import os
import sys

parameters_csv = sys.argv[1]
samples_csv = sys.argv[2]


yaml_fname = samples_csv.replace(".csv",".yaml")
print("yaml file out: ",yaml_fname)

settings = pd.read_csv(parameters_csv,names = ['parameter','setting']).set_index('parameter').to_dict()['setting']
manifest = pd.read_csv(samples_csv)


with open(yaml_fname,'w') as yf:
    for k,v in settings.items():
        yf.write(k+': ' + str(v).strip() + '\n')
    yf.write('samples:\n')
    #for samplename in manifest.loc[manifest['parameter'] == 'name', 'setting']:
    num_samples = len(manifest['sequencing_sample_name'])
    fq_dir = settings['raw_fastq_folder']
    fq_files = [x for x in os.listdir(fq_dir) if "fastq" in x]

    if num_samples == 0:
        print("ERROR; NO SAMPLES DETECTED IN MANIFEST")
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


