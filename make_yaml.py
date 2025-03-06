import pandas as pd
import os
import yaml
import log
import argparse

logger = log.createCustomLogger('root')


def make_yaml_manifest(analysis_folder,fq_dir,manifest,settings,yaml_fname):
    with open(yaml_fname,'w') as yf:
        for k,v in settings.items():
            yf.write(k+': ' + str(v).strip() + '\n')
        yf.write('samples:\n')

        num_samples = len(manifest['sequencing_sample_name'])
        fq_files = [x for x in os.listdir(fq_dir) if "fastq" in x]

        if num_samples == 0:
            logger.error("ERROR; NO SAMPLES DETECTED IN MANIFEST")
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

def parse_args():
    mainParser = argparse.ArgumentParser()
    mainParser.add_argument('--sample_manifest', '-m', help='samples manifest(csv)', required=True)
    mainParser.add_argument('--parameters', '-p', help='output directory', required=True)


    return mainParser.parse_args()


def main():
    args = parse_args()
    yaml_fname = args.sample_manifest.replace(".csv", ".yaml")
    print("yaml file out: ", yaml_fname)

    settings = pd.read_csv(args.parameters, names=['parameter', 'setting']).set_index('parameter').to_dict()['setting']
    manifest = pd.read_csv(args.sample_manifest)



if __name__ == '__main__':
    main()
