import string
import re
import gzip
import validation
import yaml
import os
import log
"""
FASTQ generator function from umi package
"""

def write_default_yaml(param_name,param):
    default_yaml = os.path.dirname(os.path.realpath(validation.__file__)) + "/default.yaml"
    with open(default_yaml, 'r') as f:
        default = yaml.load(f, Loader=yaml.FullLoader)
    default[param_name] = param
    print(f"Storing {param_name} as {param}")
    with open(default_yaml, "w") as f:
        yaml.dump(default,f,default_flow_style=False)

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


def get_parameters(manifest_data):
    default_yaml = os.path.dirname(os.path.realpath(validation.__file__)) + "/default.yaml"
    with open(default_yaml, 'r') as f:
        default = yaml.load(f, Loader=yaml.FullLoader)
    with open(manifest_data, 'r') as f:
        return_dict = yaml.load(f, Loader=yaml.FullLoader)
    default['analysis_folder'] = os.getcwd()
    validation.validateManifest(return_dict)
    for p in default:
        if not p in return_dict:
            if p == "samples":
                logger.error("No samples are found in the yaml file, please provide samples (fastq) to start with!")
                exit()
            return_dict[p] = default[p]
    return return_dict



