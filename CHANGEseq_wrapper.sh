#!/bin/bash
#SBATCH -p standard
#SBATCH --job-name IGU_CHANGEseq
#SBATCH -o %j.out
#SBATCH -e %j.err

## ACTIVATE CONDA
eval "$(conda shell.bash hook)"
conda activate cseq3

pwd=$(pwd)
## PERMISSIONS
chmod -R 775 $pwd

echo "parameters manifest(csv): $1"
echo "samples manifest(csv): $2"
echo "command: $3"
echo "sample: $4"
echo "making yaml file"
python /home/thudson/projects/IGU_CHANGEseq/make_yaml.py $1 $2



infile=$2
YAML=$(echo ${infile} | sed -e 's/\.csv$/.yaml/')
echo $YAML

python /home/thudson/projects/IGU_CHANGEseq/changeseq/changeseq.py $3 --manifest $YAML --sample $4


