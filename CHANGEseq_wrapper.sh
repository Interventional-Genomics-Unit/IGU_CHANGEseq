#!/bin/bash
#SBATCH --job-name IGU_CHANGEseq
#SBATCH -n 32
#SBATCH -o %j.out
#SBATCH -e %j.err


## ACTIVATE CONDA
eval "$(conda shell.bash hook)"
conda activate changeseqbe
SCRIPT_PATH=/home/thudson/projects/IGU_CHANGEseq/CHANGEseq_wrapper.sh
# Get the directory where the script is located
SCRIPT_DIR="$(dirname "$(realpath "$SCRIPT_PATH")")"

# Now you can reference this directory in your script
echo "parameters manifest(csv): $1"
echo "samples manifest(csv): $2"
echo "command: $3"
echo "sample: $4"
echo "making yaml file"

python $SCRIPT_DIR/make_yaml.py $1 $2



infile=$2
YAML=$(echo ${infile} | sed -e 's/\.csv$/.yaml/')
echo $YAML

python $SCRIPT_DIR/changeseq/changeseq.py $3 --manifest $YAML --sample $4


