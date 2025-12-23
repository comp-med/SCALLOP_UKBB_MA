#!/bin/sh

## script to run pQTL analysis across FinnGen stats

#! select partition
#SBATCH --partition=compute

#! select account
#SBATCH --account=sc-users

#! Specify required run time
#SBATCH --time=24:00:00

#! how many nodes
#SBATCH --nodes=1

#! how many tasks
#SBATCH --ntasks=1

#! how much memory per cpu
#SBATCH --mem-per-cpu=3G

#! how many cpus per task
#SBATCH --cpus-per-task=10

#! no other jobs can be run on the same node
##SBATCH --exclusive

#! run as job array
#SBATCH --array=1-1114%15

#! What types of email messages do you wish to receive?
#SBATCH --mail-type=FAIL

#! set name
#SBATCH --output=slurm-%x-%j.out

## change directory
cd <path to dir>

## get the exposure
echo "Job ID: $SLURM_ARRAY_TASK_ID"
expo="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F '\t' 'NR == var {print $1}' input/input.FinnGen.cis.trans.MR.pipeline.txt)"

echo "Node ID: $SLURM_NODELIST"

echo "Exposure: ${expo}"

## set up R environment
eval "$(/opt/conda/bin/conda shell.bash hook)"
source activate Renv

## run the R script
scripts/03_MR_pipeline.R ${expo}

## deactivate conda env
conda deactivate