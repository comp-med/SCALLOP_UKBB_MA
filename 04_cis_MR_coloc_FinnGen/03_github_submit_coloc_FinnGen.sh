#!/bin/sh

## script to run cis-pQTL coloc analysis across FinnGen stats

#! select partition
#SBATCH --partition=compute

#! select account
#SBATCH --account=sc-users

#! Specify required run time
#SBATCH --time=12:00:00

#! how many nodes
#SBATCH --nodes=1

#! how many tasks
#SBATCH --ntasks=1

#! how much memory per cpu
#SBATCH --mem-per-cpu=5G

#! how many cpus per task
#SBATCH --cpus-per-task=10

#! no other jobs can be run on the same node
##SBATCH --exclusive

#! run as job array
#SBATCH --array=1%20

#! What types of email messages do you wish to receive?
#SBATCH --mail-type=FAIL

#! set name
#SBATCH --output=slurm-%x-%j.out

## change directory
cd <path to file>

## file with results
export FL='input/Olink.coloc.targets.rerun.txt'

## get the exposure
echo "Job ID: $SLURM_ARRAY_TASK_ID"
olink="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F '\t' 'NR == var {print $1}' ${FL})"
chr="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F '\t' 'NR == var {print $2}' ${FL})"
poss="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F '\t' 'NR == var {print $3}' ${FL})"
pose="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F '\t' 'NR == var {print $4}' ${FL})"

echo "Node ID: $SLURM_NODELIST"

echo "Phenotype ${olink} : Chromosome ${chr} : Locus start ${poss} : Locus end ${pose}"

## set up R environment
eval "$(/opt/conda/bin/conda shell.bash hook)"
source activate Renv

## run the R script
scripts/05_github_run_coloc.R "${olink}" "${chr}" "${poss}" "${pose}"

## deactivate conda env
conda deactivate