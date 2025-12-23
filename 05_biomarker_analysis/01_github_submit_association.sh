#!/bin/sh

## script to run association testing

#! select partition
#SBATCH --partition=compute

#! select account
#SBATCH --account=sc-users

#! Specify required run time
#SBATCH --time=72:00:00

#! how many nodes
#SBATCH --nodes=1

#! how many tasks
#SBATCH --ntasks=1

#! how many cpus per task
#SBATCH --cpus-per-task=13

#! define how much memory for each node
#SBATCH --mem-per-cpu=2G

#! no other jobs can be run on the same node
##SBATCH --exclusive

#! run as job array
#SBATCH --array=1-1447%20

#! What types of email messages do you wish to receive?
#SBATCH --mail-type=FAIL

#! set name
#SBATCH --output=slurm-%x-%j.out

## change directory
cd <path to dir>

## get the exposure (add one to avoid the header)
echo "Job ID: $SLURM_ARRAY_TASK_ID"
outc="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F '\t' 'NR == var+1 {print $1}' input/label.phecodes.update.20230524.txt)"
sex="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F '\t' 'NR == var+1 {print $4}' input/label.phecodes.update.20230524.txt)"

echo "Test outcome: ${outc} in ${sex}"

## set up R environment
source activate Renv

## run the R script
scripts/02_github_association_testing.R ${outc} ${sex}

## deactivate conda env
conda deactivate
