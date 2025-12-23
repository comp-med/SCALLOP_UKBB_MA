#!/bin/sh

## extract proxies

#SBATCH --partition=compute
#SBATCH --job-name=olink_prevalent_disease
#SBATCH --account=sc-users
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=2G
#SBATCH --array=1-845%50
#SBATCH --output=%x-%A-%2a.out

## change directory
cd <path to file>
  
## Use Array Index to select features
echo "Job ID: $SLURM_ARRAY_TASK_ID"
outc="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F '\t' 'NR == var+1 {print $3}' input/label.phecodes.20240624.txt)"
sex="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F '\t' 'NR == var+1 {print $9}' input/label.phecodes.20240624.txt)"

echo "Node ID: $SLURM_NODELIST"

## Do some logging
echo "Test outcome: ${outc} in ${sex}"
printf -v date '%(%Y-%m-%d %H:%M:%S)T\n' -1
echo "Starting script at: $date"

## This is the container to be used
R_CONTAINER='<path to file>/all_inclusive_rstudio_4.3.2.sif'

# This is the script that is executed
# Get with rstudioapi::getSourceEditorContext()$path
R_SCRIPT='<path to file>/04_github_prevalent_association_testing.R'

# Enter all directories you need, simply in a comma-separated list
BIND_DIR="<path to file>"

## The container 
singularity exec \
--bind $BIND_DIR \
$R_CONTAINER Rscript $R_SCRIPT ${outc} ${sex}

## Some more logging
printf -v date '%(%Y-%m-%d %H:%M:%S)T\n' -1
echo "Finishing script at: $date"
echo "Done!"