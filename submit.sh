#!/bin/bash
#SBATCH --job-name=phyto_test
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --output=slurm/pipeline_%j.log
#SBATCH --error=slurm/pipeline_%j.err

# Load the environment
module load conda/latest
conda activate smk

# Execute the pipeline
snakemake --use-conda --cores 8 run_assembly

