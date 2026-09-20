#!/bin/bash
#SBATCH --job-name=get_data
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4

module load sratoolkit || conda activate smk

cd data/hifi

# Download the REAL full dataset (replace SRR_NUMBER with your chosen accession)
fasterq-dump --threads 4 SRR28741072
