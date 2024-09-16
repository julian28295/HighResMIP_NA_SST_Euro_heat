#!/bin/bash
#SBATCH --account=bm1398
#SBATCH --job-name=HighResMIP_bootstrap
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=08:00:00
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --constraint="512G|1024G"
#SBATCH --mem=0

# Load modules
module load python3
which python
# Execute python script
python HighResMIP_bootstrap_map_composites_pr.py
