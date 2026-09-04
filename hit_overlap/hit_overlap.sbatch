#!/bin/bash
# Author: Anna Korda
# Usage: sbatch hit_overlap/hit_overlap.sbatch [cutoff]
#   e.g. sbatch hit_overlap/hit_overlap.sbatch
#   e.g. sbatch hit_overlap/hit_overlap.sbatch -6.0
#SBATCH --job-name=hit_overlap
#SBATCH --partition=normal
#SBATCH --cpus-per-task=2
#SBATCH --output=logs/hit_overlap_%j.out
#SBATCH --error=logs/hit_overlap_%j.err

set -e
cd "$SLURM_SUBMIT_DIR"

CUTOFF="${1:--7.0}"

module load Miniconda3/202601
eval "$(conda shell.bash hook)"
conda activate /users/gpcr/annak/ultra-large/envs/luna-chem

python hit_overlap/hit_overlap.py \
    --vs_results_dir "/users/gpcr/annak/ultra-large/vs_results" \
    --out_dir "/users/gpcr/annak/ultra-large/analysis_results/hit_overlap" \
    --conformations c1,ref1,c5 \
    --cutoff "$CUTOFF"
