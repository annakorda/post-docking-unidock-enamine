#!/bin/bash
# Author: Anna Korda
# Usage: sbatch snn_similarity/snn_similarity.sbatch [cutoff]
#   e.g. sbatch snn_similarity/snn_similarity.sbatch
#   e.g. sbatch snn_similarity/snn_similarity.sbatch -6.0
#SBATCH --job-name=snn_similarity
#SBATCH --partition=normal
#SBATCH --cpus-per-task=32
#SBATCH --output=logs/snn_similarity_%j.out
#SBATCH --error=logs/snn_similarity_%j.err

set -e
cd "$SLURM_SUBMIT_DIR"

CUTOFF="${1:--7.0}"

module load Miniconda3/202601
eval "$(conda shell.bash hook)"
conda activate /users/gpcr/annak/ultra-large/envs/luna-chem

python snn_similarity/snn_similarity.py \
    --vs_results_dir "/users/gpcr/annak/ultra-large/vs_results" \
    --out_dir "/users/gpcr/annak/ultra-large/analysis_results/snn_similarity" \
    --conformations c1,ref1,c5 \
    --cutoff "$CUTOFF" \
    --workers "$SLURM_CPUS_PER_TASK"
