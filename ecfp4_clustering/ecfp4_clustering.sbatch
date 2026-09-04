#!/bin/bash
# Author: Anna Korda
# Usage: sbatch ecfp4_clustering/ecfp4_clustering.sbatch <conformation> [cutoff]
#   e.g. sbatch ecfp4_clustering/ecfp4_clustering.sbatch c1
#   e.g. sbatch ecfp4_clustering/ecfp4_clustering.sbatch c1 -6.0
#SBATCH --job-name=ecfp4_clustering
#SBATCH --partition=normal
#SBATCH --cpus-per-task=32
#SBATCH --output=logs/ecfp4_clustering_%j.out
#SBATCH --error=logs/ecfp4_clustering_%j.err

set -e
cd "$SLURM_SUBMIT_DIR"

CONFORMATION="$1"
CUTOFF="${2:--7.0}"
if [ -z "$CONFORMATION" ]; then
    echo "Usage: sbatch ecfp4_clustering/ecfp4_clustering.sbatch <conformation> [cutoff]"
    exit 1
fi

module load Miniconda3/202601
eval "$(conda shell.bash hook)"
conda activate /users/gpcr/annak/ultra-large/envs/luna-chem

python ecfp4_clustering/ecfp4_clustering.py \
    --conformation "$CONFORMATION" \
    --vs_results_dir "/users/gpcr/annak/ultra-large/vs_results" \
    --out_dir "/users/gpcr/annak/ultra-large/analysis_results/ecfp4_clustering" \
    --cutoff "$CUTOFF" \
    --thresholds 0.90,0.85,0.80,0.75,0.70 \
    --workers "$SLURM_CPUS_PER_TASK"
