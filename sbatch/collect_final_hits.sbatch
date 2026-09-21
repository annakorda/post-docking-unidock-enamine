#!/bin/bash
# Author: Anna Korda
# Usage: sbatch sbatch/collect_final_hits.sbatch <conformation> [chunk_size]
#   e.g. sbatch sbatch/collect_final_hits.sbatch c1
#   e.g. sbatch sbatch/collect_final_hits.sbatch c1 10000
#SBATCH --job-name=collect_final_hits
#SBATCH --partition=normal
#SBATCH --cpus-per-task=32
#SBATCH --time=02:00:00
#SBATCH --output=logs/collect_final_hits_%j.out
#SBATCH --error=logs/collect_final_hits_%j.err

set -e
cd "$SLURM_SUBMIT_DIR"

CONFORMATION="$1"
CHUNK_SIZE="${2:-50000}"
if [ -z "$CONFORMATION" ]; then
    echo "Usage: sbatch sbatch/collect_final_hits.sbatch <conformation> [chunk_size]"
    exit 1
fi

module load Miniconda3/202601
eval "$(conda shell.bash hook)"
conda activate "$SLURM_SUBMIT_DIR/envs/post-dock"

python collect_final_hits.py \
    --conformation "$CONFORMATION" \
    --vs_results_dir "vs_results" \
    --chunk_size "$CHUNK_SIZE" \
    --workers "$SLURM_CPUS_PER_TASK"
