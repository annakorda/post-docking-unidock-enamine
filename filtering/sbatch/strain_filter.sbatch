#!/bin/bash
# Author: Anna Korda
# Usage: sbatch filtering/sbatch/strain_filter.sbatch <conformation> [--dry_run]
#   e.g. sbatch filtering/sbatch/strain_filter.sbatch c5 --dry_run
#   e.g. sbatch filtering/sbatch/strain_filter.sbatch c5
#SBATCH --job-name=strain_filter
#SBATCH --partition=normal
#SBATCH --cpus-per-task=32
#SBATCH --time=04:00:00
#SBATCH --output=filtering/logs/strain_filter_%j.out
#SBATCH --error=filtering/logs/strain_filter_%j.err

set -e
cd "$SLURM_SUBMIT_DIR"

CONFORMATION="$1"
shift || true
EXTRA_ARGS="$@"   # e.g. --dry_run
if [ -z "$CONFORMATION" ]; then
    echo "Usage: sbatch filtering/sbatch/strain_filter.sbatch <conformation> [--dry_run]"
    exit 1
fi

module load Miniconda3/202601
eval "$(conda shell.bash hook)"
conda activate "$SLURM_SUBMIT_DIR/envs/post-dock"

python filtering/strain_filter.py \
    --conformation "$CONFORMATION" \
    --vs_results_dir "vs_results" \
    --torsion_lib_xml "filtering/data/TL_2.1_VERSION_6.xml" \
    --workers "$SLURM_CPUS_PER_TASK" \
    $EXTRA_ARGS
