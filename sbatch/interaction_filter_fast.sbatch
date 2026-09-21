#!/bin/bash
# Author: Anna Korda
# Usage: sbatch sbatch/interaction_filter_fast.sbatch <conformation> <receptor_pdb_path>
#   e.g. sbatch sbatch/interaction_filter_fast.sbatch c5 receptors/c5.pdb
#SBATCH --job-name=interaction_filter_fast
#SBATCH --partition=normal
#SBATCH --cpus-per-task=32
#SBATCH --time=04:00:00
#SBATCH --output=logs/interaction_filter_fast_%j.out
#SBATCH --error=logs/interaction_filter_fast_%j.err

set -e
cd "$SLURM_SUBMIT_DIR"

CONFORMATION="$1"
RECEPTOR_PDB="$2"
if [ -z "$CONFORMATION" ] || [ -z "$RECEPTOR_PDB" ]; then
    echo "Usage: sbatch sbatch/interaction_filter_fast.sbatch <conformation> <receptor_pdb_path>"
    exit 1
fi

module load Miniconda3/202601
eval "$(conda shell.bash hook)"
conda activate "$SLURM_SUBMIT_DIR/envs/post-dock"

python interaction_filter_fast.py \
    --conformation "$CONFORMATION" \
    --poses_dir "vs_results/ad4_redock_${CONFORMATION}" \
    --receptor_pdb "$RECEPTOR_PDB" \
    --out_dir "vs_results" \
    --workers "$SLURM_CPUS_PER_TASK"
