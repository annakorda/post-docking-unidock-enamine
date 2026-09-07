#!/bin/bash
# Author: Anna Korda
# Run AFTER mmgbsa_rescore.sbatch has finished for this conformation.
# Usage: sbatch mmgbsa_rescore/mmgbsa_save_top_complexes.sbatch <conformation> [top_pct] [cutoff]
#   e.g. sbatch mmgbsa_rescore/mmgbsa_save_top_complexes.sbatch c1
#   e.g. sbatch mmgbsa_rescore/mmgbsa_save_top_complexes.sbatch c1 1.0 -7.0
#SBATCH --job-name=mmgbsa_save_top_complexes
#SBATCH --partition=long
#SBATCH --cpus-per-task=64
#SBATCH --output=logs/mmgbsa_save_top_complexes_%j.out
#SBATCH --error=logs/mmgbsa_save_top_complexes_%j.err

set -e
cd "$SLURM_SUBMIT_DIR"

CONFORMATION="$1"
TOP_PCT="${2:-1.0}"
CUTOFF="${3:--7.0}"
if [ -z "$CONFORMATION" ]; then
    echo "Usage: sbatch mmgbsa_rescore/mmgbsa_save_top_complexes.sbatch <conformation> [top_pct] [cutoff]"
    exit 1
fi

module load Miniconda3/202601
eval "$(conda shell.bash hook)"
conda activate /users/gpcr/annak/ultra-large/envs/gbsa

SCRATCH_ROOT="${TMPDIR:-/tmp}/mmgbsa_top_${CONFORMATION}_${SLURM_JOB_ID}"
mkdir -p "$SCRATCH_ROOT"

python mmgbsa_rescore/mmgbsa_save_top_complexes.py \
    --conformation "$CONFORMATION" \
    --vs_results_dir "/users/gpcr/annak/ultra-large/vs_results" \
    --receptor "receptor_prep/${CONFORMATION}_gmxready.pdb" \
    --mmgbsa_out_dir "/users/gpcr/annak/ultra-large/analysis_results/mmgbsa_rescore" \
    --cutoff "$CUTOFF" \
    --top_pct "$TOP_PCT" \
    --workers "$SLURM_CPUS_PER_TASK" \
    --scratch_root "$SCRATCH_ROOT"

rm -rf "$SCRATCH_ROOT"
