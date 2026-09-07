#!/bin/bash
# Author: Anna Korda
# Usage: sbatch mmgbsa_rescore/mmgbsa_rescore.sbatch <conformation> [cutoff]
#   e.g. sbatch mmgbsa_rescore/mmgbsa_rescore.sbatch c1
#SBATCH --job-name=mmgbsa_rescore
#SBATCH --partition=long
#SBATCH --cpus-per-task=64
#SBATCH --output=logs/mmgbsa_rescore_%j.out
#SBATCH --error=logs/mmgbsa_rescore_%j.err

set -e
cd "$SLURM_SUBMIT_DIR"

CONFORMATION="$1"
CUTOFF="${2:--7.0}"
if [ -z "$CONFORMATION" ]; then
    echo "Usage: sbatch mmgbsa_rescore/mmgbsa_rescore.sbatch <conformation> [cutoff]"
    exit 1
fi

module load Miniconda3/202601
eval "$(conda shell.bash hook)"
conda activate /users/gpcr/annak/ultra-large/envs/gbsa

# node-local scratch, not shared GPFS -- this script churns ~7 small files per
# compound across ~50-100k compounds, real inode/IO pressure on a shared fs
SCRATCH_ROOT="${TMPDIR:-/tmp}/mmgbsa_${CONFORMATION}_${SLURM_JOB_ID}"
mkdir -p "$SCRATCH_ROOT"

python mmgbsa_rescore/mmgbsa_rescore.py \
    --conformation "$CONFORMATION" \
    --vs_results_dir "/users/gpcr/annak/ultra-large/vs_results" \
    --receptor "receptor_prep/${CONFORMATION}_gmxready.pdb" \
    --out_dir "/users/gpcr/annak/ultra-large/analysis_results/mmgbsa_rescore" \
    --cutoff "$CUTOFF" \
    --workers "$SLURM_CPUS_PER_TASK" \
    --scratch_root "$SCRATCH_ROOT"

rm -rf "$SCRATCH_ROOT"
