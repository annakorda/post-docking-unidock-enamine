#!/bin/bash
# Author: Anna Korda
# Run AFTER mmgbsa_rescore.sbatch has finished for this conformation.
# Usage: sbatch mmgbsa_rescore/mmgbsa_save_top_complexes.sbatch <conformation> [top_pct] [cutoff] [workers]
#   e.g. sbatch mmgbsa_rescore/mmgbsa_save_top_complexes.sbatch c1
#   e.g. sbatch mmgbsa_rescore/mmgbsa_save_top_complexes.sbatch c1 1.0 -7.0 32
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
# see mmgbsa_rescore.sbatch -- not tied 1:1 to --cpus-per-task=64, real MPI
# singleton-init contention at high concurrency
WORKERS="${4:-16}"
if [ -z "$CONFORMATION" ]; then
    echo "Usage: sbatch mmgbsa_rescore/mmgbsa_save_top_complexes.sbatch <conformation> [top_pct] [cutoff] [workers]"
    exit 1
fi

module load Miniconda3/202601
eval "$(conda shell.bash hook)"
conda activate /users/gpcr/annak/ultra-large/envs/gbsa

export OMPI_MCA_btl=self
export OMPI_MCA_orte_tmpdir_base="/dev/shm/ompi_${SLURM_JOB_ID}"
mkdir -p "$OMPI_MCA_orte_tmpdir_base"

# /dev/shm, not /tmp -- see mmgbsa_rescore.sbatch for why (real run showed
# severe stalling; /dev/shm is RAM-backed and genuinely node-local)
SCRATCH_ROOT="/dev/shm/mmgbsa_top_${CONFORMATION}_${SLURM_JOB_ID}"
mkdir -p "$SCRATCH_ROOT"

python mmgbsa_rescore/mmgbsa_save_top_complexes.py \
    --conformation "$CONFORMATION" \
    --vs_results_dir "/users/gpcr/annak/ultra-large/vs_results" \
    --receptor "receptor_prep/${CONFORMATION}_gmxready.pdb" \
    --mmgbsa_out_dir "/users/gpcr/annak/ultra-large/analysis_results/mmgbsa_rescore" \
    --cutoff "$CUTOFF" \
    --top_pct "$TOP_PCT" \
    --workers "$WORKERS" \
    --scratch_root "$SCRATCH_ROOT"

rm -rf "$SCRATCH_ROOT" "$OMPI_MCA_orte_tmpdir_base"
