#!/bin/bash
# Author: Anna Korda
# Only scores the top --top_pct of each conformation (default 1%) -- real
# measured baseline on shiva is ~65-72s/compound (vs ~25-32s on a fast local
# machine) plus real, gradual concurrency-scaling cost (measured: 1.00x at
# 1 worker, 1.06x at 4, 1.29x at 8, 1.72x at 16), so the full ~106k-230k
# compound sets aren't practical here. Top 1% is a real, small, bounded set
# per conformation (c1 ~1,064, ref1 ~512, c5 ~713) that's actually the most
# interesting compounds anyway (best AD4-ranked already).
#
# Usage: sbatch mmgbsa_rescore/mmgbsa_rescore.sbatch <conformation> [top_pct] [cutoff] [workers]
#   e.g. sbatch mmgbsa_rescore/mmgbsa_rescore.sbatch c1
#   e.g. sbatch mmgbsa_rescore/mmgbsa_rescore.sbatch c1 1.0 -7.0 16
#SBATCH --job-name=mmgbsa_rescore
#SBATCH --partition=long
#SBATCH --cpus-per-task=32
#SBATCH --output=logs/mmgbsa_rescore_%j.out
#SBATCH --error=logs/mmgbsa_rescore_%j.err

set -e
cd "$SLURM_SUBMIT_DIR"

CONFORMATION="$1"
TOP_PCT="${2:-1.0}"
CUTOFF="${3:--7.0}"
# 16 workers is the real, measured sweet spot (1.72x slowdown, not a wall) --
# see the header comment above for the actual data this is based on
WORKERS="${4:-16}"
if [ -z "$CONFORMATION" ]; then
    echo "Usage: sbatch mmgbsa_rescore/mmgbsa_rescore.sbatch <conformation> [top_pct] [cutoff] [workers]"
    exit 1
fi

module load Miniconda3/202601
eval "$(conda shell.bash hook)"
conda activate /users/gpcr/annak/ultra-large/envs/gbsa

export OMPI_MCA_btl=self
export OMPI_MCA_orte_tmpdir_base="/dev/shm/ompi_${SLURM_JOB_ID}"
mkdir -p "$OMPI_MCA_orte_tmpdir_base"

SCRATCH_ROOT="/dev/shm/mmgbsa_${CONFORMATION}_${SLURM_JOB_ID}"
mkdir -p "$SCRATCH_ROOT"

python mmgbsa_rescore/mmgbsa_rescore.py \
    --conformation "$CONFORMATION" \
    --vs_results_dir "/users/gpcr/annak/ultra-large/vs_results" \
    --receptor "receptor_prep/${CONFORMATION}_gmxready.pdb" \
    --out_dir "/users/gpcr/annak/ultra-large/analysis_results/mmgbsa_rescore" \
    --cutoff "$CUTOFF" \
    --top_pct "$TOP_PCT" \
    --workers "$WORKERS" \
    --scratch_root "$SCRATCH_ROOT"

rm -rf "$SCRATCH_ROOT" "$OMPI_MCA_orte_tmpdir_base"
