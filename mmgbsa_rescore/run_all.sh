#!/bin/bash
# Author: Anna Korda
# Runs MM-GBSA on all 3 conformations sequentially (not in parallel --
# each one already uses --workers processes internally; running the 3
# conformations concurrently on top of that would oversubscribe the
# machine). Resumable: safe to re-run after a crash/interrupt, already-
# done chunks are skipped (see run_mmgbsa_q5.py's own resume logic).
#
# Usage: bash run_all.sh [workers] [scratch_root]
#   workers:      passed through as --workers (default: all cores)
#   scratch_root: passed through as --scratch_root (default: system tmp --
#                 point this at local fast storage if available, e.g. /dev/shm)
# unigbsa does a per-process OpenMPI singleton-init internally (via mpi4py)
# -- flagged in the original sbatch-based run as a real, possibly-
# contributing factor to concurrency degradation, same class of issue as
# disk contention. Same two mitigations here: OMPI_MCA_btl=self (skip
# network transport probing for a single-process singleton), and moving
# OpenMPI's own session dir to /dev/shm (fast, and out of the way of
# --scratch_root's own per-compound dirs). Harmless either way, not
# confirmed as the primary fix on its own, but cheap and real.
export OMPI_MCA_btl=self
export OMPI_MCA_orte_tmpdir_base="/dev/shm/ompi_$$"
mkdir -p "$OMPI_MCA_orte_tmpdir_base"
trap 'rm -rf "$OMPI_MCA_orte_tmpdir_base"' EXIT

set -e
cd "$(dirname "$0")"

WORKERS="${1:-}"
SCRATCH="${2:-}"
WORKERS_ARG=()
SCRATCH_ARG=()
[ -n "$WORKERS" ] && WORKERS_ARG=(--workers "$WORKERS")
[ -n "$SCRATCH" ] && SCRATCH_ARG=(--scratch_root "$SCRATCH")

mkdir -p logs results

for conf in c1 ref1 c5; do
    echo "=== $conf: starting $(date) ==="
    python run_mmgbsa_q5.py \
        --conformation "$conf" \
        --sdf "../nested_bm_clustering/q5/q5_mmgbsa_${conf}.sdf" \
        --receptor "../receptor_prep/${conf}_gmxready.pdb" \
        --out_dir results \
        "${WORKERS_ARG[@]}" "${SCRATCH_ARG[@]}" \
        2>&1 | tee -a "logs/mmgbsa_${conf}.log"
    echo "=== $conf: done $(date) ==="
done

echo "All 3 conformations done. Results under results/mmgbsa_<conf>/mmgbsa_results_<conf>.csv"
