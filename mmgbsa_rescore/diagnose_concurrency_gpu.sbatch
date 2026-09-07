#!/bin/bash
# Author: Anna Korda
# Same diagnostic as diagnose_concurrency.sbatch, targeting gpu_new/gloin
# instead of the normal/long partition -- CPU-only workload (no CUDA used
# anywhere in this pipeline), just using gloin's newer CPUs. If the scheduler
# refuses this without a GPU resource request even though we don't use one,
# uncomment the --gres line below.
# Usage: sbatch mmgbsa_rescore/diagnose_concurrency_gpu.sbatch <conformation>
#   e.g. sbatch mmgbsa_rescore/diagnose_concurrency_gpu.sbatch c1
#SBATCH --job-name=diagnose_concurrency_gpu
#SBATCH --partition=gpu_new
#SBATCH --nodelist=gloin
# #SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=64
#SBATCH --time=01:00:00
#SBATCH --output=logs/diagnose_concurrency_gpu_%j.out
#SBATCH --error=logs/diagnose_concurrency_gpu_%j.err

set -e
cd "$SLURM_SUBMIT_DIR"

CONFORMATION="$1"
if [ -z "$CONFORMATION" ]; then
    echo "Usage: sbatch mmgbsa_rescore/diagnose_concurrency_gpu.sbatch <conformation>"
    exit 1
fi

module load Miniconda3/202601
eval "$(conda shell.bash hook)"
conda activate /users/gpcr/annak/ultra-large/envs/gbsa

python mmgbsa_rescore/diagnose_concurrency.py \
    --receptor "receptor_prep/${CONFORMATION}_gmxready.pdb" \
    --vs_results_dir "/users/gpcr/annak/ultra-large/vs_results" \
    --conformation "$CONFORMATION" --cutoff -7.0 \
    --levels 1,4,8,16,32,64 --per_worker 2 --timeout 180 \
    --out_csv "logs/concurrency_report_gpu_${CONFORMATION}_${SLURM_JOB_ID}.csv" \
    --scratch_root /dev/shm
