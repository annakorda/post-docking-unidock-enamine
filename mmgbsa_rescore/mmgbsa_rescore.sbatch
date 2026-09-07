#!/bin/bash
# Author: Anna Korda
# Usage: sbatch mmgbsa_rescore/mmgbsa_rescore.sbatch <conformation> [cutoff] [workers]
#   e.g. sbatch mmgbsa_rescore/mmgbsa_rescore.sbatch c1
#   e.g. sbatch mmgbsa_rescore/mmgbsa_rescore.sbatch c1 -7.0 32
#SBATCH --job-name=mmgbsa_rescore
#SBATCH --partition=long
#SBATCH --cpus-per-task=64
#SBATCH --output=logs/mmgbsa_rescore_%j.out
#SBATCH --error=logs/mmgbsa_rescore_%j.err

set -e
cd "$SLURM_SUBMIT_DIR"

CONFORMATION="$1"
CUTOFF="${2:--7.0}"
# --workers is deliberately NOT tied 1:1 to --cpus-per-task=64 above: each
# compound's gmx_MMPBSA step does a real per-process OpenMPI singleton init
# (confirmed in unigbsa's own source -- it calls mpi4py's MPI.COMM_WORLD.Barrier()
# even without mpirun), and a real production run at --workers=64 stalled hard
# (processes stuck 45-70+ min accumulating only ~20-30s of actual CPU time, many
# in D/disk-wait state) -- almost certainly contention between 64 simultaneous
# MPI inits on one node, not a real compute bottleneck. Keeping the full 64-core
# allocation (headroom) but running fewer concurrent compound-jobs at once cuts
# that contention directly; raise --cpus-per-task above 64 too if the cluster
# has bigger nodes and this is still not enough headroom.
WORKERS="${3:-16}"
if [ -z "$CONFORMATION" ]; then
    echo "Usage: sbatch mmgbsa_rescore/mmgbsa_rescore.sbatch <conformation> [cutoff] [workers]"
    exit 1
fi

module load Miniconda3/202601
eval "$(conda shell.bash hook)"
conda activate /users/gpcr/annak/ultra-large/envs/gbsa

# tell OpenMPI not to bother with real inter-process transport setup for what's
# always a 1-rank job here, and keep its own session/rendezvous files off
# whatever /tmp actually is on this node
export OMPI_MCA_btl=self
export OMPI_MCA_orte_tmpdir_base="/dev/shm/ompi_${SLURM_JOB_ID}"
mkdir -p "$OMPI_MCA_orte_tmpdir_base"

# /dev/shm for our own per-compound throwaway dirs too -- RAM-backed, always
# node-local regardless of how /tmp is mounted; real footprint here (~3.7MB x
# up to $WORKERS concurrent) is trivial for RAM
SCRATCH_ROOT="/dev/shm/mmgbsa_${CONFORMATION}_${SLURM_JOB_ID}"
mkdir -p "$SCRATCH_ROOT"

python mmgbsa_rescore/mmgbsa_rescore.py \
    --conformation "$CONFORMATION" \
    --vs_results_dir "/users/gpcr/annak/ultra-large/vs_results" \
    --receptor "receptor_prep/${CONFORMATION}_gmxready.pdb" \
    --out_dir "/users/gpcr/annak/ultra-large/analysis_results/mmgbsa_rescore" \
    --cutoff "$CUTOFF" \
    --workers "$WORKERS" \
    --scratch_root "$SCRATCH_ROOT"

rm -rf "$SCRATCH_ROOT" "$OMPI_MCA_orte_tmpdir_base"
