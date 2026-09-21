# mmgbsa_rescore

Author: Anna Korda

MM-GBSA rescoring of `nested_bm_clustering/`'s curated Q5 hit list (Uni-GBSA,
single-point energy-minimized pose plus GB solvation: `mode=em`, Uni-GBSA's
own validated default, not MD). Runs on `receptor_prep/`'s gmx-ready
receptors and `nested_bm_clustering/q5/q5_mmgbsa_<conf>.sdf`.

## Environment

**Load-bearing gotcha**: a plain `pip install unigbsa` pulls in
`openbabel-wheel` as a transitive dependency, which silently shadows the
conda-installed `openbabel` and is ABI-incompatible with it. `obabel` then
fails on every ligand with `undefined symbol:
_ZN9OpenBabel8OBPlugin7DisplayERSsPKcS3_`. `install_pip.sh` installs with
`--no-deps` and repairs it. Don't improvise a plain `pip install` here.

```
conda env create -f environment.yml -p /path/to/envs/gbsa
conda activate /path/to/envs/gbsa
bash install_pip.sh
```

## Run

```
bash run_all.sh <workers> <scratch_root>
```

Runs c1/ref1/c5 sequentially, each parallelized internally across
`--workers` processes. Chunked and resumable: safe to interrupt/rerun,
already-completed chunks skip instantly. `<scratch_root>` should be fast
local storage (`/dev/shm` if there's enough RAM); each compound writes
~7 files/~3.7MB to its own throwaway dir, deleted right after scoring.

Or run one conformation directly:
```
python run_mmgbsa_q5.py --conformation c1 \
    --sdf ../nested_bm_clustering/q5/q5_mmgbsa_c1.sdf \
    --receptor ../receptor_prep/c1_gmxready.pdb \
    --out_dir results --workers 10 --scratch_root /dev/shm
```

Measured per-compound cost: ~25-32s on a fast local machine, ~65-72s on a
shared cluster at low concurrency, with gradual (not catastrophic) further
slowdown as concurrency rises:

| workers | mean s/compound | slowdown |
|---|---|---|
| 1 | 72.2 | 1.00x |
| 4 | 75.9 | 1.06x |
| 8 | 90.7 | 1.29x |
| 16 | 120.8 | 1.72x |

Two cheap mitigations for that degradation, both wired into `run_all.sh`:
point `--scratch_root` at genuinely fast/local storage (disk contention was
the main suspected cause), and `OMPI_MCA_btl=self` plus a `/dev/shm`-based
OpenMPI session dir (unigbsa does a per-process OpenMPI singleton-init
internally via `mpi4py`, a second, smaller contributing factor).

## Output

Under `results/mmgbsa_<conf>/`:
- `mmgbsa_results_<conf>.csv`: every compound that succeeded, ranked by
  `mmgbsa_dg` (most negative first), with `ad4_score`, `cluster_id`,
  `role` (medoid/best_scorer/medoid+best_scorer/singleton_outlier), and
  `source` (cluster/singleton) carried through from step 7.
- `mmgbsa_failed_<conf>.csv`: anything that errored or timed out, same
  columns, real error message.
- `mmgbsa_top20_vs_ad4_<conf>.csv`: top 20 MM-GBSA hits with their AD4
  rank within the Q5 set, a quick sanity check that MM-GBSA isn't
  reordering toward compounds AD4/clustering already considered weak.

Live progress is logged per compound, not just per chunk, with a running
ok/fail count and an explicit flag if the failure rate exceeds 30% after
the first 20+ compounds. A systemic problem (wrong receptor, broken env)
shows up within minutes, not only after the full run finishes.

## `clean_failed_chunks.py`

Deletes any `mmgbsa_part_*.csv` chunk file that contains a real failure (a
row with a non-empty `error` column). Run before resubmitting so those
compounds get genuinely recomputed instead of resumability trusting stale
failures/timeouts from an earlier run. Chunks with no failures are left
alone.

## Open items

- **Saving top complexes** (`complex.pdb` for the best final scorers, not
  just the delta-G) isn't wired up for the curated Q5 list yet. Each
  compound's minimized complex is currently deleted right after scoring.
- **Combining AD4 + MM-GBSA into one score**: not yet implemented,
  deliberately. Run the real thing first and check results before
  deciding. Combining raw score *values* (even Z-scored) across a docking
  score and a free energy is a documented pitfall (incompatible
  scales/units). The literature-backed alternative is rank-based
  consensus: plain rank-averaging, or Exponential Consensus Ranking
  (Palacio-Rodriguez, Lans, Cavasotto, Cossio, *Sci. Rep.* 2019,
  10.1038/s41598-019-41594-3), which converts each method's score to a
  rank first, then combines via `p(r) = (1/sigma) * exp(-r/sigma)` summed
  across methods. Foundational consensus-scoring concept: Charifson et
  al., *J. Med. Chem.* 1999, PMID 10602695.
