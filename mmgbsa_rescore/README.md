# mmgbsa_rescore

Author: Anna Korda

MM-GBSA rescoring of each conformation's final hit set (Uni-GBSA,
single-point energy-minimized pose + GB solvation -- `mode=em`, not MD,
matching Uni-GBSA's own validated default). Runs on top of
`receptor_prep/`'s gmx-ready receptors and `apply_score_cutoff.py`'s already-
built `final_hits_<conf>_cutoff<C>.csv` + `_partNN.sdf`.

## Real benchmark, and why this only scores the top 1%

Per-compound cost is genuinely higher on shiva than a fast local machine --
measured directly (not estimated): **~65-72s/compound on shiva** at
low concurrency vs ~25-32s locally, for the identical code and receptor.
On top of that, real concurrency-scaling cost was measured directly on shiva
(same node, same code, only worker count changed):

| workers | mean s/compound | slowdown |
|---|---|---|
| 1 | 72.2 | 1.00x |
| 4 | 75.9 | 1.06x |
| 8 | 90.7 | 1.29x |
| 16 | 120.8 | 1.72x |

That's a real, gradual, measured curve, not a guess -- it also means an
earlier attempt to run this at 64 workers correctly stalled hard (many
`gmx_MMPBSA`/`unigbsa-pipeline` processes stuck in D/disk-wait state,
accumulating only ~20-30s of CPU time after 45-70+ min of wall-clock
existence). At this real cost, the full ~106k-230k compound sets aren't
practical here.

**Fix: only score the top `--top_pct` by AD4 rank** (default 1%, a
real contiguous prefix of the already-rank-ordered
`final_hits_<conf>_cutoff<C>.csv` -- no re-sorting needed). That's a small,
bounded, real set per conformation: c1 ~1,064, ref1 ~512, c5 ~713 -- already
the most interesting compounds (best AD4-ranked), and at even the slower
16-worker rate (1.72x), c1's ~1,064 compounds is ~1,064 x 72s x 1.72 / 16
=~ 2.3 hours, not days.

`--workers` defaults to 16 (not tied to `--cpus-per-task=32`) as a real,
measured middle ground -- tune via the 4th CLI arg if a given run still
looks slow (check with `ps aux | grep unigbsa-pipeline` on the compute node
and watch whether CPU time tracks wall-clock time). `OMPI_MCA_btl=self` and
moving OpenMPI's own session directory to `/dev/shm` are also set (harmless,
possibly-contributing mitigations for the real per-process OpenMPI
singleton-init `unigbsa` does internally via `mpi4py`), though the measured
curve above shows the degradation is gradual rather than a hard wall, so
these aren't confirmed as the primary fix on their own.

Each compound writes ~7 files / ~3.7MB (without `--verbose`) to its own
throwaway working directory, deleted immediately after its result is parsed.
Use `--scratch_root` to point this at genuinely node-local storage (`/dev/shm`
in the sbatch scripts).

## Environment

**Real, load-bearing gotcha, found the hard way**: a plain
`pip install unigbsa` pulls in `openbabel-wheel` as a transitive dependency,
which silently shadows the conda-installed `openbabel` above it and is
ABI-incompatible with it -- `obabel` then fails with `undefined symbol:
_ZN9OpenBabel8OBPlugin7DisplayERSsPKcS3_` on every single ligand (100% failure
rate, not a flaky/occasional thing). Root cause: conda's `openbabel` package
here provides the CLI binary + C++ library but **not** an importable Python
`openbabel` module (needed by Uni-GBSA's net-charge step) -- that only comes
from the pip wheel, which also overwrites the conda CLI binary on install.
`install_pip.sh` installs `unigbsa`/`lickit` with `--no-deps`, adds
`openbabel-wheel` explicitly for the Python module, then force-reinstalls the
conda `openbabel` package to restore its native CLI binary (which lives at a
different file than the Python module, so this repair doesn't undo the
Python import). Verified both `obabel -V` (native, conda) and
`python -c "from openbabel import openbabel"` (pip wheel) work together
afterward, and re-ran real compounds successfully.

```
conda env create -f environment.yml -p /users/gpcr/annak/ultra-large/envs/gbsa
conda activate /users/gpcr/annak/ultra-large/envs/gbsa
bash install_pip.sh
```

## Run

```
sbatch mmgbsa_rescore/mmgbsa_rescore.sbatch c1
sbatch mmgbsa_rescore/mmgbsa_rescore.sbatch ref1
sbatch mmgbsa_rescore/mmgbsa_rescore.sbatch c5
```

Top 1% by AD4 rank (default, override with a 2nd CLI arg), 16 workers
(default, override with a 4th CLI arg), `--partition=long`, resumable
(skips chunks whose output already exists), logs real periodic progress
with ETA.

Output per conformation, under `analysis_results/mmgbsa_rescore/mmgbsa_<conf>/`:
- `mmgbsa_results_<conf>.csv` -- every compound that succeeded, ranked by
  `mmgbsa_dg` (most negative first), **with everything already known about
  it**: `ad4_score`, `total_TEU`, `single_TEU` alongside the MM-GBSA number --
  one CSV, not several. Only the final delta G is kept from Uni-GBSA's own
  output; its full per-term energy decomposition (van der Waals,
  electrostatic, polar/non-polar solvation, ...) is discarded per compound,
  not accumulated anywhere.
- `mmgbsa_failed_<conf>.csv` -- same columns, for anything that errored or
  timed out, with the real error message.
- `mmgbsa_top20_vs_ad4_<conf>.csv` -- the top 20 MM-GBSA hits with their AD4
  score and AD4-only rank (within the full set), for a quick sanity check:
  do the MM-GBSA-best compounds also look reasonable by AD4, or is MM-GBSA
  picking things AD4 considered bad? (Real small-scale test: top-10 MM-GBSA
  hits all had AD4 ranks within 1-10 too -- a reordering within the already-
  good set, not noise.)

## `mmgbsa_save_top_complexes.py` (run after the above, per conformation)

`mmgbsa_rescore.py` keeps only the delta G per compound -- its minimized
complex is deleted right after scoring. This script re-runs unigbsa-pipeline
for just the top `--top_pct` percent (default 1%) of the already-ranked
`mmgbsa_results_<conf>.csv` and this time keeps `complex.pdb`. **Note**: now
that `mmgbsa_rescore.py` itself already only scores the top 1% by AD4 rank,
this script's default top 1% is 1% *of that* (~10-11 compounds/conformation)
-- pass a higher `--top_pct` (e.g. 50 or 100) if you want structures for
more of the already-small MM-GBSA-scored set.

```
sbatch mmgbsa_rescore/mmgbsa_save_top_complexes.sbatch c1
sbatch mmgbsa_rescore/mmgbsa_save_top_complexes.sbatch ref1
sbatch mmgbsa_rescore/mmgbsa_save_top_complexes.sbatch c5
```

Output under `analysis_results/mmgbsa_rescore/mmgbsa_<conf>/top_complexes/`:
`complexes/<compound_id>_complex.pdb` (one per selected compound) and
`top_complexes_manifest_<conf>.csv` (rank, compound_id, mmgbsa_dg, ad4_score,
path to its complex.pdb).

## Combining AD4 + MM-GBSA into one score (open item)

Not yet implemented -- deliberately, since Anna asked to run the real thing
first and check the results before deciding this. Real literature finding
for when that's next: combining raw score *values* (even Z-scored) is a
documented pitfall (incompatible scales/units between a docking score and a
free energy). The literature-backed alternative is rank-based consensus --
either plain rank-averaging, or Exponential Consensus Ranking (ECR; Palacio-
Rodriguez, Lans, Cavasotto, Cossio, *Sci. Rep.* 2019,
10.1038/s41598-019-41594-3), which converts each method's score to a rank
first, then combines via `p(r) = (1/sigma) * exp(-r/sigma)` summed across
methods -- shown to outperform Z-score/raw averaging and to be robust to its
one free parameter (sigma). Foundational consensus-scoring concept: Charifson
et al., *J. Med. Chem.* 1999, PMID 10602695. No published example found doing
exactly "AD4 + MM-GBSA on an ultra-large REAL-space screen" as one blended
score -- the pattern in that literature (e.g. Lyu et al. 2019-style funnels)
is MM-GBSA as a downstream re-ranking/filter stage on an already-shortlisted
set, which is also what this script's own place in the pipeline already is.
