# mmgbsa_rescore

Author: Anna Korda

MM-GBSA rescoring of each conformation's final hit set (Uni-GBSA,
single-point energy-minimized pose + GB solvation -- `mode=em`, not MD,
matching Uni-GBSA's own validated default). Runs on top of
`receptor_prep/`'s gmx-ready receptors and `apply_score_cutoff.py`'s already-
built `final_hits_<conf>_cutoff<C>.csv` + `_partNN.sdf`.

## Real benchmark (not estimated)

Timed 12 real compounds across all 3 receptors (c1, ref1, c5): consistently
**~25 seconds wall-clock per compound**, and essentially single-threaded
internally (`-nt 1`, `-nt 2`, and `-nt 8` all measured ~110-120% CPU) -- more
threads per compound don't help, so the real lever is running many compounds
as **separate concurrent processes**, which is what this script does
(`ProcessPoolExecutor` across `--workers`, not one job with many threads).

Real compound counts (cutoff -7.0): c1 = 106,380, ref1 = 51,195, c5 = 71,264,
**228,839 total**. At ~25s/compound on 64 concurrent workers:
`228,839 x 25s / 64 ~= 24.8 hours` -- about 1 day, matching the original
target.

Each compound writes ~7 files / ~3.7MB (without `--verbose`) to its own
throwaway working directory, deleted immediately after its result is parsed
-- at 228,839 compounds that's ~850GB / ~1.6M inodes if ever kept, so nothing
survives except the one real number (`TOTAL`, the MM-GBSA delta G) pulled out
of each compound's result CSV. Use `--scratch_root` to point this at
node-local storage, not shared GPFS (the sbatch script already does this via
`$TMPDIR`).

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

64 cores, `--partition=long`, no time limit, resumable (skips chunks whose
output already exists), logs real periodic progress with ETA (throttled to
~100 log lines total regardless of compound count, so c1's ~5,300 chunks
don't spam the log).

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
complex is deleted right after scoring, since keeping all ~230k would be
~850GB. This script re-runs unigbsa-pipeline for just the top `--top_pct`
percent (default 1%) of the already-ranked `mmgbsa_results_<conf>.csv` and
this time keeps `complex.pdb` (the real minimized protein+ligand structure,
confirmed by atom count -- all standard residues + ACE/NME caps + the ligand
as `MOL`). A real, small, bounded re-run: 1% of c1/ref1/c5's full sets is
~1,064 / ~512 / ~713 compounds, not a second full-scale pass.

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
