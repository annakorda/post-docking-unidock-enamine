# post-docking-unidock-enamine

Post-docking analysis for the 5-HT1A ultra-large screen (Enamine REAL Space,
~192M compounds vs. conformations c1/ref1/c5, docked+AD4-redocked on CINECA
and MareNostrum). This repo runs on Barcelona (`shiva`) only, against the
results copied over from both clusters.

## Where the data lives

This repo is cloned inside `~/ultra-large/`, alongside the data it operates
on (not tracked in git -- too large):

```
~/ultra-large/
├── post-docking-unidock-enamine/   <- this repo
├── envs/luna-chem/                 <- conda env (see below)
├── receptors/                      <- c1.pdb / c5.pdb / ref1.pdb + per-conformation *_build/
│                                       (AD4 PDBQT + autogrid4 maps -- not what this script reads)
├── top1_undocked/                  <- pre-docking top-1% compound tars, one per conformation
├── vs_results/                     <- AD4-redocked poses (scores.csv + *.sdf), one dir per conformation
└── scores_vina/                    <- full 192M-compound vina-pass scores (scores only, no poses)
    ├── c1/  ref1/  c5/
```

## Environment

```
conda env create -f environment.yml -p ~/ultra-large/envs/luna-chem
conda activate ~/ultra-large/envs/luna-chem
bash install_pip.sh
```

The env is named `luna-chem` and `install_pip.sh` installs LUNA for
historical reasons (an earlier version of `interaction_filter.py` used it --
see git history), but the current filter (`interaction_filter_fast.py`)
only needs RDKit, which is already a dependency in this same env. No need
to rebuild the env because of this.

## interaction_filter_fast.py

Filters AD4-redocked poses by a direct geometric criterion: is the ligand's
**charged nitrogen** atom within **5.0A** of **Asp116's CG carbon** (the
carbon bonded to both carboxylate oxygens)? Same methodology and cutoff as
this project's own `analyze_asp_contact.py` (protonated N via formal
charge, distance to Asp116's CG, 5.0A). Pure RDKit, no receptor-wide
interaction perception of any kind -- the whole check is one distance
calculation per charged N atom found.

Molecules are parsed with `sanitize=False` and are **never modified**:
no kekulization, no valence/aromaticity perception, no coordinate changes.
Only the raw connection table (atoms, bonds, coordinates, formal charges
from `M CHG` lines) is read. Compounds that pass are copied byte-for-byte
(`shutil.copy2`) into the output tar -- verified locally with `diff` against
the original file.

Run via:
```
sbatch sbatch/interaction_filter_fast.sbatch <conformation> <receptor_pdb_path>
```

Output, per conformation, under `vs_results/interactions_passed_<conf>/`:
- `interactions_passed_<conf>.tar.gz` -- flat SDFs of the compounds that
  passed, split across equal-sized internal parts (concatenated gzip --
  extract with `tar -ixzf`, plain `-xzf` will silently only give you the
  first part). No nested folders inside.
- `interaction_filter_report_<conf>.csv` -- one row per compound
  (`compound_id, passed, error`), the real pass/fail/error record for
  every compound processed, not just the ones that passed.
- `run_summary_<conf>.txt` -- human-readable summary written inside this
  same folder: starting compound count, pass/fail/error counts, the exact
  criterion and receptor coordinate used, output file sizes.

Resumable and parallel (chunked `ProcessPoolExecutor`, same pattern as the
rest of this project's large-file-count scripts) -- safe to re-run after a
partial failure or timeout, already-done chunks are skipped. The `sbatch`
job's own log is written with frequent per-chunk progress lines (real
running pass/fail/error counts, not just start/end), including periodic
updates during the initial file-listing walk, which can itself take a
while over ~1.9M files on GPFS.

### Validated locally

Real ligand (`benchmarking/build/prepped_udt/`, `M CHG 1 6 1` -- an actual
formal charge, not synthetic) + real `c5_capsfixed.pdb` receptor. Two
controls: the ligand left untouched (charged N nowhere near the receptor)
vs. rigid-translated so its charged N sits 2.8A from the real Asp116 CG
coordinate (99.500, 110.803, 111.008, read directly from the PDB and
confirmed to be the only match -- the script errors out rather than
guessing if that residue is ambiguous or missing). Result: `negative_control`
-> fail, `positive_control` -> pass, correctly. Scaled to 200 real compounds
(100 duplicated positive + 100 duplicated negative): **100 pass, 100 fail,
exactly right**, in ~0.2s total (~0.5-1ms/compound with real per-file I/O).
The saved output SDF was diffed byte-for-byte against the original input --
identical, confirming nothing gets modified.

An earlier LUNA-based version of this filter (`interaction_filter.py`,
removed from this repo -- see git history if needed) checked the same
Asp116 contact via LUNA's general-purpose interaction-perception machinery.
Real production throughput was ~15-19s/compound even after four targeted
fixes (receptor-perception caching, restricting which interaction types
LUNA computes, disabling its per-operation GPFS logging, and searching only
the ligand's own atom groups instead of the whole cached neighborhood) --
at ~1.9M compounds that doesn't finish in a reasonable time. This script
replaces it.

## Open items

- `c1`'s AD4-redock tar is still mid-transfer as of this writing.
