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

Checks **every pose in the SDF, not just the first**. unidock's real
defaults (`unidock --help`, confirmed, not assumed) are `--num_modes 9` and
`--energy_range 3` (kcal/mol) -- AD4 redock output can hold up to 9 poses
per compound, and every pose actually written is already guaranteed by
unidock itself to be within 3 kcal/mol of the best-scored pose. A compound
passes if *any* of its poses makes the contact, not just its top-ranked-by-
score one -- verified with a real multi-pose test file where pose 1 alone
would fail but pose 2 hits; the compound correctly passes, and the full
multi-pose file is copied through untouched (both poses still present,
byte-identical to the input).

Molecules are parsed with `sanitize=False` and are **never modified**:
no kekulization, no valence/aromaticity perception, no coordinate changes.
Only the raw connection table (atoms, bonds, coordinates, formal charges
from `M CHG` lines) is read. Compounds that pass are copied byte-for-byte
(`shutil.copy2`) into the output tar -- verified locally with `diff` against
the original file.

**When a compound's winning pose isn't the first one in its SDF, that pose
gets moved to the front** (pure text-level split on the `$$$$` delimiter and
reassembly, no chemistry library involved -- if the split/rejoin doesn't
reproduce the original file exactly, it falls back to a plain untouched
copy rather than risk writing something wrong). Every individual pose's
bytes stay identical; only their order changes, and only for compounds
where it matters. Verified with a real multi-pose file (pose 1 alone would
fail, pose 2 hits) -- the winning pose correctly ends up first, both poses
still present and byte-identical to the original.

**Pose selection is salt-bridge criterion first, best score second**: among
a compound's poses, only those within 5.0A even qualify; the winner among
*those* is whichever has the best (most negative) AD4 score (`ENERGY=` tag,
same regex as `run_benchmark.py`/`extract_top_poses.py`, read from that
specific pose's own text block so distance and score always belong to the
same pose) -- not necessarily the geometrically closest one. Verified with
a real 3-pose case: pose 1 fails geometry despite a deceptively good score,
pose 2 is closer (2.0A) but scores worse (-5.0), pose 3 is farther (4.5A,
still under the cutoff) but scores best (-9.0) -- the compound correctly
picks pose 3, both in the ranking and in the reordered output file.

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
- `reordered_compounds_<conf>.csv` -- one row per compound whose winning
  pose wasn't already first (`compound_id, winning_pose_index_0based,
  distance_A`), so it's possible to check exactly which compounds and
  which pose got picked.
- `ranked_hits_<conf>.csv` -- every passing compound's winning-pose score,
  **sorted best (most negative) first** (`rank, compound_id,
  winning_pose_score`) -- the actual hit-priority list. Compounds whose
  winning pose has no `ENERGY=` tag are excluded here (logged as a warning
  with a real count) but still present in the full report.
- `run_summary_<conf>.txt` -- human-readable summary written inside this
  same folder: starting compound count, pass/fail/error/reordered counts,
  the exact criterion and receptor coordinate used, output file sizes.

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

## strain_filter.py (second filter stage, runs after interaction_filter_fast.py)

Filters salt-bridge-passing compounds by conformational strain, using the
real method from Gu, Smith, Yang, Irwin, Shoichet, "Ligand Strain Energy in
Large Library Docking," JCIM 2021 (full citation +
[interaction_filter_references.txt](interaction_filter_references.txt))
-- independently reimplemented (`torsion_lib.py`) from their paper and
public source since their repo has no license, and verified against their
own published reference output on 3 real molecules to floating-point
precision before use. Uses their real, openly-distributed torsion-library
data file (`data/TL_2.1_VERSION_6.xml`, unmodified) -- this is the
precomputed statistical result of their Cambridge Structural Database
analysis, so no CSD access of our own is needed.

Checks the exact same winning pose interaction_filter_fast.py already
selected per compound (salt bridge first, best AD4 score among qualifying
poses second) -- never re-opens pose selection at this stage. A compound
is "strained" (filtered out) if `total_TEU >= 7.0 OR single_TEU >= 1.8`,
the paper's own general-purpose default (from their 40-target DUD-E
benchmark average) -- no serotonin receptor is in DUD-E, so this isn't
calibrated to our target specifically; their real dopamine D4 case study
(the closer aminergic-GPCR analog, same conserved Asp3.32 architecture)
used Total>=6.0/Single>=1.8 from real experimental hit-rate data, which
nearly doubled hit rate (24.2% -> 34.0%) at the cost of discarding 60% of
compounds -- worth revisiting if the DUD-E-general default doesn't behave
well on our own real numbers.

Molecules are never modified -- the torsion library's SMARTS use lowercase
aromatic atoms (`[cX3]`, `[n]`), which need RDKit's aromaticity perception
to match; our SDFs are Kekulized (no file-level aromatic bond type, unlike
the paper's own mol2 input) so a *partial* sanitization (ring perception +
aromaticity only, no valence/H changes) is applied to an in-memory copy --
verified coordinates/charges/atom-count stay byte-identical, and that this
step is a no-op on already-aromatic-flagged input (their own mol2 test
data reproduced identical numbers with or without it).

Run via:
```
sbatch sbatch/strain_filter.sbatch <conformation> --dry_run   # counts only, no survivors tar
sbatch sbatch/strain_filter.sbatch <conformation>              # writes strain_filtered_<conf>.tar.gz
```
Output under `vs_results/strain_filtered_<conf>/`: `strain_filter_report_<conf>.csv`
(compound_id, total_TEU, single_TEU, flagged, n_torsions, strained, error),
`run_summary_<conf>.txt`, and (non-dry-run) `strain_filtered_<conf>.tar.gz`
-- the compounds that passed *both* filters, same flat/untouched-file
convention as interaction_filter_fast.py's own output.

**Real result on c1's top-10 ranked hits**: 9 of 10 flagged strained at the
default thresholds -- a real, substantial cut, not noise. Matches the
paper's own documented pattern: a high-scoring docked pose can achieve
that score *by* adopting a strained conformation (their Figures 3/5/8 show
real examples of exactly this).

## Open items

- Real full-scale counts (via `--dry_run`) not yet run for any conformation
  as of this writing -- the numbers above are from a 10-compound sample.
