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
│                                       (AD4 PDBQT + autogrid4 maps -- not what LUNA reads, see below)
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

`environment.yml` is conda-only (pymol, openbabel, biopython, the usual
scientific stack). `install_pip.sh` installs LUNA and its pip-only deps in
the specific order that actually works -- see the comments in that file for
why the order matters, don't reorder it.

## interaction_filter.py

Filters AD4-redocked poses down to the ones that actually make a real
receptor contact, instead of trusting AD4 score alone: does the ligand form
an **Ionic** or **Salt bridge** interaction (via
[LUNA](https://luna-toolkit.readthedocs.io)) with **Asp116**? This is the
same salt-bridge contact this project's DUDE benchmark used to validate AD4
as the more reliable scoring function for this target.

LUNA perceives charged/ionizable atoms itself (via RDKit's standard
pharmacophore feature typing) -- no manual charge-parsing needed.

`--receptor_pdb` must point at the plain `receptors/<conf>.pdb`, not the
`_build/receptor_ad4.pdbqt` -- LUNA's `pdb_path` is parsed with Biopython's
standard PDB parser, which doesn't understand PDBQT's AutoDock-specific atom
types (`OA`, `HD`, ...) or its merged nonpolar hydrogens.

Run via:
```
sbatch sbatch/interaction_filter.sbatch <conformation> <receptor_pdb_path>
```

Output, per conformation, under `vs_results/interactions_filtered_<conf>/`:
- `interactions_filtered_<conf>.tar.gz` -- flat SDFs of the compounds that
  passed, split across internal parts (concatenated gzip -- extract with
  `tar -ixzf`, plain `-xzf` will silently only give you the first part).
  No nested folders inside.
- `interaction_filter_report_<conf>.csv` -- one row per compound
  (`compound_id, passed, error`), the real pass/fail/error count for every
  compound processed, not just the ones that passed.

Resumable and parallel (chunked `ProcessPoolExecutor`, same pattern as the
rest of this project's large-file-count scripts) -- safe to re-run after a
partial failure or timeout, already-done chunks are skipped.

## Open items

- **`--chunk_size` (default 100) hasn't been benchmarked for real** -- unlike
  this project's pure file-copy chunking elsewhere, this is a real
  per-compound LUNA interaction-calculation cost. Run a small chunk first
  (e.g. the 15 sample ligands already pulled locally during AD4 spot-checks)
  and time it before trusting the default at full scale.
- `c1`'s AD4-redock tar is still mid-transfer as of this writing.
