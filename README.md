# Post-docking-unidock-enamine

Post-docking analysis for the 5-HT1A ultra-large screen 

**Initial Library:** Enamine REAL Space 69B

**Docked:** 192M compounds vs. conformations c1/ref1/c5, vina+AD4 on CINECA and MareNostrum). 

## The filtering funnel, in order

Each step runs on the previous step's own output.

Steps 0-2 run in different repos (for now)
Steps 3-8 are this repo.

| # | step | script | c1 | ref1 | c5 |
|---|------|--------|----|------|-----|
| 0 | Enamine REAL Space, full | (separate repo) | 69,000,000,000 total | | |
| 1 | Property/PAINS/Similarity pre-screen | (separate repo) | **191,904,689** total | | |
| 2 | Top 1% by vina score -> AD4 redock | (separate repo) | **1,919,047** each | | |
| 3 | Salt-bridge interaction filter (Asp116, 5.0A) | `filtering/interaction_filter_fast.py` | **1,780,753** | **1,042,591** | **1,539,369** |
| 4 | Strain filter (Total>=7.0 OR Single>=1.8 TEU) | `filtering/strain_filter.py` -> `filtering/collect_final_hits.py` | **117,404** | **56,904** | **96,478** |
| 5 | AD4 score cutoff (<=-7.0) | `filtering/apply_score_cutoff.py` | **106,380** | **51,195** | **71,264** |
| 6 | Dedup stereoisomers/tautomers | `dedup_stereoisomers.py` | **95,468** | **47,790** | **65,364** |
| 7 | BM & Tanimoto clustering, Q1-Q5 selection | `nested_bm_clustering/` | **2,770** | **1,502** | **1,794** |
| 8 | MM-GBSA rescoring | `mmgbsa_rescore/run_mmgbsa_q5.py` | - | - | - |

### Step 1's pre-screen (69B -> 191,904,689): 

MW 200-500 Da, LogP 0-5, 1-2 positive charge centers, 0 negative charge centers, <=2 chiral centers, <=10 rotatable bonds, TPSA <=140, 20-40 heavy atoms, 2-6 rings, 2-12 heteroatoms, Enamine class S only (M and U dropped), PAINS filtered out, Tanimoto <0.35 to 5,250 known GPCRdb actives.

## Data

Everything lives inside this one cloned directory. `vs_results/`, `receptors/`, and `envs/` are real data/output, gitignored.  Clone the repo, then those three get created by running the pipeline:

```
post-docking-unidock-enamine/      
├── envs/post-dock/                 <- conda env (gitignored, you create it)
├── receptors/                      <- c1.pdb / c5.pdb / ref1.pdb + AD4 PDBQT builds
│                                       (gitignored, you provide these)
├── receptor_prep/                  <- gmx-ready receptors, tracked 
├── vs_results/                     <- steps 3-6's real output (gitignored, created by
│                                       running the pipeline)
├── nested_bm_clustering/           <- step 7, tracked (Q1-Q5 scripts + real output)
└── mmgbsa_rescore/                 <- step 8, tracked (scripts); results/ gitignored
```

Exact filenames at every step: **`DIRECTORY_STRUCTURE.md`**.

## Environment
```
conda env create -f environment.yml -p envs/post-dock
conda activate envs/post-dock
bash install_pip.sh
```

## Step 3: `filtering/interaction_filter_fast.py`

**Does a compound's docked pose put a charged nitrogen within 5.0A of Asp116's CG carbon?** Checks every pose in the SDF (unidock's real default is up to 9 poses/compound), picks the winner among salt-bridge-qualifying poses by best AD4 score, reorders it to the front if it wasn't already there. Molecules are never modified (`sanitize=False`, byte-identical pass-through). 

```
sbatch filtering/sbatch/interaction_filter_fast.sbatch <conformation> <receptor_pdb_path>
```

Output under `vs_results/interactions_passed_<conf>/`: the passing SDFs
(tar.gz), `interaction_filter_report_<conf>.csv` (every compound's
pass/fail/error), `reordered_compounds_<conf>.csv`, `ranked_hits_<conf>.csv`
(best-score-first), `run_summary_<conf>.txt`.

## Step 4: `filtering/strain_filter.py`

**Filters the interaction filter's winning poses by conformational strain**,
using Gu, Smith, Yang, Irwin, Shoichet, "Ligand Strain Energy in Large
Library Docking," JCIM 2021 (citation + DOI in
[filtering/interaction_filter_references.txt](filtering/interaction_filter_references.txt);
independently reimplemented in `filtering/torsion_lib.py`, verified
against the paper's own reference output before use). A pose is
"strained" if `total_TEU >= 7.0 OR single_TEU >= 1.8` (the paper's own
general default).

```
sbatch filtering/sbatch/strain_filter.sbatch <conformation> --dry_run   # counts only
sbatch filtering/sbatch/strain_filter.sbatch <conformation>              # writes survivors tar
```

Output under `vs_results/strain_filtered_<conf>/`:
`strain_filter_report_<conf>.csv`, `run_summary_<conf>.txt`,
`strain_filtered_<conf>.tar.gz`.

## `filtering/collect_final_hits.py`

**Joins the interaction & strain filter reports, ranks by AD4 score**, splits
into multi-SDF files (default 50,000/file).

```
sbatch filtering/sbatch/collect_final_hits.sbatch <conformation> [chunk_size]
```

Output under `vs_results/final_hits_<conf>/`: `final_hits_<conf>.csv`
(rank, compound_id, ad4_score, total_TEU, single_TEU),
`score_range_report_<conf>.txt`, `final_hits_<conf>_partNN.sdf`.

## Step 5: `filtering/apply_score_cutoff.py`

Cuts `collect_final_hits.py`'s already-ranked output to `ad4_score <=
--cutoff` (-7.0 used for the run above). 
```
python filtering/apply_score_cutoff.py --conformation <conf> \
    --vs_results_dir vs_results --cutoff -7.0
```

Output under `vs_results/final_hits_<conf>_cutoff<C>/`: same columns as
above, re-split into `_partNN.sdf` files.

## Step 6: `dedup_stereoisomers/`

Many rows in step 5's output are the same 2D compound, enumerated as
different stereoisomers/tautomers by the original library build. **Keeps
only the best-AD4-scoring version of each.**

```
python dedup_stereoisomers/dedup_stereoisomers.py --conformation <conf> \
    --vs_results_dir vs_results \
    --out_dir vs_results/dedup --cutoff -7.0
```

Output under `vs_results/dedup/final_hits_<conf>_cutoff<C>_dedup/`. Step 7
(below) needs its own `--vs_results_dir` pointed at this `--out_dir`
value (`vs_results/dedup`), not at the top-level `vs_results/`.

## Step 7: `nested_bm_clustering/`

**Two-pass clustering (Bemis-Murcko scaffold groups, then Tanimoto/Butina
within large groups) plus a 5-question pipeline** (Q1-Q5, full reasoning in
`nested_bm_clustering/README.md`) that **picks a validated, statistically
grounded set of compounds for MM-GBSA**, not just "top N by AD4 score."

![Q5](nested_bm_clustering/q5/q5_visualization.png)

Output: `q5/q5_mmgbsa_<conf>.sdf` + `q5/q5_mmgbsa_input_list.csv`.

## Step 8: `mmgbsa_rescore/`

**MM-GBSA rescoring** (Uni-GBSA, single-point energy-minimized pose + GB
solvation) of step 7's curated list. Needs `receptor_prep/`'s gmx-ready
receptors. Full usage, environment setup, and real benchmark numbers in
`mmgbsa_rescore/README.md`.

```
bash mmgbsa_rescore/run_all.sh <workers> <scratch_root>
```

Output under `mmgbsa_rescore/results/mmgbsa_<conf>/`:
`mmgbsa_results_<conf>.csv` (ranked by `mmgbsa_dg`, with `ad4_score`,
`cluster_id`, `role`, `source` carried through from step 7).
