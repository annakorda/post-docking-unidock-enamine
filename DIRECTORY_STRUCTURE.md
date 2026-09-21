# Directory structure: exact names, in and out, per step

Everything lives inside the cloned repo directory. Two kinds of content:

1. **Tracked** (in git): code, small real deliverables (`nested_bm_clustering/`'s
   Q1-Q5 output, `receptor_prep/`).
2. **Gitignored** (real, but too large/instance-specific to commit):
   `envs/`, `receptors/`, `vs_results/`, `mmgbsa_rescore/results/`,
   `mmgbsa_rescore/logs/`.

```
post-docking-unidock-enamine/                 <- clone this, cd into it, run everything from here
├── filtering/                                tracked -- steps 3-5's code
│   ├── interaction_filter_fast.py, strain_filter.py, collect_final_hits.py,
│   │   apply_score_cutoff.py, torsion_lib.py, extract_review_set.py,
│   │   interaction_filter_references.txt
│   ├── data/TL_2.1_VERSION_6.xml             torsion library, used by strain_filter.py
│   ├── sbatch/*.sbatch                       submission wrappers for steps 3-4-collect
│   └── logs/                                 gitignored contents (*.out/*.err), dir tracked
├── envs/post-dock/                           gitignored, you create this (conda env)
├── receptors/
│   ├── c1.pdb  ref1.pdb  c5.pdb               gitignored, YOU provide these (raw receptors)
│   └── ...AD4 PDBQT builds
├── receptor_prep/                             tracked, already here
│   └── c1_gmxready.pdb  ref1_gmxready.pdb  c5_gmxready.pdb
└── vs_results/                                gitignored, created by running steps 3-6
    ├── ad4_redock_<conf>/                     YOU provide this (steps 0-2's output,
    │   └── <compound_id>.sdf                  separate repo). One SDF per compound,
    │                                           filename = compound_id, up to 9 poses inside.
    ├── interactions_passed_<conf>/            step 3 output
    │   ├── interactions_passed_<conf>.tar.gz
    │   ├── interaction_filter_report_<conf>.csv
    │   ├── reordered_compounds_<conf>.csv
    │   ├── ranked_hits_<conf>.csv
    │   └── run_summary_<conf>.txt
    ├── strain_filtered_<conf>/                step 4 output
    │   ├── strain_filtered_<conf>.tar.gz
    │   ├── strain_filter_report_<conf>.csv
    │   └── run_summary_<conf>.txt
    ├── final_hits_<conf>/                     filtering/collect_final_hits.py output
    │   ├── final_hits_<conf>.csv
    │   ├── final_hits_<conf>_part01.sdf, _part02.sdf, ...
    │   └── score_range_report_<conf>.txt
    ├── final_hits_<conf>_cutoff<C>/           step 5 output (C = e.g. -7.0 -> "-7")
    │   ├── final_hits_<conf>_cutoff<C>.csv
    │   └── final_hits_<conf>_cutoff<C>_part01.sdf, ...
    └── dedup/                                 step 6's SEPARATE --out_dir
        └── final_hits_<conf>_cutoff<C>_dedup/
            ├── final_hits_<conf>_cutoff<C>_dedup.csv
            └── final_hits_<conf>_cutoff<C>_dedup.sdf   (single file, not split)
```

Steps 3-5 all take one `--vs_results_dir vs_results` and find/create their
own named subfolder there. Step 6 (`dedup_stereoisomers.py`) is different:
it takes `--vs_results_dir` (reads step 5's output from there) AND a
separate `--out_dir` (writes here instead, conventionally
`vs_results/dedup`). Step 7 (`nested_bm_clustering.py`)'s own
`--vs_results_dir` must then point at step 6's `--out_dir` value
(`vs_results/dedup`), not at the top-level `vs_results/` -- it's looking
for `final_hits_<conf>_cutoff<C>_dedup/` as a direct child of whatever
path you pass it.

## Steps 7-8, tracked

```
post-docking-unidock-enamine/
├── environment.yml, install_pip.sh          the post-dock env
├── nested_bm_clustering/
│   ├── nested_bm_clustering.py              step 7a: reads vs_results/dedup/
│   │                                        final_hits_<conf>_cutoff<C>_dedup/,
│   │                                        writes here:
│   ├── bm_group_summary_<conf>.csv
│   ├── bm_scaffold_assignment_<conf>.csv
│   ├── nested_cluster_assignment_<conf>_<threshold>.csv   (x12 per conf)
│   ├── nested_cluster_summary_<conf>.csv
│   ├── q1/ ... q4/, q_singletons/           step 7b: Q1-Q5, each script's
│   │                                        own CSV + visualization, and its
│   │                                        own DEDUP_DIR derived automatically
│   │                                        as vs_results/dedup relative to the repo
│   └── q5/
│       ├── q5_mmgbsa_c1.sdf, _ref1.sdf, _c5.sdf   <- what step 8 reads
│       ├── q5_mmgbsa_input_list.csv
│       └── q5_kept_clusters_summary.csv
└── mmgbsa_rescore/
    ├── run_mmgbsa_q5.py, run_all.sh
    ├── logs/mmgbsa_<conf>.log               gitignored, created on run
    └── results/                             gitignored, created on run
        └── mmgbsa_<conf>/
            ├── mmgbsa_part_00000.csv, ...    per-chunk, intermediate
            ├── mmgbsa_results_<conf>.csv     <- THE deliverable
            ├── mmgbsa_failed_<conf>.csv
            └── mmgbsa_top20_vs_ad4_<conf>.csv
```

## Cloning fresh and dropping in your MM-GBSA results

`mmgbsa_rescore/results/` is gitignored, generated fresh by
`run_mmgbsa_q5.py`/`run_all.sh`. To bring results made on another machine
(e.g. faramir) into a fresh clone:

```bash
git clone <this repo>
cd post-docking-unidock-enamine/mmgbsa_rescore
mkdir -p results
scp -r faramir:/path/to/mmgbsa_rescore/results/* results/
```
