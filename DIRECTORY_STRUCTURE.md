# Directory structure: exact names, in and out, per step

Two kinds of storage in this pipeline:

1. **External** (`~/ultra-large/`, not tracked in git): steps 3-6's real
   output. Too large to commit, one flat dir per step/conformation.
2. **Repo-relative** (tracked structure, gitignored contents): steps 7-8's
   output. Small enough to sit next to the code that made it.

## External: `~/ultra-large/`

```
~/ultra-large/
├── envs/post-dock/                          conda env (see root README.md)
├── receptors/
│   ├── c1.pdb  ref1.pdb  c5.pdb              raw receptors, YOU provide these
│   └── ...AD4 PDBQT builds
├── receptor_prep/
│   └── c1_gmxready.pdb  ref1_gmxready.pdb  c5_gmxready.pdb
│                                              gmx-ready receptors -- already in
│                                              this repo's receptor_prep/, copy
│                                              or symlink from there
└── vs_results/
    ├── ad4_redock_<conf>/                    YOU provide this (steps 0-2's
    │   └── <compound_id>.sdf                 output, separate repo). One SDF
    │                                          per compound, filename = compound_id,
    │                                          up to 9 poses inside each.
    ├── interactions_passed_<conf>/           step 3 output
    │   ├── interactions_passed_<conf>.tar.gz
    │   ├── interaction_filter_report_<conf>.csv
    │   ├── reordered_compounds_<conf>.csv
    │   ├── ranked_hits_<conf>.csv
    │   └── run_summary_<conf>.txt
    ├── strain_filtered_<conf>/               step 4 output
    │   ├── strain_filtered_<conf>.tar.gz
    │   ├── strain_filter_report_<conf>.csv
    │   └── run_summary_<conf>.txt
    ├── final_hits_<conf>/                    collect_final_hits.py output
    │   ├── final_hits_<conf>.csv
    │   ├── final_hits_<conf>_part01.sdf, _part02.sdf, ...
    │   └── score_range_report_<conf>.txt
    ├── final_hits_<conf>_cutoff<C>/          step 5 output (C = e.g. -7.0 -> "-7")
    │   ├── final_hits_<conf>_cutoff<C>.csv
    │   └── final_hits_<conf>_cutoff<C>_part01.sdf, ...
    └── final_hits_<conf>_cutoff<C>_dedup/    step 6 output
        ├── final_hits_<conf>_cutoff<C>_dedup.csv
        └── final_hits_<conf>_cutoff<C>_dedup.sdf   (single file, not split)
```

Every step 3-6 script takes `--vs_results_dir ~/ultra-large/vs_results`
and finds/creates its own named subfolder there. Nothing to configure
beyond that one path.

## Repo-relative: inside the cloned repo

```
post-docking-unidock-enamine/
├── environment.yml, install_pip.sh          the post-dock env
├── nested_bm_clustering/
│   ├── nested_bm_clustering.py              step 7a: reads
│   │                                        ~/ultra-large/vs_results/
│   │                                        final_hits_<conf>_cutoff<C>_dedup/,
│   │                                        writes here:
│   ├── bm_group_summary_<conf>.csv
│   ├── bm_scaffold_assignment_<conf>.csv
│   ├── nested_cluster_assignment_<conf>_<threshold>.csv   (x12 per conf)
│   ├── nested_cluster_summary_<conf>.csv
│   ├── q1/ ... q4/, q_singletons/           step 7b: Q1-Q5, each script's
│   │                                        own CSV + visualization
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


