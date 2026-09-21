# hit_selection

Author: Anna Korda

Downstream of `mmgbsa_rescore/`: narrows the MM-GBSA-scored set down to a
final, ranked pick list. Runs separately per conformation throughout --
c1/ref1/c5 never get merged into one shared list, each keeps its own
survivors and its own final ranking.

## Step 9: `01_cluster_rescore/`

Join `mmgbsa_rescore/results/mmgbsa_<conf>/mmgbsa_results_<conf>.csv`'s
`mmgbsa_dg` onto the existing `cluster_id` column (from the AD4 stage,
already in that same file) -- no re-clustering. Within each cluster
(including size-2), keep only the best MM-GBSA scorer. Singletons pass
through as-is.

Output: `cluster_best_<conf>.csv`.

## Step 10: `02_nitrile_check/`

For nitrile-bearing survivors of step 9: align to PDB 8FYL (vilazodone-
bound 5-HT1A), compare nitrile position/orientation and nearest
contacting residue, judge real contact vs. non-specific.

Put the downloaded 8FYL structure at `02_nitrile_check/structures/8FYL.pdb`.

Output: `nitrile_geometry_<conf>.csv`.

## Step 11: `03_novelty_chembl/`

Known-actives novelty check against CHEMBL214 (5-HT1A): pull Ki/IC50
<=1000 nM actives, dedupe, ECFP4 both sets, max Tanimoto per candidate
against all actives. Flag >=0.35 to nearest known active. Separate
sub-check restricted to nitrile-bearing CHEMBL actives (vilazodone,
analogs) for step 10's nitrile candidates specifically.

Output: `novelty_<conf>.csv` -- full max-Tc distribution per candidate,
not just pass/fail, with borderline cases (near 0.35) flagged separately
from clearly novel ones.

## Step 12: `04_admet/`

SwissADME (BOILED-Egg, CNS MPO) + hERG (admetSAR/pkCSM) on whatever
survives step 11. Removes poor CNS fit or cardiotox liability.

Output: `admet_<conf>.csv`.

## Step 13: `05_final_rank/`

Whatever survives steps 9-12, ranked by MM-GBSA `mmgbsa_dg`. Full ranking
kept, not just the cut -- if a top-30 pick turns out bad on inspection,
the next-ranked one is right there.

Output: `final_ranked_<conf>.csv` (full ranked list, `rank` column,
`in_top30` boolean) -- 30 picks/conformation, 90 total, but the full
list stays available as fallback.

## Step 14: `06_visual_inspection/`

Manual review of the final picks (poses/complexes) in Maestro/PyMOL. No
script -- output/notes land here as they're produced.
