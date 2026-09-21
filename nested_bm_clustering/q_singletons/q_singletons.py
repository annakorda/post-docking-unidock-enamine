#!/usr/bin/env python3
"""
Author: Anna Korda

Singleton compounds (no cluster to validate their score against) get a
plain z-score cutoff instead of a significance test -- there's no second
member to test agreement with, so this isn't "trustworthy," it's just
"among the single best individual scorers." z<=-3.0 (standard 3-sigma
outlier convention) at Tc=0.60, matching the Tc used for Q3/Q4.
"""
import csv
from pathlib import Path

HERE = Path(__file__).parent
NESTED_DIR = HERE.parent
DEDUP_DIR = Path("/home/annie/Desktop/alphavs/revision_2026/ultra-large-results/dedup")
CONFS = ["c1", "ref1", "c5"]
TC = "0.60"
Z_CUTOFF = -3.0


def load_scores_and_zscores(conf):
    path = DEDUP_DIR / f"final_hits_{conf}_cutoff-7_dedup" / f"final_hits_{conf}_cutoff-7_dedup.csv"
    with open(path) as f:
        rows = list(csv.DictReader(f))
    scores = {r["compound_id"]: float(r["ad4_score"]) for r in rows}
    vals = list(scores.values())
    mean = sum(vals) / len(vals)
    std = (sum((v - mean) ** 2 for v in vals) / len(vals)) ** 0.5
    z = {cid: (s - mean) / std for cid, s in scores.items()}
    return scores, z


def main():
    all_rows = []
    for conf in CONFS:
        scores, z = load_scores_and_zscores(conf)
        path = NESTED_DIR / f"nested_cluster_assignment_{conf}_{TC}.csv"
        with open(path) as f:
            assignment = {r["compound_id"]: r["cluster_id"] for r in csv.DictReader(f)}
        by_cluster = {}
        for cid, clid in assignment.items():
            by_cluster.setdefault(clid, []).append(cid)

        singletons = [ids[0] for ids in by_cluster.values() if len(ids) == 1]
        qualifying = [cid for cid in singletons if z[cid] <= Z_CUTOFF]
        qualifying.sort(key=lambda cid: z[cid])

        for cid in qualifying:
            all_rows.append({"conformation": conf, "compound_id": cid,
                              "ad4_score": f"{scores[cid]:.3f}", "z": f"{z[cid]:.4f}"})
        print(f"{conf}: {len(singletons):,} singletons -> {len(qualifying):,} with z<={Z_CUTOFF}")

    out_path = HERE / "qualifying_singletons.csv"
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["conformation", "compound_id", "ad4_score", "z"])
        w.writeheader()
        w.writerows(all_rows)
    print(f"\nWrote {out_path} ({len(all_rows):,} rows)")


if __name__ == "__main__":
    main()
