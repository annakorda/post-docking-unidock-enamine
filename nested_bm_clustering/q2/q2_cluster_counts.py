#!/usr/bin/env python3
"""
Author: Anna Korda

Q2: for each Tc in the Q1-derived range, per conformation, how many
clusters (one number per conformation-Tc pair) -- straight from the
already-written nested_cluster_summary_<conf>.csv, no recompute.
"""
import csv
from pathlib import Path

HERE = Path(__file__).parent
NESTED_DIR = HERE.parent
CONFS = ["c1", "ref1", "c5"]


def load_tc_range():
    with open(NESTED_DIR / "q1" / "tc_range_result.csv") as f:
        rows = list(csv.DictReader(f))
    lo, hi = None, None
    for r in rows:
        if r["conformation"] == "COMBINED_ROUNDED":
            lo, hi = float(r["lower_bound_raw"]), float(r["upper_bound_raw"])
    all_thresholds = [0.90, 0.85, 0.80, 0.75, 0.70, 0.65, 0.60, 0.55, 0.50, 0.45, 0.40, 0.35]
    return [t for t in sorted(all_thresholds) if lo <= t <= hi]


def main():
    tc_range = load_tc_range()
    print(f"Tc range from Q1: {tc_range}\n")

    rows_out = []
    for conf in CONFS:
        with open(NESTED_DIR / f"nested_cluster_summary_{conf}.csv") as f:
            n_clusters_by_t = {float(r["threshold"]): int(r["n_clusters"]) for r in csv.DictReader(f)}
        for t in tc_range:
            rows_out.append({"conformation": conf, "tc": f"{t:.2f}", "n_clusters": n_clusters_by_t[t]})

    out_path = HERE / "q2_cluster_counts.csv"
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["conformation", "tc", "n_clusters"])
        w.writeheader()
        w.writerows(rows_out)
    print(f"Wrote {out_path}\n")

    print(f"{'conf':<6}" + "".join(f"{t:.2f}".rjust(10) for t in tc_range))
    for conf in CONFS:
        vals = [r["n_clusters"] for r in rows_out if r["conformation"] == conf]
        print(f"{conf:<6}" + "".join(f"{v:,}".rjust(10) for v in vals))


if __name__ == "__main__":
    main()
