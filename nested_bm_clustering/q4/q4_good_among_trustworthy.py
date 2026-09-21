#!/usr/bin/env python3
"""
Author: Anna Korda

Q4, final: of Q3's trustworthy (std_z<0.5) clusters, which are also
significantly better than the pool? One-sided one-sample z-test on
mean_z, RAW p<0.05 -- no BH correction, deliberately, for consistency
with Q3 (both use the same threshold philosophy now).

Requiring a cluster to pass two independent p<0.05-equivalent criteria
at once (Q3's std cutoff acts as a hard descriptive gate, Q4's z-test is
the one real significance test) is not formally equivalent to an FDR
correction -- Q3 isn't a p-value, so there's no "joint alpha" to compute.
The actual justification is simpler: Q3 removes clusters whose average
isn't representative of the whole cluster (large internal disagreement),
and Q4 only asks the significance question of clusters that already
passed that bar -- so the significance test is answering "is the mean
of an already-validated-as-representative sample better than the pool,"
which is the question it was built for.
"""
import csv
import math
from pathlib import Path

HERE = Path(__file__).parent
Q3_PATH = HERE.parent / "q3" / "q3_per_cluster.csv"
CONFS = ["c1", "ref1", "c5"]
ALPHA = 0.05


def norm_cdf(x):
    return 0.5 * (1 + math.erf(x / math.sqrt(2)))


def main():
    with open(Q3_PATH) as f:
        rows = [r for r in csv.DictReader(f) if r["trustworthy"] == "True"]
    for r in rows:
        r["count"] = int(r["count"])
        r["mean_z"] = float(r["mean_z"])

    tcs = sorted({r["tc"] for r in rows}, key=float)
    out_rows = []
    summary_rows = []

    for conf in CONFS:
        for tc in tcs:
            sub = [r for r in rows if r["conformation"] == conf and r["tc"] == tc]
            n_good = 0
            for r in sub:
                t = r["mean_z"] * math.sqrt(r["count"])
                p = norm_cdf(t)
                good = p <= ALPHA
                if good:
                    n_good += 1
                out_rows.append({
                    "conformation": conf, "tc": tc, "cluster_id": r["cluster_id"], "count": r["count"],
                    "mean_z": f"{r['mean_z']:.4f}", "min_z": r["min_z"], "std_z": r["std_z"],
                    "t_stat": f"{t:.4f}", "p_value": f"{p:.6g}", "good": good,
                })
            summary_rows.append({"conformation": conf, "tc": tc, "n_trustworthy": len(sub),
                                  "n_good": n_good, "pct_good": f"{100*n_good/max(1,len(sub)):.1f}"})
            print(f"{conf} Tc={tc}: {len(sub):,} trustworthy -> {n_good:,} also good "
                  f"({100*n_good/max(1,len(sub)):.1f}%)")

    per_cluster_path = HERE / "q4_per_cluster.csv"
    with open(per_cluster_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["conformation", "tc", "cluster_id", "count", "mean_z", "min_z",
                                           "std_z", "t_stat", "p_value", "good"])
        w.writeheader()
        w.writerows(out_rows)
    print(f"\nWrote {per_cluster_path} ({len(out_rows):,} rows)")

    summary_path = HERE / "q4_summary.csv"
    with open(summary_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["conformation", "tc", "n_trustworthy", "n_good", "pct_good"])
        w.writeheader()
        w.writerows(summary_rows)
    print(f"Wrote {summary_path}")


if __name__ == "__main__":
    main()
