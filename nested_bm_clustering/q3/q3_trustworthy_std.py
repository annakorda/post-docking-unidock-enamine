#!/usr/bin/env python3
"""
Author: Anna Korda

Q3, final: trustworthy = std_z < STD_CUTOFF (sample std, ddof=1, of a
cluster's own member z-scores; STD_CUTOFF=1.0, ~0.6-0.7 kcal/mol real AD4
spread). Deliberately NOT a significance test -- a chi-squared variance
test was tried and rejected: verified on real data that it favors large
clusters for the same reason a mean-based significance test does (more
members = more power to prove a given spread is "significantly" low),
regardless of the cluster's actual spread. Confirmed at an earlier
candidate cutoff (0.5): 10-15% of chi2-passing clusters had std_z>=0.5
(loose by eye, passed on n alone, counts 7-101) while ~70% of genuinely
tight (std_z<0.3) 2-3-member clusters failed the test outright, purely
from lack of power.

Trustworthiness is a descriptive property (how much do these members
actually disagree), not an estimation problem -- unlike "good" (Q4), where
more data genuinely earning more confidence in a mean IS the right idea,
so a significance test belongs there instead.

Singleton clusters (count=1) have no spread to measure -- excluded here,
handled separately (see q_singletons/).
"""
import csv
from pathlib import Path

HERE = Path(__file__).parent
NESTED_DIR = HERE.parent
DEDUP_DIR = NESTED_DIR.parent / "vs_results" / "dedup"  # repo_root/vs_results/dedup, dedup_stereoisomers.py's --out_dir
CONFS = ["c1", "ref1", "c5"]
STD_CUTOFF = 1.0  # ~0.6-0.7 kcal/mol, matched to AD4's commonly-cited noise floor


def load_tc_range():
    with open(NESTED_DIR / "q1" / "tc_range_result.csv") as f:
        rows = list(csv.DictReader(f))
    for r in rows:
        if r["conformation"] == "COMBINED_ROUNDED":
            lo, hi = float(r["lower_bound_raw"]), float(r["upper_bound_raw"])
    all_thresholds = [0.90, 0.85, 0.80, 0.75, 0.70, 0.65, 0.60, 0.55, 0.50, 0.45, 0.40, 0.35]
    return [t for t in sorted(all_thresholds) if lo <= t <= hi]


def load_zscores(conf):
    path = DEDUP_DIR / f"final_hits_{conf}_cutoff-7_dedup" / f"final_hits_{conf}_cutoff-7_dedup.csv"
    with open(path) as f:
        rows = list(csv.DictReader(f))
    scores = {r["compound_id"]: float(r["ad4_score"]) for r in rows}
    vals = list(scores.values())
    mean = sum(vals) / len(vals)
    std = (sum((v - mean) ** 2 for v in vals) / len(vals)) ** 0.5
    return {cid: (s - mean) / std for cid, s in scores.items()}


def sample_std(vals):
    n = len(vals)
    if n < 2:
        return 0.0
    m = sum(vals) / n
    return (sum((v - m) ** 2 for v in vals) / (n - 1)) ** 0.5


def main():
    tc_range = load_tc_range()
    per_cluster_rows = []
    summary_rows = []

    for conf in CONFS:
        z = load_zscores(conf)
        for t in tc_range:
            path = NESTED_DIR / f"nested_cluster_assignment_{conf}_{t:.2f}.csv"
            with open(path) as f:
                assignment = {r["compound_id"]: r["cluster_id"] for r in csv.DictReader(f)}
            by_cluster = {}
            for cid, cluster_id in assignment.items():
                by_cluster.setdefault(cluster_id, []).append(z[cid])

            n_singleton, n_testable, n_trustworthy = 0, 0, 0
            for cluster_id, zs in by_cluster.items():
                count = len(zs)
                if count == 1:
                    n_singleton += 1
                    continue
                n_testable += 1
                mean_z = sum(zs) / count
                min_z = min(zs)
                std_z = sample_std(zs)
                trustworthy = std_z < STD_CUTOFF
                if trustworthy:
                    n_trustworthy += 1
                per_cluster_rows.append({
                    "conformation": conf, "tc": f"{t:.2f}", "cluster_id": cluster_id, "count": count,
                    "mean_z": f"{mean_z:.4f}", "min_z": f"{min_z:.4f}", "std_z": f"{std_z:.4f}",
                    "trustworthy": trustworthy,
                })

            summary_rows.append({
                "conformation": conf, "tc": f"{t:.2f}", "n_clusters_total": len(by_cluster),
                "n_singleton": n_singleton, "n_testable": n_testable, "n_trustworthy": n_trustworthy,
                "pct_trustworthy_of_testable": f"{100 * n_trustworthy / max(1,n_testable):.2f}",
            })
            print(f"{conf} Tc={t:.2f}: {len(by_cluster):,} clusters ({n_singleton:,} singleton, "
                  f"{n_testable:,} testable) -> {n_trustworthy:,} trustworthy "
                  f"({100 * n_trustworthy / max(1,n_testable):.1f}%)")

    per_cluster_path = HERE / "q3_per_cluster.csv"
    with open(per_cluster_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["conformation", "tc", "cluster_id", "count", "mean_z", "min_z",
                                           "std_z", "trustworthy"])
        w.writeheader()
        w.writerows(per_cluster_rows)
    print(f"\nWrote {per_cluster_path} ({len(per_cluster_rows):,} rows, non-singleton clusters only)")

    summary_path = HERE / "q3_summary.csv"
    with open(summary_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["conformation", "tc", "n_clusters_total", "n_singleton",
                                           "n_testable", "n_trustworthy", "pct_trustworthy_of_testable"])
        w.writeheader()
        w.writerows(summary_rows)
    print(f"Wrote {summary_path}")


if __name__ == "__main__":
    main()
