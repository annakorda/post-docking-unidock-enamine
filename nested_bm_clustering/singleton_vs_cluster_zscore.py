#!/usr/bin/env python3
"""
Author: Anna Korda

Dot plot: cluster size (x, log scale) vs. mean AD4 z-score (y), one panel
per conformation, at a single Tc. Every cluster is a dot -- singletons and
multi-member clusters are different colors, and among multi-member
clusters, Q3-trustworthy and Q4-good are further highlighted, so it's
visible at a glance whether singletons cluster at extreme z (a real
selection-bias risk: a single lucky compound looks "great" with nothing
to average it against) compared to real multi-member clusters.

Tc=0.60 used -- the peak of the flat trustworthy/good stretch in the Q3/Q4
tables. Rerun with a different TC below if a different threshold is wanted.
"""
import csv
from pathlib import Path

import matplotlib.pyplot as plt

HERE = Path(__file__).parent
DEDUP_DIR = Path("/home/annie/Desktop/alphavs/revision_2026/ultra-large-results/dedup")
CONFS = ["c1", "ref1", "c5"]
TC = "0.60"

CATEGORY_COLORS = {
    "singleton": "#B0B0B0",
    "multi, not trustworthy": "#F2C14E",
    "trustworthy, not good (tier2)": "#3D7A8C",
    "trustworthy AND good (tier1)": "#C0392B",
}
CATEGORY_ORDER = ["singleton", "multi, not trustworthy", "trustworthy, not good (tier2)", "trustworthy AND good (tier1)"]


def load_zscores(conf):
    path = DEDUP_DIR / f"final_hits_{conf}_cutoff-7_dedup" / f"final_hits_{conf}_cutoff-7_dedup.csv"
    with open(path) as f:
        rows = list(csv.DictReader(f))
    scores = {r["compound_id"]: float(r["ad4_score"]) for r in rows}
    vals = list(scores.values())
    mean = sum(vals) / len(vals)
    std = (sum((v - mean) ** 2 for v in vals) / len(vals)) ** 0.5
    return {cid: (s - mean) / std for cid, s in scores.items()}, mean, std


def load_flags(conf):
    """cluster_id -> 'trustworthy' or 'good' flag, from the already-computed Q3/Q4 tables."""
    trustworthy_ids, good_ids = set(), set()
    with open(HERE / "q3" / "q3_per_cluster.csv") as f:
        for r in csv.DictReader(f):
            if r["conformation"] == conf and r["tc"] == TC and r["trustworthy"] == "True":
                trustworthy_ids.add(r["cluster_id"])
    with open(HERE / "q4" / "q4_per_cluster.csv") as f:
        for r in csv.DictReader(f):
            if r["conformation"] == conf and r["tc"] == TC and r["good"] == "True":
                good_ids.add(r["cluster_id"])
    return trustworthy_ids, good_ids


def main():
    fig, axes = plt.subplots(1, 3, figsize=(17, 5.5), sharey=True)

    for ax, conf in zip(axes, CONFS):
        z, pool_mean, pool_std = load_zscores(conf)
        trustworthy_ids, good_ids = load_flags(conf)

        path = HERE / f"nested_cluster_assignment_{conf}_{TC}.csv"
        with open(path) as f:
            assignment = {r["compound_id"]: r["cluster_id"] for r in csv.DictReader(f)}
        by_cluster = {}
        for cid, cluster_id in assignment.items():
            by_cluster.setdefault(cluster_id, []).append(z[cid])

        points = {cat: {"x": [], "y": []} for cat in CATEGORY_ORDER}
        for cluster_id, zs in by_cluster.items():
            count = len(zs)
            mean_z = sum(zs) / count
            if count == 1:
                cat = "singleton"
            elif cluster_id in good_ids:
                cat = "trustworthy AND good (tier1)"
            elif cluster_id in trustworthy_ids:
                cat = "trustworthy, not good (tier2)"
            else:
                cat = "multi, not trustworthy"
            points[cat]["x"].append(count)
            points[cat]["y"].append(mean_z)

        # dashed lines at the average z of each trustworthy tier -- makes the
        # "opposite sides of the pool mean" split visible, not just implied by dots
        for cat in ("trustworthy, not good (tier2)", "trustworthy AND good (tier1)"):
            ys = points[cat]["y"]
            if ys:
                avg = sum(ys) / len(ys)
                ax.axhline(avg, color=CATEGORY_COLORS[cat], linewidth=1.3, linestyle="--", alpha=0.8)

        # draw least-prominent categories first so the two trustworthy tiers sit on top
        prominent = ("trustworthy, not good (tier2)", "trustworthy AND good (tier1)")
        for cat in CATEGORY_ORDER:
            p = points[cat]
            if not p["x"]:
                continue
            size = 22 if cat in prominent else 8
            alpha = 0.85 if cat in prominent else 0.35
            ax.scatter(p["x"], p["y"], s=size, alpha=alpha, color=CATEGORY_COLORS[cat],
                       edgecolors="black" if cat == "trustworthy AND good (tier1)" else "none", linewidths=0.4,
                       label=f"{cat} (n={len(p['x']):,})")

        ax.set_xscale("log")
        ax.set_xlabel("cluster size (compounds)", fontsize=11, fontweight="bold")
        if conf == CONFS[0]:
            ax.set_ylabel("mean AD4 z-score (pool-normalized)", fontsize=11, fontweight="bold")
        ax.axhline(0, color="black", linewidth=0.8, linestyle=":", alpha=0.6)
        ax.set_title(f"{conf}  (Tc={TC})", fontsize=13, fontweight="bold")
        ax.grid(alpha=0.25, linestyle="--", linewidth=0.5)
        ax.set_axisbelow(True)
        ax.spines["top"].set_visible(False)
        ax.legend(frameon=True, fontsize=8, edgecolor="black", loc="lower right")

        # real AD4 score on the right axis -- z and AD4 are a linear transform of
        # each other per conformation (z = (ad4-mean)/std), pool mean/std computed
        # fresh per conformation, not shared across panels
        def z_to_ad4(zz, m=pool_mean, s=pool_std):
            return zz * s + m

        def ad4_to_z(aa, m=pool_mean, s=pool_std):
            return (aa - m) / s

        secax = ax.secondary_yaxis("right", functions=(z_to_ad4, ad4_to_z))
        secax.set_ylabel("mean AD4 score (kcal/mol)", fontsize=11, fontweight="bold")

    fig.suptitle("Singletons vs. multi-member clusters: cluster size vs. mean z-score",
                 fontsize=14, fontweight="bold", y=1.02)
    plt.tight_layout()
    png_path = HERE / "singleton_vs_cluster_zscore.png"
    svg_path = HERE / "singleton_vs_cluster_zscore.svg"
    plt.savefig(png_path, bbox_inches="tight", facecolor="white", dpi=200)
    plt.savefig(svg_path, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Wrote {png_path}")
    print(f"Wrote {svg_path}")


if __name__ == "__main__":
    main()
