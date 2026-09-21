#!/usr/bin/env python3
"""
Author: Anna Korda

Compares the old flat Tanimoto/Butina pass (ecfp4_clustering_dedup/) against
the new two-pass Bemis-Murcko-then-Tanimoto approach (nested_bm_clustering/),
both run on the SAME deduplicated input so the comparison isolates the
algorithm change, not the dedup step. One subplot per conformation: total
cluster count vs. Tanimoto threshold, old (solid) vs. new (dashed).
"""
import csv
from pathlib import Path

import matplotlib.pyplot as plt

HERE = Path(__file__).parent
FLAT_DIR = HERE.parent / "ecfp4_clustering_dedup"
COLORS = {"c1": "#3D7A8C", "ref1": "#C97B3D", "c5": "#7B5EA7"}
CONFS = ["c1", "ref1", "c5"]


def load_summary(path_dir, prefix, conf):
    with open(path_dir / f"{prefix}_{conf}.csv") as f:
        rows = list(csv.DictReader(f))
    rows.sort(key=lambda r: float(r["threshold"]))
    thresholds = [float(r["threshold"]) for r in rows]
    n_clusters = [int(r["n_clusters"]) for r in rows]
    largest = [int(r["largest_cluster_size"]) for r in rows]
    n_compounds = int(rows[0]["n_compounds"])
    return thresholds, n_clusters, largest, n_compounds


def main():
    fig, axes = plt.subplots(1, 3, figsize=(16, 5), sharey=False)
    for ax, conf in zip(axes, CONFS):
        t_flat, n_flat, largest_flat, n_flat_compounds = load_summary(FLAT_DIR, "cluster_summary", conf)
        t_nested, n_nested, largest_nested, n_nested_compounds = load_summary(HERE, "nested_cluster_summary", conf)
        # thresholds is ascending after load_summary's own sort, so index 0 is the
        # loosest (minimum) threshold actually tested -- not assumed to be exactly
        # 0.35, in case --thresholds ever changes
        min_t_flat, min_t_nested = t_flat[0], t_nested[0]

        ax.plot(t_flat, n_flat, marker="o", markersize=6, linewidth=2, linestyle="-",
                color=COLORS[conf], label=f"flat (largest@{min_t_flat:.2f}={largest_flat[0]:,})")
        ax.plot(t_nested, n_nested, marker="s", markersize=6, linewidth=2, linestyle="--",
                color=COLORS[conf], alpha=0.7, label=f"BM+Tanimoto (largest@{min_t_nested:.2f}={largest_nested[0]:,})")

        ax.set_xlabel("ECFP4 Tanimoto similarity threshold", fontsize=11, fontweight="bold")
        if conf == CONFS[0]:
            ax.set_ylabel("total cluster count", fontsize=11, fontweight="bold")
        ax.set_title(f"{conf} (n={n_flat_compounds:,}, deduplicated)", fontsize=12, fontweight="bold")
        ax.invert_xaxis()
        ax.grid(alpha=0.3, linestyle="--", linewidth=0.5)
        ax.set_axisbelow(True)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.legend(frameon=True, fontsize=9, edgecolor="black", loc="upper right")

    fig.suptitle("Flat Tanimoto/Butina vs. Bemis-Murcko + Tanimoto/Butina (same deduplicated input)",
                 fontsize=14, fontweight="bold", y=1.02)
    plt.tight_layout()
    plt.savefig(HERE / "nested_vs_flat_cluster_count.svg", bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("Wrote nested_vs_flat_cluster_count.svg")


if __name__ == "__main__":
    main()
