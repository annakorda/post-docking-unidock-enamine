#!/usr/bin/env python3
"""
Author: Anna Korda

Visualizes how Q1 derived the Tc range: cluster-count curve (left axis)
and top-10-%-of-pool curve (right axis) vs Tc, one panel per conformation,
with the derived lower/upper bounds marked as vertical lines and the
accepted range shaded.
"""
import csv
from pathlib import Path

import matplotlib.pyplot as plt

HERE = Path(__file__).parent
NESTED_DIR = HERE.parent
CONFS = ["c1", "ref1", "c5"]
COLORS = {"c1": "#3D7A8C", "ref1": "#C97B3D", "c5": "#7B5EA7"}
THRESHOLDS = [0.90, 0.85, 0.80, 0.75, 0.70, 0.65, 0.60, 0.55, 0.50, 0.45, 0.40, 0.35]


def load_n_clusters(conf):
    with open(NESTED_DIR / f"nested_cluster_summary_{conf}.csv") as f:
        rows = {float(r["threshold"]): int(r["n_clusters"]) for r in csv.DictReader(f)}
    return [rows[t] for t in THRESHOLDS]


def load_top10_pct(conf):
    from collections import Counter
    vals = []
    for t in THRESHOLDS:
        path = NESTED_DIR / f"nested_cluster_assignment_{conf}_{t:.2f}.csv"
        with open(path) as f:
            cluster_ids = [row["cluster_id"] for row in csv.DictReader(f)]
        counts = Counter(cluster_ids)
        top10_sum = sum(sorted(counts.values(), reverse=True)[:10])
        vals.append(100 * top10_sum / len(cluster_ids))
    return vals


def load_bounds():
    with open(HERE / "tc_range_result.csv") as f:
        for r in csv.DictReader(f):
            if r["conformation"] == "COMBINED_ROUNDED":
                return float(r["lower_bound_raw"]), float(r["upper_bound_raw"])


def main():
    lo, hi = load_bounds()
    fig, axes = plt.subplots(1, 3, figsize=(17, 5.5))

    for ax, conf in zip(axes, CONFS):
        n_clusters = load_n_clusters(conf)
        top10 = load_top10_pct(conf)
        color = COLORS[conf]

        l1, = ax.plot(THRESHOLDS, n_clusters, marker="o", markersize=6, linewidth=2.2,
                       color=color, label="total cluster count")
        ax.set_xlabel("ECFP4 Tanimoto threshold (Tc)", fontsize=11, fontweight="bold")
        ax.set_ylabel("total cluster count", fontsize=11, fontweight="bold", color=color)
        ax.tick_params(axis="y", labelcolor=color)
        ax.invert_xaxis()

        ax2 = ax.twinx()
        l2, = ax2.plot(THRESHOLDS, top10, marker="s", markersize=6, linewidth=2.2, linestyle="--",
                        color="#555555", alpha=0.85, label="top-10 clusters % of pool")
        ax2.set_ylabel("top-10 clusters % of pool", fontsize=11, fontweight="bold", color="#555555")
        ax2.tick_params(axis="y", labelcolor="#555555")

        ax.axvspan(lo, hi, color=color, alpha=0.10, zorder=0)
        ax.axvline(hi, color=color, linewidth=1.3, linestyle=":", alpha=0.9)
        ax.axvline(lo, color=color, linewidth=1.3, linestyle=":", alpha=0.9)
        ax.text(hi, ax.get_ylim()[1]*0.97, f" upper={hi:.2f}\n (slope->half-max)", fontsize=8,
                va="top", ha="left" if hi < (THRESHOLDS[0]+THRESHOLDS[-1])/2 else "right")
        ax.text(lo, ax.get_ylim()[1]*0.97, f" lower={lo:.2f}\n (climb-rate local min)", fontsize=8,
                va="top", ha="right")

        ax.set_title(f"{conf}", fontsize=13, fontweight="bold", color=color)
        ax.grid(alpha=0.25, linestyle="--", linewidth=0.5)
        ax.set_axisbelow(True)
        if conf == CONFS[0]:
            ax.legend(handles=[l1, l2], loc="lower left", fontsize=9, frameon=True, edgecolor="black")

    fig.suptitle("Q1: Tc range derived from curve shape alone (shaded = accepted range)",
                 fontsize=14, fontweight="bold", y=1.03)
    plt.tight_layout()
    for ext in ("png", "svg"):
        kw = dict(dpi=200) if ext == "png" else {}
        plt.savefig(HERE / f"q1_visualization.{ext}", bbox_inches="tight", facecolor="white", **kw)
    plt.close(fig)
    print(f"Wrote q1_visualization.png / .svg")


if __name__ == "__main__":
    main()
