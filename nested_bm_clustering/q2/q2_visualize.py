#!/usr/bin/env python3
"""
Author: Anna Korda

Visualizes Q2: cluster count per Tc, all 3 conformations overlaid.
"""
import csv
from pathlib import Path

import matplotlib.pyplot as plt

HERE = Path(__file__).parent
CONFS = ["c1", "ref1", "c5"]
COLORS = {"c1": "#3D7A8C", "ref1": "#C97B3D", "c5": "#7B5EA7"}


def main():
    with open(HERE / "q2_cluster_counts.csv") as f:
        rows = list(csv.DictReader(f))

    fig, ax = plt.subplots(figsize=(8, 5.5))
    for conf in CONFS:
        sub = [r for r in rows if r["conformation"] == conf]
        sub.sort(key=lambda r: float(r["tc"]))
        ts = [float(r["tc"]) for r in sub]
        ns = [int(r["n_clusters"]) for r in sub]
        ax.plot(ts, ns, marker="o", markersize=7, linewidth=2.2, color=COLORS[conf], label=conf)
        ax.annotate(f"{ns[-1]:,}", (ts[-1], ns[-1]), textcoords="offset points", xytext=(6, 4),
                    fontsize=8, color=COLORS[conf])

    ax.set_xlabel("ECFP4 Tanimoto threshold (Tc)", fontsize=12, fontweight="bold")
    ax.set_ylabel("total cluster count", fontsize=12, fontweight="bold")
    ax.set_title("Q2: cluster count across the Q1-derived Tc range", fontsize=13, fontweight="bold", pad=12)
    ax.invert_xaxis()
    ax.grid(alpha=0.3, linestyle="--", linewidth=0.5)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(frameon=True, fontsize=10, edgecolor="black")
    plt.tight_layout()
    for ext in ("png", "svg"):
        kw = dict(dpi=200) if ext == "png" else {}
        plt.savefig(HERE / f"q2_visualization.{ext}", bbox_inches="tight", facecolor="white", **kw)
    plt.close(fig)
    print("Wrote q2_visualization.png / .svg")


if __name__ == "__main__":
    main()
