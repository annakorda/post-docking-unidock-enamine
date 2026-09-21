#!/usr/bin/env python3
"""
Author: Anna Korda

Visualizes Q5: final MM-GBSA input list composition, per conformation --
stacked bar, cluster-derived picks (medoid/best-scorer, deduplicated)
vs. singleton outliers, with the grand total annotated above each bar.
"""
import csv
from pathlib import Path

import matplotlib.pyplot as plt

HERE = Path(__file__).parent
CONFS = ["c1", "ref1", "c5"]
CAT_COLORS = {"cluster-derived": "#2E86AB", "singleton outlier": "#B0B0B0"}


def main():
    with open(HERE / "q5_mmgbsa_input_list.csv") as f:
        rows = list(csv.DictReader(f))

    counts = {conf: {"cluster-derived": 0, "singleton outlier": 0} for conf in CONFS}
    for r in rows:
        cat = "singleton outlier" if r["source"] == "singleton" else "cluster-derived"
        counts[r["conformation"]][cat] += 1

    totals = {conf: sum(counts[conf].values()) for conf in CONFS}
    grand_total = sum(totals.values())

    fig, ax = plt.subplots(figsize=(8, 6))
    x = range(len(CONFS))
    bottoms = [0] * len(CONFS)
    for cat in ("cluster-derived", "singleton outlier"):
        heights = [counts[conf][cat] for conf in CONFS]
        bars = ax.bar(x, heights, bottom=bottoms, color=CAT_COLORS[cat], label=cat,
                      edgecolor="white", linewidth=0.8)
        for i, (h, b) in enumerate(zip(heights, bottoms)):
            if h > 150:
                ax.text(i, b + h / 2, f"{h:,}", ha="center", va="center", fontsize=10,
                        color="white", fontweight="bold")
        bottoms = [b + h for b, h in zip(bottoms, heights)]

    for i, conf in enumerate(CONFS):
        ax.text(i, totals[conf] + max(totals.values()) * 0.02, f"{totals[conf]:,}",
                ha="center", fontsize=12, fontweight="bold")

    ax.set_xticks(list(x))
    ax.set_xticklabels(CONFS, fontsize=13, fontweight="bold")
    ax.set_ylabel("compounds in final MM-GBSA list", fontsize=12, fontweight="bold")
    ax.set_ylim(0, max(totals.values()) * 1.15)
    ax.set_title(f"Q5: final MM-GBSA input list (grand total = {grand_total:,})",
                 fontsize=13, fontweight="bold", pad=14)
    ax.grid(axis="y", alpha=0.3, linestyle="--", linewidth=0.5)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(frameon=True, fontsize=10, edgecolor="black", loc="upper right")
    plt.tight_layout()
    for ext in ("png", "svg"):
        kw = dict(dpi=200) if ext == "png" else {}
        plt.savefig(HERE / f"q5_visualization.{ext}", bbox_inches="tight", facecolor="white", **kw)
    plt.close(fig)
    print("Wrote q5_visualization.png / .svg")


if __name__ == "__main__":
    main()
