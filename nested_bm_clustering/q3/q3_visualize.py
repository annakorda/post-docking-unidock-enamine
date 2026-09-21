#!/usr/bin/env python3
"""
Author: Anna Korda

Visualizes Q3: distribution of std_z (internal spread) at Tc=0.60, one
panel per conformation, cutoff (0.5) marked -- trustworthy (std_z<0.5)
shaded differently from not.
"""
import csv
from pathlib import Path

import matplotlib.pyplot as plt

HERE = Path(__file__).parent
CONFS = ["c1", "ref1", "c5"]
COLORS = {"c1": "#3D7A8C", "ref1": "#C97B3D", "c5": "#7B5EA7"}
TC = "0.60"
STD_CUTOFF = 1.0


def main():
    with open(HERE / "q3_per_cluster.csv") as f:
        rows = [r for r in csv.DictReader(f) if r["tc"] == TC]
    for r in rows:
        r["std_z"] = float(r["std_z"])

    fig, axes = plt.subplots(1, 3, figsize=(16, 5), sharey=False)
    for ax, conf in zip(axes, CONFS):
        vals = [r["std_z"] for r in rows if r["conformation"] == conf]
        n_trust = sum(1 for v in vals if v < STD_CUTOFF)
        bins = [i * 0.05 for i in range(0, 41)]
        ax.hist([v for v in vals if v < STD_CUTOFF], bins=bins, color=COLORS[conf], alpha=0.85,
                label=f"trustworthy (n={n_trust:,})")
        ax.hist([v for v in vals if v >= STD_CUTOFF], bins=bins, color="#B0B0B0", alpha=0.6,
                label=f"not trustworthy (n={len(vals)-n_trust:,})")
        ax.axvline(STD_CUTOFF, color="black", linewidth=1.5, linestyle="--")
        ax.text(STD_CUTOFF, ax.get_ylim()[1] if ax.get_ylim()[1] else 1, f" cutoff={STD_CUTOFF}",
                fontsize=9, va="top")
        ax.set_xlabel("std_z (internal spread, sample std of member z-scores)", fontsize=10, fontweight="bold")
        if conf == CONFS[0]:
            ax.set_ylabel("number of clusters", fontsize=11, fontweight="bold")
        ax.set_title(f"{conf}  (Tc={TC})", fontsize=13, fontweight="bold", color=COLORS[conf])
        ax.grid(alpha=0.25, linestyle="--", linewidth=0.5)
        ax.set_axisbelow(True)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.legend(frameon=True, fontsize=9, edgecolor="black")

    fig.suptitle("Q3: trustworthy = std_z < 0.5 (descriptive, size-independent)",
                 fontsize=14, fontweight="bold", y=1.03)
    plt.tight_layout()
    for ext in ("png", "svg"):
        kw = dict(dpi=200) if ext == "png" else {}
        plt.savefig(HERE / f"q3_visualization.{ext}", bbox_inches="tight", facecolor="white", **kw)
    plt.close(fig)
    print("Wrote q3_visualization.png / .svg")


if __name__ == "__main__":
    main()
