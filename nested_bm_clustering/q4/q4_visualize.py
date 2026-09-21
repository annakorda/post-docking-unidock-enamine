#!/usr/bin/env python3
"""
Author: Anna Korda

Visualizes Q4: among trustworthy clusters, cluster size (x, log) vs
mean_z (y) at Tc=0.60, with the significance boundary curve overlaid
(mean_z = -1.645/sqrt(count), the p=0.05 one-tailed critical line for
t_stat = mean_z*sqrt(count)) -- makes visible why "good" needs a more
extreme mean_z for small clusters than for large ones (this is the same
count-dependence Q3 deliberately avoided, but it's the CORRECT behavior
here, since Q4 answers an estimation question, not a descriptive one).
"""
import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).parent
CONFS = ["c1", "ref1", "c5"]
COLORS = {"c1": "#3D7A8C", "ref1": "#C97B3D", "c5": "#7B5EA7"}
TC = "0.60"
Z_CRIT = -1.645  # one-tailed p=0.05


def main():
    with open(HERE / "q4_per_cluster.csv") as f:
        rows = [r for r in csv.DictReader(f) if r["tc"] == TC]
    for r in rows:
        r["count"] = int(r["count"])
        r["mean_z"] = float(r["mean_z"])
        r["good"] = r["good"] == "True"

    fig, axes = plt.subplots(1, 3, figsize=(16, 5.5), sharey=True)
    counts_range = np.linspace(2, 200, 300)
    boundary = Z_CRIT / np.sqrt(counts_range)

    for ax, conf in zip(axes, CONFS):
        sub = [r for r in rows if r["conformation"] == conf]
        good = [r for r in sub if r["good"]]
        not_good = [r for r in sub if not r["good"]]

        ax.scatter([r["count"] for r in not_good], [r["mean_z"] for r in not_good],
                   s=6, alpha=0.25, color="#B0B0B0", label=f"trustworthy only (n={len(not_good):,})")
        ax.scatter([r["count"] for r in good], [r["mean_z"] for r in good],
                   s=20, alpha=0.85, color=COLORS[conf], edgecolors="black", linewidths=0.3,
                   label=f"trustworthy AND good (n={len(good):,})")
        ax.plot(counts_range, boundary, color="black", linewidth=1.3, linestyle="--",
                label="p=0.05 boundary (mean_z=-1.645/√n)")

        ax.set_xscale("log")
        ax.set_xlabel("cluster size (compounds)", fontsize=11, fontweight="bold")
        if conf == CONFS[0]:
            ax.set_ylabel("mean AD4 z-score", fontsize=11, fontweight="bold")
        ax.axhline(0, color="black", linewidth=0.7, linestyle=":", alpha=0.5)
        ax.set_title(f"{conf}  (Tc={TC})", fontsize=13, fontweight="bold", color=COLORS[conf])
        ax.grid(alpha=0.25, linestyle="--", linewidth=0.5)
        ax.set_axisbelow(True)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.legend(frameon=True, fontsize=8, edgecolor="black", loc="lower left")

    fig.suptitle("Q4: good = significantly better than pool (raw p<0.05) among the trustworthy set",
                 fontsize=14, fontweight="bold", y=1.03)
    plt.tight_layout()
    for ext in ("png", "svg"):
        kw = dict(dpi=200) if ext == "png" else {}
        plt.savefig(HERE / f"q4_visualization.{ext}", bbox_inches="tight", facecolor="white", **kw)
    plt.close(fig)
    print("Wrote q4_visualization.png / .svg")


if __name__ == "__main__":
    main()
