#!/usr/bin/env python3
"""
Author: Anna Korda

Visualizes the singleton z-cutoff choice: qualifying-singleton count vs.
z threshold, per conformation, showing how fast the count collapses as
the bar tightens, with the chosen cutoff (z<=-3.0) marked. This is the
real basis for -3.0 -- not a round number, the point where singleton
count drops to a scale comparable to the cluster tier instead of
swamping it.
"""
import csv
from pathlib import Path

import matplotlib.pyplot as plt

HERE = Path(__file__).parent
DEDUP_DIR = Path("/home/annie/Desktop/alphavs/revision_2026/ultra-large-results/dedup")
NESTED_DIR = HERE.parent
CONFS = ["c1", "ref1", "c5"]
COLORS = {"c1": "#3D7A8C", "ref1": "#C97B3D", "c5": "#7B5EA7"}
TC = "0.60"
CUTOFFS = [-1.645, -1.96, -2.326, -2.576, -3.0, -3.29, -3.5, -3.72, -4.0]
CHOSEN = -3.0


def load_zscores(conf):
    path = DEDUP_DIR / f"final_hits_{conf}_cutoff-7_dedup" / f"final_hits_{conf}_cutoff-7_dedup.csv"
    with open(path) as f:
        rows = list(csv.DictReader(f))
    scores = {r["compound_id"]: float(r["ad4_score"]) for r in rows}
    vals = list(scores.values())
    mean = sum(vals) / len(vals)
    std = (sum((v - mean) ** 2 for v in vals) / len(vals)) ** 0.5
    return {cid: (s - mean) / std for cid, s in scores.items()}


def main():
    fig, ax = plt.subplots(figsize=(9, 6))
    cluster_totals = {"c1": 1735, "ref1": 959, "c5": 1053}  # trustworthy+good cluster counts, Tc=0.60

    for conf in CONFS:
        z = load_zscores(conf)
        path = NESTED_DIR / f"nested_cluster_assignment_{conf}_{TC}.csv"
        with open(path) as f:
            assignment = {r["compound_id"]: r["cluster_id"] for r in csv.DictReader(f)}
        by_cluster = {}
        for cid, clid in assignment.items():
            by_cluster.setdefault(clid, []).append(cid)
        singleton_z = [z[ids[0]] for ids in by_cluster.values() if len(ids) == 1]

        counts = [sum(1 for s in singleton_z if s <= c) for c in CUTOFFS]
        ax.plot(CUTOFFS, counts, marker="o", markersize=7, linewidth=2.2, color=COLORS[conf], label=conf)
        chosen_count = sum(1 for s in singleton_z if s <= CHOSEN)
        ax.annotate(f"{chosen_count}", (CHOSEN, chosen_count), textcoords="offset points",
                    xytext=(8, 4), fontsize=9, color=COLORS[conf], fontweight="bold")
        ax.axhline(cluster_totals[conf], color=COLORS[conf], linewidth=1, linestyle=":", alpha=0.5)

    ax.axvline(CHOSEN, color="black", linewidth=1.5, linestyle="--", label=f"chosen cutoff ({CHOSEN})")
    ax.set_yscale("log")
    ax.set_xlabel("z-score cutoff (singleton qualifies if z ≤ cutoff)", fontsize=12, fontweight="bold")
    ax.set_ylabel("qualifying singleton count (log scale)", fontsize=12, fontweight="bold")
    ax.set_title("Singleton z-cutoff: how fast the pool shrinks as the bar tightens\n"
                 "(dotted lines = that conformation's trustworthy+good cluster count, for scale)",
                 fontsize=12, fontweight="bold", pad=12)
    ax.invert_xaxis()
    ax.grid(alpha=0.3, linestyle="--", linewidth=0.5)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(frameon=True, fontsize=9, edgecolor="black")
    plt.tight_layout()
    for ext in ("png", "svg"):
        kw = dict(dpi=200) if ext == "png" else {}
        plt.savefig(HERE / f"singleton_zcutoff_choice.{ext}", bbox_inches="tight", facecolor="white", **kw)
    plt.close(fig)
    print("Wrote singleton_zcutoff_choice.png / .svg")


if __name__ == "__main__":
    main()
