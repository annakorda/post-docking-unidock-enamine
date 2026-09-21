#!/usr/bin/env python3
"""
Author: Anna Korda

Per-conformation trustworthy/good/tier2 table (Q3 = trustworthy, internal
consistency chi-squared test; Q4 = good, significantly better than pool,
tested only on Q3's survivors) across the full Q1 Tc range -- one table
per conformation, not summed, plus a rendered image (3 panels).
"""
import csv
from pathlib import Path

import matplotlib.pyplot as plt

HERE = Path(__file__).parent
CONFS = ["c1", "ref1", "c5"]
COLORS = {"c1": "#3D7A8C", "ref1": "#C97B3D", "c5": "#7B5EA7"}


def main():
    with open(HERE / "q3" / "q3_summary.csv") as f:
        q3 = list(csv.DictReader(f))
    with open(HERE / "q4" / "q4_summary.csv") as f:
        q4 = list(csv.DictReader(f))

    tcs = sorted({r["tc"] for r in q3}, key=float)

    per_conf_rows = {}
    for conf in CONFS:
        rows = []
        for tc in tcs:
            q3r = next(r for r in q3 if r["conformation"] == conf and r["tc"] == tc)
            q4r = next((r for r in q4 if r["conformation"] == conf and r["tc"] == tc), None)
            trust = int(q3r["n_trustworthy"])
            good = int(q4r["n_good"]) if q4r else 0
            rows.append({"tc": tc, "n_full": int(q3r["n_clusters_total"]),
                         "n_non_singleton": int(q3r["n_testable"]), "trustworthy": trust,
                         "good": good, "tier2": trust - good})
        per_conf_rows[conf] = rows

        out_path = HERE / f"q3q4_summary_{conf}.csv"
        with open(out_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["tc", "n_full", "n_non_singleton", "trustworthy", "good", "tier2"])
            w.writeheader()
            w.writerows(rows)
        print(f"Wrote {out_path}")

        print(f"\n{conf}")
        print(f"{'tc':<6}{'full':>10}{'non_singleton':>15}{'trustworthy':>13}{'good':>8}{'tier2':>8}")
        for r in rows:
            print(f"{r['tc']:<6}{r['n_full']:>10,}{r['n_non_singleton']:>15,}{r['trustworthy']:>13,}"
                  f"{r['good']:>8,}{r['tier2']:>8,}")

    # --- render image: one table panel per conformation ---
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    for ax, conf in zip(axes, CONFS):
        ax.axis("off")
        rows = per_conf_rows[conf]
        col_labels = ["Tc", "full", "non-singleton", "trustworthy", "good", "tier2"]
        cell_text = [[r["tc"], f"{r['n_full']:,}", f"{r['n_non_singleton']:,}", f"{r['trustworthy']:,}",
                      f"{r['good']:,}", f"{r['tier2']:,}"] for r in rows]
        table = ax.table(cellText=cell_text, colLabels=col_labels, loc="center", cellLoc="center")
        table.auto_set_font_size(False)
        table.set_fontsize(9)
        table.scale(1, 1.5)
        for (row, col), cell in table.get_celld().items():
            if row == 0:
                cell.set_facecolor(COLORS[conf])
                cell.set_text_props(color="white", fontweight="bold")
            elif row % 2 == 0:
                cell.set_facecolor("#f2f2f2")
        ax.set_title(conf, fontsize=13, fontweight="bold", color=COLORS[conf], pad=14)

    fig.suptitle("Trustworthy (internal-consistency chi² test, p<0.05) / good (significantly "
                  "better than pool, of the trustworthy) / tier2 = trustworthy - good",
                  fontsize=11, y=1.04)
    plt.tight_layout()
    png_path = HERE / "q3q4_summary_table.png"
    svg_path = HERE / "q3q4_summary_table.svg"
    plt.savefig(png_path, bbox_inches="tight", facecolor="white", dpi=200)
    plt.savefig(svg_path, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"\nWrote {png_path}")
    print(f"Wrote {svg_path}")


if __name__ == "__main__":
    main()
