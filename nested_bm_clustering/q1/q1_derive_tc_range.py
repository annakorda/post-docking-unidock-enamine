#!/usr/bin/env python3
"""
Author: Anna Korda

Q1: derive a Tc range from the BM+Tanimoto curves already on hand (no new
clustering -- just aggregating nested_cluster_summary_<conf>.csv and the
already-written nested_cluster_assignment_<conf>_<threshold>.csv files).

Upper bound: the Tc where the total-cluster-count curve's slope (dN/dTc)
falls back to <=50% of its own peak value, scanning from the tightest
threshold downward. Below-peak-on-the-tight-side means BM+Tanimoto is
still doing substantial extra splitting; above the crossing, it's mostly
flat/near the BM-only ceiling.

Lower bound: the top-10-clusters-% curve's climb rate (as Tc loosens) isn't
monotonic -- real data shows it rises toward the tight end, dips to a local
minimum, then surges again toward the loose end. That surge (climb rate
itself increasing again past the minimum) is the onset of a chaining/step
event, not a real chemistry signal -- so the lower bound is the Tc at that
local minimum, the last point before the curve starts accelerating.

Both derivatives use central differences on the actual tested threshold
grid (spacing 0.05); reported at each midpoint. All curves computed fresh
per conformation, then reconciled into ONE range applied consistently to
all three (same convention as pick_tanimoto_cutoff.py used before).
"""
import csv
from collections import Counter
from pathlib import Path

HERE = Path(__file__).parent
NESTED_DIR = HERE.parent
CONFS = ["c1", "ref1", "c5"]
THRESHOLDS = [0.90, 0.85, 0.80, 0.75, 0.70, 0.65, 0.60, 0.55, 0.50, 0.45, 0.40, 0.35]


def load_n_clusters(conf):
    with open(NESTED_DIR / f"nested_cluster_summary_{conf}.csv") as f:
        rows = {float(r["threshold"]): int(r["n_clusters"]) for r in csv.DictReader(f)}
    return [rows[t] for t in THRESHOLDS]


def load_top10_pct(conf):
    vals = []
    for t in THRESHOLDS:
        path = NESTED_DIR / f"nested_cluster_assignment_{conf}_{t:.2f}.csv"
        with open(path) as f:
            cluster_ids = [row["cluster_id"] for row in csv.DictReader(f)]
        counts = Counter(cluster_ids)
        pool_size = len(cluster_ids)
        top10_sum = sum(sorted(counts.values(), reverse=True)[:10])
        vals.append(100 * top10_sum / pool_size)
    return vals


def first_diff(ts, ys):
    """dy/dt at each midpoint between consecutive ts."""
    mids, slopes = [], []
    for i in range(len(ts) - 1):
        dt = ts[i + 1] - ts[i]
        mids.append((ts[i] + ts[i + 1]) / 2)
        slopes.append((ys[i + 1] - ys[i]) / dt)
    return mids, slopes


def main():
    per_conf = {}
    for conf in CONFS:
        ts_desc = THRESHOLDS
        ts_asc = list(reversed(ts_desc))
        n_clusters_desc = load_n_clusters(conf)
        n_clusters_asc = list(reversed(n_clusters_desc))
        top10_desc = load_top10_pct(conf)
        top10_asc = list(reversed(top10_desc))

        mids1_n, slope1_n = first_diff(ts_asc, n_clusters_asc)
        mids1_t10, slope1_t10 = first_diff(ts_asc, top10_asc)
        mids2_t10, slope2_t10 = first_diff(mids1_t10, slope1_t10)

        per_conf[conf] = dict(ts_asc=ts_asc, n_clusters_asc=n_clusters_asc, top10_asc=top10_asc,
                               mids1_n=mids1_n, slope1_n=slope1_n,
                               mids1_t10=mids1_t10, slope1_t10=slope1_t10,
                               mids2_t10=mids2_t10, slope2_t10=slope2_t10)

    # --- write per-conformation derivative tables ---
    deriv_path = HERE / "curve_derivatives.csv"
    with open(deriv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["conformation", "tc_low", "tc_high", "midpoint",
                     "n_clusters_slope1", "top10_pct_slope1"])
        for conf in CONFS:
            d = per_conf[conf]
            for i, mid in enumerate(d["mids1_n"]):
                w.writerow([conf, f"{d['ts_asc'][i]:.2f}", f"{d['ts_asc'][i+1]:.2f}", f"{mid:.3f}",
                            f"{d['slope1_n'][i]:.1f}", f"{d['slope1_t10'][i]:.4f}"])
    print(f"Wrote {deriv_path}")

    deriv2_path = HERE / "top10_second_derivative.csv"
    with open(deriv2_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["conformation", "midpoint_low", "midpoint_high", "midpoint2", "top10_pct_slope2"])
        for conf in CONFS:
            d = per_conf[conf]
            for i, mid2 in enumerate(d["mids2_t10"]):
                w.writerow([conf, f"{d['mids1_t10'][i]:.3f}", f"{d['mids1_t10'][i+1]:.3f}", f"{mid2:.3f}",
                            f"{d['slope2_t10'][i]:.5f}"])
    print(f"Wrote {deriv2_path}")

    # --- apply rules per conformation ---
    bounds = {}
    print(f"\n{'conf':<6}{'peak_slope1_n':>15}{'upper_bound_Tc':>17}{'lower_bound_Tc':>17}")
    for conf in CONFS:
        d = per_conf[conf]
        mids1_n, slope1_n = d["mids1_n"], d["slope1_n"]
        peak = max(slope1_n)
        half = 0.5 * peak
        peak_idx = slope1_n.index(peak)
        # upper bound: scan from tight end (last midpoint, ~0.875) down toward peak,
        # find first midpoint (in descending Tc) whose slope1_n exceeds half-max
        upper = None
        for i in range(len(mids1_n) - 1, -1, -1):
            if slope1_n[i] >= half:
                upper = mids1_n[i]
                break
        if upper is None:
            upper = mids1_n[-1]

        # lower bound: the top10-% "climb rate as Tc loosens" (-slope1_t10) isn't
        # monotonic -- it rises from the tight end to a local peak, dips to a local
        # minimum, then surges again toward the loose end. That surge past the
        # minimum is the real chaining-event signal (second derivative flips from
        # decelerating to reaccelerating); the flat local peak near the tight end is
        # a different, uninteresting kind of flatness (near-singleton ceiling), not
        # what we want. So: find the first interior local minimum of the climb-rate
        # curve scanning from the loose end inward -- the point right before the
        # curve starts accelerating again.
        mids1_t10, slope1_t10 = d["mids1_t10"], d["slope1_t10"]
        climb_rate = [-s for s in slope1_t10]  # ascending Tc order (loose -> tight)
        lower = None
        for i in range(1, len(climb_rate) - 1):
            if climb_rate[i] < climb_rate[i - 1] and climb_rate[i] < climb_rate[i + 1]:
                lower = mids1_t10[i]
                break
        if lower is None:
            lower = mids1_t10[0]

        bounds[conf] = (lower, upper)
        print(f"{conf:<6}{peak:>15.1f}{upper:>17.3f}{lower:>17.3f}")

    lowers = [b[0] for b in bounds.values()]
    uppers = [b[1] for b in bounds.values()]
    combined_lower = max(lowers)   # tightest (highest) of the per-conf lower bounds -> safe for all
    combined_upper = min(uppers)   # loosest (lowest) of the per-conf upper bounds -> safe for all
    print(f"\nPer-conformation bounds: {bounds}")
    print(f"Combined range (intersection, safe for all 3 conformations): "
          f"[{combined_lower:.3f}, {combined_upper:.3f}]")

    # round to the nearest tested threshold grid point for practical use
    def nearest_tested(x):
        return min(THRESHOLDS, key=lambda t: abs(t - x))

    lo_t, hi_t = nearest_tested(combined_lower), nearest_tested(combined_upper)
    if lo_t > hi_t:
        lo_t, hi_t = hi_t, lo_t
    tested_range = [t for t in sorted(THRESHOLDS) if lo_t <= t <= hi_t]
    print(f"Rounded to tested grid: [{lo_t:.2f}, {hi_t:.2f}] -> thresholds to carry forward: {tested_range}")

    with open(HERE / "tc_range_result.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["conformation", "lower_bound_raw", "upper_bound_raw"])
        for conf in CONFS:
            w.writerow([conf, f"{bounds[conf][0]:.3f}", f"{bounds[conf][1]:.3f}"])
        w.writerow(["COMBINED", f"{combined_lower:.3f}", f"{combined_upper:.3f}"])
        w.writerow(["COMBINED_ROUNDED", f"{lo_t:.2f}", f"{hi_t:.2f}"])
    print(f"Wrote {HERE / 'tc_range_result.csv'}")


if __name__ == "__main__":
    main()
