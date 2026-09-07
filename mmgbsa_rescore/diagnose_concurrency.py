#!/usr/bin/env python3
"""
Author: Anna Korda

Finds the real, safe concurrency ceiling on this specific node instead of
guessing worker counts one sbatch job at a time. Runs the EXACT same
per-compound function mmgbsa_rescore.py uses (copied verbatim, not
reimplemented) at increasing concurrency levels (1, 4, 8, 16, 32, 64) inside
ONE job, each level using a small fixed batch (2 compounds per worker, so
every level's *expected* wall time is the same ~50s if nothing degrades --
any real per-worker slowdown as concurrency rises shows up directly as a
rising wall time / slowdown_factor, not something you have to infer from a
stalled sbatch log).

A per-compound --timeout keeps a genuinely stuck level from hanging the
whole diagnostic -- a level that times out gets recorded as such and the
script moves on to the next one rather than blocking forever.

Usage:
    python diagnose_concurrency.py \\
        --receptor /users/gpcr/annak/ultra-large/post-docking-unidock-enamine/receptor_prep/c1_gmxready.pdb \\
        --vs_results_dir /users/gpcr/annak/ultra-large/vs_results \\
        --conformation c1 --cutoff -7.0 \\
        --levels 1,4,8,16,32,64 --out_csv concurrency_report.csv --scratch_root /dev/shm
"""
import argparse
import csv
import os
import shutil
import subprocess
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path


def log(msg):
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}", flush=True)


def _split_sdf_blocks(text):
    blocks = []
    start = 0
    idx = text.find("$$$$")
    while idx != -1:
        end = idx + 4
        if text[end:end + 2] == "\r\n":
            end += 2
        elif text[end:end + 1] == "\n":
            end += 1
        blocks.append(text[start:end])
        start = end
        idx = text.find("$$$$", start)
    if text[start:].strip():
        blocks.append(text[start:])
    return blocks


def load_hits(vs_results_dir, conf, cutoff, n_needed):
    hit_dir = Path(vs_results_dir) / f"final_hits_{conf}_cutoff{cutoff:g}"
    csv_path = hit_dir / f"final_hits_{conf}_cutoff{cutoff:g}.csv"
    with open(csv_path) as f:
        ids = [r["compound_id"] for r in csv.DictReader(f)]
    part_files = sorted(hit_dir.glob(f"final_hits_{conf}_cutoff{cutoff:g}_part*.sdf"))
    blocks = []
    for pf in part_files:
        blocks.extend(_split_sdf_blocks(pf.read_text()))
        if len(blocks) >= n_needed:
            break
    return list(zip(ids[:len(blocks)], blocks))[:n_needed]


# same function as mmgbsa_rescore.py's _run_one_compound, copied verbatim so
# this diagnostic tests the exact real code path, not a reimplementation
def _run_one_compound(compound_id, block, receptor_path, scratch_root, timeout_s):
    work_dir = tempfile.mkdtemp(prefix=f"gbsa_{compound_id[:40]}_", dir=scratch_root)
    t0 = time.time()
    try:
        ligand_path = Path(work_dir) / f"{compound_id}.sdf"
        ligand_path.write_text(block)
        result_csv = Path(work_dir) / "result.csv"
        proc = subprocess.run(
            ["unigbsa-pipeline", "-i", str(receptor_path), "-l", str(ligand_path),
             "-o", str(result_csv), "-nt", "1"],
            cwd=work_dir, capture_output=True, text=True, timeout=timeout_s,
        )
        wall = time.time() - t0
        if not result_csv.is_file():
            return compound_id, None, f"no_output_rc{proc.returncode}", wall
        with open(result_csv) as f:
            row = next(csv.DictReader(f), None)
        if row is None or row.get("status") != "S":
            return compound_id, None, f"bad_status", wall
        return compound_id, float(row["TOTAL"]), None, wall
    except subprocess.TimeoutExpired:
        return compound_id, None, "timeout", time.time() - t0
    except Exception as e:
        return compound_id, None, f"exception: {e}", time.time() - t0
    finally:
        shutil.rmtree(work_dir, ignore_errors=True)


def run_level(level, items, receptor_path, scratch_root, timeout_s):
    t0 = time.time()
    per_compound_walls = []
    n_ok, n_err = 0, 0
    with ProcessPoolExecutor(max_workers=level) as pool:
        futures = [pool.submit(_run_one_compound, cid, block, receptor_path, scratch_root, timeout_s)
                   for cid, block in items]
        for fut in as_completed(futures):
            _, dg, err, wall = fut.result()
            per_compound_walls.append(wall)
            if dg is not None:
                n_ok += 1
            else:
                n_err += 1
    total_wall = time.time() - t0
    return total_wall, per_compound_walls, n_ok, n_err


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--receptor", required=True)
    ap.add_argument("--vs_results_dir", required=True)
    ap.add_argument("--conformation", required=True)
    ap.add_argument("--cutoff", type=float, default=-7.0)
    ap.add_argument("--levels", default="1,4,8,16,32,64")
    ap.add_argument("--per_worker", type=int, default=2,
                     help="compounds per worker at each level -- kept fixed so every level's "
                          "expected wall time is the same if nothing degrades")
    ap.add_argument("--timeout", type=int, default=180,
                     help="per-compound timeout -- keeps one stuck level from hanging the whole run")
    ap.add_argument("--out_csv", required=True)
    ap.add_argument("--scratch_root", default="/dev/shm")
    args = ap.parse_args()

    receptor_path = str(Path(args.receptor).resolve())
    levels = [int(x) for x in args.levels.split(",")]
    max_needed = sum(level * args.per_worker for level in levels)  # non-overlapping compounds per level

    log(f"Loading {max_needed} real compounds from {args.conformation} (cutoff {args.cutoff:g}) ...")
    all_items = load_hits(args.vs_results_dir, args.conformation, args.cutoff, max_needed)
    log(f"{len(all_items)} loaded, receptor={receptor_path}")

    baseline_per_compound = None
    rows = []
    cursor = 0
    for level in levels:
        n_items = level * args.per_worker
        items = all_items[cursor:cursor + n_items]
        cursor += n_items
        if len(items) < n_items:
            log(f"  level {level}: not enough real compounds left ({len(items)}/{n_items}), skipping")
            continue

        log(f"Level {level}: running {len(items)} real compounds ({args.per_worker}/worker) ...")
        total_wall, per_compound_walls, n_ok, n_err = run_level(
            level, items, receptor_path, args.scratch_root, args.timeout)
        mean_wall = sum(per_compound_walls) / len(per_compound_walls) if per_compound_walls else 0
        max_wall = max(per_compound_walls) if per_compound_walls else 0

        if baseline_per_compound is None and level == 1:
            baseline_per_compound = mean_wall

        expected_total = args.per_worker * (baseline_per_compound or mean_wall)
        slowdown = total_wall / expected_total if expected_total else float("nan")

        row = {"level": level, "n_compounds": len(items), "n_ok": n_ok, "n_err": n_err,
               "total_wall_s": f"{total_wall:.1f}", "mean_per_compound_s": f"{mean_wall:.1f}",
               "max_per_compound_s": f"{max_wall:.1f}", "expected_total_s": f"{expected_total:.1f}",
               "slowdown_factor": f"{slowdown:.2f}"}
        rows.append(row)
        log(f"  level {level}: total={total_wall:.1f}s, mean/compound={mean_wall:.1f}s, "
            f"max/compound={max_wall:.1f}s, slowdown={slowdown:.2f}x, ok={n_ok}, err={n_err}")

        # write incrementally -- crash-safe, and gives a partial answer even
        # if a later (higher) level hangs past its own per-compound timeouts
        with open(args.out_csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(row.keys()))
            w.writeheader()
            w.writerows(rows)

    log(f"Wrote {args.out_csv}")
    log("Real safe concurrency ceiling = the highest level where slowdown_factor stays close to 1.0 "
        "(e.g. <2x) -- above that, the level itself is the problem, not any single compound.")


if __name__ == "__main__":
    main()
