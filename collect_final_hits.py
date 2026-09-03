#!/usr/bin/env python3
"""
Author: Anna Korda

Final collection step, run after strain_filter.py (dry-run is enough --
this script never recomputes strain, it only reads that run's report):
takes every compound that passed BOTH filters (salt bridge + not strained),
orders them by AD4 score (best/most negative first, same convention as
ranked_hits_<conf>.csv), and writes:
  - final_hits_<conf>.csv        -- rank, compound_id, ad4_score, total_TEU,
                                     single_TEU for every survivor
  - score_range_report_<conf>.txt -- real min/max/mean/median across
                                     ad4_score, total_TEU, single_TEU
  - final_hits_<conf>_partNN.sdf  -- the actual structures, split into
                                     --chunk_size-sized multi-molecule SDFs
                                     (default 50000), rank-ordered within
                                     and across files (part01 has the best
                                     scores, part02 the next best, etc.)

Only the winning pose already selected earlier in the pipeline is used --
same pose, same file bytes, never re-derived or modified. Currently only
includes compounds the strain filter could confidently score (strained ==
False in the report) -- compounds it couldn't score (rigid/flagged) are
NOT included here yet; that's a real, deliberate open item (Anna: "if we
cant score them well we should keep them... lets keep for the future").

Same chunked/resumable/parallel pattern as the rest of this project, with
a fetch stage (parallel, I/O-bound, same real per-file GPFS cost as
strain_filter.py) followed by a fast ordered merge (work-chunks are
processed in the same rank-ascending order they were split in, so
concatenating them in chunk-index order reproduces the correct global
rank order with no extra sorting step).

Usage:
    python collect_final_hits.py --conformation c1 \\
        --vs_results_dir /users/gpcr/annak/ultra-large/vs_results \\
        --chunk_size 50000 --workers 32
"""
import argparse
import csv
import os
import statistics
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


def build_ranked_survivor_list(vs_dir, conf):
    """Real join of three already-computed CSVs -- no recomputation, just
    lookup. Returns a list of dicts, sorted best-score-first."""
    passed_dir = vs_dir / f"interactions_passed_{conf}"
    strained_dir = vs_dir / f"strain_filtered_{conf}"
    poses_dir = vs_dir / f"ad4_redock_{conf}"

    strain_info = {}
    with open(next(strained_dir.glob("strain_filter_report_*.csv"))) as f:
        for row in csv.DictReader(f):
            if row["strained"] == "False":
                strain_info[row["compound_id"]] = (row["total_TEU"], row["single_TEU"])
    log(f"{len(strain_info):,} compounds passed the strain filter")

    scores = {}
    with open(next(passed_dir.glob("ranked_hits_*.csv"))) as f:
        for row in csv.DictReader(f):
            if row["compound_id"] in strain_info:
                scores[row["compound_id"]] = float(row["winning_pose_score"])

    id_to_path = {}
    with open(poses_dir / "scores.csv") as f:
        for row in csv.DictReader(f):
            if row["compound_id"] in scores:
                id_to_path[row["compound_id"]] = row["path"]

    winning_idx = {}
    reordered_csv = next(passed_dir.glob("reordered_compounds_*.csv"), None)
    if reordered_csv:
        with open(reordered_csv) as f:
            for row in csv.DictReader(f):
                if row["compound_id"] in scores:
                    winning_idx[row["compound_id"]] = int(row["winning_pose_index_0based"])

    survivors = []
    missing = 0
    for compound_id, score in scores.items():
        if compound_id not in id_to_path:
            missing += 1
            continue
        total_teu, single_teu = strain_info[compound_id]
        survivors.append({
            "compound_id": compound_id,
            "ad4_score": score,
            "total_TEU": total_teu,
            "single_TEU": single_teu,
            "path": poses_dir / id_to_path[compound_id],
            "pose_idx": winning_idx.get(compound_id, 0),
        })
    if missing:
        log(f"WARNING: {missing} strain-passing compound(s) had no path in scores.csv -- skipped")

    survivors.sort(key=lambda d: d["ad4_score"])
    for rank, d in enumerate(survivors, start=1):
        d["rank"] = rank
    return survivors


def _fetch_chunk(chunk, chunk_idx, parts_dir_str):
    parts_dir = Path(parts_dir_str)
    part_path = parts_dir / f"fetch_{chunk_idx:05d}.txt"
    if part_path.exists():
        return chunk_idx, None, len(chunk)

    out_texts = []
    n_missing = 0
    for row in chunk:
        try:
            text = row["path"].read_text()
            blocks = _split_sdf_blocks(text)
            pose_idx = row["pose_idx"]
            if not (0 <= pose_idx < len(blocks)):
                n_missing += 1
                continue
            out_texts.append(blocks[pose_idx])
        except Exception:
            n_missing += 1
    part_path.write_text("".join(out_texts))
    return chunk_idx, None, len(chunk) - n_missing


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--conformation", required=True)
    ap.add_argument("--vs_results_dir", required=True)
    ap.add_argument("--chunk_size", type=int, default=50000,
                     help="compounds per output multi-SDF file")
    ap.add_argument("--fetch_chunk_size", type=int, default=1000,
                     help="work-unit size for the parallel fetch stage -- unrelated to "
                          "--chunk_size, just controls parallel granularity")
    ap.add_argument("--workers", type=int, default=None)
    args = ap.parse_args()

    conf = args.conformation
    vs_dir = Path(args.vs_results_dir)
    out_root = vs_dir / f"final_hits_{conf}"
    out_root.mkdir(parents=True, exist_ok=True)

    log(f"Building ranked survivor list for {conf} ...")
    survivors = build_ranked_survivor_list(vs_dir, conf)
    n_total = len(survivors)
    log(f"{n_total:,} compounds passed both filters, ranked by AD4 score")
    if n_total == 0:
        raise SystemExit("No survivors -- nothing to write.")

    # Real score-range report, computed directly from the already-joined
    # data (no extra I/O needed for this part).
    def _stats(vals):
        vals = [float(v) for v in vals]
        return {"min": min(vals), "max": max(vals), "mean": statistics.mean(vals),
                "median": statistics.median(vals)}

    ad4_stats = _stats([s["ad4_score"] for s in survivors])
    total_stats = _stats([s["total_TEU"] for s in survivors])
    single_stats = _stats([s["single_TEU"] for s in survivors])

    report_path = out_root / f"score_range_report_{conf}.txt"
    with open(report_path, "w") as f:
        f.write(f"Final hits score range -- conformation {conf}\n")
        f.write(f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Total compounds (passed salt bridge + strain filters): {n_total:,}\n\n")
        for name, s in [("AD4 score (kcal/mol, winning pose)", ad4_stats),
                        ("Total strain (TEU)", total_stats),
                        ("Single strain (TEU)", single_stats)]:
            f.write(f"{name}:\n")
            f.write(f"  min={s['min']:.4f}  max={s['max']:.4f}  "
                    f"mean={s['mean']:.4f}  median={s['median']:.4f}\n")
    log(f"Wrote {report_path}")
    log(f"  AD4 score range: {ad4_stats['min']:.2f} to {ad4_stats['max']:.2f} kcal/mol "
        f"(median {ad4_stats['median']:.2f})")

    # Ranked CSV -- every survivor, in final order.
    final_csv_path = out_root / f"final_hits_{conf}.csv"
    with open(final_csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["rank", "compound_id", "ad4_score", "total_TEU", "single_TEU"])
        for s in survivors:
            w.writerow([s["rank"], s["compound_id"], f"{s['ad4_score']:.4f}",
                       s["total_TEU"], s["single_TEU"]])
    log(f"Wrote {final_csv_path}")

    # Parallel fetch stage -- work-chunks split in rank order, so
    # concatenating fetch_NNNNN.txt files in index order reproduces the
    # correct global rank order with no extra sorting.
    parts_dir = out_root / ".parts"
    parts_dir.mkdir(parents=True, exist_ok=True)
    fetch_chunks = [survivors[i:i + args.fetch_chunk_size]
                    for i in range(0, n_total, args.fetch_chunk_size)]
    n_fetch_chunks = len(fetch_chunks)
    workers = args.workers or len(os.sched_getaffinity(0))

    have = {i for i in range(n_fetch_chunks) if (parts_dir / f"fetch_{i:05d}.txt").exists()}
    todo = [i for i in range(n_fetch_chunks) if i not in have]
    log(f"Fetching real pose text: {n_fetch_chunks} work-chunks of up to "
        f"{args.fetch_chunk_size}, {len(have)} already done, {workers} workers")

    t0 = time.time()
    if todo:
        done_count = 0
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = {pool.submit(_fetch_chunk, fetch_chunks[i], i, str(parts_dir)): i
                       for i in todo}
            for fut in as_completed(futures):
                idx, err, n_fetched = fut.result()
                done_count += 1
                rate = done_count / (time.time() - t0) if time.time() > t0 else 0
                eta = (len(todo) - done_count) / rate if rate else float("inf")
                log(f"  fetch chunk {done_count}/{len(todo)} done "
                    f"({time.time() - t0:.0f}s elapsed, ETA {eta:.0f}s)")

    log("all fetches done -- writing final rank-ordered multi-SDF files ...")
    part_num = 1
    count_in_part = 0
    out_f = open(out_root / f"final_hits_{conf}_part{part_num:02d}.sdf", "w")
    written = 0
    for i in range(n_fetch_chunks):
        text = (parts_dir / f"fetch_{i:05d}.txt").read_text()
        if not text:
            continue
        for block in _split_sdf_blocks(text):
            if count_in_part >= args.chunk_size:
                out_f.close()
                part_num += 1
                count_in_part = 0
                out_f = open(out_root / f"final_hits_{conf}_part{part_num:02d}.sdf", "w")
            out_f.write(block)
            count_in_part += 1
            written += 1
    out_f.close()

    import shutil
    shutil.rmtree(parts_dir, ignore_errors=True)

    log(f"Wrote {written:,} compounds across {part_num} multi-SDF file(s) "
        f"(up to {args.chunk_size:,} compounds each) in "
        f"{out_root}/final_hits_{conf}_partNN.sdf")
    if written != n_total:
        log(f"WARNING: wrote {written:,} but expected {n_total:,} -- "
            f"{n_total - written} compound(s) missing (check for fetch errors)")


if __name__ == "__main__":
    main()
