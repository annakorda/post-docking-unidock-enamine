#!/usr/bin/env python3
"""
Author: Anna Korda

MM-GBSA rescoring of a conformation's final hit set (apply_score_cutoff.py's
already-built final_hits_<conf>_cutoff<C>.csv + _partNN.sdf), via Uni-GBSA
(single-point energy-minimized pose, GB solvation -- unigbsa-pipeline's own
"mode=em" default, not MD; real per-compound wall time benchmarked locally
at ~25s across 4 real compounds spanning all 3 receptors, essentially
single-threaded internally -- -nt 1/2/8 all gave ~110-120% CPU, so the real
lever is running many compounds AS SEPARATE CONCURRENT PROCESSES, not more
threads per compound).

Each compound gets its own throwaway working directory (unigbsa-pipeline
writes ~7 files / ~3.7MB per run without --verbose) that is deleted right
after its result is parsed -- at ~230k compounds total across all 3
conformations, keeping these would be ~850GB and ~1.6M inodes, so nothing
from a compound's run is kept except the one real number extracted from its
result CSV.

Chunked/resumable/parallel, same pattern as the rest of this repo: split the
hit list into chunks, ProcessPoolExecutor across --workers, each worker runs
its chunk's compounds one at a time (so total concurrency = --workers, not
--workers x threads-per-job), write a per-chunk CSV, skip chunks whose output
already exists on rerun, periodic progress logging.

Usage:
    python mmgbsa_rescore.py --conformation c1 \\
        --vs_results_dir /users/gpcr/annak/ultra-large/vs_results \\
        --receptor /users/gpcr/annak/ultra-large/receptor_prep/c1_gmxready.pdb \\
        --out_dir /users/gpcr/annak/ultra-large/analysis_results/mmgbsa_rescore \\
        --cutoff -7.0 --workers 64
"""
import argparse
import csv
import math
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


def load_hits(vs_results_dir, conf, cutoff, top_pct=None):
    hit_dir = Path(vs_results_dir) / f"final_hits_{conf}_cutoff{cutoff:g}"
    csv_path = hit_dir / f"final_hits_{conf}_cutoff{cutoff:g}.csv"
    if not csv_path.is_file():
        raise SystemExit(f"{csv_path} not found -- run apply_score_cutoff.py for {conf} first")
    with open(csv_path) as f:
        rows = list(csv.DictReader(f))
    ids = [r["compound_id"] for r in rows]
    # carry through every column collect_final_hits.py already computed (ad4_score,
    # total_TEU, single_TEU) so the final MM-GBSA output is one CSV with everything,
    # not just the mmgbsa number on its own
    other_info = {r["compound_id"]: {"ad4_score": r["ad4_score"], "total_TEU": r["total_TEU"],
                                      "single_TEU": r["single_TEU"]} for r in rows}

    part_files = sorted(hit_dir.glob(f"final_hits_{conf}_cutoff{cutoff:g}_part*.sdf"))
    if not part_files:
        raise SystemExit(f"No SDF parts found in {hit_dir}")
    blocks = []
    for pf in part_files:
        blocks.extend(_split_sdf_blocks(pf.read_text()))
    if len(blocks) != len(ids):
        raise SystemExit(f"{conf}: {len(ids)} CSV rows but {len(blocks)} SDF blocks -- "
                          f"expected the same rank order, refusing to guess a correspondence")

    if top_pct is not None:
        # rows are already rank-ordered best AD4 score first (apply_score_cutoff.py's own
        # guarantee), so the top N% is just a contiguous prefix -- no re-sorting needed
        n_top = max(1, math.ceil(len(ids) * top_pct / 100.0))
        ids, blocks = ids[:n_top], blocks[:n_top]

    return ids, blocks, other_info


def _run_one_compound(compound_id, block, receptor_path, scratch_root, timeout_s):
    work_dir = tempfile.mkdtemp(prefix=f"gbsa_{compound_id[:40]}_", dir=scratch_root)
    try:
        ligand_path = Path(work_dir) / f"{compound_id}.sdf"
        ligand_path.write_text(block)
        result_csv = Path(work_dir) / "result.csv"
        proc = subprocess.run(
            ["unigbsa-pipeline", "-i", str(receptor_path), "-l", str(ligand_path),
             "-o", str(result_csv), "-nt", "1"],
            cwd=work_dir, capture_output=True, text=True, timeout=timeout_s,
        )
        if not result_csv.is_file():
            return compound_id, None, f"no_output_rc{proc.returncode}: {proc.stderr[-300:]}"
        with open(result_csv) as f:
            row = next(csv.DictReader(f), None)
        if row is None or row.get("status") != "S":
            return compound_id, None, f"bad_status: {row}"
        return compound_id, float(row["TOTAL"]), None
    except subprocess.TimeoutExpired:
        return compound_id, None, "timeout"
    except Exception as e:
        return compound_id, None, f"exception: {e}"
    finally:
        shutil.rmtree(work_dir, ignore_errors=True)


FIELDNAMES = ["compound_id", "ad4_score", "total_TEU", "single_TEU", "mmgbsa_dg", "error"]


def _process_chunk(args):
    chunk, chunk_idx, receptor_path, other_info, out_dir_str, scratch_root, timeout_s = args
    out_dir = Path(out_dir_str)
    part_path = out_dir / f"mmgbsa_part_{chunk_idx:05d}.csv"
    if part_path.exists():
        return chunk_idx, len(chunk), True  # resumed from a prior run -- not real new work

    rows = []
    for compound_id, block in chunk:
        _, dg, err = _run_one_compound(compound_id, block, receptor_path, scratch_root, timeout_s)
        info = other_info[compound_id]
        rows.append({"compound_id": compound_id, "ad4_score": info["ad4_score"],
                     "total_TEU": info["total_TEU"], "single_TEU": info["single_TEU"],
                     "mmgbsa_dg": dg if dg is not None else "", "error": err or ""})

    with open(part_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDNAMES)
        w.writeheader()
        w.writerows(rows)
    return chunk_idx, len(chunk), False


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--conformation", required=True)
    ap.add_argument("--vs_results_dir", required=True)
    ap.add_argument("--receptor", required=True, help="gmx-ready receptor PDB from receptor_prep/")
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--cutoff", type=float, default=-7.0)
    ap.add_argument("--top_pct", type=float, default=None,
                     help="only score the top N%% by AD4 rank (contiguous prefix, already "
                          "rank-ordered) -- omit to score the full set")
    ap.add_argument("--chunk_size", type=int, default=20,
                     help="compounds per work-chunk -- small, since each compound is its own "
                          "~25s subprocess and chunks are the resume granularity")
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--timeout", type=int, default=300, help="per-compound subprocess timeout, seconds")
    ap.add_argument("--scratch_root", default=None,
                     help="fast local scratch for per-compound throwaway dirs (default: system tmp -- "
                          "use node-local storage on a cluster, not shared GPFS, if available)")
    args = ap.parse_args()

    conf = args.conformation
    receptor_path = str(Path(args.receptor).resolve())
    if not Path(receptor_path).is_file():
        raise SystemExit(f"{receptor_path} not found")
    if args.scratch_root:
        Path(args.scratch_root).mkdir(parents=True, exist_ok=True)

    out_dir = Path(args.out_dir) / f"mmgbsa_{conf}"
    out_dir.mkdir(parents=True, exist_ok=True)
    workers = args.workers or len(os.sched_getaffinity(0))

    top_label = f", top {args.top_pct:g}%" if args.top_pct is not None else ""
    log(f"Loading {conf} hits (cutoff {args.cutoff:g}{top_label}) ...")
    ids, blocks, other_info = load_hits(args.vs_results_dir, conf, args.cutoff, args.top_pct)
    n_total = len(ids)
    log(f"{n_total:,} compounds loaded for {conf}, receptor={receptor_path}, workers={workers}")

    items = list(zip(ids, blocks))
    chunks = [items[i:i + args.chunk_size] for i in range(0, n_total, args.chunk_size)]
    n_chunks = len(chunks)
    have = sum(1 for i in range(n_chunks) if (out_dir / f"mmgbsa_part_{i:05d}.csv").exists())
    log(f"{n_chunks:,} work-chunks of up to {args.chunk_size}, {have:,} already done, {workers} workers")

    tasks = [(chunks[i], i, receptor_path, other_info, str(out_dir), args.scratch_root, args.timeout)
              for i in range(n_chunks)]

    # log roughly every ~1% of chunks (at least every chunk if there are few) --
    # at c1's real scale (~106k compounds / 20 per chunk = ~5,320 chunks) logging
    # every single chunk would be thousands of near-identical lines; this still
    # gives real periodic progress + a live ETA, just not spammed
    log_every = max(1, n_chunks // 100)

    t0 = time.time()
    done_count = 0
    new_done_count = 0  # excludes resumed (instant, not real work) chunks -- rate/ETA use only this
    done_compounds = 0  # every chunk (resumed or new) passes through the loop below exactly
    # once, so this is the sole accumulator -- no separate have-based pre-seed (that
    # double-counted resumed chunks, since they're also counted as they "complete" below)
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(_process_chunk, t): t[1] for t in tasks}
        for fut in as_completed(futures):
            chunk_idx, n_in_chunk, was_resumed = fut.result()
            done_count += 1
            done_compounds += n_in_chunk
            if not was_resumed:
                new_done_count += 1
            if done_count % log_every == 0 or done_count == n_chunks:
                elapsed = time.time() - t0
                rate = new_done_count / elapsed if elapsed > 0 else 0
                eta_s = (n_chunks - done_count) / rate if rate else float("inf")
                pct = 100 * done_count / n_chunks
                log(f"  chunk {done_count:,}/{n_chunks:,} ({pct:.1f}%) done -- "
                    f"~{done_compounds:,}/{n_total:,} compounds -- "
                    f"{elapsed / 60:.1f} min elapsed ({new_done_count:,} chunks actually computed "
                    f"this run, {done_count - new_done_count:,} resumed instantly) -- "
                    f"ETA {eta_s / 60:.1f} min ({eta_s / 3600:.1f} h)")

    log("all chunks done -- merging ...")
    all_rows = []
    for i in range(n_chunks):
        with open(out_dir / f"mmgbsa_part_{i:05d}.csv") as f:
            all_rows.extend(csv.DictReader(f))

    n_ok = sum(1 for r in all_rows if r["mmgbsa_dg"])
    n_err = len(all_rows) - n_ok
    log(f"{n_ok:,} succeeded, {n_err:,} failed/errored")

    # one combined CSV: everything already known about each compound (ad4_score,
    # total_TEU, single_TEU) plus the mmgbsa_dg -- only the final delta G is kept
    # from unigbsa-pipeline's own output, not its full energy-component breakdown
    scored = [r for r in all_rows if r["mmgbsa_dg"]]
    scored.sort(key=lambda r: float(r["mmgbsa_dg"]))
    final_path = out_dir / f"mmgbsa_results_{conf}.csv"
    with open(final_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["rank", "compound_id", "mmgbsa_dg", "ad4_score", "total_TEU", "single_TEU"])
        for rank, r in enumerate(scored, start=1):
            w.writerow([rank, r["compound_id"], r["mmgbsa_dg"], r["ad4_score"],
                       r["total_TEU"], r["single_TEU"]])
    log(f"Wrote {final_path}")

    failed_path = out_dir / f"mmgbsa_failed_{conf}.csv"
    with open(failed_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDNAMES)
        w.writeheader()
        w.writerows(r for r in all_rows if not r["mmgbsa_dg"])
    log(f"Wrote {failed_path}")

    # sanity check: do the top MM-GBSA hits also look reasonable by AD4?
    top_n = min(20, len(scored))
    sanity_path = out_dir / f"mmgbsa_top{top_n}_vs_ad4_{conf}.csv"
    with open(sanity_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["mmgbsa_rank", "compound_id", "mmgbsa_dg", "ad4_score",
                   "total_TEU", "single_TEU", "ad4_rank_within_full_set"])
        ad4_sorted_ids = sorted(other_info, key=lambda cid: float(other_info[cid]["ad4_score"]))
        ad4_rank = {cid: i + 1 for i, cid in enumerate(ad4_sorted_ids)}
        for rank, r in enumerate(scored[:top_n], start=1):
            w.writerow([rank, r["compound_id"], r["mmgbsa_dg"], r["ad4_score"],
                       r["total_TEU"], r["single_TEU"], ad4_rank[r["compound_id"]]])
    log(f"Wrote {sanity_path} (top {top_n} MM-GBSA hits with their AD4 score and AD4-only rank, "
        f"for a quick sanity cross-check)")


if __name__ == "__main__":
    main()
