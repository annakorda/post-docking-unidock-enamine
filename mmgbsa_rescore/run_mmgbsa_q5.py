#!/usr/bin/env python3
"""
Author: Anna Korda

MM-GBSA rescoring of the Q5 curated hit list (q5_mmgbsa_<conf>.sdf --
medoid + best-scorer per trustworthy+good cluster, plus qualifying
singletons; see ALL_QUESTIONS.md), via Uni-GBSA (single-point energy-
minimized pose, GB solvation -- unigbsa-pipeline's own "mode=em" default,
not MD).

Same core scoring/chunking/resume/progress-logging machinery as the
original mmgbsa_rescore.py this was adapted from -- the only real
difference is where the input comes from: that script reads
apply_score_cutoff.py's `final_hits_<conf>_cutoff<C>.csv` + `_partNN.sdf`
convention (the full, uncurated hit set, sliced by --top_pct); this one
reads a single q5_mmgbsa_<conf>.sdf directly (already the fully curated
set -- no further top-N slicing here, we score all of it) and carries
CLUSTER_ID/ROLE/SOURCE through into the final output alongside AD4_SCORE
and the new mmgbsa_dg.

Each compound gets its own throwaway working directory (unigbsa-pipeline
writes ~7 files / ~3.7MB per run without --verbose), deleted right after
its result is parsed -- nothing is kept from a compound's run except the
one real number extracted from its result CSV.

Usage:
    python run_mmgbsa_q5.py --conformation c1 \\
        --sdf ../nested_bm_clustering/q5/q5_mmgbsa_c1.sdf \\
        --receptor ../receptor_prep/c1_gmxready.pdb \\
        --out_dir results \\
        --workers 8
"""
import argparse
import csv
import os
import re
import shutil
import subprocess
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

PROPERTY_RE = re.compile(r">\s*<(\w+)>\s*\n(.*?)\n\n", re.DOTALL)


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


def load_hits(sdf_path):
    text = Path(sdf_path).read_text()
    blocks = _split_sdf_blocks(text)
    ids, meta = [], {}
    for block in blocks:
        props = dict(PROPERTY_RE.findall(block))
        cid = props.get("COMPOUND_ID")
        if cid is None:
            raise SystemExit(f"{sdf_path}: a block is missing the COMPOUND_ID property tag -- "
                              f"expected q5_mmgbsa_<conf>.sdf's own format, refusing to guess an id")
        ids.append(cid)
        meta[cid] = {"cluster_id": props.get("CLUSTER_ID", ""), "role": props.get("ROLE", ""),
                     "source": props.get("SOURCE", ""), "ad4_score": props.get("AD4_SCORE", "")}
    return ids, blocks, meta


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


FIELDNAMES = ["compound_id", "conformation", "cluster_id", "role", "source", "ad4_score", "mmgbsa_dg", "error"]


def _process_chunk(args):
    chunk, chunk_idx, conf, receptor_path, meta, out_dir_str, scratch_root, timeout_s = args
    out_dir = Path(out_dir_str)
    part_path = out_dir / f"mmgbsa_part_{chunk_idx:05d}.csv"
    if part_path.exists():
        # resumed from a prior run -- not real new work, but still read its real
        # ok/fail counts so the running totals stay accurate across a resume
        with open(part_path) as f:
            prior_rows = list(csv.DictReader(f))
        n_ok = sum(1 for r in prior_rows if r["mmgbsa_dg"])
        return chunk_idx, len(chunk), True, n_ok, len(prior_rows) - n_ok

    rows = []
    for compound_id, block in chunk:
        _, dg, err = _run_one_compound(compound_id, block, receptor_path, scratch_root, timeout_s)
        info = meta[compound_id]
        rows.append({"compound_id": compound_id, "conformation": conf, "cluster_id": info["cluster_id"],
                     "role": info["role"], "source": info["source"], "ad4_score": info["ad4_score"],
                     "mmgbsa_dg": dg if dg is not None else "", "error": err or ""})

    with open(part_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDNAMES)
        w.writeheader()
        w.writerows(rows)
    n_ok = sum(1 for r in rows if r["mmgbsa_dg"] != "")
    return chunk_idx, len(chunk), False, n_ok, len(chunk) - n_ok


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--conformation", required=True)
    ap.add_argument("--sdf", required=True, help="q5_mmgbsa_<conf>.sdf")
    ap.add_argument("--receptor", required=True, help="gmx-ready receptor PDB")
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--chunk_size", type=int, default=20,
                     help="compounds per work-chunk -- small, since each compound is its own "
                          "long-running subprocess and chunks are the resume granularity")
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--timeout", type=int, default=300, help="per-compound subprocess timeout, seconds")
    ap.add_argument("--scratch_root", default=None,
                     help="fast local scratch for per-compound throwaway dirs (default: system tmp -- "
                          "point this at genuinely local/fast storage if available, e.g. /dev/shm)")
    args = ap.parse_args()

    conf = args.conformation
    receptor_path = str(Path(args.receptor).resolve())
    if not Path(receptor_path).is_file():
        raise SystemExit(f"{receptor_path} not found")
    if not Path(args.sdf).is_file():
        raise SystemExit(f"{args.sdf} not found")
    if args.scratch_root:
        Path(args.scratch_root).mkdir(parents=True, exist_ok=True)

    out_dir = Path(args.out_dir) / f"mmgbsa_{conf}"
    out_dir.mkdir(parents=True, exist_ok=True)
    workers = args.workers or len(os.sched_getaffinity(0))

    log(f"Loading {conf} hits from {args.sdf} ...")
    ids, blocks, meta = load_hits(args.sdf)
    n_total = len(ids)
    log(f"{n_total:,} compounds loaded for {conf}, receptor={receptor_path}, workers={workers}")

    items = list(zip(ids, blocks))
    chunks = [items[i:i + args.chunk_size] for i in range(0, n_total, args.chunk_size)]
    n_chunks = len(chunks)
    have = sum(1 for i in range(n_chunks) if (out_dir / f"mmgbsa_part_{i:05d}.csv").exists())
    log(f"{n_chunks:,} work-chunks of up to {args.chunk_size}, {have:,} already done, {workers} workers")

    tasks = [(chunks[i], i, conf, receptor_path, meta, str(out_dir), args.scratch_root, args.timeout)
              for i in range(n_chunks)]

    # real periodic progress + ETA, not spammed -- roughly every ~1% of chunks
    log_every = max(1, n_chunks // 100)

    t0 = time.time()
    done_count = 0
    new_done_count = 0  # excludes resumed (instant, not real work) chunks -- rate/ETA use only this
    done_compounds = 0
    running_ok = 0
    running_fail = 0
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(_process_chunk, t): t[1] for t in tasks}
        for fut in as_completed(futures):
            chunk_idx, n_in_chunk, was_resumed, n_ok, n_fail = fut.result()
            done_count += 1
            done_compounds += n_in_chunk
            running_ok += n_ok
            running_fail += n_fail
            if not was_resumed:
                new_done_count += 1
            if done_count % log_every == 0 or done_count == n_chunks:
                elapsed = time.time() - t0
                rate = new_done_count / elapsed if elapsed > 0 else 0
                eta_s = (n_chunks - done_count) / rate if rate else float("inf")
                pct = 100 * done_count / n_chunks
                fail_pct = 100 * running_fail / max(1, done_compounds)
                # a high failure rate this early is worth a loud, explicit flag --
                # otherwise a systemic problem (wrong receptor, broken env, ...)
                # would silently run for hours before the final summary shows it
                flag = "  <<< HIGH FAILURE RATE, CHECK mmgbsa_failed_*.csv / logs" if (
                    done_compounds >= 20 and fail_pct > 30) else ""
                log(f"  chunk {done_count:,}/{n_chunks:,} ({pct:.1f}%) done -- "
                    f"~{done_compounds:,}/{n_total:,} compounds ({running_ok:,} ok, "
                    f"{running_fail:,} failed = {fail_pct:.1f}%) -- "
                    f"{elapsed / 60:.1f} min elapsed ({new_done_count:,} chunks actually computed "
                    f"this run, {done_count - new_done_count:,} resumed instantly) -- "
                    f"ETA {eta_s / 60:.1f} min ({eta_s / 3600:.1f} h){flag}")

    log("all chunks done -- merging ...")
    all_rows = []
    for i in range(n_chunks):
        with open(out_dir / f"mmgbsa_part_{i:05d}.csv") as f:
            all_rows.extend(csv.DictReader(f))

    n_ok = sum(1 for r in all_rows if r["mmgbsa_dg"])
    n_err = len(all_rows) - n_ok
    log(f"{n_ok:,} succeeded, {n_err:,} failed/errored")

    scored = [r for r in all_rows if r["mmgbsa_dg"]]
    scored.sort(key=lambda r: float(r["mmgbsa_dg"]))
    final_path = out_dir / f"mmgbsa_results_{conf}.csv"
    with open(final_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["rank", "compound_id", "mmgbsa_dg", "ad4_score", "cluster_id", "role", "source"])
        for rank, r in enumerate(scored, start=1):
            w.writerow([rank, r["compound_id"], r["mmgbsa_dg"], r["ad4_score"],
                       r["cluster_id"], r["role"], r["source"]])
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
        w.writerow(["mmgbsa_rank", "compound_id", "mmgbsa_dg", "ad4_score", "cluster_id", "role",
                   "source", "ad4_rank_within_q5_set"])
        ad4_sorted_ids = sorted(meta, key=lambda cid: float(meta[cid]["ad4_score"]))
        ad4_rank = {cid: i + 1 for i, cid in enumerate(ad4_sorted_ids)}
        for rank, r in enumerate(scored[:top_n], start=1):
            w.writerow([rank, r["compound_id"], r["mmgbsa_dg"], r["ad4_score"], r["cluster_id"],
                       r["role"], r["source"], ad4_rank[r["compound_id"]]])
    log(f"Wrote {sanity_path} (top {top_n} MM-GBSA hits with their AD4 score and AD4-only rank "
        f"within this Q5 set, for a quick sanity cross-check)")


if __name__ == "__main__":
    main()
