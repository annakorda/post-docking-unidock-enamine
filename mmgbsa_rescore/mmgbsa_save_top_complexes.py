#!/usr/bin/env python3
"""
Author: Anna Korda

Run AFTER mmgbsa_rescore.py has scored a conformation's full hit set (that
script deliberately keeps only the delta G, deleting every compound's
minimized complex to avoid ~850GB / ~1.6M inodes across ~230k compounds).
This script re-runs unigbsa-pipeline for just the top --top_pct percent of
that already-ranked result (default 1%) and this time keeps the minimized
complex.pdb -- a real, small, bounded re-run (e.g. ~1,064 compounds for c1's
106,380 -- ~15 min on 64 workers), not a second full-scale pass.

Same chunked/resumable/parallel/progress-logging pattern as mmgbsa_rescore.py.

Usage:
    python mmgbsa_save_top_complexes.py --conformation c1 \\
        --vs_results_dir /users/gpcr/annak/ultra-large/vs_results \\
        --receptor /users/gpcr/annak/ultra-large/receptor_prep/c1_gmxready.pdb \\
        --mmgbsa_out_dir /users/gpcr/annak/ultra-large/analysis_results/mmgbsa_rescore \\
        --cutoff -7.0 --top_pct 1.0 --workers 64
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


def load_top_slice(vs_results_dir, mmgbsa_out_dir, conf, cutoff, top_pct):
    results_csv = Path(mmgbsa_out_dir) / f"mmgbsa_{conf}" / f"mmgbsa_results_{conf}.csv"
    if not results_csv.is_file():
        raise SystemExit(f"{results_csv} not found -- run mmgbsa_rescore.py for {conf} first")
    with open(results_csv) as f:
        ranked_rows = list(csv.DictReader(f))
    n_top = max(1, math.ceil(len(ranked_rows) * top_pct / 100.0))
    top_rows = ranked_rows[:n_top]
    top_ids = {r["compound_id"] for r in top_rows}
    log(f"Top {top_pct:g}% of {len(ranked_rows):,} scored compounds = {n_top:,} compounds")

    hit_dir = Path(vs_results_dir) / f"final_hits_{conf}_cutoff{cutoff:g}"
    csv_path = hit_dir / f"final_hits_{conf}_cutoff{cutoff:g}.csv"
    if not csv_path.is_file():
        raise SystemExit(f"{csv_path} not found")
    with open(csv_path) as f:
        ids = [r["compound_id"] for r in csv.DictReader(f)]
    part_files = sorted(hit_dir.glob(f"final_hits_{conf}_cutoff{cutoff:g}_part*.sdf"))
    if not part_files:
        raise SystemExit(f"No SDF parts found in {hit_dir}")
    blocks = []
    for pf in part_files:
        blocks.extend(_split_sdf_blocks(pf.read_text()))
    if len(blocks) != len(ids):
        raise SystemExit(f"{conf}: {len(ids)} CSV rows but {len(blocks)} SDF blocks -- "
                          f"expected the same rank order, refusing to guess a correspondence")
    id_to_block = {cid: block for cid, block in zip(ids, blocks) if cid in top_ids}
    missing = top_ids - set(id_to_block)
    if missing:
        raise SystemExit(f"{len(missing)} top-ranked compound(s) not found in {csv_path} -- "
                          f"mmgbsa_results and final_hits_cutoff appear out of sync, refusing to guess")

    return [(r["compound_id"], id_to_block[r["compound_id"]], r["mmgbsa_dg"], r["ad4_score"])
            for r in top_rows]


def _run_and_save_complex(compound_id, block, receptor_path, scratch_root, timeout_s, complex_out_path):
    work_dir = tempfile.mkdtemp(prefix=f"gbsa_{compound_id[:40]}_", dir=scratch_root)
    try:
        ligand_path = Path(work_dir) / f"{compound_id}.sdf"
        ligand_path.write_text(block)
        result_csv = Path(work_dir) / "result.csv"
        # no --verbose: complex.pdb is already written by default -- --verbose
        # only adds extra large intermediate files (COM.prmtop etc.) not needed here
        proc = subprocess.run(
            ["unigbsa-pipeline", "-i", str(receptor_path), "-l", str(ligand_path),
             "-o", str(result_csv), "-nt", "1"],
            cwd=work_dir, capture_output=True, text=True, timeout=timeout_s,
        )
        if not result_csv.is_file():
            return compound_id, None, f"no_output_rc{proc.returncode}: {proc.stderr[-300:]}"
        complex_pdb = Path(work_dir) / compound_id / "complex.pdb"
        if not complex_pdb.is_file():
            return compound_id, None, "complex.pdb not found in run output"
        shutil.copy2(complex_pdb, complex_out_path)
        return compound_id, str(complex_out_path), None
    except subprocess.TimeoutExpired:
        return compound_id, None, "timeout"
    except Exception as e:
        return compound_id, None, f"exception: {e}"
    finally:
        shutil.rmtree(work_dir, ignore_errors=True)


def _process_chunk(args):
    chunk, chunk_idx, receptor_path, out_dir_str, scratch_root, timeout_s = args
    out_dir = Path(out_dir_str)
    complexes_dir = out_dir / "complexes"
    part_path = out_dir / f"complexes_part_{chunk_idx:05d}.csv"
    if part_path.exists():
        return chunk_idx, len(chunk)

    rows = []
    for compound_id, block, mmgbsa_dg, ad4_score in chunk:
        complex_out_path = complexes_dir / f"{compound_id}_complex.pdb"
        _, saved_path, err = _run_and_save_complex(compound_id, block, receptor_path, scratch_root,
                                                     timeout_s, complex_out_path)
        rows.append({"compound_id": compound_id, "mmgbsa_dg": mmgbsa_dg, "ad4_score": ad4_score,
                     "complex_pdb": saved_path or "", "error": err or ""})

    with open(part_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["compound_id", "mmgbsa_dg", "ad4_score", "complex_pdb", "error"])
        w.writeheader()
        w.writerows(rows)
    return chunk_idx, len(chunk)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--conformation", required=True)
    ap.add_argument("--vs_results_dir", required=True)
    ap.add_argument("--receptor", required=True)
    ap.add_argument("--mmgbsa_out_dir", required=True, help="mmgbsa_rescore.py's --out_dir")
    ap.add_argument("--out_dir", default=None,
                     help="default: <mmgbsa_out_dir>/mmgbsa_<conf>/top_complexes")
    ap.add_argument("--cutoff", type=float, default=-7.0)
    ap.add_argument("--top_pct", type=float, default=1.0)
    ap.add_argument("--chunk_size", type=int, default=10)
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--timeout", type=int, default=300)
    ap.add_argument("--scratch_root", default=None)
    args = ap.parse_args()

    conf = args.conformation
    receptor_path = str(Path(args.receptor).resolve())
    if not Path(receptor_path).is_file():
        raise SystemExit(f"{receptor_path} not found")
    if args.scratch_root:
        Path(args.scratch_root).mkdir(parents=True, exist_ok=True)

    out_dir = Path(args.out_dir) if args.out_dir else \
        Path(args.mmgbsa_out_dir) / f"mmgbsa_{conf}" / "top_complexes"
    (out_dir / "complexes").mkdir(parents=True, exist_ok=True)
    workers = args.workers or len(os.sched_getaffinity(0))

    log(f"Loading top {args.top_pct:g}% of {conf} (cutoff {args.cutoff:g}) ...")
    items = load_top_slice(args.vs_results_dir, args.mmgbsa_out_dir, conf, args.cutoff, args.top_pct)
    n_total = len(items)
    log(f"{n_total:,} compounds selected, receptor={receptor_path}, workers={workers}")

    chunks = [items[i:i + args.chunk_size] for i in range(0, n_total, args.chunk_size)]
    n_chunks = len(chunks)
    have = sum(1 for i in range(n_chunks) if (out_dir / f"complexes_part_{i:05d}.csv").exists())
    log(f"{n_chunks:,} work-chunks of up to {args.chunk_size}, {have:,} already done, {workers} workers")

    tasks = [(chunks[i], i, receptor_path, str(out_dir), args.scratch_root, args.timeout)
              for i in range(n_chunks)]
    log_every = max(1, n_chunks // 100)

    t0 = time.time()
    done_count = 0
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(_process_chunk, t): t[1] for t in tasks}
        for fut in as_completed(futures):
            chunk_idx, n_in_chunk = fut.result()
            done_count += 1
            if done_count % log_every == 0 or done_count == n_chunks:
                elapsed = time.time() - t0
                rate = done_count / elapsed if elapsed > 0 else 0
                eta_s = (n_chunks - done_count) / rate if rate else float("inf")
                pct = 100 * done_count / n_chunks
                log(f"  chunk {done_count:,}/{n_chunks:,} ({pct:.1f}%) done -- "
                    f"{elapsed / 60:.1f} min elapsed, ETA {eta_s / 60:.1f} min")

    log("all chunks done -- merging manifest ...")
    all_rows = []
    for i in range(n_chunks):
        with open(out_dir / f"complexes_part_{i:05d}.csv") as f:
            all_rows.extend(csv.DictReader(f))

    n_ok = sum(1 for r in all_rows if r["complex_pdb"])
    log(f"{n_ok:,} complex PDBs saved, {len(all_rows) - n_ok:,} failed")

    manifest_path = out_dir / f"top_complexes_manifest_{conf}.csv"
    with open(manifest_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["rank", "compound_id", "mmgbsa_dg", "ad4_score", "complex_pdb", "error"])
        for rank, r in enumerate(all_rows, start=1):
            w.writerow([rank, r["compound_id"], r["mmgbsa_dg"], r["ad4_score"], r["complex_pdb"], r["error"]])
    log(f"Wrote {manifest_path}")


if __name__ == "__main__":
    main()
