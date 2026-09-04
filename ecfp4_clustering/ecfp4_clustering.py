#!/usr/bin/env python3
"""
Author: Anna Korda

Taylor-Butina clustering of each conformation's final hit set (ECFP4 / Morgan
r=2 fingerprints, Tanimoto similarity) at several thresholds, reporting how many
clusters (and how big) each threshold produces. Runs on apply_score_cutoff.py's
already-built final_hits_<conf>_cutoff<C>.csv + _partNN.sdf files -- reads
structures only to fingerprint them, nothing is ever rewritten.

A dense N x N similarity matrix isn't stored anywhere (for the largest set,
~106k compounds, that would be ~45GB) -- neighbor pairs above the loosest
requested threshold are found via one upper-triangle BulkTanimotoSimilarity
pass (chunked/parallel/resumable, same pattern as the rest of this repo),
written to disk as compact (i, j, sim) edges, then Taylor-Butina itself is run
per threshold from a CSR-style neighbor index built with numpy (fast, no
python-object-per-edge overhead).

Usage:
    python ecfp4_clustering.py --conformation c1 \\
        --vs_results_dir /users/gpcr/annak/ultra-large/vs_results \\
        --out_dir /users/gpcr/annak/ultra-large/analysis_results/ecfp4_clustering \\
        --cutoff -7.0 --thresholds 0.90,0.85,0.80,0.75,0.70 --workers 32
"""
import argparse
import csv
import os
import shutil
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from rdkit.DataStructs import BulkTanimotoSimilarity

RDLogger.DisableLog("rdApp.*")


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


def load_hit_ids_and_blocks(vs_results_dir, conf, cutoff):
    hit_dir = Path(vs_results_dir) / f"final_hits_{conf}_cutoff{cutoff:g}"
    csv_path = hit_dir / f"final_hits_{conf}_cutoff{cutoff:g}.csv"
    if not csv_path.is_file():
        raise SystemExit(f"{csv_path} not found -- run apply_score_cutoff.py for {conf} first")
    with open(csv_path) as f:
        ids = [row["compound_id"] for row in csv.DictReader(f)]
    part_files = sorted(hit_dir.glob(f"final_hits_{conf}_cutoff{cutoff:g}_part*.sdf"))
    if not part_files:
        raise SystemExit(f"No SDF parts found in {hit_dir}")
    blocks = []
    for pf in part_files:
        blocks.extend(_split_sdf_blocks(pf.read_text()))
    if len(blocks) != len(ids):
        raise SystemExit(f"{conf}: {len(ids)} CSV rows but {len(blocks)} SDF blocks -- "
                          f"expected the same rank order (apply_score_cutoff.py's own guarantee), "
                          f"refusing to guess a correspondence")
    return ids, blocks


def _fingerprint_chunk(args):
    blocks_chunk, nbits = args
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=nbits)
    fps, n_errors = [], 0
    for block in blocks_chunk:
        mol = Chem.MolFromMolBlock(block)
        if mol is None:
            fps.append(None)
            n_errors += 1
            continue
        fps.append(gen.GetFingerprint(mol))
    return fps, n_errors


def compute_fingerprints(blocks, nbits, workers, chunk_size=500):
    chunks = [blocks[i:i + chunk_size] for i in range(0, len(blocks), chunk_size)]
    fps = [None] * len(blocks)
    n_errors = 0
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(_fingerprint_chunk, (chunk, nbits)): i for i, chunk in enumerate(chunks)}
        done = 0
        for fut in as_completed(futures):
            i = futures[fut]
            chunk_fps, errs = fut.result()
            fps[i * chunk_size:i * chunk_size + len(chunk_fps)] = chunk_fps
            n_errors += errs
            done += 1
            if done % 20 == 0 or done == len(chunks):
                log(f"  fingerprinting: {done}/{len(chunks)} chunks done")
    return fps, n_errors


def _find_neighbors_chunk(args):
    row_start, row_end, parts_dir_str, min_threshold = args
    part_path = Path(parts_dir_str) / f"edges_{row_start:07d}.npz"
    if part_path.exists():
        return row_start, row_end, -1

    global _FPS, _HAS_NONE
    ii, jj, ss = [], [], []
    n_fps = len(_FPS)
    for i in range(row_start, row_end):
        fp_i = _FPS[i]
        if fp_i is None or i + 1 >= n_fps:
            continue
        others = _FPS[i + 1:]
        if _HAS_NONE:
            valid_idx = [k for k, fp in enumerate(others) if fp is not None]
            if not valid_idx:
                continue
            others = [others[k] for k in valid_idx]
        else:
            valid_idx = None
        sims = BulkTanimotoSimilarity(fp_i, others)
        for k, sim in enumerate(sims):
            if sim >= min_threshold:
                j = (i + 1 + valid_idx[k]) if _HAS_NONE else (i + 1 + k)
                ii.append(i)
                jj.append(j)
                ss.append(sim)

    np.savez(part_path, i=np.array(ii, dtype=np.int32), j=np.array(jj, dtype=np.int32),
              s=np.array(ss, dtype=np.float32))
    return row_start, row_end, len(ii)


def _init_neighbor_worker(fps):
    global _FPS, _HAS_NONE
    _FPS = fps
    _HAS_NONE = any(fp is None for fp in fps)


def find_all_neighbor_edges(fps, min_threshold, out_dir, conf, workers, row_chunk=300):
    parts_dir = out_dir / f".edges_{conf}"
    parts_dir.mkdir(parents=True, exist_ok=True)
    n = len(fps)
    row_starts = list(range(0, n, row_chunk))
    tasks = [(rs, min(rs + row_chunk, n), str(parts_dir), min_threshold) for rs in row_starts]

    have = sum(1 for rs in row_starts if (parts_dir / f"edges_{rs:07d}.npz").exists())
    log(f"  neighbor search: {len(tasks)} row-chunks of up to {row_chunk} rows, "
        f"{have} already done, {workers} workers")

    t0 = time.time()
    with ProcessPoolExecutor(max_workers=workers, initializer=_init_neighbor_worker,
                              initargs=(fps,)) as pool:
        futures = [pool.submit(_find_neighbors_chunk, t) for t in tasks]
        done = 0
        for fut in as_completed(futures):
            fut.result()
            done += 1
            if done % 50 == 0 or done == len(tasks):
                elapsed = time.time() - t0
                rate = done / elapsed if elapsed > 0 else 0
                eta = (len(tasks) - done) / rate if rate else float("inf")
                log(f"  neighbor search: {done}/{len(tasks)} row-chunks done "
                    f"({elapsed:.0f}s elapsed, ETA {eta:.0f}s)")

    all_i, all_j, all_s = [], [], []
    for rs in row_starts:
        z = np.load(parts_dir / f"edges_{rs:07d}.npz")
        all_i.append(z["i"])
        all_j.append(z["j"])
        all_s.append(z["s"])
    edge_i = np.concatenate(all_i) if all_i else np.array([], dtype=np.int32)
    edge_j = np.concatenate(all_j) if all_j else np.array([], dtype=np.int32)
    edge_s = np.concatenate(all_s) if all_s else np.array([], dtype=np.float32)
    log(f"  {len(edge_i):,} total edges with similarity >= {min_threshold}")
    return edge_i, edge_j, edge_s


def build_csr_adjacency(n, edge_i, edge_j):
    """Symmetric CSR neighbor index from one-directional (i<j) edges."""
    src = np.concatenate([edge_i, edge_j])
    dst = np.concatenate([edge_j, edge_i])
    order = np.argsort(src, kind="stable")
    src, dst = src[order], dst[order]
    counts = np.bincount(src, minlength=n)
    indptr = np.zeros(n + 1, dtype=np.int64)
    np.cumsum(counts, out=indptr[1:])
    return indptr, dst


def taylor_butina(n, indptr, indices):
    neighbor_counts = np.diff(indptr)
    order = np.argsort(-neighbor_counts, kind="stable")
    assigned = np.zeros(n, dtype=bool)
    cluster_id = np.full(n, -1, dtype=np.int64)
    next_cluster = 0

    for node in order:
        if assigned[node]:
            continue
        neighbors = indices[indptr[node]:indptr[node + 1]]
        unassigned_neighbors = neighbors[~assigned[neighbors]]
        cluster_id[node] = next_cluster
        assigned[node] = True
        if len(unassigned_neighbors) > 0:
            cluster_id[unassigned_neighbors] = next_cluster
            assigned[unassigned_neighbors] = True
        next_cluster += 1

    return cluster_id, next_cluster


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--conformation", required=True)
    ap.add_argument("--vs_results_dir", required=True)
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--cutoff", type=float, default=-7.0)
    ap.add_argument("--thresholds", default="0.90,0.85,0.80,0.75,0.70")
    ap.add_argument("--nbits", type=int, default=2048)
    ap.add_argument("--workers", type=int, default=None)
    args = ap.parse_args()

    workers = args.workers or len(os.sched_getaffinity(0))
    thresholds = sorted((float(t) for t in args.thresholds.split(",")), reverse=True)
    min_threshold = min(thresholds)

    conf = args.conformation
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    log(f"Loading {conf} hits (cutoff {args.cutoff:g}) ...")
    ids, blocks = load_hit_ids_and_blocks(args.vs_results_dir, conf, args.cutoff)
    n = len(ids)
    log(f"{n:,} compounds loaded for {conf}")

    log(f"Computing ECFP4 (Morgan r=2, {args.nbits} bits) fingerprints, {workers} workers ...")
    fps, n_fp_errors = compute_fingerprints(blocks, args.nbits, workers)
    if n_fp_errors:
        log(f"  WARNING: {n_fp_errors} molecule(s) failed to parse -- excluded from clustering")

    edge_i, edge_j, edge_s = find_all_neighbor_edges(fps, min_threshold, out_dir, conf, workers)

    summary_rows = []
    for threshold in thresholds:
        mask = edge_s >= threshold
        t_i, t_j = edge_i[mask], edge_j[mask]
        indptr, indices = build_csr_adjacency(n, t_i, t_j)
        cluster_id, n_clusters = taylor_butina(n, indptr, indices)

        sizes = np.bincount(cluster_id[cluster_id >= 0])
        n_singletons = int((sizes == 1).sum())
        log(f"  threshold {threshold:.2f}: {n_clusters:,} clusters "
            f"({n_singletons:,} singletons, largest={sizes.max()})")

        sizes_path = out_dir / f"cluster_sizes_{conf}_{threshold:.2f}.csv"
        with open(sizes_path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["cluster_id", "size"])
            for cid, size in enumerate(sizes):
                w.writerow([cid, int(size)])

        assign_path = out_dir / f"cluster_assignment_{conf}_{threshold:.2f}.csv"
        with open(assign_path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["compound_id", "cluster_id"])
            for compound_id, cid in zip(ids, cluster_id):
                w.writerow([compound_id, int(cid)])

        summary_rows.append({
            "conformation": conf, "threshold": f"{threshold:.2f}", "n_compounds": n,
            "n_clusters": n_clusters, "n_singletons": n_singletons,
            "largest_cluster_size": int(sizes.max()), "mean_cluster_size": f"{sizes.mean():.3f}",
            "median_cluster_size": f"{float(np.median(sizes)):.1f}",
        })

    summary_path = out_dir / f"cluster_summary_{conf}.csv"
    with open(summary_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["conformation", "threshold", "n_compounds", "n_clusters",
                                           "n_singletons", "largest_cluster_size",
                                           "mean_cluster_size", "median_cluster_size"])
        w.writeheader()
        w.writerows(summary_rows)
    log(f"Wrote {summary_path}")

    shutil.rmtree(out_dir / f".edges_{conf}", ignore_errors=True)

    run_summary_path = out_dir / f"run_summary_{conf}.txt"
    with open(run_summary_path, "w") as f:
        f.write(f"ECFP4/Taylor-Butina clustering -- conformation {conf}, cutoff {args.cutoff:g}\n")
        f.write(f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"{n:,} compounds ({n_fp_errors} fingerprint errors)\n")
        f.write(f"Fingerprints: Morgan r=2 (ECFP4), {args.nbits} bits\n")
        f.write(f"Thresholds tested: {', '.join(f'{t:.2f}' for t in thresholds)}\n\n")
        for row in summary_rows:
            f.write(f"  threshold {row['threshold']}: {row['n_clusters']:,} clusters, "
                    f"{row['n_singletons']:,} singletons, largest={row['largest_cluster_size']}, "
                    f"mean size={row['mean_cluster_size']}\n")
    log(f"Wrote {run_summary_path}")


if __name__ == "__main__":
    main()
