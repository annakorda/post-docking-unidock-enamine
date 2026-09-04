#!/usr/bin/env python3
"""
Author: Anna Korda

Cross-conformation "SNN" (Similarity to Nearest Neighbor) between the three
final hit sets -- the standard metric from the MOSES molecular-generation
benchmark (Polykovskiy et al. 2020, "Molecular Sets (MOSES): A Benchmarking
Platform for Molecular Generation Models," Front. Pharmacol.):

    SNN(A, B) = mean over m in A of [ max over m' in B of Tanimoto(m, m') ]

i.e. for every compound in set A, its best (nearest-neighbor) ECFP4 Tanimoto
match anywhere in set B, averaged over all of A. Non-symmetric by construction
(same convention as hit_overlap.py's matrix) since it's an average over A's own
compounds, not B's.

This needs no dense matrix: for each compound in A we only ever keep the max
similarity against B (BulkTanimotoSimilarity against all of B, then discard
the row) -- memory cost is O(|B|) fingerprints, not O(|A|*|B|) similarities.
Runs all ordered (A, B) pairs among the requested conformations, including the
diagonal (A, A) with self-matches excluded -- a free bonus internal-diversity
number, not part of the cross-set comparison Anna asked for, clearly labeled.

Usage:
    python snn_similarity.py --vs_results_dir /users/gpcr/annak/ultra-large/vs_results \\
        --out_dir /users/gpcr/annak/ultra-large/analysis_results/snn_similarity \\
        --conformations c1,ref1,c5 --cutoff -7.0 --workers 32
"""
import argparse
import csv
import os
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
            if done % 50 == 0 or done == len(chunks):
                log(f"    fingerprinting: {done}/{len(chunks)} chunks done")
    return fps, n_errors


def _init_snn_worker(fps_b, ids_b, self_pair):
    global _FPS_B, _IDS_B, _SELF_PAIR
    _FPS_B = fps_b
    _IDS_B = ids_b
    _SELF_PAIR = self_pair


def _nn_chunk(args):
    row_start, row_end, ids_a_chunk, fps_a_chunk = args
    global _FPS_B, _IDS_B, _SELF_PAIR
    results = []
    for offset, (compound_id_a, fp_a) in enumerate(zip(ids_a_chunk, fps_a_chunk)):
        global_i = row_start + offset
        if fp_a is None:
            results.append((compound_id_a, None, 0.0))
            continue
        sims = np.array(BulkTanimotoSimilarity(fp_a, _FPS_B), dtype=np.float32)
        if _SELF_PAIR:
            sims[global_i] = -1.0
        best_j = int(np.argmax(sims))
        results.append((compound_id_a, _IDS_B[best_j], float(sims[best_j])))
    return results


def compute_snn_pair(name_a, ids_a, fps_a, name_b, ids_b, fps_b, workers, out_dir, chunk_size=500):
    self_pair = (name_a == name_b)
    chunks = []
    for i in range(0, len(ids_a), chunk_size):
        chunks.append((i, min(i + chunk_size, len(ids_a)), ids_a[i:i + chunk_size], fps_a[i:i + chunk_size]))

    all_results = []
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=workers, initializer=_init_snn_worker,
                              initargs=(fps_b, ids_b, self_pair)) as pool:
        futures = {pool.submit(_nn_chunk, c): idx for idx, c in enumerate(chunks)}
        ordered = [None] * len(chunks)
        done = 0
        for fut in as_completed(futures):
            idx = futures[fut]
            ordered[idx] = fut.result()
            done += 1
            if done % 20 == 0 or done == len(chunks):
                elapsed = time.time() - t0
                rate = done / elapsed if elapsed > 0 else 0
                eta = (len(chunks) - done) / rate if rate else float("inf")
                log(f"    {name_a} -> {name_b}: {done}/{len(chunks)} chunks done "
                    f"({elapsed:.0f}s elapsed, ETA {eta:.0f}s)")
    for chunk_results in ordered:
        all_results.extend(chunk_results)

    pairs_path = out_dir / f"snn_pairs_{name_a}_to_{name_b}.csv"
    with open(pairs_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["compound_id", "nearest_neighbor_compound_id", "tanimoto_similarity"])
        valid_sims = []
        for compound_id, nn_id, sim in all_results:
            w.writerow([compound_id, nn_id if nn_id is not None else "", f"{sim:.4f}" if nn_id else ""])
            if nn_id is not None:
                valid_sims.append(sim)

    snn_score = sum(valid_sims) / len(valid_sims) if valid_sims else 0.0
    log(f"  SNN({name_a} -> {name_b}) = {snn_score:.4f}  (wrote {pairs_path.name})")
    return snn_score


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vs_results_dir", required=True)
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--conformations", default="c1,ref1,c5")
    ap.add_argument("--cutoff", type=float, default=-7.0)
    ap.add_argument("--nbits", type=int, default=2048)
    ap.add_argument("--workers", type=int, default=None)
    args = ap.parse_args()

    workers = args.workers or len(os.sched_getaffinity(0))
    confs = args.conformations.split(",")
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    ids, fps = {}, {}
    for conf in confs:
        log(f"Loading {conf} hits (cutoff {args.cutoff:g}) ...")
        conf_ids, blocks = load_hit_ids_and_blocks(args.vs_results_dir, conf, args.cutoff)
        log(f"  {len(conf_ids):,} compounds -- computing ECFP4 fingerprints, {workers} workers ...")
        conf_fps, n_errors = compute_fingerprints(blocks, args.nbits, workers)
        if n_errors:
            log(f"  WARNING: {n_errors} molecule(s) failed to parse")
        ids[conf], fps[conf] = conf_ids, conf_fps

    matrix = {a: {} for a in confs}
    for a in confs:
        for b in confs:
            log(f"Computing SNN({a} -> {b}) ...")
            matrix[a][b] = compute_snn_pair(a, ids[a], fps[a], b, ids[b], fps[b], workers, out_dir)

    matrix_path = out_dir / "snn_matrix.csv"
    with open(matrix_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([""] + confs)
        for a in confs:
            w.writerow([a] + [f"{matrix[a][b]:.4f}" for b in confs])
    log(f"Wrote {matrix_path}")

    summary_path = out_dir / "run_summary.txt"
    with open(summary_path, "w") as f:
        f.write(f"SNN (Similarity to Nearest Neighbor) -- cutoff {args.cutoff:g}\n")
        f.write(f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("SNN(A,B) = mean over compounds in A of their best ECFP4 Tanimoto match in B\n")
        f.write("Non-symmetric: SNN(A,B) != SNN(B,A) in general (see MOSES benchmark, "
                "Polykovskiy et al. 2020)\n")
        f.write("Diagonal entries (A,A) exclude self-matches -- internal diversity, not "
                "part of the cross-set comparison\n\n")
        for conf in confs:
            f.write(f"{conf}: {len(ids[conf]):,} compounds\n")
        f.write("\n")
        for a in confs:
            for b in confs:
                f.write(f"  SNN({a} -> {b}) = {matrix[a][b]:.4f}\n")
    log(f"Wrote {summary_path}")


if __name__ == "__main__":
    main()
