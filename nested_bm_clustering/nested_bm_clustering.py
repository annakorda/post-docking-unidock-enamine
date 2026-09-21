#!/usr/bin/env python3
"""
Author: Anna Korda

Two-pass clustering, replacing ecfp4_clustering.py's single flat Butina pass:

Pass 1 -- Bemis-Murcko, deterministic, no threshold. Every compound's generic
scaffold (ring systems + linkers only, atom types abstracted away via RDKit's
MakeScaffoldGeneric -- verified: aromatic vs aliphatic ring, N/O/F/Cl
substitution etc. don't matter, only ring/linker topology does) is computed
once. Two compounds land in the same BM group only if that skeleton's
canonical SMILES is identical.

Pass 2 -- Tanimoto/Butina, but only INSIDE each BM group, and only where a
group is large enough that sub-clustering could matter (--min_group_size,
default 10 -- below that, the whole group is already one coherent scaffold
family, treated as a single cluster with no Butina run on it at all). For
each large group, real all-pairs Tanimoto is computed once (in memory --
even the largest realistic group is far smaller than the full compound set,
so this doesn't need the disk-chunked approach the full flat pass needed),
then Taylor-Butina is run per --thresholds value using that group's own
similarity data. Total cluster count at a given Tc = the number of untouched
small BM groups, plus the sum of Butina sub-clusters produced inside every
large BM group at that Tc.

Runs on apply_score_cutoff.py's cutoff CSV/SDF (point --vs_results_dir at
the deduped set from dedup_stereoisomers.py -- that's the real intended
input here, redundant stereoisomers/tautomers should already be removed
before this runs, not re-discovered as "BM group size 1 vs many").

Usage:
    python nested_bm_clustering.py --conformation c1 \\
        --vs_results_dir /home/annie/Desktop/alphavs/revision_2026/ultra-large-results/dedup \\
        --out_dir /home/annie/Desktop/alphavs/revision_2026/analysis/.../nested_bm_clustering \\
        --cutoff -7.0 --thresholds 0.90,...,0.35 --min_group_size 10 --workers 8
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
from rdkit.Chem.Scaffolds import MurckoScaffold
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


def load_hits(vs_results_dir, conf, cutoff, dedup=True):
    stem = f"final_hits_{conf}_cutoff{cutoff:g}" + ("_dedup" if dedup else "")
    hit_dir = Path(vs_results_dir) / stem
    csv_path = hit_dir / f"{stem}.csv"
    if not csv_path.is_file():
        raise SystemExit(f"{csv_path} not found")
    with open(csv_path) as f:
        ids = [r["compound_id"] for r in csv.DictReader(f)]
    if dedup:
        # dedup_stereoisomers.py writes one plain .sdf, not _partNN.sdf
        sdf_path = hit_dir / f"{stem}.sdf"
        if not sdf_path.is_file():
            raise SystemExit(f"{sdf_path} not found")
        blocks = _split_sdf_blocks(sdf_path.read_text())
    else:
        part_files = sorted(hit_dir.glob(f"{stem}_part*.sdf"))
        blocks = []
        for pf in part_files:
            blocks.extend(_split_sdf_blocks(pf.read_text()))
    if len(blocks) != len(ids):
        raise SystemExit(f"{conf}: {len(ids)} CSV rows but {len(blocks)} SDF blocks -- "
                          f"expected the same order, refusing to guess a correspondence")
    return ids, blocks


def _fp_and_scaffold_chunk(args):
    blocks_chunk, nbits = args
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=nbits)
    fps, scaffolds, n_errors = [], [], 0
    for block in blocks_chunk:
        mol = Chem.MolFromMolBlock(block)
        if mol is None:
            fps.append(None)
            scaffolds.append(None)
            n_errors += 1
            continue
        fps.append(gen.GetFingerprint(mol))
        try:
            scaffold = MurckoScaffold.GetScaffoldForMol(mol)
            generic = MurckoScaffold.MakeScaffoldGeneric(scaffold)
            scaffolds.append(Chem.MolToSmiles(generic))
        except Exception:
            scaffolds.append(None)
            n_errors += 1
    return fps, scaffolds, n_errors


def compute_fps_and_scaffolds(blocks, nbits, workers, chunk_size=500):
    chunks = [blocks[i:i + chunk_size] for i in range(0, len(blocks), chunk_size)]
    fps = [None] * len(blocks)
    scaffolds = [None] * len(blocks)
    n_errors = 0
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(_fp_and_scaffold_chunk, (chunk, nbits)): i for i, chunk in enumerate(chunks)}
        done = 0
        for fut in as_completed(futures):
            i = futures[fut]
            chunk_fps, chunk_scaffolds, errs = fut.result()
            fps[i * chunk_size:i * chunk_size + len(chunk_fps)] = chunk_fps
            scaffolds[i * chunk_size:i * chunk_size + len(chunk_scaffolds)] = chunk_scaffolds
            n_errors += errs
            done += 1
            if done % 20 == 0 or done == len(chunks):
                log(f"  fingerprints + BM scaffolds: {done}/{len(chunks)} chunks done")
    return fps, scaffolds, n_errors


def group_by_scaffold(ids, scaffolds):
    groups = {}
    for i, (cid, scaf) in enumerate(zip(ids, scaffolds)):
        if scaf is None:
            scaf = f"__ERROR__{cid}"  # each error compound is its own singleton group
        groups.setdefault(scaf, []).append(i)
    return groups


def build_csr_adjacency(n, edge_i, edge_j):
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


def _process_large_group(args):
    """One large BM group's full in-memory pairwise Tanimoto + per-threshold
    Butina. Returns {threshold: local_cluster_id_array} for this group's
    members (local indices 0..G-1, caller maps back to global compound_ids)."""
    group_fps, thresholds = args
    g = len(group_fps)
    min_threshold = min(thresholds)

    edge_i, edge_j, edge_s = [], [], []
    for i in range(g - 1):
        sims = BulkTanimotoSimilarity(group_fps[i], group_fps[i + 1:])
        for k, sim in enumerate(sims):
            if sim >= min_threshold:
                edge_i.append(i)
                edge_j.append(i + 1 + k)
                edge_s.append(sim)
    edge_i = np.array(edge_i, dtype=np.int64)
    edge_j = np.array(edge_j, dtype=np.int64)
    edge_s = np.array(edge_s, dtype=np.float32)

    result = {}
    for t in thresholds:
        mask = edge_s >= t
        indptr, indices = build_csr_adjacency(g, edge_i[mask], edge_j[mask])
        cluster_id, _ = taylor_butina(g, indptr, indices)
        result[t] = cluster_id
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--conformation", required=True)
    ap.add_argument("--vs_results_dir", required=True)
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--cutoff", type=float, default=-7.0)
    ap.add_argument("--no_dedup_input", action="store_true",
                     help="input is the plain (non-deduplicated) final_hits set, not dedup_stereoisomers.py's output")
    ap.add_argument("--thresholds", default="0.90,0.85,0.80,0.75,0.70,0.65,0.60,0.55,0.50,0.45,0.40,0.35")
    ap.add_argument("--min_group_size", type=int, default=10,
                     help="BM groups smaller than this are treated as one cluster with no Butina run on them")
    ap.add_argument("--nbits", type=int, default=2048)
    ap.add_argument("--workers", type=int, default=None)
    args = ap.parse_args()

    conf = args.conformation
    workers = args.workers or len(os.sched_getaffinity(0))
    thresholds = sorted((float(t) for t in args.thresholds.split(",")), reverse=True)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    log(f"Loading {conf} hits (cutoff {args.cutoff:g}) ...")
    ids, blocks = load_hits(args.vs_results_dir, conf, args.cutoff, dedup=not args.no_dedup_input)
    n = len(ids)
    log(f"{n:,} compounds loaded for {conf}")

    log(f"Computing ECFP4 fingerprints + Bemis-Murcko generic scaffolds, {workers} workers ...")
    fps, scaffolds, n_errors = compute_fps_and_scaffolds(blocks, args.nbits, workers)
    if n_errors:
        log(f"  WARNING: {n_errors} molecule(s) failed to parse/scaffold -- each treated as its own singleton group")

    groups = group_by_scaffold(ids, scaffolds)
    group_sizes = sorted((len(v) for v in groups.values()), reverse=True)
    n_groups = len(groups)
    n_large = sum(1 for s in group_sizes if s >= args.min_group_size)
    log(f"Pass 1 (Bemis-Murcko): {n_groups:,} distinct scaffold groups "
        f"({n_large:,} with >= {args.min_group_size} members, largest={group_sizes[0]})")

    bm_summary_path = out_dir / f"bm_group_summary_{conf}.csv"
    with open(bm_summary_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["bm_group_id", "size", "scaffold_smiles"])
        for gid, (scaf, idxs) in enumerate(sorted(groups.items(), key=lambda kv: -len(kv[1]))):
            w.writerow([gid, len(idxs), scaf])
    log(f"Wrote {bm_summary_path}")

    assignment_path = out_dir / f"bm_scaffold_assignment_{conf}.csv"
    bm_group_id_of = {}
    with open(assignment_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["compound_id", "bm_group_id", "scaffold_smiles"])
        for gid, (scaf, idxs) in enumerate(sorted(groups.items(), key=lambda kv: -len(kv[1]))):
            for i in idxs:
                bm_group_id_of[i] = gid
                w.writerow([ids[i], gid, scaf])
    log(f"Wrote {assignment_path}")

    # Pass 2: only large groups get Butina; small groups are 1 cluster each,
    # computed instantly with no similarity math at all
    large_groups = [(scaf, idxs) for scaf, idxs in groups.items() if len(idxs) >= args.min_group_size]
    log(f"Pass 2 (Tanimoto/Butina): running on {len(large_groups):,} large group(s) "
        f"covering {sum(len(idxs) for _, idxs in large_groups):,} compounds, {workers} workers")

    tasks = [([fps[i] for i in idxs], thresholds) for _, idxs in large_groups]
    large_group_results = [None] * len(large_groups)
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(_process_large_group, t): gi for gi, t in enumerate(tasks)}
        done = 0
        for fut in as_completed(futures):
            gi = futures[fut]
            large_group_results[gi] = fut.result()
            done += 1
            if done % max(1, len(large_groups) // 20) == 0 or done == len(large_groups):
                elapsed = time.time() - t0
                log(f"  pass 2: {done}/{len(large_groups)} large groups done ({elapsed:.0f}s elapsed)")

    # assemble final nested cluster assignment + summary per threshold
    for t in thresholds:
        final_cluster_id = {}  # compound index -> global cluster id string
        next_id = 0
        # small groups: one cluster each
        for scaf, idxs in groups.items():
            if len(idxs) < args.min_group_size:
                for i in idxs:
                    final_cluster_id[i] = next_id
                next_id += 1
        # large groups: however many sub-clusters Butina found at this threshold
        for (scaf, idxs), result in zip(large_groups, large_group_results):
            local_cluster_id = result[t]
            local_to_global = {}
            for local_i, i in enumerate(idxs):
                lc = int(local_cluster_id[local_i])
                if lc not in local_to_global:
                    local_to_global[lc] = next_id
                    next_id += 1
                final_cluster_id[i] = local_to_global[lc]

        sizes = np.bincount(list(final_cluster_id.values()))
        n_clusters = len(sizes)
        n_singletons = int((sizes == 1).sum())

        assign_path = out_dir / f"nested_cluster_assignment_{conf}_{t:.2f}.csv"
        with open(assign_path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["compound_id", "cluster_id"])
            for i, cid in enumerate(ids):
                w.writerow([cid, final_cluster_id[i]])

        log(f"  threshold {t:.2f}: {n_clusters:,} total clusters "
            f"({n_singletons:,} singletons, largest={int(sizes.max())})")

        summary_row = {"conformation": conf, "threshold": f"{t:.2f}", "n_compounds": n,
                       "n_clusters": n_clusters, "n_singletons": n_singletons,
                       "largest_cluster_size": int(sizes.max()), "mean_cluster_size": f"{sizes.mean():.3f}",
                       "median_cluster_size": f"{float(np.median(sizes)):.1f}"}
        summary_path = out_dir / f"nested_cluster_summary_{conf}.csv"
        write_header = not summary_path.exists() or t == thresholds[0]
        with open(summary_path, "a" if not write_header else "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(summary_row.keys()))
            if write_header:
                w.writeheader()
            w.writerow(summary_row)

    log(f"Wrote nested_cluster_summary_{conf}.csv and nested_cluster_assignment_{conf}_<threshold>.csv")


if __name__ == "__main__":
    main()
