#!/usr/bin/env python3
"""
Author: Anna Korda

Second filter stage, run AFTER interaction_filter_fast.py: checks each
salt-bridge-passing compound's already-selected winning pose (same pose
already reordered to the front by interaction_filter_fast.py, never a
different or re-picked one) for conformational strain, using the real
torsion-strain method from Gu, Smith, Yang, Irwin, Shoichet, "Ligand Strain
Energy in Large Library Docking", J Chem Inf Model 2021 (see
interaction_filter_references.txt and torsion_lib.py for the real citation
and implementation notes -- verified against the paper authors' own
published reference output before use, not just "should work").

A pose is "strained" (filtered out) if:
    total_TEU  >= --total_cutoff   (default 7.0, the paper's own DUD-E
                                     40-target general-purpose default)
 OR single_TEU >= --single_cutoff  (default 1.8, same source)

Molecules are never modified. Ring/aromaticity perception (needed for the
torsion library's aromatic SMARTS to match -- our SDFs are Kekulized, no
file-level aromatic bond type, unlike the paper's own mol2 input) is applied
via RDKit's partial sanitization (SANITIZE_SYMMRINGS | SANITIZE_SETAROMATICITY
only -- no valence correction, no H adjustment, no coordinate changes,
verified locally: coordinates/charges/atom-count stay byte-identical) to an
in-memory copy only; nothing is ever written back to the original files.

--dry_run: computes and reports real pass/fail/error counts without writing
the survivors tar (still writes the full per-compound CSV report, since
that's the "how many pass" answer, not the expensive part).

Same chunked/resumable/parallel pattern as interaction_filter_fast.py.

Usage:
    python strain_filter.py --conformation c5 \\
        --vs_results_dir vs_results \\
        --torsion_lib_xml data/TL_2.1_VERSION_6.xml \\
        --workers 32 --dry_run
"""
import argparse
import csv
import os
import shutil
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from rdkit import Chem

from torsion_lib import load_torsion_library, compute_strain

COPY_BUF = 1024 * 1024 * 4
# Only ring perception + aromaticity flagging -- everything else about the
# molecule (coordinates, formal charges, atom count, bond orders) is left
# exactly as parsed. Needed because our SDFs are Kekulized (no file-level
# aromatic bond type), so the torsion library's lowercase-aromatic SMARTS
# (e.g. [cX3], [n]) would otherwise silently match nothing -- verified: 0
# matches without this, correct matches with it, coordinates/charges
# confirmed byte-identical before/after.
SANITIZE_OPS = Chem.SANITIZE_SYMMRINGS | Chem.SANITIZE_SETAROMATICITY

# Module-level, populated once in main() before the ProcessPoolExecutor is
# created. On Linux, ProcessPoolExecutor forks (the platform default), so
# every worker inherits this already-parsed XML via copy-on-write memory --
# no need to re-parse the same 2.6MB file separately in each of N workers.
_TORSION_LIB_ROOT = None


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


def load_job_list(passed_dir, poses_dir):
    """One entry per salt-bridge-passing compound: (compound_id, real SDF
    path, winning pose index). Real path comes from scores.csv (same
    technique as extract_review_set.py) -- avoids scanning the large
    interactions_passed tar. Winning pose index comes from
    reordered_compounds_<conf>.csv when present, else 0 (pose 1 was
    already the winner, per interaction_filter_fast.py's own logic)."""
    passing_ids = set()
    with open(next(passed_dir.glob("interaction_filter_report_*.csv"))) as f:
        for row in csv.DictReader(f):
            if row["passed"] == "True":
                passing_ids.add(row["compound_id"])

    id_to_path = {}
    with open(poses_dir / "scores.csv") as f:
        for row in csv.DictReader(f):
            if row["compound_id"] in passing_ids:
                id_to_path[row["compound_id"]] = row["path"]

    winning_idx = {}
    reordered_csv = next(passed_dir.glob("reordered_compounds_*.csv"), None)
    if reordered_csv:
        with open(reordered_csv) as f:
            for row in csv.DictReader(f):
                if row["compound_id"] in passing_ids:
                    winning_idx[row["compound_id"]] = int(row["winning_pose_index_0based"])

    jobs = []
    missing = 0
    for compound_id in sorted(passing_ids):
        if compound_id not in id_to_path:
            missing += 1
            continue
        jobs.append((compound_id, poses_dir / id_to_path[compound_id],
                     winning_idx.get(compound_id, 0)))
    if missing:
        log(f"WARNING: {missing} passing compound(s) had no path in scores.csv -- skipped")
    return jobs


def prepare_mol_for_strain(sdf_path, pose_idx):
    """Reads the real file, extracts the winning pose's own text block
    (never a different one), builds an RDKit mol from it (sanitize=False,
    removeHs=False -- same as the reference tool's own choice), applies
    the targeted partial sanitization above. Returns None if the pose
    can't be parsed."""
    text = sdf_path.read_text()
    blocks = _split_sdf_blocks(text)
    if not (0 <= pose_idx < len(blocks)):
        return None
    mol = Chem.MolFromMolBlock(blocks[pose_idx], sanitize=False, removeHs=False)
    if mol is None or mol.GetNumConformers() == 0:
        return None
    Chem.SanitizeMol(mol, sanitizeOps=SANITIZE_OPS, catchErrors=True)
    return mol


def _process_chunk(chunk_jobs, chunk_idx, total_cutoff, single_cutoff,
                    parts_dir_str, dry_run):
    parts_dir = Path(parts_dir_str)
    report_path = parts_dir / f"report_{chunk_idx:05d}.csv"
    part_path = parts_dir / f"part_{chunk_idx:05d}.tar.gz"
    done_marker = part_path if not dry_run else report_path
    if report_path.exists() and (dry_run or part_path.exists()):
        return chunk_idx, None, None

    stage_dir = parts_dir / f".stage_{chunk_idx}"
    if not dry_run:
        stage_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for compound_id, sdf_path, pose_idx in chunk_jobs:
        try:
            mol = prepare_mol_for_strain(sdf_path, pose_idx)
            if mol is None:
                rows.append({"compound_id": compound_id, "total_TEU": "", "single_TEU": "",
                             "flagged": "", "n_torsions": "", "strained": "",
                             "error": "could not parse winning pose"})
                continue
            total, single, flagged, n = compute_strain(mol, _TORSION_LIB_ROOT)
            if total is None or flagged:
                err = "no torsions matched (rigid)" if total is None else "flagged (unreliable estimate)"
                rows.append({"compound_id": compound_id, "total_TEU": "" if total is None else f"{total:.4f}",
                             "single_TEU": "" if single is None else f"{single:.4f}",
                             "flagged": flagged, "n_torsions": n, "strained": "", "error": err})
                continue
            strained = (total >= total_cutoff) or (single >= single_cutoff)
            rows.append({"compound_id": compound_id, "total_TEU": f"{total:.4f}",
                         "single_TEU": f"{single:.4f}", "flagged": False, "n_torsions": n,
                         "strained": strained, "error": ""})
            if not strained and not dry_run:
                shutil.copy2(sdf_path, stage_dir / f"{compound_id}.sdf")
        except Exception as exc:
            rows.append({"compound_id": compound_id, "total_TEU": "", "single_TEU": "",
                         "flagged": "", "n_torsions": "", "strained": "", "error": str(exc)})

    with open(report_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["compound_id", "total_TEU", "single_TEU",
                                          "flagged", "n_torsions", "strained", "error"])
        w.writeheader()
        w.writerows(rows)

    if not dry_run:
        passing = [p.name for p in stage_dir.glob("*.sdf")]
        tmp_part = parts_dir / f"part_{chunk_idx:05d}.tar.gz.partial-{os.getpid()}"
        if passing:
            list_file = parts_dir / f"list_{chunk_idx:05d}.txt"
            list_file.write_text("\n".join(passing) + "\n")
            r = subprocess.run(["tar", "-czf", str(tmp_part), "-C", str(stage_dir), "-T", str(list_file)],
                                capture_output=True, text=True)
            list_file.unlink(missing_ok=True)
        else:
            r = subprocess.run(["tar", "-czf", str(tmp_part), "--files-from", os.devnull],
                                capture_output=True, text=True)
        shutil.rmtree(stage_dir, ignore_errors=True)
        if r.returncode != 0:
            tmp_part.unlink(missing_ok=True)
            return chunk_idx, r.stderr, None
        tmp_part.rename(part_path)

    n_pass = sum(1 for row in rows if row["strained"] is False)
    n_strained = sum(1 for row in rows if row["strained"] is True)
    n_err = sum(1 for row in rows if row["error"])
    return chunk_idx, None, (n_pass, n_strained, n_err)


def _init_worker(xml_path):
    global _TORSION_LIB_ROOT
    if _TORSION_LIB_ROOT is None:
        _TORSION_LIB_ROOT = load_torsion_library(xml_path)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--conformation", required=True)
    ap.add_argument("--vs_results_dir", required=True)
    ap.add_argument("--torsion_lib_xml", required=True)
    ap.add_argument("--total_cutoff", type=float, default=7.0)
    ap.add_argument("--single_cutoff", type=float, default=1.8)
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--chunk_size", type=int, default=500)
    ap.add_argument("--dry_run", action="store_true",
                     help="report real pass/strained/error counts and write the per-compound "
                          "CSV report, but skip writing the (expensive) survivors tar")
    args = ap.parse_args()

    conf = args.conformation
    vs_dir = Path(args.vs_results_dir)
    passed_dir = vs_dir / f"interactions_passed_{conf}"
    poses_dir = vs_dir / f"ad4_redock_{conf}"
    if not passed_dir.is_dir():
        raise SystemExit(f"{passed_dir} not found -- run interaction_filter_fast.py for this conformation first")

    global _TORSION_LIB_ROOT
    log(f"Loading torsion library from {args.torsion_lib_xml} ...")
    _TORSION_LIB_ROOT = load_torsion_library(args.torsion_lib_xml)
    log(f"Torsion library loaded ({len(_TORSION_LIB_ROOT.findall('hierarchyClass'))} hierarchy classes)")

    log(f"Building job list from {passed_dir} ...")
    jobs = load_job_list(passed_dir, poses_dir)
    n_total = len(jobs)
    log(f"{n_total:,} salt-bridge-passing compounds to check for strain "
        f"(cutoffs: total>={args.total_cutoff} OR single>={args.single_cutoff} => strained)")
    if args.dry_run:
        log("DRY RUN: will report real counts, will NOT write the survivors tar")

    out_root = vs_dir / f"strain_filtered_{conf}"
    parts_dir = out_root / ".parts"
    parts_dir.mkdir(parents=True, exist_ok=True)

    chunks = [jobs[i:i + args.chunk_size] for i in range(0, n_total, args.chunk_size)]
    n_chunks = len(chunks)
    workers = args.workers or len(os.sched_getaffinity(0))

    def _chunk_done(i):
        cid = f"{i:05d}"
        report_ok = (parts_dir / f"report_{cid}.csv").exists()
        if args.dry_run:
            return report_ok
        return report_ok and (parts_dir / f"part_{cid}.tar.gz").exists()

    have = {i for i in range(n_chunks) if _chunk_done(i)}
    todo = [i for i in range(n_chunks) if i not in have]
    log(f"{n_chunks} chunks of up to {args.chunk_size}, {len(have)} already done, {workers} workers")

    total_pass = total_strained = total_err = 0
    t0 = time.time()
    if todo:
        done_count = 0
        with ProcessPoolExecutor(max_workers=workers, initializer=_init_worker,
                                  initargs=(args.torsion_lib_xml,)) as pool:
            futures = {pool.submit(_process_chunk, chunks[i], i, args.total_cutoff,
                                    args.single_cutoff, str(parts_dir), args.dry_run): i
                       for i in todo}
            for fut in as_completed(futures):
                idx, err, counts = fut.result()
                if err:
                    print(err, file=sys.stderr)
                    raise SystemExit(f"chunk {idx} failed")
                done_count += 1
                if counts:
                    n_pass, n_strained, n_err = counts
                    total_pass += n_pass
                    total_strained += n_strained
                    total_err += n_err
                rate = done_count / (time.time() - t0) if time.time() > t0 else 0
                eta = (len(todo) - done_count) / rate if rate else float("inf")
                log(f"  chunk {done_count}/{len(todo)} done ({time.time() - t0:.0f}s elapsed, "
                    f"ETA {eta:.0f}s) -- pass={total_pass} strained={total_strained} err={total_err} "
                    f"(of {total_pass + total_strained + total_err:,}/{n_total:,})")

    log("all chunks done -- merging per-compound reports ...")
    report_path = out_root / f"strain_filter_report_{conf}.csv"
    grand_pass = grand_strained = grand_err = 0
    with open(report_path, "w", newline="") as out_f:
        w = csv.writer(out_f)
        w.writerow(["compound_id", "total_TEU", "single_TEU", "flagged", "n_torsions", "strained", "error"])
        for i in range(n_chunks):
            with open(parts_dir / f"report_{i:05d}.csv") as in_f:
                for row in csv.DictReader(in_f):
                    w.writerow([row["compound_id"], row["total_TEU"], row["single_TEU"],
                               row["flagged"], row["n_torsions"], row["strained"], row["error"]])
                    if row["strained"] == "False":
                        grand_pass += 1
                    elif row["strained"] == "True":
                        grand_strained += 1
                    else:
                        grand_err += 1

    summary_path = out_root / f"run_summary_{conf}.txt"
    with open(summary_path, "w") as f:
        f.write(f"Strain filter run summary -- conformation {conf}{' (DRY RUN)' if args.dry_run else ''}\n")
        f.write(f"Finished: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Criterion: strained if total_TEU>={args.total_cutoff} OR single_TEU>={args.single_cutoff}\n")
        f.write(f"Torsion library: {args.torsion_lib_xml}\n\n")
        f.write(f"Salt-bridge-passing compounds checked: {n_total:,}\n")
        f.write(f"Passed strain filter (kept):  {grand_pass:,}\n")
        f.write(f"Strained (filtered out):      {grand_strained:,}\n")
        f.write(f"Errored/unscoreable:          {grand_err:,}\n")
        f.write(f"\nFull per-compound report: {report_path.name}\n")
    log(f"Wrote {report_path}")
    log(f"Wrote {summary_path}")

    if args.dry_run:
        shutil.rmtree(parts_dir, ignore_errors=True)
        log(f"DRY RUN TOTAL: pass={grand_pass:,} strained={grand_strained:,} error={grand_err:,} "
            f"out of {n_total:,} -- no survivors tar written")
        return

    tar_path = out_root / f"strain_filtered_{conf}.tar.gz"
    tmp_out = f"{tar_path}.partial"
    with open(tmp_out, "wb") as out_f:
        for i in range(n_chunks):
            with open(parts_dir / f"part_{i:05d}.tar.gz", "rb") as in_f:
                shutil.copyfileobj(in_f, out_f, COPY_BUF)
    os.rename(tmp_out, tar_path)
    shutil.rmtree(parts_dir, ignore_errors=True)
    size_gb = os.path.getsize(tar_path) / (1024 ** 3)
    log(f"Wrote {tar_path} ({size_gb:.2f} GB)")
    log(f"TOTAL: pass={grand_pass:,} strained={grand_strained:,} error={grand_err:,} "
        f"out of {n_total:,} compounds ({time.time() - t0:.0f}s)")
    log(f"extract poses with:  tar -ixzf {tar_path.name}   (-i required: concatenated multi-part gzip)")


if __name__ == "__main__":
    main()
